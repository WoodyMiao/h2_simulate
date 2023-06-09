#!/usr/bin/env python

import os
import logging
import argparse
from itertools import product

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm, chi2_contingency

from rpy2 import robjects
from rpy2.robjects.packages import importr
from pysnptools.snpreader import Bed, SnpData

parser = argparse.ArgumentParser(description='Simulate genotypes by HapSim and estimate h2 by HESS and EHE')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--gene-bed', type=str, required=True, help='A BED file without header; the fourth column is IDs')
parser.add_argument('--gene', type=str, required=True, help='A gene in the BED file to be simulated')
parser.add_argument('--vcf-ref', type=str, required=True, help='A phased VCF file of a reference population')
parser.add_argument('--hm3-snp', type=str, required=True, help='A file of a list of HapMap3 SNP ids')

args = parser.parse_args()
gene = args.gene
h2g = 1e-3
alpha = -0.25
maf_min = 0.01
prevalence = [0.1, 0.01]
pqtl_lst = [0.01, 0.1, 1]
n_gwa = 20000
n_case = int(n_gwa / 2)
n_rep = 100
n_ld_str = ['5e2', '2e3']
n_ld_int = [int(float(n)) for n in n_ld_str]
n_pop0 = int(1e7)  # the sample size used for realizing beta
n_sample_once = int(1e5)  # the sample size for each subprocess to sample
n_ld_total = np.sum(n_ld_int) * n_rep  # 2.5e5
tot_n_case = n_case * n_rep  # 1e6
gwa_prefix = 'plink_gwa_prvl0.01_pqtl0.1_alpha-0.25'


def sample_once(_):
    X_ = np.random.multivariate_normal(np.zeros(m), C, (2, n_sample_once))  # 2 * n_sample_once * m
    X_ = np.int8(X_ < percent_point)  # 2 * n_sample_once * m
    return X_[0] + X_[1]  # n_sample_once * m; int8


logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info(f'Start filtering the VCF file of {gene}')
# Read the BED file and the VCF file
gene_bed = pd.read_csv(args.gene_bed, sep='\t', header=None, index_col=[3]).loc[gene]
vcf = pd.read_csv(args.vcf_ref, sep='\t', comment='#', header=None)
vcf_header = vcf.columns.to_list()
vcf_header[:9] = ['CHR', 'BP', 'SNP', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
vcf.columns = vcf_header

# Extract SNPs inside genes
vcf = vcf.iloc[np.where((vcf.loc[:, 'CHR'] == gene_bed[0]) &
                        (vcf.loc[:, 'BP'] >= gene_bed[1]) &
                        (vcf.loc[:, 'BP'] <= gene_bed[2]))[0]]

# Extract reference haplotypes and SNP information from the VCF file
vcf = vcf[(~vcf[['CHR', 'BP']].duplicated(keep=False)) &
          (~vcf['SNP'].duplicated(keep=False)) & vcf['SNP'].str.match(r'^rs\d+$') &
          vcf['REF'].isin(['A', 'C', 'G', 'T']) & vcf['ALT'].isin(['A', 'C', 'G', 'T']) &
          np.all(vcf.loc[:, 9:].isin(['0|0', '0|1', '1|0', '1|1']), axis=1)]
ref_haplo = np.concatenate((vcf.loc[:, 9:].applymap(lambda x: x[0]).values.astype(np.int8).T,
                            vcf.loc[:, 9:].applymap(lambda x: x[2]).values.astype(np.int8).T))

# Filter by MAF and make a PLINK bim file
allele_frq = ref_haplo.mean(axis=0)
extract = (allele_frq > maf_min) & (allele_frq < 1 - maf_min)
ref_haplo = ref_haplo[:, extract]
n_haplo, m = ref_haplo.shape
allele_frq = allele_frq[extract]
percent_point = norm.ppf(allele_frq)  # m

vcf = vcf.loc[extract, ['CHR', 'SNP', 'QUAL', 'BP', 'ALT', 'REF']]
plink_bim = vcf.rename({'QUAL': 'CM', 'ALT': 'A1', 'REF': 'A2'}, axis=1)
plink_bim['CM'] = 0
plink_bim['hm3SNP'] = plink_bim['SNP'].isin(np.loadtxt(args.hm3_snp, dtype=str))
del vcf, vcf_header, extract

# Write a one-column HapMap3 SNP list
dir_gene = f'{args.out_dir}/{gene}'
os.system(f'mkdir -p {dir_gene}')
plink_bim.loc[plink_bim['hm3SNP'], 'SNP'].to_csv(f'{dir_gene}/plink.hm3snp', header=False, index=False)

# Write SNP files for LDER
plink_bim[['SNP', 'A1', 'A2']].to_csv(f'{dir_gene}/lder.snp', sep='\t', index=False)
plink_bim.loc[plink_bim['hm3SNP'], ['SNP', 'A1', 'A2']].to_csv(f'{dir_gene}/lder.hm3snp', sep='\t', index=False)

# Write a plink BIM file
plink_bim[['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']].to_csv(
    f'{dir_gene}/plink.bim', sep='\t', index=False, header=False)

# Write a region file for KGGSEE
with open(f'{dir_gene}/region.kggsee', 'w') as o:
    print(f'{gene_bed[0]}\t{gene_bed[1]}\t{gene_bed[2]}\t{gene}', file=o)

# Write a region file for HESS and LDER
with open(f'{dir_gene}/region.hess', 'w') as o:
    print('chr\tstart\tstop', file=o)
    print(f'chr{gene_bed[0]}\t{gene_bed[1]}\t{gene_bed[2]}', file=o)
os.system(f'cp {dir_gene}/region.hess {dir_gene}/region.lder')

# Write the region file for LDAK
with open(f'{dir_gene}/region.ldak', 'w') as o:
    print(f'{gene}\t{gene_bed[0]}\t{gene_bed[1]}\t{gene_bed[2]}', file=o)

logging.info(f'Start sampling pop0 of {gene}')
# Calculate an MVN covariance matrix using HapSim
importr('hapsim')
haplodata = robjects.r('haplodata')
haplodata = haplodata(robjects.r.matrix(robjects.IntVector(ref_haplo.T.reshape(-1)), nrow=n_haplo))
C = np.array(dict(zip(haplodata.names, list(haplodata)))['cor'])  # m * m

X_pop0 = np.concatenate(list(map(sample_once, range(int(n_pop0 / n_sample_once)))))
mean_pop0 = X_pop0.mean(axis=0)  # m
std_pop0 = X_pop0.std(axis=0)  # m
X_pop0_stdz = (X_pop0 - mean_pop0) / std_pop0

frq_pop0 = mean_pop0 / 2
beta_var = (frq_pop0 * (1 - frq_pop0)) ** (1 + alpha)
thresholds = norm.isf(prevalence)

logging.info(f'Start realizing effect sizes of {gene}')
X_idx_used = list()
X_case_dict = dict()
X_ctrl_dict = dict()
realized_beta = dict()
n_collected = pd.DataFrame(dtype=int, index=prevalence, columns=pqtl_lst)
for (th, prvl), pqtl in product(zip(thresholds, prevalence), pqtl_lst):
    # Realize beta: sample; then, scale to fit target h2g
    m_qtl = int(np.ceil(m * pqtl))
    idx_qtl = np.random.choice(m, size=m_qtl, replace=False)
    beta = np.zeros(m)
    beta[idx_qtl] = np.random.multivariate_normal(np.zeros(m_qtl), np.diag(beta_var[idx_qtl]))  # m_qtl

    genetic_eff = X_pop0_stdz @ beta
    genetic_eff_std = genetic_eff.std()
    scale_factor = h2g ** 0.5 / genetic_eff_std
    realized_beta[(prvl, pqtl)] = beta * scale_factor

    # Collect individuals from pop0 for GWAS samples
    y_norm = genetic_eff * scale_factor + np.random.normal(0, np.sqrt(1 - h2g), n_pop0)
    y = y_norm > th
    i_case = np.where(y)[0]
    n_collected.loc[prvl, pqtl] = i_case.shape[0]
    X_case_dict[(prvl, pqtl)] = [X_pop0[i_case]]
    X_idx_used.append(i_case)

    i_ctrl = np.where(~y)[0][:i_case.shape[0]]
    X_ctrl_dict[(prvl, pqtl)] = [X_pop0[i_ctrl]]
    X_idx_used.append(i_ctrl)

del X_pop0_stdz

logging.info(f'Start writing PLINK bfiles of reference LD samples of {gene}')
# Collect remaining individuals for LD samples
X_idx_used = np.unique(np.concatenate(X_idx_used))
X_pop0 = np.delete(X_pop0, X_idx_used, axis=0)[:n_ld_total]

iid = np.arange(1, 1 + tot_n_case * 2 + n_ld_total).astype(str)[:, None]
iid = np.concatenate((iid, iid), axis=1)  # n_pop * 2

X_ld_lst = list()
iid_ld_lst = list()
ld_prefixes = list()
start = 0
for a, n_ld in zip(n_ld_str, n_ld_int):
    end = start + n_rep * n_ld
    X_ld_lst.append(X_pop0[start:end].reshape(n_rep, n_ld, m))
    iid_ld_lst.append(iid[start:end].reshape(n_rep, n_ld, 2))
    ld_prefixes.append(f'plink_ld{a}')
    start = end
del X_pop0

for j in range(n_rep):
    dir_rep = f'{dir_gene}/rep{j}'
    os.system(f'mkdir -p {dir_rep}')
    for X_ld, iid_ld, ld_prefix in zip(X_ld_lst, iid_ld_lst, ld_prefixes):
        Bed.write(f'{dir_rep}/{ld_prefix}.bed', count_A1=True, _require_float32_64=False,
                  snpdata=SnpData(val=X_ld[j], iid=iid_ld[j], sid=plink_bim['SNP'],
                                  pos=plink_bim[['CHR', 'CM', 'BP']], _require_float32_64=False))
        os.remove(f'{dir_rep}/{ld_prefix}.bim')
        os.symlink('../plink.bim', f'{dir_rep}/{ld_prefix}.bim')
del X_ld_lst

logging.info(f'Continue sampling for association tests of {gene}')
n_round = 1
n_each_round = int(1e6)
while np.any(n_collected < tot_n_case):
    logging.info(f'Start round {n_round} of sampling for {gene}')
    X = np.concatenate(list(map(sample_once, range(int(n_each_round / n_sample_once)))))
    X_stdz = (X - mean_pop0) / std_pop0

    for (th, prvl), pqtl in product(zip(thresholds, prevalence), pqtl_lst):
        if n_collected.loc[prvl, pqtl] < tot_n_case:
            y_norm = X_stdz @ realized_beta[(prvl, pqtl)] + np.random.normal(0, np.sqrt(1 - h2g), n_each_round)
            y = y_norm > th
            i_case = np.where(y)[0]
            n_collected.loc[prvl, pqtl] += i_case.shape[0]
            X_case_dict[(prvl, pqtl)].append(X[i_case])
            i_ctrl = np.where(~y)[0][:i_case.shape[0]]
            X_ctrl_dict[(prvl, pqtl)].append(X[i_ctrl])
            logging.info(f'prevalence={prvl}, pqtl={pqtl}, {n_collected.loc[prvl, pqtl]} cases sampled.')

    n_round += 1

logging.info(f'Start writing PLINK bfiles of GWAS samples of {gene}')
# Divide individuals into n_rep samples for each parameter set
X_gwa_dict = dict()
for prvl, pqtl in product(prevalence, pqtl_lst):
    X_case_dict[(prvl, pqtl)] = np.concatenate(X_case_dict[(prvl, pqtl)])[: tot_n_case].reshape(n_rep, n_case, m)
    X_ctrl_dict[(prvl, pqtl)] = np.concatenate(X_ctrl_dict[(prvl, pqtl)])[: tot_n_case].reshape(n_rep, n_case, m)
    X_gwa_dict[(prvl, pqtl)] = np.concatenate((X_case_dict[(prvl, pqtl)], X_ctrl_dict[(prvl, pqtl)]), axis=1)

iid_gwa = iid[:n_gwa]
fam_gwa = pd.DataFrame(0, index=np.arange(n_gwa), columns=np.arange(6))
fam_gwa.iloc[:, :2] = iid_gwa
fam_gwa.iloc[:n_case, 5] = 2
fam_gwa.iloc[n_case:, 5] = 1
fam_gwa.to_csv(f'{dir_gene}/plink_gwa.fam', sep=' ', index=False, header=False)
for j in range(n_rep):
    dir_rep = f'{dir_gene}/rep{j}'
    Bed.write(f'{dir_rep}/{gwa_prefix}.bed', count_A1=True, _require_float32_64=False,
              snpdata=SnpData(val=X_gwa_dict[(0.01, 0.1)][j], iid=iid_gwa, sid=plink_bim['SNP'],
                              pos=plink_bim[['CHR', 'CM', 'BP']], _require_float32_64=False))
    os.remove(f'{dir_rep}/{gwa_prefix}.bim')
    os.remove(f'{dir_rep}/{gwa_prefix}.fam')
    os.symlink('../plink.bim', f'{dir_rep}/{gwa_prefix}.bim')
    os.symlink('../plink_gwa.fam', f'{dir_rep}/{gwa_prefix}.fam')


def write_sumstats(z_, stem):
    p_ = norm.sf(np.abs(z_)) * 2

    for j_ in range(n_rep):
        prefix = f'{dir_gene}/rep{j_}/{stem}'

        # Write a summary statistic file for KGGSEE, HESS, and LDSC
        sumstat = plink_bim.copy()
        sumstat['Z'] = z_[j_]
        sumstat['N'] = n_gwa
        sumstat['P'] = p_[j_]
        sumstat[['CHR', 'BP', 'P', 'SNP', 'A1', 'A2', 'Z', 'N']] \
            .to_csv(f'{prefix}.sumstat.gz', sep='\t', index=False)

        # Write a summary statistic file for LDER
        sumstat.rename({'SNP': 'snp', 'CHR': 'chr', 'A1': 'a0', 'A2': 'a1', 'Z': 'z'}, axis=1)[
            ['snp', 'chr', 'a0', 'a1', 'z']].to_csv(
            f'{prefix}.lder.sumstat.gz', sep='\t', index=False)

        # Write a summary statistic file for LDAK
        sumstat.rename({'SNP': 'Predictor', 'N': 'n'}, axis=1)[['Predictor', 'A1', 'A2', 'n', 'Z']] \
            .to_csv(f'{prefix}.ldak.sumstat', sep='\t', index=False)


y_stdz = np.concatenate((np.ones(n_case), -np.ones(n_case)))[:, None]  # n_gwa * 1
y_unit = np.concatenate((np.ones(n_case), np.zeros(n_case)))[:, None]  # n_gwa * 1
for prvl, pqtl in product(prevalence, pqtl_lst):
    logging.info(f'Performing association tests of {gene} for prvl={prvl} pqtl={pqtl} ...')
    X_stdz = (X_gwa_dict[(prvl, pqtl)] - mean_pop0) / std_pop0  # n_rep * n_gwa * m

    logging.info(f'Performing linear regressions of {gene} for prvl={prvl} pqtl={pqtl} ...')
    z = np.swapaxes(X_stdz, 1, 2) @ y_stdz / np.sqrt(n_gwa)  # n_rep * m * 1
    write_sumstats(z, f'prvl{prvl}_pqtl{pqtl}_alpha{alpha}.linear')

    logging.info(f'Performing 1df chi-square allelic test of {gene} for prvl={prvl} pqtl={pqtl} ...')
    #       A1  A2
    # case   a   b
    # ctrl   c   d
    a = X_case_dict[(prvl, pqtl)].sum(axis=1)  # n_rep * m
    b = n_gwa * 2 - a
    c = X_ctrl_dict[(prvl, pqtl)].sum(axis=1)  # n_rep * m
    d = n_gwa * 2 - c
    z = np.empty((n_rep, m, 1))
    for i in range(n_rep):
        for j in range(m):
            table = np.array([[a[i, j], b[i, j]], [c[i, j], d[i, j]]])
            z_abs = chi2_contingency(table)[0] ** 0.5
            o_r = (table[0, 0] / table[1, 0]) / (table[0, 1] / table[1, 1])
            if o_r > 1:
                z[i, j, 0] = z_abs
            else:
                z[i, j, 0] = -z_abs
    write_sumstats(z, f'prvl{prvl}_pqtl{pqtl}_alpha{alpha}.chi2')

    logging.info(f'Performing logistic regression of {gene} with prvl={prvl}, pqtl={pqtl} ...')
    z = np.empty((n_rep, m, 1))
    for i in range(n_rep):
        for j in range(m):
            logit_res = sm.Logit(y_unit, X_stdz[i, :, j]).fit()
            z[i, j, 0] = logit_res.tvalues[0]
    write_sumstats(z, f'prvl{prvl}_pqtl{pqtl}_alpha{alpha}.logit')

logging.info(f'Done simulation of {m} SNPs in {gene}')
