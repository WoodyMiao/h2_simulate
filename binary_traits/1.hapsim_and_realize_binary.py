#!/usr/bin/env python

import os
import logging
import argparse
from itertools import product
from multiprocessing import Pool

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm, chi2_contingency

from rpy2 import robjects
from rpy2.robjects.packages import importr
from pysnptools.snpreader import Bed, SnpData

parser = argparse.ArgumentParser(description='Simulate genotypes by HapSim and estimate h2 by HESS and EHE')
parser.add_argument('--nt', type=int, default=1, help='Number of processes running in parallel')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--gene-bed', type=str, required=True, help='A BED file without header; the fourth column is IDs')
parser.add_argument('--vcf-ref', type=str, required=True, help='A phased VCF file of a reference population')
parser.add_argument('--hm3-snp', type=str, required=True, help='A file of a list of HapMap3 SNP ids')
parser.add_argument('--n-gwa', type=float, default=2e4, help='Size of a sample for association tests')
parser.add_argument('--n-ld', type=str, default='5e2,2e3', help='Sizes of samples for LD panels')
parser.add_argument('--n-rep', type=float, default=100, help='Size of a sample for association tests')
parser.add_argument('--n-pop-beta', type=float, default=1e7,
                    help='The sample size used for realizing beta; must be a multiple of 1e5')
parser.add_argument('--h2g', type=float, default=1e-3, help='Target h2g for the simulation')
parser.add_argument('--maf-min', type=float, default=0.01, help='Ignore SNPs with MAF < MAF_MIN')
parser.add_argument('--neg-alpha', type=float, default=0.25, help='Power value in the LDAK-Thin Model')
parser.add_argument('--pqtl-vals', type=str, default='1,0.1,0.01', help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--prevalence-vals', type=str, default='0.1,0.01',
                    help='prevalence values for the liability threshold modeling of dichotamous phenotypes')

args = parser.parse_args()
h2g = args.h2g
alpha = -args.neg_alpha
maf_min = args.maf_min
n_gwa = int(args.n_gwa)
n_case = int(n_gwa / 2)
n_rep = int(args.n_rep)
n_ld_str = args.n_ld.split(',')
n_ld_int = [int(float(n)) for n in n_ld_str]
pqtl_lst = [float(a) for a in args.pqtl_vals.split(',')]
prevalence_lst = [float(a) for a in args.prevalence_vals.split(',')]
n_pop_beta = int(args.n_pop_beta)
n_ld_total = np.sum(n_ld_int) * n_rep  # 2.5e5
tot_n_case = n_case * n_rep  # 1e6


def hapsim_one_gene(gene):
    # Simulate for each region defined in the gene BED file
    importr('hapsim')

    # Get reference haplotypes of the i-th gene
    snp_idx = np.where((plink_bim.iloc[:, 0] == gene_bed.loc[gene, 0]) &
                       (plink_bim.iloc[:, 3] >= gene_bed.loc[gene, 1]) &
                       (plink_bim.iloc[:, 3] <= gene_bed.loc[gene, 2]))[0]
    ref_haplo_i = ref_haplo[:, snp_idx]
    n_haplo, m = ref_haplo_i.shape
    allele_frq_i = allele_frq[snp_idx]
    percent_point = norm.ppf(allele_frq_i)  # m
    plink_bim_i = plink_bim.iloc[snp_idx]

    if m == 0:
        logging.info(f'Skipped {m} SNPs in {gene}.')
        return 0, 0

    logging.info(f'Start simulating {m} SNPs in {gene} ...')
    # Write a one-column HapMap3 SNP list
    dir_gene = f'{args.out_dir}/{gene}'
    os.system(f'mkdir -p {dir_gene}')
    plink_bim_i.loc[plink_bim_i['hm3SNP'], 'SNP'].to_csv(f'{dir_gene}/plink.hm3snp', header=False, index=False)

    # Write SNP files for LDER
    plink_bim_i[['SNP', 'A1', 'A2']].to_csv(f'{dir_gene}/lder.snp', sep='\t', index=False)
    plink_bim_i.loc[plink_bim_i['hm3SNP'], ['SNP', 'A1', 'A2']].to_csv(f'{dir_gene}/lder.hm3snp', sep='\t', index=False)

    # Write a plink BIM file
    plink_bim_i[['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']].to_csv(
        f'{dir_gene}/plink.bim', sep='\t', index=False, header=False)

    # Write a region file for KGGSEE
    with open(f'{dir_gene}/region.kggsee', 'w') as o:
        print(f'{gene_bed.loc[gene, 0]}\t{gene_bed.loc[gene, 1]}\t{gene_bed.loc[gene, 2]}\t{gene}', file=o)

    # Write a region file for HESS and LDER
    with open(f'{dir_gene}/region.hess', 'w') as o:
        print('chr\tstart\tstop', file=o)
        print(f'chr{gene_bed.loc[gene, 0]}\t{gene_bed.loc[gene, 1]}\t{gene_bed.loc[gene, 2]}', file=o)
    os.system(f'cp {dir_gene}/region.hess {dir_gene}/region.lder')

    # Write the region file for LDAK
    with open(f'{dir_gene}/region.ldak', 'w') as o:
        print(f'{gene}\t{gene_bed.loc[gene, 0]}\t{gene_bed.loc[gene, 1]}\t{gene_bed.loc[gene, 2]}', file=o)

    logging.info(f'Start sampling pop0 of {gene} ...')
    # Calculate an MVN covariance matrix using HapSim
    haplodata = robjects.r('haplodata')
    haplodata = haplodata(robjects.r.matrix(robjects.IntVector(ref_haplo_i.T.reshape(-1)), nrow=n_haplo))
    C = np.array(dict(zip(haplodata.names, list(haplodata)))['cor'])  # m * m

    def sample_once(_):
        X_ = np.random.multivariate_normal(np.zeros(m), C, (2, n_sample_once))  # 2 * n_sample_once * m
        X_ = np.int8(X_ < percent_point)  # 2 * n_sample_once * m
        return X_[0] + X_[1]  # n_sample_once * m; int8

    n_sample_once = int(1e5)
    X_pop0 = np.concatenate(list(map(sample_once, range(int(n_pop_beta / n_sample_once)))))
    mean_pop0 = X_pop0.mean(axis=0)  # m
    std_pop0 = X_pop0.std(axis=0)  # m
    X_pop0_stdz = (X_pop0 - mean_pop0) / std_pop0

    frq_pop0 = mean_pop0 / 2
    beta_var = (frq_pop0 * (1 - frq_pop0)) ** (1 + alpha)
    thresholds = norm.isf(prevalence_lst)

    logging.info(f'Start realizing effect sizes of {gene} ...')
    X_idx_used = list()
    X_case_dict = dict()
    X_ctrl_dict = dict()
    realized_beta = dict()
    n_collected = pd.DataFrame(dtype=int, index=prevalence_lst, columns=pqtl_lst)
    for (th, prvl), pqtl in product(zip(thresholds, prevalence_lst), pqtl_lst):
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
        y_norm = genetic_eff * scale_factor + np.random.normal(0, np.sqrt(1 - h2g), n_pop_beta)
        y = y_norm > th
        i_case = np.where(y)[0]
        n_collected.loc[prvl, pqtl] = i_case.shape[0]
        X_case_dict[(prvl, pqtl)] = [X_pop0[i_case]]
        X_idx_used.append(i_case)

        i_ctrl = np.where(~y)[0][:i_case.shape[0]]
        X_ctrl_dict[(prvl, pqtl)] = [X_pop0[i_ctrl]]
        X_idx_used.append(i_ctrl)
    del X_pop0_stdz

    logging.info(f'Start writing PLINK bfiles of reference LD samples for {gene} ...')
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
                      snpdata=SnpData(val=X_ld[j], iid=iid_ld[j], sid=plink_bim_i['SNP'],
                                      pos=plink_bim_i[['CHR', 'CM', 'BP']], _require_float32_64=False))
            os.remove(f'{dir_rep}/{ld_prefix}.bim')
            os.symlink('../plink.bim', f'{dir_rep}/{ld_prefix}.bim')
    del X_ld_lst

    logging.info(f'Continue sampling for association tests of {gene}.')
    n_each_round = int(1e6)
    while np.any(n_collected < tot_n_case):
        X = np.concatenate(list(map(sample_once, range(int(n_each_round / n_sample_once)))))
        X_stdz = (X - mean_pop0) / std_pop0

        for (th, prvl), pqtl in product(zip(thresholds, prevalence_lst), pqtl_lst):
            if n_collected.loc[prvl, pqtl] < tot_n_case:
                y_norm = X_stdz @ realized_beta[(prvl, pqtl)] + np.random.normal(0, np.sqrt(1 - h2g), n_each_round)
                y = y_norm > th
                i_case = np.where(y)[0]
                n_collected.loc[prvl, pqtl] += i_case.shape[0]
                X_case_dict[(prvl, pqtl)].append(X[i_case])
                i_ctrl = np.where(~y)[0][:i_case.shape[0]]
                X_ctrl_dict[(prvl, pqtl)].append(X[i_ctrl])

    logging.info(f'Writing PLINK bfiles of GWAS samples for {gene} ...')
    # Divide individuals into n_rep samples for each parameter set
    X_gwa_dict = dict()
    for prvl, pqtl in product(prevalence_lst, pqtl_lst):
        X_case_dict[(prvl, pqtl)] = np.concatenate(X_case_dict[(prvl, pqtl)])[: tot_n_case].reshape(n_rep, n_case, m)
        X_ctrl_dict[(prvl, pqtl)] = np.concatenate(X_ctrl_dict[(prvl, pqtl)])[: tot_n_case].reshape(n_rep, n_case, m)
        X_gwa_dict[(prvl, pqtl)] = np.concatenate((X_case_dict[(prvl, pqtl)], X_ctrl_dict[(prvl, pqtl)]), axis=1)

    iid_gwa = iid[:n_gwa]
    fam_gwa = pd.DataFrame(0, index=np.arange(n_gwa), columns=np.arange(6))
    fam_gwa.iloc[:, :2] = iid_gwa
    fam_gwa.iloc[:n_case, 5] = 2
    fam_gwa.iloc[n_case:, 5] = 1
    gwa_prefix = 'plink_gwa_prvl0.01_pqtl0.1_alpha-0.25'
    fam_gwa.to_csv(f'{dir_gene}/plink_gwa.fam', sep=' ', index=False, header=False)
    for j in range(n_rep):
        dir_rep = f'{dir_gene}/rep{j}'
        Bed.write(f'{dir_rep}/{gwa_prefix}.bed', count_A1=True, _require_float32_64=False,
                  snpdata=SnpData(val=X_gwa_dict[(prevalence_lst[0], pqtl_lst[0])][j], _require_float32_64=False,
                                  iid=iid_gwa, sid=plink_bim_i['SNP'], pos=plink_bim_i[['CHR', 'CM', 'BP']]))
        os.remove(f'{dir_rep}/{gwa_prefix}.bim')
        os.remove(f'{dir_rep}/{gwa_prefix}.fam')
        os.symlink('../plink.bim', f'{dir_rep}/{gwa_prefix}.bim')
        os.symlink('../plink_gwa.fam', f'{dir_rep}/{gwa_prefix}.fam')

    def write_sumstats(z_, stem):
        p_ = norm.sf(np.abs(z_)) * 2

        for j_ in range(n_rep):
            prefix = f'{dir_gene}/rep{j_}/{stem}'

            # Write a summary statistic file for KGGSEE, HESS, and LDSC
            sumstat = plink_bim_i.copy()
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

        return None

    logging.info(f'Start performing association tests of {gene}.')
    y_stdz = np.concatenate((np.ones(n_case), -np.ones(n_case)))[:, None]  # n_gwa * 1
    y_unit = np.concatenate((np.ones(n_case), np.zeros(n_case)))[:, None]  # n_gwa * 1
    for prvl, pqtl in product(prevalence_lst, pqtl_lst):
        X_stdz = (X_gwa_dict[(prvl, pqtl)] - mean_pop0) / std_pop0  # n_rep * n_gwa * m
        z = np.swapaxes(X_stdz, 1, 2) @ y_stdz / np.sqrt(n_gwa)  # n_rep * m * 1
        write_sumstats(z, f'prvl{prvl}_pqtl{pqtl}_alpha{alpha}.linear')

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

        z = np.empty((n_rep, m, 1))
        for i in range(n_rep):
            for j in range(m):
                logit_res = sm.Logit(y_unit, X_stdz[i, :, j]).fit()
                z[i, j, 0] = logit_res.tvalues[0]
        write_sumstats(z, f'prvl{prvl}_pqtl{pqtl}_alpha{alpha}.logit')

    logging.info(f'Done simulation of {m} SNPs in {gene}.')
    return m, plink_bim_i['hm3SNP'].sum()


def get_idx_inside_gene(gene_):
    # return of indices of SNPs inside genes
    return np.where((vcf.loc[:, 'CHR'] == gene_bed.loc[gene_, 0]) &
                    (vcf.loc[:, 'BP'] >= gene_bed.loc[gene_, 1]) &
                    (vcf.loc[:, 'BP'] <= gene_bed.loc[gene_, 2]))[0]


logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info(f'Read the VCF file and extract SNPs ...')
# Read the BED file and the VCF file
gene_bed = pd.read_csv(args.gene_bed, sep='\t', header=None, index_col=[3])
vcf = pd.read_csv(args.vcf_ref, sep='\t', comment='#', header=None)
vcf_header = vcf.columns.to_list()
vcf_header[:9] = ['CHR', 'BP', 'SNP', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
vcf.columns = vcf_header
vcf = vcf.iloc[np.sort(np.unique(np.concatenate(list(map(get_idx_inside_gene, gene_bed.index)))))]

logging.info(f'Filter SNPs and make a DataFrame for PLINK BIM files ...')
# Extract reference haplotypes and SNP information from the VCF file
vcf = vcf[(~vcf[['CHR', 'BP']].duplicated(keep=False)) &
          (~vcf['SNP'].duplicated(keep=False)) & vcf['SNP'].str.match(r'^rs\d+$') &
          vcf['REF'].isin(['A', 'C', 'G', 'T']) & vcf['ALT'].isin(['A', 'C', 'G', 'T']) &
          np.all(vcf.loc[:, 9:].isin(['0|0', '0|1', '1|0', '1|1']), axis=1)]
ref_haplo = np.concatenate((vcf.loc[:, 9:].applymap(lambda x: x[0]).values.astype(np.int8).T,
                            vcf.loc[:, 9:].applymap(lambda x: x[2]).values.astype(np.int8).T))

# Filter by MAF
allele_frq = ref_haplo.mean(axis=0)
extract = (allele_frq > maf_min) & (allele_frq < 1 - maf_min)
ref_haplo = ref_haplo[:, extract]
allele_frq = allele_frq[extract]

# Make a PLINK bim file
plink_bim = vcf.loc[extract, ['CHR', 'SNP', 'QUAL', 'BP', 'ALT', 'REF']].rename(
    {'QUAL': 'CM', 'ALT': 'A1', 'REF': 'A2'}, axis=1)
plink_bim['CM'] = 0
plink_bim['hm3SNP'] = plink_bim['SNP'].isin(np.loadtxt(args.hm3_snp, dtype=str))
del vcf, vcf_header, extract

# Perform simulations
logging.info(f'Start simulating {gene_bed.shape[0]} genes in {args.nt} threads')
snp_counts = Pool(args.nt).map(hapsim_one_gene, gene_bed.index)
snp_counts = pd.DataFrame(snp_counts, index=gene_bed.index, columns=['allSNP', 'hm3SNP'])
snp_counts.index.name = None
snp_counts.to_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t')
logging.info(f'Done')
