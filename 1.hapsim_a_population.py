#!/usr/bin/env python

import os
import time
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm
from rpy2 import robjects
from rpy2.robjects.packages import importr
from pysnptools.snpreader import Bed, SnpData
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Simulate genotypes by HapSim and estimate h2 by HESS and EHE')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--gene-bed', type=str, required=True, help='A BED file without header; the fourth column is IDs')
parser.add_argument('--vcf-ref', type=str, required=True, help='A phased VCF file of a reference population')
parser.add_argument('--hm3-snp', type=str, required=True, help='File of a list of HapMap3 SNP ids')
parser.add_argument('--nt', type=int, default=1, help='Number of processes running in parallel')
parser.add_argument('--n-gwa', type=float, default=5e3, help='Size of a sample for association tests')
parser.add_argument('--n-ld', type=str, default='5e2', help='Sizes of samples for LD panels')
parser.add_argument('--n-rep', type=float, default=100, help='Size of a sample for association tests')
parser.add_argument('--maf-min', type=float, default=0.01, help='Ignore SNPs with MAF < MAF_MIN')

args = parser.parse_args()
n_gwa = int(args.n_gwa)
n_ld_str = args.n_ld.split(',')
n_ld_int = [int(float(n)) for n in n_ld_str]
n_rep = int(args.n_rep)
n_pop = n_rep * (n_gwa + np.sum(n_ld_int))  # population size


def hapsim_one_gene(gene_):
    # Simulate for each region defined in the gene BED file
    importr('hapsim')

    # Get reference haplotypes of the i-th gene
    snp_idx = np.where((plink_bim.iloc[:, 0] == gene_bed.loc[gene_, 0]) &
                       (plink_bim.iloc[:, 3] >= gene_bed.loc[gene_, 1]) &
                       (plink_bim.iloc[:, 3] <= gene_bed.loc[gene_, 2]))[0]
    ref_haplo_i = ref_haplo[:, snp_idx]
    allele_frq_i = allele_frq[snp_idx]
    plink_bim_i = plink_bim.iloc[snp_idx]
    n_haplo, m = ref_haplo_i.shape
    if m == 0:
        logging.info(f'Skipped {m} SNPs in {gene_}.')
        return 0, 0

    # Write a one-column HapMap3 SNP list
    dir_gene = f'{args.out_dir}/{gene_}'
    os.system(f'mkdir -p {dir_gene}')
    plink_bim_i.loc[plink_bim_i['hm3SNP'], 'SNP'].to_csv(f'{dir_gene}/plink.hm3snp', header=False, index=False)

    # Write SNP files for LDER
    plink_bim_i[['SNP', 'A1', 'A2']].to_csv(f'{dir_gene}/lder.snp', sep='\t', index=False)
    plink_bim_i.loc[plink_bim_i['hm3SNP'], ['SNP', 'A1', 'A2']].to_csv(f'{dir_gene}/lder.hm3snp', sep='\t', index=False)

    # Calculate an MVN covariance matrix using HapSim
    haplodata = robjects.r('haplodata')
    haplodata = haplodata(robjects.r.matrix(robjects.IntVector(ref_haplo_i.T.reshape(-1)), nrow=n_haplo))
    C = np.array(dict(zip(haplodata.names, list(haplodata)))['cor'])  # m * m

    # Sample genotypes
    percent_point = norm.ppf(allele_frq_i)  # m
    X012 = np.random.multivariate_normal(np.zeros(m), C, (2, n_pop))  # 2 * n_pop * m
    X012 = np.int8(X012 < percent_point)  # 2 * n_pop * m
    X012 = X012[0] + X012[1]  # n_pop * m; int8

    # Write PLINK files
    iid = np.arange(1, n_pop + 1).astype(str)[:, None]
    iid = np.concatenate((iid, iid), axis=1)  # n_pop * 2
    Bed.write(f'{dir_gene}/plink.bed', count_A1=True, _require_float32_64=False,
              snpdata=SnpData(val=X012, iid=iid, sid=plink_bim_i['SNP'],
                              pos=plink_bim_i[['CHR', 'CM', 'BP']], _require_float32_64=False))

    # Since the alleles in .bim file written by "Bed.write" are artificial, rewrite the .bim file with actual alleles.
    os.remove(f'{dir_gene}/plink.bim')
    plink_bim_i[['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']].to_csv(
        f'{dir_gene}/plink.bim', sep='\t', index=False, header=False)

    # Write a region file for KGGSEE
    with open(f'{dir_gene}/region.kggsee', 'w') as o:
        print(f'{gene_bed.loc[gene_, 0]}\t{gene_bed.loc[gene_, 1]}\t{gene_bed.loc[gene_, 2]}\t{gene_}', file=o)

    # Write a region file for HESS and LDER
    with open(f'{dir_gene}/region.hess', 'w') as o:
        print('chr\tstart\tstop', file=o)
        print(f'chr{gene_bed.loc[gene_, 0]}\t{gene_bed.loc[gene_, 1]}\t{gene_bed.loc[gene_, 2]}', file=o)
    os.system(f'cp {dir_gene}/region.hess {dir_gene}/region.lder')

    # Write the region file for LDAK
    with open(f'{dir_gene}/region.ldak', 'w') as o:
        print(f'{gene_}\t{gene_bed.loc[gene_, 0]}\t{gene_bed.loc[gene_, 1]}\t{gene_bed.loc[gene_, 2]}', file=o)

    # Take GWAS samples and LD samples and write PINK files
    tot_n_gwa = n_rep * n_gwa

    X012_gwa = X012[:tot_n_gwa].reshape(n_rep, n_gwa, m)
    iid_gwa = iid[:tot_n_gwa].reshape(n_rep, n_gwa, 2)

    X012_ld_lst = list()
    iid_ld_lst = list()
    prefixes = list()
    start = tot_n_gwa

    for a, n_ld in zip(n_ld_str, n_ld_int):
        end = start + n_rep * n_ld
        X012_ld_lst.append(X012[start:end].reshape(n_rep, n_ld, m))
        iid_ld_lst.append(iid[start:end].reshape(n_rep, n_ld, 2))
        prefixes.append(f'plink_ld{a}')
        start = end

    for j in range(n_rep):
        dir_ij = f'{dir_gene}/rep{j}'
        os.system(f'mkdir -p {dir_ij}')

        Bed.write(f'{dir_ij}/plink_gwa.bed', count_A1=True, _require_float32_64=False,
                  snpdata=SnpData(val=X012_gwa[j], iid=iid_gwa[j], sid=plink_bim_i['SNP'],
                                  pos=plink_bim_i[['CHR', 'CM', 'BP']], _require_float32_64=False))
        os.remove(f'{dir_ij}/plink_gwa.bim')
        os.symlink('../plink.bim', f'{dir_ij}/plink_gwa.bim')

        for X012_ld, iid_ld, prefix in zip(X012_ld_lst, iid_ld_lst, prefixes):
            Bed.write(f'{dir_ij}/{prefix}.bed', count_A1=True, _require_float32_64=False,
                      snpdata=SnpData(val=X012_ld[j], iid=iid_ld[j], sid=plink_bim_i['SNP'],
                                      pos=plink_bim_i[['CHR', 'CM', 'BP']], _require_float32_64=False))
            os.remove(f'{dir_ij}/{prefix}.bim')
            os.symlink('../plink.bim', f'{dir_ij}/{prefix}.bim')

    logging.info(f'Simulated {m} SNPs in {gene_}.')
    return m, plink_bim_i['hm3SNP'].sum()


logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
logging.info(f'Reading the VCF file and filter SNPs ...')

# Read the BED file and the VCF file
gene_bed = pd.read_csv(args.gene_bed, sep='\t', header=None, index_col=[3])
vcf_full = pd.read_csv(args.vcf_ref, sep='\t', comment='#', header=None)
vcf_header = vcf_full.columns.to_list()
vcf_header[:9] = ['CHR', 'BP', 'SNP', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
vcf_full.columns = vcf_header


# Extract SNPs inside genes
def get_idx_inside_gene(gene_):
    return np.where((vcf_full.loc[:, 'CHR'] == gene_bed.loc[gene_, 0]) &
                    (vcf_full.loc[:, 'BP'] >= gene_bed.loc[gene_, 1]) &
                    (vcf_full.loc[:, 'BP'] <= gene_bed.loc[gene_, 2]))[0]


idx_arrays = Pool(args.nt).map(get_idx_inside_gene, gene_bed.index)
vcf = vcf_full.iloc[np.sort(np.unique(np.concatenate(idx_arrays)))]
del vcf_full

# Extract reference haplotypes and SNP information from the VCF file
vcf = vcf[(~vcf[['CHR', 'BP']].duplicated(keep=False)) &
          (~vcf['SNP'].duplicated(keep=False)) & vcf['SNP'].str.match(r'^rs\d+$') &
          vcf['REF'].isin(['A', 'C', 'G', 'T']) & vcf['ALT'].isin(['A', 'C', 'G', 'T']) &
          np.all(vcf.loc[:, 9:].isin(['0|0', '0|1', '1|0', '1|1']), axis=1)]
ref_haplo = np.concatenate((vcf.loc[:, 9:].applymap(lambda x: x[0]).values.astype(np.int8).T,
                            vcf.loc[:, 9:].applymap(lambda x: x[2]).values.astype(np.int8).T))
# Filter by MAF
allele_frq = ref_haplo.mean(axis=0)
extract = (allele_frq > args.maf_min) & (allele_frq < 1 - args.maf_min)
ref_haplo = ref_haplo[:, extract]
allele_frq = allele_frq[extract]

# Make a PLINK bim file
plink_bim = vcf.loc[extract, ['CHR', 'SNP', 'QUAL', 'BP', 'ALT', 'REF']].rename(
    {'QUAL': 'CM', 'ALT': 'A1', 'REF': 'A2'}, axis=1)
plink_bim['CM'] = 0
plink_bim['hm3SNP'] = plink_bim['SNP'].isin(np.loadtxt(args.hm3_snp, dtype=str))
del vcf, vcf_header, extract

# Perform simulations
logging.info(f'Start the simulations.')
snp_counts = Pool(args.nt).map(hapsim_one_gene, gene_bed.index)
snp_counts = pd.DataFrame(snp_counts, index=gene_bed.index, columns=['allSNP', 'hm3SNP'])
snp_counts.index.name = None
snp_counts.to_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t')
logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
