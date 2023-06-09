#!/usr/bin/env python

import os
import time
import logging
import argparse
import numpy as np
import pandas as pd

# Set the executables
gcta_exe = 'gcta --threads 1'
hess_exe = '/app/user/ml/hess-0.5-chr1/hess.py'
ldak_exe = 'ldak'
lder_r = 'R'
lder_plinkLD = 'plinkLD.py --thread 1'
ldsc_exe = 'ldsc.py'
plink_exe = 'plink --threads 1'
kggsee_exe = 'java -Xmx4g -jar /app/pmglab/kggsee/kggsee.jar --nt 1'

parser = argparse.ArgumentParser(description='Make shell scripts of h2 estamating programs')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--gene-list', type=str, default=None, help='File of a list of genes to be included')
parser.add_argument('--n-ld', type=str, default='5e2', help='Sizes of samples for LD panels')
parser.add_argument('--n-rep', type=float, default=100, help='Number of repititions')
parser.add_argument('--extract-hm3', action='store_true', default=False,
                    help='Extract hm3 SNPs from "plink.hm3snp" and insert "hm3_" to the output file names')
parser.add_argument('--skip-mkdir', action='store_true', default=False, help='Skip commands of mkdir')
parser.add_argument('--skip-gcta', action='store_true', default=False, help='Skip commands for GCTA')
parser.add_argument('--skip-assoc-ld', action='store_true', default=False,
                    help='Skip commands on samples of association tests')
parser.add_argument('--old-suffix', action='store_true', default=False, help='To be compatible with old runs')

args = parser.parse_args()
n_rep = int(args.n_rep)
if args.old_suffix:
    ld_sfx = ['ld1', 'ld2']
else:
    ld_sfx = [f'ld{n}' for n in args.n_ld.split(',')]

if not args.skip_assoc_ld:
    ld_sfx = ['gwa'] + ld_sfx

if not args.gene_list:
    snp_counts = pd.read_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t', index_col=0)
    if args.extract_hm3:
        gene_list = snp_counts.loc[snp_counts.hm3SNP >= 3].index.values
    else:
        gene_list = snp_counts.loc[snp_counts.allSNP >= 3].index.values
else:
    gene_list = np.loadtxt(args.gene_list, dtype=str)

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')

targets = list()
targets_gcta = list()
cmd_dict = {method: [] for method in ['kggsee', 'ldak_sumher', 'ldsc', 'lder']}
if (not args.skip_gcta) and (not args.extract_hm3):
    cmd_dict['gcta'] = []

for j in range(n_rep):
    logging.info(f'Processing commands for rep{j} ...')

    for gene in gene_list:
        dir_gene = f'{args.out_dir}/{gene}'
        dir_rep = f'{args.out_dir}/{gene}/rep{j}'

        if not args.skip_mkdir:
            os.system(f'mkdir -p {dir_rep}/plinkLD/ldetect-data')
            os.system(f'ln -fs ../../../region.lder {dir_rep}/plinkLD/ldetect-data/fourier_ls-all.bed')

        if (not args.skip_gcta) and (not args.extract_hm3):
            targets_gcta.append(f' {gene}rep{j}')

            # Calculate a GRM for GATK
            cmd_dict['gcta'].append(
                f'{gene}rep{j}:\n'
                f'\t{gcta_exe} --bfile {dir_rep}/plink_gwa --make-grm --out {dir_rep}/gcta >/dev/null 2>&1\n'
            )

        for k in ld_sfx:
            targets.append(f' {gene}rep{j}{k}')

            if not args.extract_hm3:
                hm3k = k

                # Make a VCF file for KGGSEE
                cmd_dict['kggsee'].append(
                    f'{gene}rep{j}{k}:\n'
                    f'\t{plink_exe} --recode vcf-fid bgz --real-ref-alleles '
                    f'--out {dir_rep}/plink_{k} --bfile {dir_rep}/plink_{k} >/dev/null 2>&1\n'
                    f'\trm {dir_rep}/plink_{k}.log {dir_rep}/plink_{k}.nosex\n'
                )
            else:
                hm3k = f'hm3_{k}'
                # Make a VCF file for KGGSEE and make a set of PLINK files with only HM3 SNPs.
                # Since the PLINK files are needed by all methods, make "makefile.hm3_step1.kggsee" first.
                cmd_dict['kggsee'].append(
                    f'{gene}rep{j}{k}:\n'
                    f'\t{plink_exe} --recode vcf-fid bgz --real-ref-alleles '
                    f'--extract {dir_gene}/plink.hm3snp --make-bed '
                    f'--out {dir_rep}/plink_{hm3k} --bfile {dir_rep}/plink_{k} >/dev/null 2>&1\n'
                    f'\trm {dir_rep}/plink_{hm3k}.log {dir_rep}/plink_{hm3k}.nosex\n'
                )

            # Calculate LD scores for LDAK-SumHer
            cmd_dict['ldak_sumher'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{ldak_exe} --thin {dir_rep}/ldak_sumher_{hm3k} --bfile {dir_rep}/plink_{hm3k} '
                f'--window-prune .98 --window-kb 100 >/dev/null 2>&1\n'
                "\tawk '{print $$1, 1}' "
                f'{dir_rep}/ldak_sumher_{hm3k}.in >{dir_rep}/ldak_sumher_{hm3k}.weights.thin\n'
                f'\t{ldak_exe} --calc-tagging {dir_rep}/ldak_sumher_{hm3k}.ldak_thin --bfile {dir_rep}/plink_{hm3k} '
                f'--weights {dir_rep}/ldak_sumher_{hm3k}.weights.thin --power -0.25 --window-kb 100 >/dev/null 2>&1\n'
            )

            # Calculate LD scores for LDSC
            cmd_dict['ldsc'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{ldsc_exe} --yes-really --bfile {dir_rep}/plink_{hm3k} --l2 --ld-wind-kb 100 '
                f'--out {dir_rep}/ldsc_{hm3k} >/dev/null 2>&1\n'
            )

            # Calculate LD scores for LDER
            cmd_dict['lder'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{lder_plinkLD} --bfile {dir_rep}/plink_{hm3k} --block {dir_gene}/region.lder '
                f'--snplist {dir_gene}/lder.hm3snp --output {dir_rep}/lder_{hm3k}.LD >/dev/null 2>&1\n'
            )

# Write a makefile for each method that calls all shell scripts for every genes and repetitions.
for method, cmd in cmd_dict.items():
    if not args.extract_hm3:
        makefile = f'{args.out_dir}/makefile.step1.{method}'
    else:
        makefile = f'{args.out_dir}/makefile.hm3_step1.{method}'

    logging.info(f'Writing {makefile} ...')
    if method == 'gcta':
        first_line = f'{makefile}.done: ' + \
                     ' '.join(targets_gcta) + \
                     f'\n\ttouch {makefile}.done\n'
    else:
        first_line = f'{makefile}.done: ' + \
                     ' '.join(targets) + \
                     f'\n\ttouch {makefile}.done\n'

    with open(f'{makefile}', 'w') as o:
        o.write(first_line)
        o.writelines(cmd)

logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
