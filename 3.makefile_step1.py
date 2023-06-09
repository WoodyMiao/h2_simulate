#!/usr/bin/env python

import os
import time
import logging
import argparse
import pandas as pd
from shutil import rmtree

# Set the executables
gcta_exe = 'gcta'
hess_exe = 'hess.py'
ldak_exe = 'ldak'
lder_r = 'R'
lder_plinkLD = 'plinkLD.py'
ldsc_exe = 'ldsc.py'
plink_exe = 'plink'
kggsee_exe = 'java -Xmx4g -jar /app/pmglab/kggsee/kggsee.jar --nt 1'

parser = argparse.ArgumentParser(description='Make shell scripts of h2 estamating programs')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--n-rep', type=float, default=100, help='Number of repititions')
parser.add_argument('--skip-mkdir', action='store_true', default=False, help='Skip the commands of mkdir')

args = parser.parse_args()
n_rep = int(float(args.n_rep))
snp_counts = pd.read_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t', index_col=0)
snp_counts = snp_counts[snp_counts.allSNP > 0]

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')

targets = list()
cmds_kggsee_ldak_sumher = list()
cmds_ldsc_lder = list()

for gene in snp_counts.index:
    logging.info(f'Processing commands for {gene} ...')
    dir_gene = f'{args.out_dir}/{gene}'

    for j in range(n_rep):
        dir_rep = f'{args.out_dir}/{gene}/rep{j}'

        if not args.skip_mkdir:
            os.system(f'mkdir -p {dir_rep}/plinkLD/ldetect-data')
            os.system(f'ln -fs ../../../region.lder {dir_rep}/plinkLD/ldetect-data/fourier_ls-all.bed')

        for k in ['gwa', 'ld1', 'ld2']:
            targets.append(f' {gene}rep{j}{k}')
            cmds_kggsee_ldak_sumher.append(
                f'{gene}rep{j}{k}:\n'
                # Make a VCF file for KGGSEE
                f'\t{plink_exe} --threads 1 --recode vcf-fid bgz --real-ref-alleles '
                f'--out {dir_rep}/plink_{k} --bfile {dir_rep}/plink_{k} >/dev/null 2>&1\n'
                f'\trm {dir_rep}/plink_{k}.log {dir_rep}/plink_{k}.nosex\n'
                # Calculate LD scores for LDAK-SumHer
                f'\t{ldak_exe} --thin {dir_rep}/ldak_sumher_{k} --bfile {dir_rep}/plink_{k} '
                f'--window-prune .98 --window-kb 100 >/dev/null 2>&1\n'
                "\tawk '{print $$1, 1}' "
                f'{dir_rep}/ldak_sumher_{k}.in >{dir_rep}/ldak_sumher_{k}.weights.thin\n'
                f'\t{ldak_exe} --calc-tagging {dir_rep}/ldak_sumher_{k}.ldak_thin --bfile {dir_rep}/plink_{k} '
                f'--weights {dir_rep}/ldak_sumher_{k}.weights.thin --power -0.25 --window-kb 100 >/dev/null 2>&1\n'
            )

            cmds_ldsc_lder.append(
                f'{gene}rep{j}{k}:\n'
                # Calculate LD scores for LDSC
                f'\t{ldsc_exe} --yes-really --bfile {dir_rep}/plink_{k} --l2 --ld-wind-kb 100 '
                f'--out {dir_rep}/ldsc_{k} >/dev/null 2>&1\n'
                # Calculate LD scores for LDER
                f'\t{lder_plinkLD} --thread 1 --bfile {dir_rep}/plink_{k} --block {dir_gene}/region.lder '
                f'--snplist {dir_gene}/lder.snp --output {dir_rep}/lder_{k}.LD >/dev/null 2>&1\n'
            )

# Write a makefile for each method that calls all shell scripts for every genes and repetitions.
logging.info(f'Writing {args.out_dir}/makefile.step1.kggsee_ldak_sumher ...')
target_kggsee_ldak_sumher = f'{args.out_dir}/makefile.step1.kggsee_ldak_sumher.done: ' + \
                            ' '.join(targets) + \
                            f'\n\ttouch {args.out_dir}/makefile.step1.kggsee_ldak_sumher.done\n'

with open(f'{args.out_dir}/makefile.step1.kggsee_ldak_sumher', 'w') as o:
    o.write(target_kggsee_ldak_sumher)
    o.writelines(cmds_kggsee_ldak_sumher)

logging.info(f'Writing {args.out_dir}/makefile.step1.ldsc_lder ...')
target_ldsc_lder = f'{args.out_dir}/makefile.step1.done: ' + \
                   ' '.join(targets) + \
                   f'\n\ttouch {args.out_dir}/makefile.step1.ldsc_lder.done\n'
with open(f'{args.out_dir}/makefile.step1.ldsc_lder', 'w') as o:
    o.write(target_ldsc_lder)
    o.writelines(cmds_ldsc_lder)

logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
