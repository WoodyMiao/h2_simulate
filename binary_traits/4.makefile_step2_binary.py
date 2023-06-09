#!/usr/bin/env python

import os
import time
import logging
import argparse
import numpy as np
import pandas as pd
from itertools import product
from multiprocessing.dummy import Pool

# Set the executables
gcta_exe = 'gcta --threads 1'
hess_exe = '/app/user/ml/hess-0.5-chr1/hess.py'
ldak_exe = 'ldak'
lder_r = 'R'
ldsc_exe = 'ldsc.py'
kggsee_exe = 'java -Xmx16g -jar /app/pmglab/kggsee/kggsee.jar --nt 1'

parser = argparse.ArgumentParser(description='Make shell scripts of h2 estamating programs')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--gene-list', type=str, default=None, help='File of a list of genes to be included')
parser.add_argument('--chrom', type=str, required=True, help='The chromosome where the data from')
parser.add_argument('--n-gwa', type=float, default=2e4, help='Size of a sample for association tests')
parser.add_argument('--n-ld', type=str, default='5e2,2e3', help='Sizes of samples for LD panels')
parser.add_argument('--n-rep', type=float, default=100, help='Size of a sample for association tests')
parser.add_argument('--h2g', type=float, default=1e-3, help='Target h2g for the simulation')
parser.add_argument('--maf-min', type=float, default=0.01, help='Ignore SNPs with MAF < MAF_MIN')
parser.add_argument('--neg-alpha', type=float, default=0.25, help='Power value in the LDAK-Thin Model')
parser.add_argument('--pqtl-vals', type=str, default='1,0.1,0.01', help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--prevalence-vals', type=str, default='0.1,0.01',
                    help='prevalence values for the liability threshold modeling of dichotamous phenotypes')
parser.add_argument('--hess', action='store_true', default=False, help='Write shell scripts to run HESS')
parser.add_argument('--kggsee', action='store_true', default=False, help='Write shell scripts to run KGGSEE')
parser.add_argument('--ldak-gbat', action='store_true', default=False, help='Write shell scripts to run LDAK-GBAT')
parser.add_argument('--ldak-sumher', action='store_true', default=False, help='Write shell scripts to run LDAK-SumHer')
parser.add_argument('--lder', action='store_true', default=False, help='Write shell scripts to run LDER')
parser.add_argument('--ldsc', action='store_true', default=False, help='Write shell scripts to run LDSC')
parser.add_argument('--skip-mkdir', action='store_true', default=False, help='Skip the commands of mkdir')
parser.add_argument('--makefile-only', action='store_true', default=False, help='Only remake makefiles')
parser.add_argument('--jobs', type=int, default=100, help='The number jobs for the GNU make to run simultaneously')
parser.add_argument('--assoc-test', type=str, default='logit',
                    help='chi2|linear|logit, part of the suffix of the input files.')

args = parser.parse_args()
h2g = args.h2g
alpha = -args.neg_alpha
maf_min = args.maf_min
n_gwa = int(args.n_gwa)
n_case = int(n_gwa / 2)
n_rep = int(args.n_rep)
n_ld_str = args.n_ld.split(',')
n_ld_int = [int(float(n)) for n in n_ld_str]
ld_sfx = [f'ld{n}' for n in args.n_ld.split(',')]
pqtl_lst = args.pqtl_vals.split(',')
prevalence_lst = args.prevalence_vals.split(',')
par_tup_list = list(product(prevalence_lst, pqtl_lst))

if not args.gene_list:
    snp_counts = pd.read_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t', index_col=0)
    gene_list = snp_counts.loc[snp_counts.allSNP >= 3].index.values
else:
    gene_list = np.loadtxt(args.gene_list, dtype=str)

methods = list()
if args.kggsee:
    methods.append('kggsee')
if args.hess:
    methods.append('hess')
if args.ldsc:
    methods.append('ldsc')
if args.lder:
    methods.append('lder')
if args.ldak_sumher:
    methods.append('ldak_sumher')
if args.ldak_gbat:
    methods.append('ldak_gbat')

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')


def one_parameter_set_makefile(prvl_pqtl):
    prvl, pqtl = prvl_pqtl
    par_str = f'prvl{prvl}_pqtl{pqtl}_alpha{alpha}.{args.assoc_test}'

    if not args.makefile_only:

        for j in range(n_rep):
            logging.info(f'Processing commands for {par_str} rep{j} ...')

            for gene in gene_list:
                dir_gene = f'{args.out_dir}/{gene}'
                dir_rep = f'{args.out_dir}/{gene}/rep{j}'
                dir_out = f'{dir_rep}/{par_str}'

                if not args.skip_mkdir:
                    os.system(f'mkdir -p {dir_out}')

                if args.kggsee:
                    with open(f'{dir_out}/kggsee.sh', 'w') as o:
                        for k in ld_sfx:
                            print(f'{kggsee_exe} --filter-maf-le {maf_min} --gene-herit '
                                  f'--vcf-ref {dir_rep}/plink_{k}.vcf.gz --sum-file {dir_rep}/{par_str}.sumstat.gz '
                                  f'--nmiss-col N --regions-bed {dir_gene}/region.kggsee '
                                  f'--out {dir_out}/kggsee_{k} >/dev/null 2>&1\n', file=o)

                if args.hess:
                    with open(f'{dir_out}/hess.sh', 'w') as o:
                        for k in ld_sfx:
                            print(f'{hess_exe} --min-maf {maf_min} '
                                  f'--local-hsqg {dir_rep}/{par_str}.sumstat.gz --chrom {args.chrom} '
                                  f'--bfile {dir_rep}/plink_{k} --partition {dir_gene}/region.hess '
                                  f'--out {dir_out}/hess_{k}.step1 >/dev/null 2>&1', file=o)
                            print(f'{hess_exe} --prefix {dir_out}/hess_{k}.step1 '
                                  f'--out {dir_out}/hess_{k}.step2 >/dev/null 2>&1\n', file=o)

                if args.ldsc:
                    with open(f'{dir_out}/ldsc.sh', 'w') as o:
                        for k in ld_sfx:
                            print(f'{ldsc_exe} --h2 {dir_rep}/{par_str}.sumstat.gz --no-intercept '
                                  f'--ref-ld {dir_rep}/ldsc_{k} --w-ld {dir_rep}/ldsc_{k} --n-blocks 2 '
                                  f'--out {dir_out}/ldsc_{k} >/dev/null 2>&1\n', file=o)

                if args.lder:
                    with open(f'{dir_out}/lder.sh', 'w') as o:
                        print(f'mkdir -p {dir_out}/plinkLD/ldetect-data', file=o)
                        print(f'ln -fs ../../../../region.lder {dir_out}/plinkLD/ldetect-data/fourier_ls-all.bed\n',
                              file=o)

                        for nld, k in zip(n_ld_int, ld_sfx):
                            if k == 'gwa':
                                insample = 'T'
                            else:
                                insample = 'F'
                            print(f'ln -fs ../lder_{k}.LD {dir_out}/LD.shrink', file=o)
                            print(f'{lder_r} -q -e "library(LDER); library(data.table); '
                                  f'assoc <- fread(\'{dir_rep}/{par_str}.lder.sumstat.gz\'); '
                                  f'runLDER(assoc=assoc, n.gwas={n_gwa}, n.ld={nld}, path=\'{dir_out}\', '
                                  f'LD.insample={insample}, a=0, cores=1, type=\'boot\', n.bs=0)" '
                                  f'>{dir_out}/lder_{k}.result 2>/dev/null', file=o)
                            print(f'rm {dir_out}/LD.shrink\n', file=o)

                        print(f'rm -r {dir_out}/plinkLD', file=o)

                if args.ldak_sumher:
                    with open(f'{dir_out}/ldak_sumher.sh', 'w') as o:
                        for k in ld_sfx:
                            print(f'{ldak_exe} --sum-hers {dir_out}/ldak_sumher_{k} '
                                  f'--summary {dir_rep}/{par_str}.ldak.sumstat '
                                  f'--tagfile {dir_rep}/ldak_sumher_{k}.ldak_thin.tagging >/dev/null 2>&1\n', file=o)

                if args.ldak_gbat:
                    with open(f'{dir_out}/ldak_gbat.sh', 'w') as o:
                        for k in ld_sfx:
                            print(f'{ldak_exe} --cut-genes {dir_out}/ldak_gbat_{k} --genefile {dir_gene}/region.ldak '
                                  f'--bfile {dir_rep}/plink_{k} >/dev/null 2>&1', file=o)
                            print(f'{ldak_exe} --calc-genes-reml {dir_out}/ldak_gbat_{k} '
                                  f'--summary {dir_rep}/{par_str}.ldak.sumstat --bfile {dir_rep}/plink_{k} '
                                  f'--ignore-weights YES --power -0.25 --allow-ambiguous YES >/dev/null 2>&1\n', file=o)

    # Write a makefile for each method that calls all shell scripts for every genes and repetitions.
    makefile_dict_ = dict()
    for method_ in methods:
        makefile = f'{args.out_dir}/makefile.{par_str}.{method_}'
        logging.info(f'Writing {makefile} ...')
        target_final = f'{makefile}.done:'
        target_list = list()
        for j in range(n_rep):
            for gene in gene_list:
                target_final += f' {gene}rep{j}'
                target_list.append(
                    f'{gene}rep{j}:\n\tsh {args.out_dir}/{gene}/rep{j}/{par_str}/{method_}.sh\n')

        target_final += f'\n\ttouch {makefile}.done\n'
        with open(makefile, 'w') as o:
            o.write(target_final)
            o.writelines(target_list)

        makefile_dict_[method_] = makefile

    return makefile_dict_


list_of_dict = Pool(len(par_tup_list)).map(one_parameter_set_makefile, par_tup_list)
for method in methods:
    with open(f'{args.out_dir}/makefile.step2.{args.assoc_test}.{method}.sh', 'w') as O:
        if method in ['hess', 'lder', 'ldsc']:
            O.write('''
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
\n''')

        if method == 'kggsee':
            nj = int(args.jobs / 2.5)
        else:
            nj = args.jobs

        for makefile_dict in list_of_dict:
            O.write(f'make -kj {nj} -f {makefile_dict[method]} >{makefile_dict[method]}.log 2>&1\n')

logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
