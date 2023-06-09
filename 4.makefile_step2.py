#!/usr/bin/env python

import os
import time
import logging
import argparse
import pandas as pd
from itertools import product
from multiprocessing import Pool

# Set the executables
gcta_exe = 'gcta'
hess_exe = 'hess.py'
ldak_exe = 'ldak'
lder_r = 'R'
lder_plinkLD = 'plinkLD.py'
ldsc_exe = 'ldsc.py'
plink_exe = 'plink'
kggsee_exe = 'java -Xmx16g -jar /app/pmglab/kggsee/kggsee.jar --nt 1'

parser = argparse.ArgumentParser(description='Make shell scripts of h2 estamating programs')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--process', type=int, default=1, help='Number of processes running in parallel')
parser.add_argument('--chrom', type=str, required=True, help='The chromosome where the data from')
parser.add_argument('--n-gwa', type=float, default=2e4, help='Size of a sample for association tests')
parser.add_argument('--n-ld1', type=float, default=2e3, help='Size of a sample for LD panel 1')
parser.add_argument('--n-ld2', type=float, default=500, help='Size of a sample for LD panel 2')
parser.add_argument('--n-rep', type=float, default=100, help='Size of a sample for association tests')
parser.add_argument('--h2g-vals', type=str, default='1e-2,1e-3,1e-4', help='Target h2g values for the simulation')
parser.add_argument('--pqtl-vals', type=str, default='1,0.1,0.01', help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--neg-alpha-vals', type=str, default='0.25,1.0', help='Power values in the LDAK-Thin Model')
parser.add_argument('--maf-min', type=float, default=0.01, help='Ignore SNPs with MAF < MAF_MIN')
parser.add_argument('--gcta', action='store_true', default=False, help='Write shell scripts to run GCTA')
parser.add_argument('--hess', action='store_true', default=False, help='Write shell scripts to run HESS')
parser.add_argument('--kggsee', action='store_true', default=False, help='Write shell scripts to run KGGSEE')
parser.add_argument('--ldak-gbat', action='store_true', default=False, help='Write shell scripts to run LDAK-GBAT')
parser.add_argument('--ldak-sumher', action='store_true', default=False, help='Write shell scripts to run LDAK-SumHer')
parser.add_argument('--lder', action='store_true', default=False, help='Write shell scripts to run LDER')
parser.add_argument('--ldsc', action='store_true', default=False, help='Write shell scripts to run LDSC')

args = parser.parse_args()
n_gwa = int(float(args.n_gwa))
n_ld1 = int(float(args.n_ld1))
n_ld2 = int(float(args.n_ld2))
n_rep = int(float(args.n_rep))
h2g_list = args.h2g_vals.split(',')
pqtl_list = args.pqtl_vals.split(',')
neg_alpha_list = args.neg_alpha_vals.split(',')
snp_counts = pd.read_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t', index_col=0)
snp_counts = snp_counts[snp_counts.allSNP > 0]

methods = list()
if args.gcta:
    methods.append('gcta')
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


def one_parameter_set_makefile(pqtl_negalpha_h2g):
    pqtl, neg_alpha, h2g = pqtl_negalpha_h2g
    parameter_str = f'pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}'

    for gene in snp_counts.index:
        logging.info(f'Processing commands for pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g} {gene} ...')
        dir_gene = f'{args.out_dir}/{gene}'

        for j in range(n_rep):
            dir_rep = f'{args.out_dir}/{gene}/rep{j}'
            dir_out = f'{dir_rep}/{parameter_str}'
            os.system(f'mkdir -p {dir_out}')

            if args.gcta:
                with open(f'{dir_out}/gcta.sh', 'w') as o:
                    print(f'{gcta_exe} --bfile {dir_rep}/plink_gwa --make-grm '
                          f'--out {dir_out}/gcta >/dev/null 2>&1', file=o)
                    print(f'{gcta_exe} --threads 2 --reml --pheno {dir_rep}/{parameter_str}.pheno '
                          f'--grm {dir_out}/gcta --out {dir_out}/gcta >/dev/null 2>&1', file=o)

            if args.kggsee:
                with open(f'{dir_out}/kggsee.sh', 'w') as o:
                    for k in ['gwa', 'ld1', 'ld2']:
                        print(f'{kggsee_exe} --filter-maf-le {args.maf_min} --gene-herit '
                              f'--vcf-ref {dir_rep}/plink_{k}.vcf.gz --sum-file {dir_rep}/{parameter_str}.sumstat.gz '
                              f'--nmiss-col N --regions-bed {dir_gene}/region.kggsee '
                              f'--out {dir_out}/kggsee_{k} >/dev/null 2>&1', file=o)

            if args.hess:
                with open(f'{dir_out}/hess.sh', 'w') as o:
                    for k in ['gwa', 'ld1', 'ld2']:
                        print(f'{hess_exe} --min-maf {args.maf_min} --local-hsqg {dir_rep}/{parameter_str}.sumstat.gz '
                              f'--chrom {args.chrom} --bfile {dir_rep}/plink_{k} --partition {dir_gene}/region.hess '
                              f'--out {dir_out}/hess_{k}.step1 >/dev/null 2>&1', file=o)
                        print(f'{hess_exe} --prefix {dir_out}/hess_{k}.step1 '
                              f'--out {dir_out}/hess_{k}.step2 >/dev/null 2>&1', file=o)

            if args.ldsc:
                with open(f'{dir_out}/ldsc.sh', 'w') as o:
                    for k in ['gwa', 'ld1', 'ld2']:
                        print(f'{ldsc_exe} --h2 {dir_rep}/{parameter_str}.sumstat.gz --no-intercept '
                              f'--ref-ld {dir_rep}/ldsc_{k} --w-ld {dir_rep}/ldsc_{k} --n-blocks 2 '
                              f'--out {dir_rep}/ldsc_{k} >/dev/null 2>&1', file=o)

            if args.lder:
                with open(f'{dir_out}/lder.sh', 'w') as o:
                    for nld, k in zip([n_gwa, n_ld1, n_ld2], ['gwa', 'ld1', 'ld2']):
                        print(f'ln -s ../lder_{k}.LD {dir_out}/LD.shrink', file=o)
                        print(f'{lder_r} -q -e "library(LDER); library(data.table); '
                              f'assoc <- fread(\'{dir_rep}/{parameter_str}.lder.sumstat.gz\'); '
                              f'runLDER(assoc=assoc, n.gwas={n_gwa}, path=\'{dir_out}\', LD.insample=T, n.ld={nld}, '
                              f'a=0, cores=1, type=\'boot\', n.bs=0)" >{dir_out}/lder_{k}.result 2>/dev/null', file=o)
                        print(f'rm {dir_out}/LD.shrink', file=o)

            if args.ldak_sumher:
                with open(f'{dir_out}/ldak_sumher.sh', 'w') as o:
                    for k in ['gwa', 'ld1', 'ld2']:
                        print(f'{ldak_exe} --sum-hers {dir_out}/ldak_sumher_{k} '
                              f'--summary {dir_rep}/{parameter_str}.ldak.sumstat '
                              f'--tagfile {dir_rep}/ldak_sumher_{k}.ldak_thin.tagging >/dev/null 2>&1', file=o)

            if args.ldak_gbat:
                with open(f'{dir_out}/ldak_gbat.sh', 'w') as o:
                    for k in ['gwa', 'ld1', 'ld2']:
                        print(f'{ldak_exe} --cut-genes {dir_out}/ldak_gbat_{k} --genefile {dir_gene}/region.ldak '
                              f'--bfile {dir_rep}/plink_{k} >/dev/null 2>&1', file=o)
                        print(f'{ldak_exe} --calc-genes-reml {dir_out}/ldak_gbat_{k} '
                              f'--summary {dir_rep}/{parameter_str}.ldak.sumstat '
                              f'--bfile {dir_rep}/plink_{k} --ignore-weights YES '
                              f'--power -0.25 --allow-ambiguous YES >/dev/null 2>&1', file=o)

    # Write a makefile for each method that calls all shell scripts for every genes and repetitions.
    makefile_dict_ = dict()
    for method_ in methods:
        makefile = f'{args.out_dir}/makefile.{parameter_str}.{method_}'
        logging.info(f'Writing {makefile} ...')
        target_final = f'{makefile}.done:'
        target_list = list()
        for j in range(n_rep):
            for gene in snp_counts.index:
                target_final += f' {gene}rep{j}'
                target_list.append(f'{gene}rep{j}:\n\tsh {args.out_dir}/{gene}/rep{j}/{parameter_str}/{method_}.sh\n')

        target_final += f'\n\ttouch {makefile}.done\n'
        with open(makefile, 'w') as o:
            o.write(target_final)
            o.writelines(target_list)

        makefile_dict_[method_] = makefile

    return makefile_dict_


parameter_sets = list(product(pqtl_list, neg_alpha_list, h2g_list))
list_of_dict = Pool(args.process).map(one_parameter_set_makefile, parameter_sets)
for method in methods:
    with open(f'{args.out_dir}/makefile.step2.{method}.pbs', 'w') as O:
        O.write(f'''#PBS -N {method}
#PBS -l nodes=pan02:ppn=40
#PBS -l mem=100g
#PBS -o /dev/null
#PBS -e /dev/null

cd {args.out_dir}\n
''')
        for makefile_dict in list_of_dict:
            O.write(f'make -j 40 -f {makefile_dict[method]} >{makefile_dict[method]}.log 2>&1\n')

logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
