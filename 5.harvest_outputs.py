#!/usr/bin/env python

import re
import time
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm
from itertools import product
from multiprocessing.dummy import Pool

parser = argparse.ArgumentParser(description='Harvest results of h2 estimating programs')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--out-suffix', type=str, required=True, help='Suffix of the output files of this script')
parser.add_argument('--gene-list', type=str, default=None, help='File of a list of genes to be included')
parser.add_argument('--n-ld', type=str, default='5e2', help='Sizes of samples for LD panels')
parser.add_argument('--n-rep', type=float, default=100, help='Number of repititions')
parser.add_argument('--h2g-vals', type=str, default='1e-2,1e-3', help='Target h2g values for the simulation')
parser.add_argument('--pqtl-vals', type=str, default='1.0,0.1', help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--neg-alpha-vals', type=str, default='0.25,1.0', help='Power values in the LDAK-Thin Model')
parser.add_argument('--prevalence-vals', type=str, default='0.1,0.01', help='prevalence values of binary traits')
parser.add_argument('--nt', type=int, default=10, help='Number of threads running in parallel')
parser.add_argument('--hm3-h2', action='store_true', default=False,
                    help='Read results output by the "--hm3-h2" flag set in 4.makefile_step2.py')
parser.add_argument('--hm3-only', action='store_true', default=False,
                    help='Read results output by the "--hm3-only" flag set in 4.makefile_step2.py')
parser.add_argument('--gcta', action='store_true', default=False, help='Harvest GCTA results')
parser.add_argument('--hess', action='store_true', default=False, help='Harvest HESS results')
parser.add_argument('--kggsee', action='store_true', default=False, help='Harvest EHE results')
parser.add_argument('--ldak-gbat', action='store_true', default=False, help='Harvest LDAK-GBAT results')
parser.add_argument('--ldak-sumher', action='store_true', default=False, help='Harvest DAK-SumHer results')
parser.add_argument('--lder', action='store_true', default=False, help='Harvest LDER results')
parser.add_argument('--ldsc', action='store_true', default=False, help='Harvest LDSC results')
parser.add_argument('--binary-assoc-test', type=str, default=None,
                    help='[chi2|linear|logit], if specified, consider the results as from binary traits')
parser.add_argument('--skip-assoc-ld', action='store_true', default=False,
                    help='Skip results on LD matrices calculated by samples of association tests')
parser.add_argument('--old-suffix', action='store_true', default=False, help='To be compatible with old runs')

args = parser.parse_args()
if args.skip_assoc_ld:
    R = []
    ld_sfx = []
else:
    R = ['Rgwa']
    ld_sfx = ['gwa']
if args.n_ld:
    R += [f'R{n}' for n in args.n_ld.split(',')]
    if args.old_suffix:
        ld_sfx += ['ld1', 'ld2']
    else:
        ld_sfx += [f'ld{n}' for n in args.n_ld.split(',')]

n_rep = int(float(args.n_rep))
h2g_lst = args.h2g_vals.split(',')
pqtl_lst = args.pqtl_vals.split(',')
neg_alpha_lst = args.neg_alpha_vals.split(',')

if not args.binary_assoc_test:
    par_tup_list = [(pqtl_, neg_alpha_, h2g_)
                    for pqtl_, neg_alpha_, h2g_ in product(pqtl_lst, neg_alpha_lst, h2g_lst)]
else:
    prevalence_lst = [float(a) for a in args.prevalence_vals.split(',')]
    coeffcients = {k: (k * (1 - k) / norm.pdf(norm.ppf(k))) ** 2 * 4 for k in prevalence_lst}
    par_tup_list = [(prvl_, pqtl_, neg_alpha_, args.binary_assoc_test, h2g_)
                    for prvl_, pqtl_, neg_alpha_, h2g_ in product(prevalence_lst, pqtl_lst, neg_alpha_lst, h2g_lst)]

if not args.gene_list:
    snp_counts = pd.read_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t', index_col=0)
    gene_list = snp_counts.loc[snp_counts.allSNP >= 3].index.values
else:
    gene_list = np.loadtxt(args.gene_list, dtype=str)

if args.hm3_h2 or args.hm3_only:
    ld_sfx = [f'hm3_{a}' for a in ld_sfx]
    gene_sub_lst = np.loadtxt('a.coding_gene_bed/chr1.one4th.3more_hm3snps.lst', dtype=str)
    gene_list = np.intersect1d(gene_list, gene_sub_lst)

methods = list()
if args.gcta:
    methods.append('GCTA')
if args.kggsee:
    methods.append('EHE')
if args.hess:
    methods.append('HESS')
if args.ldak_gbat:
    methods.append('GBAT')
if args.ldak_sumher:
    methods.append('SumHer')
if args.ldsc:
    methods.append('LDSC')
if args.lder:
    methods.append('LDER')

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')


def read_one_par_set_result(gene, par_tup_):
    if not args.binary_assoc_test:
        pqtl_, neg_alpha_, h2g_ = par_tup_
        par_str_ = f'pqtl{pqtl_}_alpha-{neg_alpha_}_h2g{h2g_}'
    else:
        prvl_, pqtl_, neg_alpha_, assoc_test_, h2g_ = par_tup_
        par_str_ = f'prvl{prvl_}_pqtl{pqtl_}_alpha-{neg_alpha_}.{assoc_test_}'

    if args.hm3_only:
        par_str_ = 'hm3.' + par_str_

    h2_target = float(h2g_)
    dir_gene = f'{args.out_dir}/{gene}'
    result = pd.Series(dtype=float, index=pd.MultiIndex.from_product(
        (methods, R, ['h2_mean', 'h2_median', 'h2_mrb', 'h2_sd', 'se_mean', 'se_mrb'])))

    def empty_results():
        h2_ = np.empty(n_rep)
        h2_se_ = np.empty(n_rep)
        h2_[:] = np.nan
        h2_se_[:] = np.nan
        return h2_, h2_se_

    def fill_result(method, ldmatrix, h2_, h2_se_):
        if 0 < np.isnan(h2_).sum() < h2_.shape[0]:
            h2_ = h2_[~np.isnan(h2_)]
        if 0 < np.isnan(h2_se_).sum() < h2_se_.shape[0]:
            h2_se_ = h2_se_[~np.isnan(h2_se_)]
        if args.binary_assoc_test:
            h2_ *= coeffcients[prvl_]
            h2_se_ *= coeffcients[prvl_]

        result[(method, ldmatrix, 'h2_mean')] = np.mean(h2_)
        result[(method, ldmatrix, 'h2_median')] = np.median(h2_)
        result[(method, ldmatrix, 'h2_mrb')] = result[(method, ldmatrix, 'h2_mean')] / h2_target - 1

        result[(method, ldmatrix, 'h2_sd')] = np.std(h2_, ddof=1)
        if result[(method, ldmatrix, 'h2_sd')] == 0:
            result[(method, ldmatrix, 'h2_sd')] = np.nan

        result[(method, ldmatrix, 'se_mean')] = np.mean(h2_se_)
        result[(method, ldmatrix, 'se_mrb')] \
            = result[(method, ldmatrix, 'se_mean')] / result[(method, ldmatrix, 'h2_sd')] - 1

    if args.gcta:
        h2, h2_se = empty_results()
        for j in range(n_rep):
            with open(f'{dir_gene}/rep{j}/{par_str_}/gcta.hsq') as I:
                lines = I.readlines()
            line4 = lines[4].strip().split('\t')
            if line4[0] == 'V(G)/Vp':
                h2[j] = float(line4[1])
                h2_se[j] = float(line4[2])
        fill_result('GCTA', 'Rgwa', h2, h2_se)

    for r, k in zip(R, ld_sfx):
        if args.kggsee:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/{par_str_}/kggsee_{k}.gene.pvalue.txt', sep='\t', index_col=0)
                h2[j] = df.loc[gene, 'Herit']
                h2_se[j] = df.loc[gene, 'HeritSE']
            fill_result('EHE', r, h2, h2_se)

        if args.hess:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/{par_str_}/hess_{k}.step2.txt', sep='\t').iloc[0]
                h2[j] = df['local_h2g']
                h2_se[j] = df['se']
            fill_result('HESS', r, h2, h2_se)

        if args.ldak_gbat:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/{par_str_}/ldak_gbat_{k}/remls.1', sep=' ', index_col=1)
                h2[j] = df.loc[f'1_{gene}', 'Heritability']
                h2_se[j] = df.loc[f'1_{gene}', 'SD']
            fill_result('GBAT', r, h2, h2_se)

        if args.ldak_sumher:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/{par_str_}/ldak_sumher_{k}.hers', sep=' ', index_col=0)
                h2[j] = df.loc['Her_All', 'Heritability']
                # h2_se[j] = df.loc['Her_All', 'SD']
            fill_result('SumHer', r, h2, h2_se)

        if args.ldsc:
            regex = re.compile(r'^Total Observed scale h2: (-?\d\.?\d*(?:e-\d+)?) \((\d\.?\d*(?:e-\d+)?)\)$')
            h2, h2_se = empty_results()
            for j in range(n_rep):
                with open(f'{dir_gene}/rep{j}/{par_str_}/ldsc_{k}.log') as I:
                    for df in I:
                        match = regex.match(df.strip())
                        if match:
                            h2[j] = float(match.group(1))
                            # h2_se[j] = float(match.group(2))
                            break
            fill_result('LDSC', r, h2, h2_se)

        if args.lder:
            regex = re.compile(r'^(-?\d(?:\.\d+)?(?:e-\d+)?)$')
            h2, h2_se = empty_results()
            for j in range(n_rep):
                with open(f'{dir_gene}/rep{j}/{par_str_}/lder_{k}.result') as I:
                    for df in I:
                        match = regex.match(df.strip())
                        if match:
                            h2[j] = float(match.group(1))
                            break
            fill_result('LDER', r, h2, h2_se)

    result.name = gene
    return result


def harvest_one_gene(gene):
    logging.info(f'Harvesting {gene} ...')
    one_gene_dict_ = dict()
    for par_tup_ in par_tup_list:
        one_gene_dict_[par_tup_] = read_one_par_set_result(gene, par_tup_)

    return one_gene_dict_


all_genes_dict = dict()
for one_gene_dict in Pool(args.nt).map(harvest_one_gene, gene_list):
    for par_tup in par_tup_list:
        if par_tup not in all_genes_dict:
            all_genes_dict[par_tup] = list()
        all_genes_dict[par_tup].append(one_gene_dict[par_tup])

for par_tup, result_list in all_genes_dict.items():
    if not args.binary_assoc_test:
        pqtl, neg_alpha, h2g = par_tup
        par_str = f'pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}'
    else:
        prvl, pqtl, neg_alpha, assoc_test, h2g = par_tup
        par_str = f'prvl{prvl}_pqtl{pqtl}_alpha-{neg_alpha}.{assoc_test}'
    if args.hm3_only:
        par_str = 'hm3only.' + par_str

    logging.info(f'Merging {par_str} results ...')
    results = pd.concat(result_list, axis=1).T
    results_name = f'{args.out_dir}/results.{par_str}.{args.out_suffix}'
    results.to_csv(results_name, sep='\t')
    logging.info(f'Written {results_name}')

logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
