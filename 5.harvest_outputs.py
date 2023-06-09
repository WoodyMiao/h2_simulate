#!/usr/bin/env python

import os
import re
import time
import logging
import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Harvest results of h2 estimating programs')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--n-rep', type=float, required=True, help='Number of repititions')
parser.add_argument('--process', type=int, default=10, help='Number of processes running in parallel')
parser.add_argument('--gcta', action='store_true', default=False, help='Harvest GCTA results')
parser.add_argument('--hess', action='store_true', default=False, help='Harvest HESS results')
parser.add_argument('--kggsee', action='store_true', default=False, help='Harvest KGGSEE results')
parser.add_argument('--ldak-gbat', action='store_true', default=False, help='Harvest LDAK-GBAT results')
parser.add_argument('--ldak-sumher', action='store_true', default=False, help='Harvest DAK-SumHer results')
parser.add_argument('--lder', action='store_true', default=False, help='Harvest LDER results')
parser.add_argument('--ldsc', action='store_true', default=False, help='Harvest LDSC results')

args = parser.parse_args()
n_rep = int(float(args.n_rep))
a = os.listdir(args.out_dir)
genes = list()
for b in a:
    c = re.match(r'gene\d+$', b)
    if c:
        genes.append(c.group(0))

methods = list()
if args.gcta:
    methods.append('GCTA')
if args.kggsee:
    methods.append('KGGSEE')
if args.hess:
    methods.append('HESS')
if args.ldsc:
    methods.append('LDSC')
if args.lder:
    methods.append('LDER')
if args.ldak_sumher:
    methods.append('SumHer')
if args.ldak_gbat:
    methods.append('GBAT')

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')


def harvest(gene):
    logging.info(f'Harvesting {gene} ...')
    dir_gene = f'{args.out_dir}/{gene}'
    result = pd.Series(dtype=float, index=pd.MultiIndex.from_product(
        (methods, ['Rgwas', 'Rld'], ['h2_mean', 'h2_median', 'h2_sd', 'se_mean', 'se_mrb'])))

    def fill_result(method, ldmatrix, h2_, h2_se_):
        result[(method, ldmatrix, 'h2_mean')] = np.mean(h2_)
        result[(method, ldmatrix, 'h2_median')] = np.median(h2_)
        result[(method, ldmatrix, 'h2_sd')] = np.std(h2_, ddof=1)
        # result[(method, ldmatrix, 'h2_var')] = np.var(h2_, ddof=1)
        if 0 < np.isnan(h2_se_).sum() < h2_se_.shape[0]:
            h2_se_ = h2_se_[~np.isnan(h2_se_)]
        result[(method, ldmatrix, 'se_mean')] = np.mean(h2_se_)
        # result[(method, ldmatrix, 'var_mean')] = np.mean(h2_se_ ** 2)
        return None

    def empty_results():
        h2_ = np.empty(n_rep)
        h2_se_ = np.empty(n_rep)
        h2_[:] = np.nan
        h2_se_[:] = np.nan
        return h2_, h2_se_

    if args.gcta:
        h2, h2_se = empty_results()
        for j in range(n_rep):
            with open(f'{dir_gene}/rep{j}/gcta.hsq') as I:
                lines = I.readlines()
            line4 = lines[4].strip().split('\t')
            if line4[0] == 'V(G)/Vp':
                h2[j] = float(line4[1])
                h2_se[j] = float(line4[2])
        fill_result('GCTA', 'Rgwas', h2, h2_se)

    for k in ['gwas', 'ld']:
        if args.kggsee:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/kggsee_{k}.gene.pvalue.txt', sep='\t', index_col=0)
                h2[j] = df.loc[gene, 'Herit']
                h2_se[j] = df.loc[gene, 'HeritSE']
            fill_result('KGGSEE', f'R{k}', h2, h2_se)

        if args.hess:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/hess_{k}.step2.txt', sep='\t').iloc[0]
                h2[j] = df['local_h2g']
                h2_se[j] = df['se']
            fill_result('HESS', f'R{k}', h2, h2_se)

        if args.ldsc:
            regex = re.compile(r'^Total Observed scale h2: (-?\d\.?\d*(?:e-\d+)?) \((\d\.?\d*(?:e-\d+)?)\)$')
            h2, h2_se = empty_results()
            for j in range(n_rep):
                with open(f'{dir_gene}/rep{j}/ldsc_{k}.log') as I:
                    for df in I:
                        match = regex.match(df.strip())
                        if match:
                            h2[j] = float(match.group(1))
                            h2_se[j] = float(match.group(2))
                            break
            fill_result('LDSC', f'R{k}', h2, h2_se)

        if args.lder:
            regex = re.compile(r'^(-?\d(?:\.\d+)?(?:e-\d+)?)$')
            h2, h2_se = empty_results()
            for j in range(n_rep):
                with open(f'{dir_gene}/rep{j}/lder_{k}.result') as I:
                    for df in I:
                        match = regex.match(df.strip())
                        if match:
                            h2[j] = float(match.group(1))
                            break
            fill_result('LDER', f'R{k}', h2, h2_se)

        if args.ldak_sumher:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/ldak_sumher_{k}.hers', sep=' ', index_col=0)
                h2[j] = df.loc['Her_All', 'Heritability']
                h2_se[j] = df.loc['Her_All', 'SD']
            fill_result('SumHer', f'R{k}', h2, h2_se)

        if args.ldak_gbat:
            h2, h2_se = empty_results()
            for j in range(n_rep):
                df = pd.read_csv(f'{dir_gene}/rep{j}/ldak_gbat_{k}/remls.1', sep=' ', index_col='Gene_Name')
                h2[j] = df.loc[f'1_{gene}', 'Heritability']
                h2_se[j] = df.loc[f'1_{gene}', 'SD']
            fill_result('GBAT', f'R{k}', h2, h2_se)

    result.name = gene
    return result


results = pd.concat(Pool(args.process).map(harvest, genes), axis=1).T
results.loc[:, (slice(None), slice(None), 'se_mrb')] = results.loc[:, (slice(None), slice(None), 'se_mean')].values / \
                                                       results.loc[:, (slice(None), slice(None), 'h2_sd')].applymap(
                                                           lambda x: np.nan if x == 0 else x).values - 1
results.to_csv(f'{args.out_dir}.harvest', sep='\t')
logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
