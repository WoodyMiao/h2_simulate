#!/usr/bin/env python

import time
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm
from pysnptools.snpreader import Bed
from multiprocessing.dummy import Pool

parser = argparse.ArgumentParser(description='Realize effect sizes and perform association tests')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')
parser.add_argument('--gene-list', type=str, default=None, help='File of a list of genes to be included')
parser.add_argument('--nt', type=int, default=1, help='Number of threads running in parallel')
parser.add_argument('--h2g-vals', type=str, default='1e-2,1e-3', help='Target h2g values for the simulation')
parser.add_argument('--pqtl-vals', type=str, default='1.0,0.1', help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--neg-alpha-vals', type=str, default='0.25,1.0', help='Power values in the LDAK-Thin Model')
parser.add_argument('--prevalence-vals', type=str, default='0.1,0.01',
                    help='prevalence values for the liability threshold modeling of dichotamous phenotypes')
parser.add_argument('--n-gwa', type=float, default=2e4, help='Size of a sample for association tests')
parser.add_argument('--n-rep', type=float, default=100, help='Number of repititions')
parser.add_argument('--hm3-only', action='store_true', default=False,
                    help='Use only HapMap3 SNPs and prepend "hm3." to the output file names')

args = parser.parse_args()
n_gwa = int(args.n_gwa)
n_rep = int(args.n_rep)
out_dir = args.out_dir
h2g_list = args.h2g_vals.split(',')
pqtl_list = args.pqtl_vals.split(',')
neg_alpha_list = args.neg_alpha_vals.split(',')
prevalence_list = args.prevalence_vals.split(',')

if not args.gene_list:
    snp_counts = pd.read_csv(f'{args.out_dir}.snp_counts.tsv', sep='\t', index_col=0)
    if not args.hm3_only:
        gene_list = snp_counts.loc[snp_counts.allSNP >= 3].index.values
    else:
        gene_list = snp_counts.loc[snp_counts.hm3SNP >= 3].index.values
else:
    gene_list = np.loadtxt(args.gene_list, dtype=str)


def process_one_gene(gene):
    dir_gene = f'{out_dir}/{gene}'

    # Read plink.bim and plink.fam
    plink_bim = pd.read_csv(f'{dir_gene}/plink.bim', sep='\t', header=None)
    plink_bim.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']

    # Read genotypes and calculate MAF
    snp_on_disk = Bed(f'{dir_gene}/plink', count_A1=True)

    if args.hm3_only:
        hm3snp = np.loadtxt(f'{dir_gene}/plink.hm3snp', dtype=str)
        is_hm3snp = plink_bim.SNP.isin(hm3snp)
        plink_bim = plink_bim[is_hm3snp]
        snp_on_disk = snp_on_disk[:, is_hm3snp]

    snpdata = snp_on_disk.read()
    frq_A1 = np.mean(snpdata.val, axis=0) / 2

    # Standardize genotypes
    snpdata.standardize()
    X = snpdata.val
    n_pop, m = X.shape
    tot_n_gwa = n_rep * n_gwa
    X_gwa = X[:tot_n_gwa].reshape(n_rep, n_gwa, m)
    iid_gwa = snpdata.iid[:tot_n_gwa].reshape(n_rep, n_gwa, 2)

    columns = pd.MultiIndex.from_product((pqtl_list, neg_alpha_list, h2g_list))
    realized_beta = pd.DataFrame(dtype=float, index=np.arange(m), columns=columns)
    realized_beta.columns.names = ['prop_qtl', 'neg_alpha', 'target_h2']
    realized_h2_ = pd.DataFrame(dtype=float, index=[gene], columns=columns)
    realized_h2_.columns.names = ['prop_qtl', 'neg_alpha', 'target_h2']

    columns = pd.MultiIndex.from_product((pqtl_list, neg_alpha_list, h2g_list, ['cont'] + prevalence_list))
    phenotypes = pd.DataFrame(-9, dtype=np.int8, index=snpdata.iid[:, 0], columns=columns)
    phenotypes.columns.names = ['prop_qtl', 'neg_alpha', 'target_h2', 'prevalence']
    for pqtl in pqtl_list:
        m_qtl = int(np.ceil(m * float(pqtl)))
        idx_qtl = np.random.choice(m, size=m_qtl, replace=False)
        for neg_alpha in neg_alpha_list:
            # Realize beta
            beta = np.zeros(m)
            qtl_beta_var = (frq_A1[idx_qtl] * (1 - frq_A1[idx_qtl])) ** (1 - float(neg_alpha))
            beta[idx_qtl] = np.random.multivariate_normal(np.zeros(m_qtl), np.diag(qtl_beta_var))  # m_qtl
            genetic_eff = X @ beta
            genetic_eff_std = genetic_eff.std()
            for h2g in h2g_list:
                # Scale beta to fit target h2g
                scale_factor = float(h2g) ** 0.5 / genetic_eff_std
                genetic_eff_scale_h2g = genetic_eff * scale_factor
                # Calculate phenotypes and realized h2g
                y_norm = genetic_eff_scale_h2g + np.random.normal(0, np.sqrt(1 - float(h2g)), n_pop)
                realized_beta[(pqtl, neg_alpha, h2g)] = beta * scale_factor
                realized_h2_[(pqtl, neg_alpha, h2g)] = genetic_eff_scale_h2g.var() / y_norm.var()

                thresholds = norm.isf(np.float64(prevalence_list))
                phenotypes.loc[:, (pqtl, neg_alpha, h2g, prevalence_list)] = np.int8(y_norm[:, None] > thresholds) + 1
                phenotypes[(pqtl, neg_alpha, h2g, 'cont')] = y_norm
                # Perform association tests
                y_test = y_norm[:tot_n_gwa].reshape(n_rep, n_gwa, 1)
                z = np.swapaxes(X_gwa, 1, 2) @ y_test / np.sqrt(n_gwa)  # n_rep * m * 1
                p = norm.sf(np.abs(z)) * 2

                for j in range(n_rep):
                    if args.hm3_only:
                        out_pre = f'{dir_gene}/rep{j}/hm3.pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}'
                    else:
                        out_pre = f'{dir_gene}/rep{j}/pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}'

                    # Write phenotypes of the sample
                    pd.DataFrame(y_test[j], index=pd.MultiIndex.from_arrays(iid_gwa[j].T)).to_csv(
                        f'{out_pre}.pheno', sep=' ', header=False)

                    # Write a summary statistic file for KGGSEE, HESS, and LDSC
                    sumstat = plink_bim.copy()
                    sumstat['Z'] = z[j]
                    sumstat['N'] = n_gwa
                    sumstat['P'] = p[j]
                    sumstat[['CHR', 'BP', 'P', 'SNP', 'A1', 'A2', 'Z', 'N']] \
                        .to_csv(f'{out_pre}.sumstat.gz', sep='\t', index=False)

                    # Write a summary statistic file for LDER
                    sumstat.rename({'SNP': 'snp', 'CHR': 'chr', 'A1': 'a0', 'A2': 'a1', 'Z': 'z'}, axis=1)[
                        ['snp', 'chr', 'a0', 'a1', 'z']].to_csv(
                        f'{out_pre}.lder.sumstat.gz', sep='\t', index=False)

                    # Write a summary statistic file for LDAK
                    sumstat.rename({'SNP': 'Predictor', 'N': 'n'}, axis=1)[['Predictor', 'A1', 'A2', 'n', 'Z']] \
                        .to_csv(f'{out_pre}.ldak.sumstat', sep='\t', index=False)

    if args.hm3_only:
        realized_beta.to_csv(f'{dir_gene}/hm3.realized_effect_sizes.tsv', sep='\t')
        phenotypes.to_csv(f'{dir_gene}/hm3.realized_phenotypes.tsv', sep='\t')
    else:
        realized_beta.to_csv(f'{dir_gene}/realized_effect_sizes.tsv', sep='\t')
        phenotypes.to_csv(f'{dir_gene}/realized_phenotypes.tsv', sep='\t')

    logging.info(f'Done {gene}.')
    return realized_h2_


logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
realized_h2 = pd.concat(Pool(args.nt).map(process_one_gene, gene_list))
if args.hm3_only:
    realized_h2.to_csv(f'{out_dir}.hm3.realized_h2.tsv', sep='\t')
else:
    realized_h2.to_csv(f'{out_dir}.realized_h2.tsv', sep='\t')
logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
