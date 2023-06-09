#!/usr/bin/env python

import time
import logging
import argparse
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from multiprocessing.dummy import Pool

parser = argparse.ArgumentParser(description='Perform tests of H0:H/se(H)=0')
parser.add_argument('--nt', type=int, default=1, help='Number of threads running in parallel')

args = parser.parse_args()
out_dir = '/home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1'
gene_list = np.loadtxt('/home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.one4th.lst', dtype=str)

n_rep = 100
n_gwa = 20000
n_dist_rep = 1999

h2g = '1e-3'
pqtl = '0.1'
neg_alpha = '0.25'


def process_one_gene(gene):
    dir_gene = f'{out_dir}/{gene}'

    # Read plink.bim and plink.fam
    plink_bim = pd.read_csv(f'{dir_gene}/plink.bim', sep='\t', header=None)
    plink_bim.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']

    # Read genotypes and calculate MAF
    snp_on_disk = Bed(f'{dir_gene}/plink', count_A1=True)
    snpdata = snp_on_disk.read()
    frq_A1 = np.mean(snpdata.val, axis=0) / 2

    # Standardize genotypes
    snpdata.standardize()
    X = snpdata.val
    n_pop, m = X.shape
    tot_n_gwa = n_rep * n_gwa
    X_gwa = X[:tot_n_gwa].reshape(n_rep, n_gwa, m)

    m_qtl = int(np.ceil(m * float(pqtl)))
    idx_qtl = np.random.choice(m, size=m_qtl, replace=False)

    # Realize beta
    beta = np.zeros(m)
    qtl_beta_var = (frq_A1[idx_qtl] * (1 - frq_A1[idx_qtl])) ** (1 - float(neg_alpha))
    beta[idx_qtl] = np.random.multivariate_normal(np.zeros(m_qtl), np.diag(qtl_beta_var))  # m_qtl
    genetic_eff = X @ beta
    genetic_eff_std = genetic_eff.std()

    # Scale beta to fit target h2g
    scale_factor = float(h2g) ** 0.5 / genetic_eff_std
    genetic_eff_scale_h2g = genetic_eff * scale_factor

    # Calculate phenotypes and perform association tests
    y_norm = genetic_eff_scale_h2g + np.random.normal(0, np.sqrt(1 - float(h2g)), n_pop)
    y_test = y_norm[:tot_n_gwa].reshape(n_rep, n_gwa, 1)
    z = np.swapaxes(X_gwa, 1, 2) @ y_test / np.sqrt(n_gwa)  # n_rep * m * 1

    R = np.swapaxes(X_gwa, 1, 2) @ X_gwa / n_gwa  # n_rep * m * m
    Omega = R ** 2  # n_rep * m * m
    Omega_inv = np.linalg.pinv(Omega, rcond=1e-10, hermitian=True)  # n_rep * m * m
    # z: n_rep * m * 1; R: n_rep * m * m
    chi2cov = 4 * R * (z @ np.swapaxes(z, 1, 2)) - 2 * Omega  # n_rep * m * m
    h2 = np.sum(Omega_inv @ (z ** 2 - 1), axis=(1, 2)) / n_gwa  # n_rep
    h2var = np.sum(Omega_inv @ chi2cov @ Omega_inv, axis=(1, 2)) / n_gwa ** 2  # n_rep
    z_ehe = h2 / h2var ** 0.5  # n_rep

    p_vals = np.empty(n_rep)
    for i in range(n_rep):
        if np.isnan(z_ehe[i]):
            p_vals[i] = np.nan
            continue
        z_gwa_null = np.random.multivariate_normal(np.zeros(m), R[i], n_dist_rep)[:, :, None]
        z_ehe_null = np.sum(Omega_inv[i] @ (z_gwa_null ** 2 - 1), axis=(1, 2)) / np.sqrt(2 * np.sum(Omega_inv[i]))
        rank_in_null = np.concatenate((np.array([z_ehe[i]]), z_ehe_null)).argsort().argmin() + 1
        p_vals[i] = (min(rank_in_null, 2 + n_dist_rep - rank_in_null) / (1 + n_dist_rep) * 2)

    df = pd.DataFrame({'ehe': h2, 'z': z_ehe, 'p': p_vals})
    df.to_csv(f'{dir_gene}/ehe_pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}.tsv', sep='\t', index=False)
    logging.info(f'Done {gene}.')
    return


logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
_ = Pool(args.nt).map(process_one_gene, gene_list)
logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
