#!/usr/bin/env python

import time
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis
from pysnptools.snpreader import Bed
from multiprocessing.dummy import Pool

parser = argparse.ArgumentParser(description='Perform tests of H0:H/se(H)=0')
parser.add_argument('--nt', type=int, default=1, help='Number of threads running in parallel')

args = parser.parse_args()
out_dir = '/home/ml/work/q.EHE_paper/h.hapsim_h2/d.n2e4_2e3_5e2_rep100_chr1'
gene_list = np.loadtxt('/home/ml/work/q.EHE_paper/h.hapsim_h2/a.coding_gene_bed/chr1.one4th.lst', dtype=str)

n_rep = 100
n_gwa = 20000
tot_n_gwa = n_rep * n_gwa
n_dist_rep = int(2e5-1)
pqtl = '0.1'
neg_alpha = '0.25'


def process_one_gene(gene):
    dir_gene = f'{out_dir}/{gene}'

    # Read plink.bim and plink.fam
    plink_bim = pd.read_csv(f'{dir_gene}/plink.bim', sep='\t', header=None)
    plink_bim.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']

    # Read genotypes and calculate MAF
    snp_on_disk = Bed(f'{dir_gene}/plink', count_A1=True)
    X = snp_on_disk.read().val[:tot_n_gwa]
    m = X.shape[1]
    frq_A1 = np.mean(X, axis=0) / 2
    X -= X.mean(axis=0)
    X /= X.std(axis=0)
    R = X.T @ X / tot_n_gwa  # m * m
    Omega = R ** 2  # n_rep * m * m
    Omega_inv = np.linalg.pinv(Omega, rcond=1e-10, hermitian=True)  # m * m
    z_gwa_null = np.random.multivariate_normal(np.zeros(m), R, n_dist_rep)[:, :, None]  # n_dist_rep * m * 1
    z_ehe_null = np.sum(Omega_inv @ (z_gwa_null ** 2 - 1), axis=(1, 2)) / np.sqrt(2 * np.sum(Omega_inv))  # n_dist_rep

    # Realize beta
    beta = np.zeros(m)
    m_qtl = int(np.ceil(m * float(pqtl)))
    idx_qtl = np.random.choice(m, size=m_qtl, replace=False)
    qtl_beta_var = (frq_A1[idx_qtl] * (1 - frq_A1[idx_qtl])) ** (1 - float(neg_alpha))
    beta[idx_qtl] = np.random.multivariate_normal(np.zeros(m_qtl), np.diag(qtl_beta_var))  # m_qtl
    genetic_eff = X @ beta
    genetic_eff_std = genetic_eff.std()

    for h2g in ['1e-2', '1e-3', '1e-4']:
        # Scale beta to fit target h2g
        scale_factor = float(h2g) ** 0.5 / genetic_eff_std
        genetic_eff_scale_h2g = genetic_eff * scale_factor

        # Calculate phenotypes and perform association tests
        y = genetic_eff_scale_h2g + np.random.normal(0, np.sqrt(1 - float(h2g)), tot_n_gwa)
        z = np.swapaxes(X.reshape(n_rep, n_gwa, m), 1, 2) @ y.reshape(n_rep, n_gwa, 1) / np.sqrt(n_gwa)  # n_rep * m * 1

        # Calculate EHE z-score
        chi2cov = 4 * R * (z @ np.swapaxes(z, 1, 2)) - 2 * Omega  # n_rep * m * m
        h2 = np.sum(Omega_inv @ (z ** 2 - 1), axis=(1, 2)) / n_gwa  # n_rep
        h2var = np.sum(Omega_inv @ chi2cov @ Omega_inv, axis=(1, 2)) / n_gwa ** 2  # n_rep
        i = np.where(h2var > 0)[0]
        h2 = h2[i]
        z_ehe = h2 / h2var[i] ** 0.5  # n_rep

        # Calculate EHE p-value
        p_vals = list()
        for z_ehe_i in z_ehe:
            rank_in_null = np.concatenate((np.array([z_ehe_i]), z_ehe_null)).argsort().argmin() + 1
            p_vals.append(min(rank_in_null, 2 + n_dist_rep - rank_in_null) / (1 + n_dist_rep) * 2)

        df = pd.DataFrame({'ehe': h2, 'z': z_ehe, 'p': p_vals})
        df.to_csv(f'{dir_gene}/ehe_pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}.tsv',
                  sep='\t', index=False, float_format='%.4g')
        logging.info(f'Done {gene}.')
    return pd.DataFrame({'mean': z_ehe_null.mean(), 'variance': z_ehe_null.var(),
                         'skew': skew(z_ehe_null), 'kurtosis': kurtosis(z_ehe_null)}, index=[gene])


logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
logging.info(f'Getting started at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
moments = pd.concat(Pool(args.nt).map(process_one_gene, gene_list))
moments.to_csv('/home/ml/work/q.EHE_paper/h.hapsim_h2/h.ehe_z_null_moments.tsv', sep='\t', float_format='%.4g')
logging.info(f'Done at {time.strftime("%d %b %Y %H:%M:%S", time.localtime())}')
