#!/usr/bin/env python

import pandas as pd
from itertools import product
from settings import *


h2g_list = ['1e-2', '1e-3']
neg_alpha_list = ['1.0', '0.25']
par_tup_list = list(product(h2g_list, quantiles, neg_alpha_list, pqtl_list))

results = {h2: {q: {neg_alpha: {} for neg_alpha in neg_alpha_list} for q in quantiles} for h2 in h2g_list}
for h2, q, neg_alpha, pqtl in par_tup_list:
    results[h2][q][neg_alpha][pqtl] = pd.read_csv(f'/home/ml/work/q.EHE_paper/h.hapsim_h2/g.results/'
                                                  f'results.pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2}.one4th_{q}.tsv',
                                                  sep='\t', header=[0, 1, 2], index_col=0)

# Plot ['h2_mean', 'h2_mrb', 'h2_sd']
q_alpha_tup_list = list(product(quantiles, neg_alpha_list))
nrow = len(q_alpha_tup_list)
ncol = len(pqtl_list)
abc = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWX')).reshape(nrow, ncol)

for h2 in h2g_list:
    h2f = float(h2)
    for est, ylabel in est_ylabel.items():
        if est == 'se_mrb':
            continue

        print(f'Plotting h2={h2}, {est} ...')
        figsize = (12, 16)
        fig, ax = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
        fig.subplots_adjust(left=0.07, bottom=0.03, right=0.99, top=0.99, wspace=0, hspace=0)
        prog_plot = programs

        x = np.arange(len(prog_plot))
        bp = dict()
        for i, (q, neg_alpha) in enumerate(q_alpha_tup_list):
            for j, pqtl in enumerate(pqtl_list):
                df = results[h2][q][neg_alpha][pqtl]

                for u in [0, 1, 2]:
                    bp[u] = ax[i, j].boxplot(df.loc[:, (prog_plot, R[u], est)], positions=x + s[u],
                                             boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                    if df.loc[:, ('LDER', R[u], est)].isnull().sum() > 0:
                        ax[i, j].boxplot(df.loc[:, ('LDER', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                         boxprops=dict(facecolor=C[u], lw=lw), **bp_style)

                if i == 0 and j == 0:
                    ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}$\n'
                                    f'$α=-{neg_alpha}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                elif i == 0 and j != 0:
                    ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                elif i != 0 and j == 0:
                    ax[i, j].text(s=f'{abc[i, j]}. $α=-{neg_alpha}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                else:
                    ax[i, j].text(s=f'{abc[i, j]}', transform=ax[i, j].transAxes, **text_bbox)

                if i == nrow - 1:
                    ax[i, j].set_xticks(x)
                    ax[i, j].set_xticklabels(prog_plot)

        for i in range(nrow):
            for j in range(ncol):
                if est == 'h2_mean':
                    ax[i, j].hlines(h2f, xmin=-0.5, xmax=x[-1] + .5, colors='gray', linestyles='dashed', lw=1)
                if est == 'h2_mrb':
                    ax[i, j].hlines(0, xmin=-0.5, xmax=x[-1] + .5, colors='gray', linestyles='dashed', lw=1)
                if j == 0:
                    ax[i, j].set_ylabel(ylabel)
                    if est == 'h2_mean' or est == 'h2_sd':
                        if est == 'h2_mean':
                            ax[0, 0].set_ylim(-h2f / 3, h2f * 2.9)
                        elif est == 'h2_sd':
                            ax[0, 0].set_ylim(0, h2f * 1.99)
                        yticks = ax[i, j].get_yticks()
                        yticklabels = [f'{a * 100:.2f}%' for a in yticks]
                        ax[i, j].set_yticks(yticks)
                        ax[i, j].set_yticklabels(yticklabels)

        if est == 'h2_mean':
            ax[0, 0].set_ylim(-h2f / 3, h2f * 2.9)
        elif est == 'h2_sd':
            ax[0, 0].set_ylim(0, h2f * 1.99)
        elif est == 'h2_mrb':
            ax[0, 0].set_ylim(-1.4, 1.9)

        ax[0, 2].legend([bp[0]['boxes'][0], bp[1]['boxes'][0], bp[2]['boxes'][0]], legend_label, title=legend_title,
                        loc='upper right', labelspacing=0.2)
        prefix = f'/home/ml/work/q.EHE_paper/h.hapsim_h2/g.plots/fig.n2e4.488genes.100reps.{est}.h2g{h2}'
        fig.savefig(f'{prefix}.svg')
        fig.savefig(f'{prefix}.tif', dpi=500, pil_kwargs={'compression': 'tiff_lzw'})
        plt.close()


# Plot 'se_mrb'
p_alpha_tup_list = list(product(pqtl_list, neg_alpha_list))
nrow = len(quantiles)
ncol = len(p_alpha_tup_list)
abc = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWX')).reshape(nrow, ncol)
text_bbox['x'] = 0.04
est = 'se_mrb'
ylabel = r'MRB of SE($\hat{h}^2)$'

for h2 in h2g_list:
    h2f = float(h2)
    print(f'Plotting h2={h2}, {est} ...')

    figsize = (12, 9)
    fig, ax = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, bottom=0.03, right=0.99, top=0.99, wspace=0, hspace=0)
    prog_plot = prog_local_h2

    x = np.arange(len(prog_plot))
    bp = dict()
    for i, q in enumerate(quantiles):
        for j, (pqtl, neg_alpha) in enumerate(p_alpha_tup_list):
            df = results[h2][q][neg_alpha][pqtl]

            for u in [0, 1, 2]:
                bp[u] = ax[i, j].boxplot(df.loc[:, (prog_plot, R[u], est)], positions=x + s[u],
                                         boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                if df.loc[:, ('GBAT', R[u], est)].isnull().sum() > 0:
                    ax[i, j].boxplot(df.loc[:, ('GBAT', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                     boxprops=dict(facecolor=C[u], lw=lw), **bp_style)

            if i == 0 and j == 0:
                ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}$, $α=-{neg_alpha}$\n '
                                f'${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                              transform=ax[i, j].transAxes, **text_bbox)
            elif i == 0 and j != 0:
                ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}$, $α=-{neg_alpha}$',
                              transform=ax[i, j].transAxes, **text_bbox)
            elif i != 0 and j == 0:
                ax[i, j].text(s=f'{abc[i, j]}. ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                              transform=ax[i, j].transAxes, **text_bbox)
            else:
                ax[i, j].text(s=f'{abc[i, j]}', transform=ax[i, j].transAxes, **text_bbox)

            if i == nrow - 1:
                ax[i, j].set_xticks(x)
                ax[i, j].set_xticklabels(prog_plot)
        ax[i, 0].set_ylabel(ylabel)

    ax[0, 0].set_ylim(-1.2, 2.7)
    for i in range(nrow):
        ax[i, 0].set_ylabel(ylabel)
        for j in range(ncol):
            ax[i, j].hlines(0, xmin=-0.5, xmax=x[-1] + .5, colors='gray', linestyles='dashed', lw=1)

    ax[0, 5].legend([bp[0]['boxes'][0], bp[1]['boxes'][0], bp[2]['boxes'][0]], legend_label, title=legend_title,
                    loc=(0.01, 0.64), labelspacing=0.2)
    prefix = f'/home/ml/work/q.EHE_paper/h.hapsim_h2/g.plots/fig.n2e4.488genes.100reps.{est}.h2g{h2}'
    fig.savefig(f'{prefix}.svg')
    fig.savefig(f'{prefix}.tif', dpi=500, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
