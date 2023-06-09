#!/usr/bin/env python

import pandas as pd
from settings import *
from itertools import product

h2f = 0.001
tests = ['logit', 'chi2']

R = R[1:]
C = C[1:]
legend_label = legend_label[1:]

w = 1 / len(R) - 0.15
s = np.arange(len(R)) * w * 1.2  # shift
s -= s.mean()
bp_style['widths'] = w

results = {test: {q: {prvl: {} for prvl in prvl_list} for q in quantiles} for test in tests}
for test, q, prvl, pqtl in product(tests, quantiles, prvl_list, pqtl_list):
    results[test][q][prvl][pqtl] = pd.read_csv(f'/home/ml/work/q.EHE_paper/h.hapsim_h2/g.results/'
                                               f'results.prvl{prvl}_pqtl{pqtl}_alpha-0.25.{test}.one4th_{q}.tsv',
                                               sep='\t', header=[0, 1, 2], index_col=0)

# Plot ['h2_mean', 'h2_mrb', 'h2_sd']
q_prvl_tup_list = list(product(quantiles, prvl_list))
nrow = len(q_prvl_tup_list)
ncol = len(pqtl_list)
abc = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWX')).reshape(nrow, ncol)

for test in tests:
    for est, ylabel in est_ylabel.items():
        print(f'Plotting {test}, {est} ...')
        if est == 'se_mrb':
            figsize = (9, 16)
            prog_plot = prog_local_h2
            text_bbox['x'] = 0.03
        else:
            figsize = (12, 16)
            prog_plot = programs
            text_bbox['x'] = 0.02

        fig, ax = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
        fig.subplots_adjust(left=0.07, bottom=0.03, right=0.99, top=0.99, wspace=0, hspace=0)
        x = np.arange(len(prog_plot))
        bp = dict()
        for i, (q, prvl) in enumerate(q_prvl_tup_list):
            for j, pqtl in enumerate(pqtl_list):
                df = results[test][q][prvl][pqtl]

                for u in [0, 1]:
                    bp[u] = ax[i, j].boxplot(df.loc[:, (prog_plot, R[u], est)], positions=x + s[u],
                                             boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                    if df.loc[:, ('LDER', R[u], est)].isnull().sum() > 0:
                        ax[i, j].boxplot(df.loc[:, ('LDER', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                         boxprops=dict(facecolor=C[u], lw=lw), **bp_style)

                if i == 0 and j == 0:
                    ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}, α=-0.25$\n'
                                    f'$K={prvl}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                elif i == 0 and j != 0:
                    ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}, α=-0.25$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                elif i != 0 and j == 0:
                    ax[i, j].text(s=f'{abc[i, j]}. $K={prvl}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
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
                if est == 'h2_mrb' or est == 'se_mrb':
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
        elif est == 'se_mrb':
            ax[0, 0].set_ylim(-0.9, 2.3)

        ax[0, 2].legend([bp[0]['boxes'][0], bp[1]['boxes'][0]], legend_label, title=legend_title,
                        loc='upper right', labelspacing=0.2)
        prefix = f'/home/ml/work/q.EHE_paper/h.hapsim_h2/g.plots/fig.binary_{test}.n2e4.488genes.100reps.{est}.h2g1e-3'
        fig.savefig(f'{prefix}.tif', dpi=500, pil_kwargs={'compression': 'tiff_lzw'})
        fig.savefig(f'{prefix}.svg')
        plt.close()
