import numpy as np
import matplotlib.pyplot as plt


quantiles = ['q0-q1', 'q1-q2', 'q2-q3', 'q3-q4']
pqtl_list = ['0.01', '0.1', '1']
prvl_list = ['0.1', '0.01']


n_snps = {'q0-q1': (6, 69), 'q1-q2': (70, 120), 'q2-q3': (121, 222), 'q3-q4': (230, 2094)}
n_bt_quantiles = {'q0-q1': 122, 'q1-q2': 124, 'q2-q3': 121, 'q3-q4': 121}

n_hm3_snps = {'q0-q1': (3, 22), 'q1-q2': (3, 34), 'q2-q3': (4, 68), 'q3-q4': (13, 395)}
n_hm3_bt_quantiles = {'q0-q1': 114, 'q1-q2': 120, 'q2-q3': 119, 'q3-q4': 120}

programs = ['EHE', 'HESS', 'GBAT', 'SumHer', 'LDSC', 'LDER']
prog_local_h2 = ['EHE', 'HESS', 'GBAT']

R = ['Rgwa', 'R2e3', 'R5e2']
C = ['#555555', '#bbbbbb', '#ffffff']
legend_label = ['Same ($n=20000$)', 'Different ($n=2000$)', 'Different ($n=500$)']
legend_title = 'LD coefficient matrix'

w = 1 / len(R) - 0.15
s = np.arange(len(R)) * w * 1.2  # shift
s -= s.mean()
est_ylabel = {'h2_mean': r'Mean($\hat{h}^2$)',
              'h2_mrb': r'MRB of $\hat{h}^2$',
              'h2_sd': r'SD($\hat{h}^2$)',
              'se_mrb': r'MRB of SE($\hat{h}^2)$'}


XSMALL = 7
SMALL = 8
MEDIUM = 9

plt.rc('font', size=SMALL)  # controls default text sizes
plt.rc('axes', titlesize=SMALL)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=SMALL)  # fontsize of the x tick labels
plt.rc('legend', fontsize=XSMALL)  # legend fontsize

lw = 0.5
whiskerprops = dict(color='black', lw=lw)
medianprops = dict(color='black', lw=lw)
capprops = dict(color='black', lw=lw)
meanprops = dict(markersize=5, markeredgecolor='black', markerfacecolor='white', markeredgewidth=lw)
flierprops = dict(markersize=4, marker='o', markerfacecolor='white', markeredgecolor='black', markeredgewidth=lw)

bp_style = dict(showmeans=True, showfliers=True, patch_artist=True, whiskerprops=whiskerprops, widths=w,
                flierprops=flierprops, capprops=capprops, medianprops=medianprops, meanprops=meanprops)
text_bbox = dict(x=0.02, y=0.96, ha='left', va='top', bbox=dict(facecolor='white', alpha=0.9, edgecolor='white'))
text_only = dict(x=0.02, y=0.96, ha='left', va='top')
