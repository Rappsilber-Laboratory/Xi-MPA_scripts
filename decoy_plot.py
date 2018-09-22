import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl
from matplotlib.ticker import MaxNLocator

import matplotlib as mpl
import seaborn as sns


plot_style = {'axes.axisbelow': True,
            'axes.grid': False,
            'axes.labelcolor': '.15',
            'axes.linewidth': 1,
            'font.family': ['sans-serif'],
            'font.sans-serif': ['Arial',
                                'Liberation Sans',
                                'Bitstream Vera Sans',
                                'sans-serif'],
            'image.cmap': 'Greys',
            'legend.frameon': False,
            'legend.numpoints': 1,
            'legend.scatterpoints': 1,
              'lines.solid_capstyle': 'round',
              'xtick.direction': 'out',
              'xtick.major.size': 8,
            'xtick.minor.size': 5,
            'ytick.direction': 'out',
            'ytick.major.size': 8,
            'ytick.minor.size': 5,
            'legend.fontsize': 22,
            'font.size': 20,
            'axes.titlesize': 22,
            'xtick.labelsize': 22,
            'ytick.labelsize': 22,
            'axes.labelsize': 22,
            'svg.fonttype': 'none'}

# use seaborn default settings
sns.set_context("notebook")
sns.set_style("white")

# update the parameters in plm
mpl.rcParams.update(plot_style)
mpl.rcParams.update({'font.size': 22})

sns.set(context="poster", style="ticks", palette='deep',
        rc={"xtick.major.size": 15, "axes.labelsize": 25,
            "xtick.labelsize": 22, "ytick.labelsize": 22,
            "ytick.major.size": 15, "legend.fancybox": True})

mpl.rcParams['svg.fonttype'] = 'none'


def sort_to_td_group(row):
    if (row['Protein1decoy'] == 0) & (row['Protein2decoy'] == 0):
        return 'TT'
    elif ((row['Protein1decoy'] == 1) & (row['Protein2decoy'] == 0)) | ((row['Protein1decoy'] == 0) & (row['Protein2decoy'] == 1)):
        return 'TD'
    elif (row['Protein1decoy'] == 1) & (row['Protein2decoy'] == 1):
        return 'DD'


def sort_to_mass_groups(row):
    if row['mps_decoy_mip'] == 0:
        return '0 Da'
    elif row['mps_decoy_mip'] > 0:
        return '+1 Da to +4 Da'
    elif row['mps_decoy_mip'] < 0:
        return '-1 Da to -4 Da'


base_dir = ''
prec_df = ''

res = []
groups = ['+1 Da to +4 Da', '0 Da', '-1 Da to -4 Da']
df_all = pd.DataFrame()
for plot_file, identifier in [('decoy_4_mscon_PF_20_100_0_B160803_02', 'B160803_02'),
                              ('decoy_4_mscon_PF_20_100_0_B160803_07', 'B160803_07'),
                              ('decoy_4_mscon_PF_20_100_0_B160803_12', 'B160803_12')]:
    score_cutoff = 4.91344662865

    mps_df = pd.read_csv(prec_df, index_col=0)
    mps_df = mps_df[mps_df.index.str.contains(identifier)]
    mps_df['Scan'] = [int(x[11:]) for x in mps_df.index]

    xi_df = pd.read_csv(os.path.join(base_dir, 'Xi', 'xi_filtered', 'decoy_4', plot_file + '.csv'))
    xi_df = xi_df[xi_df['Fasta1'] == xi_df['Fasta2']]
    xi_df = xi_df.merge(mps_df[['Scan', 'mps_decoy_mip']], on='Scan')
    xi_df['group'] = xi_df.apply(sort_to_mass_groups, axis=1)
    xi_df['fdr_group'] = xi_df.apply(sort_to_td_group, axis=1)
    df_all = pd.concat([df_all, xi_df])

    if identifier == 'B160803_02':
        bins_from, bins_to = int(min(xi_df['match score']) - 0.5), int(max(xi_df['match score']))
        bins = np.arange(bins_from-1, bins_to+1, 0.5)


fig, ax = plt.subplots()
offset = 0
sns.distplot(df_all['match score'][(df_all['fdr_group'].isin(['TD', 'DD']))],
             bins=bins + offset, label='decoys', kde=False, ax=ax, norm_hist=True, color='#d95f02',
             hist_kws={'histtype': 'step', 'linewidth': 3, 'alpha': 0.4})

for group in groups:
    sns.distplot(df_all['match score'][(df_all['group'] == group) & (df_all['fdr_group'] == 'TT')], norm_hist=True,
                 kde=False, bins=bins+offset, label=group,
                 hist_kws={'histtype': 'step', 'linewidth': 3, 'alpha': 0.8}, ax=ax)
    offset+=0.1

df_all = df_all[['match score', 'psmid', 'mps_decoy_mip', 'group', 'fdr_group']]
df_all.to_csv(base_dir + '/decoy_search/decoysearch.csv')

ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
sns.despine()
ax.set_ylabel('frequency (%)')
plt.axvline(score_cutoff, color='black')
plt.legend()
plt.savefig(base_dir + '/decoy_search/decoy_histabs_norm.svg')
plt.savefig(base_dir + '/decoy_search/decoy_histabs_norm.png')
