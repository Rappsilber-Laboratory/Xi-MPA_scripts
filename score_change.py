import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
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
outdir = ''
prec_df_path = ''

prec_mps_df = pd.read_csv(prec_df_path, index_col=0)

hits_shifted = prec_mps_df[['Xi_mip', 'Xi_hit']][(prec_mps_df['Xi_mip'] != 0)]
hits_shifted['scan_id'] = hits_shifted.index

tmp = []
for xi_res in os.listdir(os.path.join(base_dir, 'Xi', 'xi', 'other_NoPrePro')):
    if xi_res not in ['other_NoPrePro_B160803_02.csv', 'other_NoPrePro_B160803_07.csv', 'other_NoPrePro_B160803_12.csv']: continue
    file_id = '_'.join(xi_res[:-4].split('_')[-2:])
    df = pd.read_csv(os.path.join(base_dir, 'Xi', 'xi', 'other_NoPrePro', xi_res))
    df['scan_id'] = df['Scan'].apply(lambda x: file_id + '_' + str(x))
    df['fdrgroup_npp'] = df.apply(sort_to_td_group, axis=1)
    tmp.append(df[['scan_id', 'match score', 'fdrgroup_npp']])

df = pd.concat(tmp)
hits_shifted = hits_shifted.merge(df, how='outer', on='scan_id')

hits_shifted = hits_shifted.rename(columns={'match score': 'score unprocessed'})

tmp = []
for xi_res in os.listdir(os.path.join(base_dir, 'Xi', 'xi', 'rel_4_mscon_PF_20_100_0')):
    if xi_res not in ['rel_4_mscon_PF_20_100_0_B160803_02.csv', 'rel_4_mscon_PF_20_100_0_B160803_07.csv', 'rel_4_mscon_PF_20_100_0_B160803_12.csv']: continue

    file_id = '_'.join(xi_res[:-4].split('_')[-2:])
    df = pd.read_csv(os.path.join(base_dir, 'Xi', 'xi', 'rel_4_mscon_PF_20_100_0', xi_res))
    df['scan_id'] = df['Scan'].apply(lambda x: file_id + '_' + str(x))
    df['fdrgroup'] = df.apply(sort_to_td_group, axis=1)
    tmp.append(df[['scan_id', 'match score', 'fdrgroup']])

df = pd.concat(tmp)
hits_shifted = hits_shifted.merge(df, how='outer', on='scan_id')

decoys = hits_shifted[(hits_shifted['fdrgroup'].isin(['TD', 'DD'])) & (hits_shifted['fdrgroup_npp'].isin(['TD', 'DD']))].copy() #
decoys.to_csv(outdir + '/scores_mps_unprocessed_decoy.csv', index=None)
decoys = decoys[~decoys['score unprocessed'].isnull()]
decoys.drop_duplicates(inplace=True)

hits_shifted = hits_shifted[hits_shifted['Xi_hit'] == True]
hits_shifted.to_csv(outdir + '/scores_mps_unprocessed.csv', index=None)
hits_shifted = hits_shifted[~hits_shifted['score unprocessed'].isnull()]
hits_shifted.drop_duplicates(inplace=True)

print stats.ttest_rel(hits_shifted['match score'], hits_shifted['score unprocessed'])[1]

bins = np.arange(min(hits_shifted['score unprocessed']), max(hits_shifted['match score']) + 0.5, 0.5)
fig, ax = plt.subplots()
sns.distplot(hits_shifted['match score'], label='Xi-MPA', kde=False, norm_hist=True,
             hist_kws={'histtype': 'step', 'linewidth': 3, 'alpha': 0.9}, ax=ax, bins=bins+0.1, color='#69afa3')
sns.distplot(hits_shifted['score unprocessed'], label='unprocessed', kde=False, norm_hist=True,
             hist_kws={'histtype': 'step', 'linewidth': 3, 'alpha': 0.8}, ax=ax, bins=bins)
sns.distplot(decoys['match score'], label='Xi-MPA decoys', kde=False, norm_hist=True,
             hist_kws={'histtype': 'step', 'linewidth': 3, 'alpha': 0.8}, ax=ax, bins=bins+0.1, color='#d95f02')
sns.distplot(decoys['score unprocessed'], label='unprocessed decoys', kde=False, norm_hist=True,
             hist_kws={'histtype': 'step', 'linewidth': 3, 'alpha': 0.4}, ax=ax, bins=bins + 0.1, color='#d95f02')

hits_shifted.to_csv(outdir + '/scores_paired.csv')
ax.set_ylabel('frequency (%)')
ax.set_xlabel('score')
ax.yaxis.set_major_locator(MaxNLocator(nbins=3))

sns.despine()
plt.legend()
plt.tight_layout()
plt.savefig(outdir + '/scores_mps_unprocessed.png')
plt.savefig(outdir + '/scores_mps_unprocessed.svg')
plt.close()
