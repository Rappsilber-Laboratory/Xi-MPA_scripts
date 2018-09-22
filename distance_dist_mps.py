import pandas as pd
import matplotlib.pyplot as plt
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


def sort_to_groups(row):
    if not row['isTT']:
        return 'decoy'
    elif (row['Xi_mip'] == 0):
        return 0
    elif (row['Xi_mip'] == -1):
        return -1
    elif (row['Xi_mip'] == -2):
        return -2
    elif (row['Xi_mip'] == -3):
        return -3
    elif (row['Xi_mip'] == -4):
        return -4


base_dir = ''
outdir = ''
distance_df_path = ''
prec_df_path = ''

distances_mps = pd.read_csv(distance_df_path)

distances_mps = distances_mps[~distances_mps.dist.isnull()]
distances_mps['scan_id'] = distances_mps.apply(lambda x: x['run'][:10] + '_' + str(x['scan']), axis=1)

prec_mps_df = pd.read_csv(prec_df_path, index_col=0)
prec_mps_df['scan_id'] = prec_mps_df.index

distances_mps = distances_mps.merge(prec_mps_df[['scan_id', 'Xi_mip', 'Xi_hit']], on='scan_id', how='left')
distances_mps['plot_group'] = distances_mps.apply(sort_to_groups, axis=1)
distances_mps = distances_mps[~distances_mps.plot_group.isnull()]

tt_df = distances_mps[distances_mps.isTT == True]


fig, ax = plt.subplots()
sns.boxplot(x='plot_group', y='dist', data=distances_mps, whis=[5, 95], order=[-4, -3, -2, -1, 0, 'decoy'])

bottom, top = ax.get_ylim()
for i, group in enumerate([-4, -3, -2, -1, 0, 'decoy']):
    ax.text(i, top, '\nn = {}'.format(len(distances_mps[distances_mps.plot_group == group])), horizontalalignment='center')

distances_mps.to_csv(outdir + 'psm_distances_by_massshift.csv')

ax.set_xlabel('mass correction (Da)')
ax.set_ylabel('distance (A)')
ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

sns.despine()
plt.tight_layout()
plt.savefig(outdir + 'distance_by_allshift.svg')
plt.savefig(outdir + 'distance_by_allshift.png')
