import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl
from scipy import stats
from matplotlib.ticker import MaxNLocator
import matplotlib_venn as venn
import distance_reader
import zipfile


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

colors_ = brewer2mpl.get_map('Set1', 'Qualitative', 5).hex_colors
frag = ['CID', 'HCD', 'ETciD', 'EThcD', 'ETD']
frag_colors = {frag[i]: colors_[i] for i in range(len(frag))}
colors_ = brewer2mpl.get_map('Accent', 'Qualitative', 5).hex_colors


def write_XiFDR_res_dict(xifdr_folder, restype='PSM'):
    settings = os.listdir(xifdr_folder)

    result_dict = {}
    for single_setting in settings:
        all_xifdr_files = os.listdir(os.path.join(xifdr_folder, single_setting))
        if restype == 'PSM':
            FDR_files = [x for x in all_xifdr_files if '_PSM_' in x]
            FDR_files = [x for x in FDR_files if '_Linear_' not in x]
            column_targets = ['scan', 'Protein1', 'Protein2', 'PepSeq1', 'PepSeq2', 'isTT', 'isTD', 'isDD',
                              'Score', 'Charge', 'LinkPos1', 'LinkPos2', 'fdrGroup']
        elif restype == 'Links':
            FDR_files = [x for x in all_xifdr_files if 'Links' in x]
            column_targets = ['Protein1', 'Protein2', 'fromSite', 'ToSite', 'isTT', 'isTD', 'isDD', 'Score', 'fdrGroup']

        elif restype == 'Pep':
            FDR_files = [x for x in all_xifdr_files if 'PeptidePairs' in x]
            column_targets = ['scan', 'Protein1', 'Protein2', 'PepSeq1', 'PepSeq2', 'isTT', 'isTD', 'isDD',
                              'Score', 'LinkPos1', 'LinkPos2'] # add charge, link pos

        for single_file in FDR_files:
            result_file = pd.read_csv(os.path.join(xifdr_folder, single_setting, single_file),
                                      index_col=False, header=0)
            result_dict[single_file.split('.csv')[0]] = result_file.loc[:, column_targets]

    return result_dict


def write_counts_to_dict(XiFDR_res_dict, filename, counts_dict):
    counts_dict[filename] = {}

    # go through result types set for fdr (from config) and write it to table
    table_sub = XiFDR_res_dict[filename]

    counts_dict[filename]['PSM'] = table_sub[table_sub['isTT'] == True].shape[0]
    counts_dict[filename]['PSM_decoy'] = table_sub[table_sub['isTT'] == False].shape[0]
    counts_dict[filename]['self'] = table_sub[(table_sub['fdrGroup'].str.contains('Within')) & (table_sub['isTT'] == True)].shape[0]
    counts_dict[filename]['between'] = table_sub[(table_sub['fdrGroup'].str.contains('Between')) & (table_sub['isTT'] == True)].shape[0]

    return counts_dict


# return csvs with counts of PSM, links or results (depending on settings in config) and decoys
def write_count_csv(XiFDR_res_dict, out_path, prefix=''):
    # return files done in XiFDR
    all_files_all_settings = [x for x in XiFDR_res_dict]

    counts_dict = {}
    # go through each file
    for i in range(len(all_files_all_settings)):
        # returns dict with counts depending on file and fdr
        counts_dict = write_counts_to_dict(XiFDR_res_dict=XiFDR_res_dict, filename=all_files_all_settings[i],
                                           counts_dict=counts_dict)

    # create dataframe
    df_counts = pd.DataFrame.from_dict(counts_dict, orient='index')
    df_counts['file'] = df_counts.apply(lambda x: '_'.join(x.name.split('_')[-11:-9]), axis=1)
    df_counts['setting'] = df_counts.apply(lambda x: '_'.join(x.name.split('_')[:-11]), axis=1)
    df_counts.index = df_counts['setting']
    df_counts.to_csv(os.path.join(out_path, prefix + 'info_psm_counts.csv'))
    df_2 = pd.DataFrame.from_dict({y: df_counts['self'][df_counts['file'] == y] for y in df_counts['file'].unique()})
    df_2.to_csv(os.path.join(out_path, prefix + 'psm_counts.csv'))
    df_betw = pd.DataFrame.from_dict({y: df_counts['between'][df_counts['file'] == y] for y in df_counts['file'].unique()})
    df_betw.to_csv(os.path.join(out_path, prefix + 'psm_counts_betw.csv'))

    return df_2


def xifdr_counts(xiFDR_folder, out_path, prefix='', restype='PSM'):
    # divide selections somehow differently
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # read in results into big dict XiFDR_res_dict[file][fdr][restype]
    xiFDR_res_dict = write_XiFDR_res_dict(xiFDR_folder, restype=restype)

    # return number of hits for each file for each fdr
    df = write_count_csv(XiFDR_res_dict=xiFDR_res_dict, out_path=out_path, prefix=prefix)

    return df


# bar plot average fold change
# in: list of ids to be plotted, list of settings, setting to calculate fold change to
def plt_fold_change(frag_dict, counts_table, out_path, name, setting_selection, labels=None, group_frag=False,
                    val_csv=False, scale_to='other_NoPrePro', ytop=None):
    def scale_values(series):
        x = series[scale_to]
        return series.apply(lambda y: float(y) / x)

    counts_table = counts_table.apply(scale_values)

    values_dict = {}
    for frag_method in frag_dict:
        score_table_sub = counts_table[[x for x in frag_dict[frag_method] if x in counts_table.columns]]
        score_table_sub = score_table_sub.loc[setting_selection]
        mean = score_table_sub.apply(np.mean, axis=1)
        sem = score_table_sub.apply(stats.sem, axis=1)
        values_dict[frag_method] = [mean, sem]

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    x_all = np.arange(len(setting_selection))
    counter = 0
    frag_method_list = [x for x in ['CID', 'HCD', 'ETD', 'ETciD', 'EThcD'] if x in frag_dict]
    if val_csv:
        df = pd.DataFrame({x: values_dict[x][0] for x in frag_method_list})
        df = pd.concat([df, pd.DataFrame({x + '_sem': values_dict[x][1] for x in frag_method_list})], axis=1)
        df.to_csv(os.path.join(out_path, name + '.csv'))
    if len(frag_method_list) == 0:
        frag_method_list = [x for x in frag_dict]
    fig, ax = plt.subplots()
    width = len(setting_selection)/15.
    for frag_method in frag_method_list:
        if group_frag:
            ax.bar(x_all + width * counter, values_dict[frag_method][0][:len(setting_selection)], width,
                   yerr=values_dict[frag_method][1][:len(setting_selection)],
                   color=frag_colors[frag_method], label=frag_method,
                   error_kw=dict(ecolor='black', alpha=0.7, lw=1, capthick=1))
        else:
            fig, ax = plt.subplots()
            width = len(x_all)/15.
            ax.bar(x_all, values_dict[frag_method][0][:len(setting_selection)], width,
                   yerr=values_dict[frag_method][1][:len(setting_selection)],
                   color=frag_colors[frag_method],
                   error_kw=dict(ecolor='black', alpha=0.7, lw=1, capthick=1))
        counter += 1
        if (frag_method == 'HCD') | (frag_method == frag_method_list[-1]) | (group_frag == False):
            ax.set_xticks(x_all)

            if labels is None:
                ax.set_xticklabels(values_dict[frag_method][0].index, rotation=20, ha='right')
            else:
                ax.set_xticklabels(labels)
            ax.yaxis.grid(True)
            ax.set_ylabel('fold change over unprocessed')
            ax.set_xlabel('preprocessing')
            plt.axhline(1, alpha=0.7, color='black', linewidth=2, linestyle='dashed')  #
            if group_frag:
                ax.legend()

            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            if ytop is not None: ax.set_ylim(top=ytop)
            sns.despine()
            plt.tight_layout()
            plt.savefig(os.path.join(out_path, frag_method + '_' + name + '.svg'), format='svg', transparent=True)
            plt.close()
            counter = 0
            fig, ax = plt.subplots()
            plt_dict = {}


def plt_fold_change_datasets(file_ls, counts_table, out_path, name, setting_selection, labels=None,
                             scale_to='other_NoPrePro', ytop=None):
    def scale_values(series):
        x = series[scale_to]
        return series.apply(lambda y: float(y) / x)

    counts_table = counts_table.apply(scale_values)

    values = []
    for files in file_ls:
        score_table_sub = counts_table[[x for x in files if x in counts_table.columns]]
        score_table_sub = score_table_sub.loc[setting_selection]
        mean = score_table_sub.apply(np.mean, axis=1)
        sem = score_table_sub.apply(stats.sem, axis=1)
        values.append([mean, sem])

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    fig, ax = plt.subplots()
    x_all = np.arange(len(file_ls))

    width = len(x_all) / 15.
    heights = np.array([values[i][0][:len(setting_selection)].tolist() for i in range(len(file_ls))])

    err = np.array([values[i][1][:len(setting_selection)].tolist() for i in range(len(file_ls))])

    for i in range(len(labels)):
        ax.bar(x_all + i*width, heights[:,i], width, yerr=err[:,i], label=labels[i],
               error_kw=dict(ecolor='black', alpha=0.7, lw=1, capthick=1))

    ax.set_xticks(x_all + width/2.*(len(labels)-1))

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.grid(True)
    ax.set_xlabel('preprocessing software')
    ax.set_ylabel('fold change over unprocessed')
    plt.axhline(1, alpha=0.7, color='black', linewidth=2, linestyle='dashed')  #

    if ytop is not None: ax.set_ylim(top=ytop)
    sns.despine()
    plt.legend(fontsize=22)
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, name + '.svg'), format='svg', transparent=True)
    plt.close()
    return ax


def return_prec_mz_xi(xi_table_dir, ids_sub=None, withins=False, mz_only=False):
    res_table = pd.DataFrame()
    for i in os.listdir(xi_table_dir):
        if (ids_sub is not None) & ('_'.join(i.split('.')[0].split('_')[-2:]) not in ids_sub):
            continue
        xi_table = pd.read_csv(os.path.join(xi_table_dir, i), index_col=None, dtype=object)
        xi_table = xi_table.sort_values('match score', ascending=False)
        xi_table.drop_duplicates('Scan', inplace=True)
        if withins:
            xi_table = xi_table[xi_table['Fasta1'] == xi_table['Fasta2']]
        xi_table_sub = xi_table[['Scan', 'PrecoursorCharge', 'PrecurserMZ', 'CalcMZ']]
        xi_table_sub.columns = ['scan', 'charge', 'mz', 'calcmz']
        xi_table_sub['scan'] = '_'.join(i.split('.')[0].split('_')[-2:]) + '_' + xi_table_sub['scan']
        res_table = pd.concat([res_table, xi_table_sub])

    res_table.index = res_table['scan']
    res_table['mz'] = res_table['mz'].astype(float)
    if not mz_only:
        res_table['charge'] = res_table['charge'].astype(int)
        res_table['calcmz'] = res_table['calcmz'].astype(float)
    else:
        del res_table['charge']
        del res_table['calcmz']

    del res_table['scan']
    return res_table


def result_scans(xifdr_folder, id_prefix=False, return_calcmz=False, withins=False):
    all_files = os.listdir(xifdr_folder)
    psm_res = [x for x in all_files if '_PSM' in x and 'Linear' not in x]

    result = []
    for i in psm_res:
        if not id_prefix:
            identifier = i.split('_')[-11:-9]
        else:
            identifier = i.split('_')[:2]
        res_table = pd.DataFrame.from_csv(os.path.join(xifdr_folder, i))
        res_table = res_table[res_table.isTT == True]
        if withins:
            res_table = res_table[res_table.fdrGroup.str.contains('Within')]
        scans = res_table.scan.tolist()
        if return_calcmz:
            calcmz = (res_table['match mass'] / res_table['match charge'] + res_table['match charge']).tolist()
            result += [('_'.join(identifier + [str(x)]), calcmz[i_scan]) for i_scan, x in enumerate(scans)]
        else:
            result += ['_'.join(identifier + [str(x)]) for x in scans]

    return result


def single_mass_hist(csv, column, outfile, cutoff=True):
    df = pd.DataFrame.from_csv(csv)
    if cutoff:
        df = df[(df.decoy == 0) & (df['match score'] > 8)]
    fig, ax = plt.subplots()
    sns.distplot(df[column])
    ax.axvline(np.median(df[column]))
    plt.savefig(outfile)


def read_intensities(mgf_file, identifier):
    spectra_file = open(mgf_file)
    spectra_dict = {}
    mz = {}
    charge = {}
    end_ident, prec_ident, charge_ident, scan_ident = 'END IONS', 'PEPMASS', 'CHARGE', 'TITLE'
    for line in spectra_file:
        if prec_ident in line:
            try:
                spectra_dict[identifier + '_' + index] = float(line.split('=')[1].split(' ')[1])
                mz[identifier + '_' + index] = float(line.split('=')[1].split(' ')[0])
            except:
                # print(identifier + '_' + index)
                spectra_dict[identifier + '_' + index] = 0
                mz[identifier + '_' + index] = float(line.split('=')[1].split(' ')[0])
        elif scan_ident in line:
            try:
                index = line.split('scan=')[1].split('_')[0]
                if '"' in index:
                    index = index[:index.find('"')]
            except IndexError:
                index = line.split('.')[-2]
        elif charge_ident in line:
            try:
                charge[identifier + '_' + index] = int(line.split('=')[1][0])
            except UnboundLocalError:
                pass
    return spectra_dict, mz, charge


def read_apl(apl_zip, identifier):
    apl_zipped = zipfile.ZipFile(os.path.join(apl_zip))
    apls = apl_zipped.namelist()
    mz = {}
    charge = {}
    for x in apls:
        print x
        apl = apl_zipped.open(x)
        end_ident, prec_ident, charge_ident, scan_ident = 'peaklist end', 'mz', 'charge', 'header'
        for line in apl:
            if prec_ident in line:
                mz_ = float(line.split('=')[1].split(' ')[0])
            elif scan_ident in line:
                index = line.split('Index: ')[1].split(' ')[0]
                mz[identifier + '_' + index] = mz_
                charge[identifier + '_' + index] = charge_
            elif charge_ident in line:
                if '.peak' in x:
                    charge_ = 0
                else:
                    charge_ = int(line.split('=')[1])

    return mz, charge


def prec_info_table(setting_list, xidir, ids, outdir, rawmgf_dir, save=True, withins=False):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    precdf = return_prec_mz_xi(xi_table_dir=os.path.join(xidir, 'xi', 'other_NoPrePro'), ids_sub=ids)
    precdf.columns = ['npp_' + x for x in precdf.columns]
    hits = result_scans(os.path.join(xidir, 'xi_xiFDR', 'other_NoPrePro'))

    precdf['npp_hit'] = [True if x in hits else False for x in precdf.index]

    for setting in setting_list:
        df = return_prec_mz_xi(xi_table_dir=os.path.join(xidir, 'xi', setting[0]), ids_sub=ids, withins=withins)
        df.columns = [setting[1] + '_' + x for x in df.columns]
        hits = result_scans(os.path.join(xidir, 'xi_xiFDR', setting[0]))
        df[setting[1] + '_hit'] = [True if x in hits else False for x in df.index]
        precdf = pd.concat([precdf, df], axis=1, join='outer')

    intensities, mz, charge = {}, {}, {}
    for single_id in ids:
        mgf_file = [x for x in os.listdir(rawmgf_dir) if single_id in x and '.mgf' in x]
        if len(mgf_file) == 0: continue
        intensities_, mz_, charge_ = read_intensities(os.path.join(rawmgf_dir, mgf_file[0]), single_id)
        charge.update(charge_)
        mz.update(mz_)
        intensities.update(intensities_)

    precdf = pd.concat([precdf, pd.Series(intensities, name='intensity')], axis=1, join='inner')
    precdf['npp_mz'] = pd.Series(mz, name='npp_mz')
    precdf['npp_charge'] = pd.Series(charge, name='npp_charge')
    precdf['mass_raw'] = precdf.apply(lambda x: x['npp_charge'] * x['npp_mz'] - x['npp_charge'] * 1.00727646658, axis=1)

    for setting in setting_list:
        precdf[setting[1] + '_mip'] = precdf.apply(lambda x: round(x['npp_charge'] * (x[setting[1] + '_mz'] - x['npp_mz'])),
                                             axis=1)
    if save:
        precdf.to_csv(os.path.join(outdir, 'prec_infos.csv'))
    return precdf


def correction_matches(prec_df, setting, outdir, ref_setting='npp', exclude_npp=False):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if len(setting)==1:
        if exclude_npp:
            hit_tables = [prec_df[(prec_df[setting[0] + '_hit'] == True) & (prec_df[ref_setting+'_hit'] != True)& (prec_df['npp_hit'] != True)]]
        else:
            hit_tables = [prec_df[(prec_df[setting[0] + '_hit'] == True) & (prec_df[ref_setting+'_hit'] != True)]]
    elif len(setting)==2:
        hit_tables = [prec_df[(prec_df[setting[0] + '_hit'] == True) & (prec_df[ref_setting+'_hit'] != True) & (prec_df[setting[1]+'_hit'] != True)],
                      prec_df[(prec_df[setting[0] + '_hit'] != True) & (prec_df[ref_setting + '_hit'] != True) & (
                      prec_df[setting[1] + '_hit'] == True)],
                      prec_df[(prec_df[setting[0] + '_hit'] == True) & (prec_df[ref_setting + '_hit'] != True) & (
                      prec_df[setting[1] + '_hit'] == True)]]

    fig, ax = plt.subplots()
    col, label = ['#008000', '#0000ff', '#ff0000'], setting + ['intersection']

    for i, matches in enumerate(hit_tables):
        if label[i] == ref_setting:
            charge = (matches[ref_setting + '_charge'] - matches[label[i] + '_charge']) != 0
            mip_any = (False) * len(matches)
            mip = [[False]*len(matches) if x!=0 else [True]*len(matches) for x in [-5, -4, -3, -2, -1, 0, 1]]
        elif label[i] == 'intersection':
            charge_1 = (matches[ref_setting + '_charge'] - matches[setting[0] + '_charge']) != 0
            charge_2 = (matches[ref_setting + '_charge'] - matches[setting[1] + '_charge']) != 0
            if (charge_1 == charge_2).all():
                charge = charge_1
            mip_any_1 = matches[setting[0] + '_mip'] != 0
            mip_any_2 = matches[setting[1] + '_mip'] != 0
            if (mip_any_1 == mip_any_2).all():
                mip_any = mip_any_2
            else:
                mip_any = mip_any_2
            mip_1 = [matches[setting[0] + '_mip'] == x for x in [-5, -4, -3, -2, -1, 0, 1]]
            mip_2 = [matches[setting[1] + '_mip'] == x for x in [-5, -4, -3, -2, -1, 0, 1]]
            if all([(mip_1[j] == mip_2[j]).all() for j in range(len(mip_1))]):
                mip = mip_2
            else:
                # ----- change this ------ #
                print 'mip', len(matches), [len(matches[mip_1[j] != mip_2[j]]) for j in range(len(mip_1))]
                mip = mip_2
        else:
            matches['MQ_charge'] = matches.apply(lambda x: x['MQ_charge'] if x['MQ_charge'] > 0 else x['npp_charge'], axis=1)
            if not ref_setting == 'npp':
                matches[label[i] + '_mip'] = matches[label[i] + '_mip'] - matches[ref_setting + '_mip']
            charge = (matches[ref_setting + '_charge'] - matches[label[i] + '_charge']) != 0
            mip_any = matches[label[i] + '_mip'] != 0
            mip = [matches[label[i] + '_mip'] == x for x in [-5, -4, -3, -2, -1, 0, 1]]

        ndiff_charge = len(matches[charge & ~mip_any])
        # diff mip
        ndiff_mip = [len(matches[x & ~charge]) for x in mip]
        ndiff_mipchar = len(matches[mip_any & charge])

        heights = ndiff_mip + [ndiff_charge, ndiff_mipchar] + [len(matches) - sum(ndiff_mip + [ndiff_charge, ndiff_mipchar])]
        print sum(heights)
        x = np.arange(len(heights))
        width = len(x)/15.

        if i==0:
            ax.bar(x, heights, width, color=col[i], lw=1.2, label=label[i], alpha=0.4
                   )
            bottom_heights = heights
        else:
            ax.bar(x, heights, width, color=col[i], lw=1.2, label=label[i], bottom=bottom_heights, alpha=0.4
                   )
            bottom_heights = [sum(y) for y in zip(heights, bottom_heights)]

    ax.set_xticks(x)
    ax.set_xticklabels(['-5 Da', '-4 Da', '-3 Da', '-2 Da', '-1 Da', '0 Da', '1 Da'] + ['charge', 'MP\n+charge', 'other'])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_ylabel('scans')
    ax.set_xlabel('precursor correction')
    plt.legend()
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '_'.join(setting+[ref_setting])+'_corrections.svg'))


def plot_properties(intensities, mz, masses, comp, labels, out, charges):
    if not os.path.exists(out):
        os.makedirs(out)

    plt_int, plt_mz, plt_mass, plt_charge = [], [], [], []
    if rel_int is not None: plt_relint = []
    for element in comp:
        plt_int.append(np.array([(intensities[x]) for x in element])) #
        plt_mz.append(np.array([mz[x] for x in element]))
        plt_mass.append(np.array([masses[x] for x in element]))
        plt_charge.append(np.array([charges[x]/np.sqrt(masses[x]) for x in element]))
        if rel_int is not None: plt_relint.append(np.array([rel_int[x]/np.sqrt(masses[x]) for x in element]))

    sizes = ['\nn = {}'.format(len(comp[i])) for i in range(len(comp))]

    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    sns.boxplot(data=plt_int, palette=sns.cubehelix_palette(len(plt_mass), dark=0.4), ax=ax)
    ax.set_xticklabels(labels)
    ax.set_ylabel('intensity')
    ax.set_xlabel('precursor mass correction')
    bottom, top = ax.get_ylim()
    ax.set_ylim(top=top+0.8)
    for i, sample_size in enumerate(sizes):
        ax.text(i, top+0.4, sample_size, horizontalalignment='center')
    print 'INTENSITY'
    for i, x in enumerate(plt_int):
        print labels[i]
        print np.median(x)
        if i != 0:
            y_ = top
            x1, x2 = i-0.01, i-1+0.01
            plt.plot([x1, x1, x2, x2], [y_, y_ + 0.2, y_ + 0.2, y_], lw=1.5, color='black')
            p = (stats.ttest_ind(x, plt_int[i-1])[1])/2
            print p
            if p < 0.0001:
                plt.text((x1 + x2) * .5, y_ + 0.2, "****", ha='center', va='bottom')
            elif (p < 0.001):
                plt.text((x1 + x2) * .5, y_ + 0.2, "***", ha='center', va='bottom')
            elif (p < 0.01):
                plt.text((x1 + x2) * .5, y_ + 0.2, "**", ha='center', va='bottom')
            elif (p < 0.05):
                plt.text((x1 + x2) * .5, y_ + 0.2, "*", ha='center', va='bottom')
            else:
                plt.text((x1 + x2) * .5, y_ + 0.2, "ns", ha='center', va='bottom')
    fig.canvas.draw()
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(out, 'intensities.svg'))

    fig, ax = plt.subplots()
    sns.boxplot(data=plt_mass, palette=sns.cubehelix_palette(len(plt_mass), dark=0.4), whis=[5, 95])
    ax.set_xticklabels(labels)
    ax.set_ylabel('precursor mass (Da)')
    ax.set_xlabel('precursor mass correction')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    bottom, top = ax.get_ylim()
    ax.set_ylim(top=top+200)
    for i, sample_size in enumerate(sizes):
        ax.text(i, top+200, sample_size, horizontalalignment='center')
    for i, x in enumerate(plt_mass):
        print labels[i]
        print np.median(x)
        if i != 0:
            y_ = top
            x1, x2 = i-0.01, i-1+0.01
            plt.plot([x1, x1, x2, x2], [y_, y_ + 100, y_ + 100, y_], lw=1.5, color='black')
            p = stats.ttest_ind(plt_int[i-1], x)[1]/2
            if p < 0.0001:
                plt.text((x1 + x2) * .5, y_ + 100, "****", ha='center', va='bottom')
            elif (p < 0.001):
                plt.text((x1 + x2) * .5, y_ + 100, "***", ha='center', va='bottom')
    plt.tight_layout()
    sns.despine()
    plt.savefig(os.path.join(out, 'mass.svg'))


def mip_dependency(prec_mz_table, setting, outdir):
    da_0 = prec_mz_table[(prec_mz_table[setting + '_mip'] == 0) & (prec_mz_table[setting + '_hit'] == True)].index
    da_1 = prec_mz_table[(prec_mz_table[setting + '_mip'] == -1) & (prec_mz_table[setting + '_hit'] == True)].index
    da_2 = prec_mz_table[(prec_mz_table[setting + '_mip'] == -2) & (prec_mz_table[setting + '_hit'] == True)].index
    da_3 = prec_mz_table[(prec_mz_table[setting + '_mip'] == -3) & (prec_mz_table[setting + '_hit'] == True)].index
    da_4 = prec_mz_table[(prec_mz_table[setting + '_mip'] == -4) & (prec_mz_table[setting + '_hit'] == True)].index

    plot_properties(comp=[da_4, da_3, da_2, da_1, da_0, ], labels=['-4 Da', '-3 Da', '-2 Da', '-1 Da', '0 Da', ] ,
                    out=outdir, intensities=prec_mz_table.intensity, mz=prec_mz_table.npp_mz, masses=prec_mz_table.mass_raw,
                    charges=prec_mz_table.npp_charge)


def delete_mod_from_seq(table, mod_list):
    for i in mod_list:
        if i == 'Mox':
            table['PepSeq1'] = table['PepSeq1'].str.replace(i, 'M')
            table['PepSeq2'] = table['PepSeq2'].str.replace(i, 'M')
        else:
            table['PepSeq1'] = table['PepSeq1'].str.replace(i, '')
            table['PepSeq2'] = table['PepSeq2'].str.replace(i, '')

    return table


def overlap_scans2(xifdr_dir, setting_ls, outpath, file_ls='', save_plot=True, withins=False):
    join_series = lambda x: '_'.join([str(y) for y in pd.Series.tolist(x)])
    xifdr_out = write_XiFDR_res_dict(xifdr_folder=xifdr_dir)

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    res_list = []
    decoys_table = pd.DataFrame()
    for single_file in xifdr_out:
        name = '_'.join(single_file.split('_')[-11:-9])
        setting = '_'.join(single_file.split('_')[:-11])

        if setting not in setting_ls or name not in file_ls:
            continue
        table = xifdr_out[single_file].copy()
        if withins:
            table = table[table.fdrGroup.str.contains('Within')]
        sub = table[['scan']]
        decoys = table[['isTT', 'isTD', 'isDD']]

        sub = sub.apply(join_series, axis=1)
        sub = [name + '_' + str(x) for x in sub]

        decoys['joined'], decoys['setting'] = sub, setting
        decoys_table = pd.concat([decoys_table, decoys])

    fdrs = {}
    for sel_setting in setting_ls:
        setting_all = decoys_table[decoys_table.setting == sel_setting]
        res_list.append(set(setting_all['joined'][setting_all.isTT == True].tolist()))

        outer_joined = np.setdiff1d(setting_all['joined'], decoys_table[decoys_table.setting != sel_setting]['joined'])
        outer = setting_all[setting_all['joined'].isin(outer_joined)]
        fdrs[sel_setting] = (sum(outer.isTD) - sum(outer.isDD)) / float(sum(outer.isTT)) * 100

    pd.Series(fdrs).to_csv(os.path.join(outpath, 'psm_fdrs.csv'))

    labels = [x.split('_')[1] for x in setting_ls]

    plt.figure()
    if len(res_list) == 2:
        diag = venn.venn2(res_list, tuple(labels))
        venn.venn2_circles(res_list)
    elif len(res_list) == 3:
        diag = venn.venn3(res_list, tuple(labels))
        venn.venn3_circles(res_list)

    if save_plot:
        plt.savefig(os.path.join(outpath, '_'.join(labels) + '_venn_scans.svg'))
    return res_list


def overlap_links(setting_ls, outpath, fig_name, xifdr_dir, file_ls='', create_plot=True, save_plot=True, withins=False):

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    def join_linkpos(x):
        a = '_'.join([str(x['fromSite']), x['Protein1']])
        b = '_'.join([str(x['ToSite']), x['Protein2']])
        res = [a, b]
        res.sort()
        return '_'.join(res)

    xifdr_out = write_XiFDR_res_dict(xifdr_folder=xifdr_dir, restype='Links')

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    all_res = pd.DataFrame()
    for single_file in xifdr_out:
        name = 'HCD'
        setting = single_file.split('_')[0]
        if setting not in setting_ls or (name not in file_ls and file_ls is not ''):
            continue
        table = xifdr_out[single_file].copy()
        if withins:
            table = table[table['fdrGroup'].str.contains('Within')]
        table['setting'] = setting
        table['joined'] = table.apply(lambda x: join_linkpos(x) + '_' + name, axis=1)
        table.to_csv(os.path.join(outpath, 'links_joined_' + setting + '.csv'))
        all_res = pd.concat([all_res, table])

    all_res = all_res.drop_duplicates(['setting', 'joined'])
    all_res.to_csv(os.path.join(outpath, 'links_joined_merged.csv'))
    res_list = []
    fdrs = {}
    for sel_setting in setting_ls:
        setting_all = all_res[all_res.setting == sel_setting]
        res_list.append(set(setting_all['joined'][setting_all.isTT == True].tolist()))

        outer_joined = np.setdiff1d(setting_all['joined'], all_res[all_res.setting != sel_setting]['joined'])
        outer = setting_all[setting_all['joined'].isin(outer_joined)]
        fdrs[sel_setting] = {'TT': sum(outer.isTT), 'TD': sum(outer.isTD), 'DD': sum(outer.isDD),
                             'FDR': (sum(outer.isTD) - sum(outer.isDD)) / float(sum(outer.isTT)) * 100}

    pd.DataFrame(fdrs).to_csv(os.path.join(outpath, 'fdrs.csv'))
    labels = setting_ls # [x.split('_')[1] for x in setting_ls]

    if create_plot:
        plt.figure()
        if len(res_list)==3:
            diag = venn.venn3(res_list, tuple(labels))
            venn.venn3_circles(res_list)
        elif len(res_list)==2:
            diag = venn.venn2(res_list, tuple(labels))
            venn.venn2_circles(res_list)
        venn_colors = sns.cubehelix_palette(3, start=1.6, rot=1, dark=.3, light=.7)

    if save_plot:
        plt.savefig(os.path.join(outpath, fig_name + '.svg'))

    for i, setting in enumerate(res_list):
        other = [x for x in range(len(res_list)) if not x == i]
        outer = np.setdiff1d(np.setdiff1d(list(res_list[i]), list(res_list[other[0]])), list(res_list[other[1]]))
        pd.DataFrame(all_res[all_res['joined'].isin(outer)]).to_csv(os.path.join(outpath, setting_ls[i] + '.csv'))

    return res_list


def plot_distogram(distances, rnd, outpath, max_distance, distances_2):
    fig, ax = plt.subplots()
    sns.violinplot(data=[rnd, distances, distances_2], bw=0.3)
    ax.set_xticklabels(['random', 'Xi-MPS', 'Xi-MPS\nunique'])
    ax.set_ylim([-7, 95])
    sns.despine()
    plt.savefig(outpath + 'box.svg')
    plt.close()
    fig, ax = plt.subplots()
    binwidth = 1
    bins = np.arange(int(min(distances) - 0.5), int(max(rnd))+2.1, 1)

    tmp_hist = np.histogram(rnd, np.arange((min(distances)), (max(rnd))+2.1, 1.5), density=True)
    ax.bar(tmp_hist[1][1:]-0.5, tmp_hist[0], alpha=1, color='#4d4d4d', label='random', width=1.45)

    hist_1 = np.histogram(distances.dropna(), bins, density=True)
    ax.bar(hist_1[1][1:]-0.5, hist_1[0], label='observed', width=binwidth, alpha=0.9, color='#cccccc')

    hist_2 = np.histogram(distances_2.dropna(), bins)

    ax.bar(hist_2[1][1:]-0.5, [x / float(len(distances)) for x in hist_2[0]], label='outer', width=binwidth,
           alpha=0.8, color='#69afa3')

    ax.axvline(x=max_distance, color='black', alpha=0.9)
    if not len(distances) == 0:
        long_dist = round(len(distances[distances > max_distance]) / float(len(distances)) * 100, 2)
        ax.text(32, 0.02, 'Long distance: ' + str(long_dist) + ' %')
    if not len(distances_2) == 0:
        long_dist_2 = round(len(distances_2.dropna()[distances_2.dropna() > max_distance]) / float(len(distances_2.dropna())) * 100, 2)
        ax.text(32, 0.01, 'Long distance outer: ' + str(long_dist_2) + ' %')
    ax.set_xlabel('distance (A)')
    ax.set_yticklabels(['{}'.format(x*100) for x in ax.get_yticks()])
    ax.set_ylabel('frequency (%)')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    sns.despine()
    ax.legend(loc=0, fontsize=22)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def get_link_distances(pdb_distances, links, outpath, name, max_distance=25, diff_color='', psm_level=False):
    if not psm_level:
        from_col, to_col = 'fromSite', 'ToSite'
    else:
        from_col, to_col = 'ProteinLinkPos1', 'ProteinLinkPos2'

    def join_linkpos(x):
        a = '_'.join([str(x[from_col]), x['Protein1']])
        b = '_'.join([str(x[to_col]), x['Protein2']])
        res = [a, b]
        res.sort()
        return '_'.join(res)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if not diff_color == '':
        diff_mark_tt = pd.read_csv(diff_color)['joined']
    else:
        diff_mark_tt = []
    links_tt = links[links.isTT == True]

    links_tt['Protein1'] = links_tt['Protein1'].apply(lambda x: x.split('DECOY:')[-1])
    links_tt['Protein2'] = links_tt['Protein2'].apply(lambda x: x.split('DECOY:')[-1])
    links_tt = links_tt[links_tt.Protein1 == links_tt.Protein2]
    links_site = links_tt[['Protein1', from_col, to_col]]
    links_tt['dist'] = None
    links_tt['joined'] = links_tt.apply(join_linkpos, axis=1)

    def return_distance(x):
        row = pdb_single_prot[pdb_single_prot['residue_1'].isin(x) &
                               pdb_single_prot['residue_2'].isin(x)]
        dist = row.distance
        if len(dist) == 1:
            if row['residue_type_1'].values[0] not in ['LYS', 'SER', 'THR', 'TYR']:
                print row
                raise StandardError
            elif row['residue_type_2'].values[0] not in ['LYS', 'SER', 'THR', 'TYR']:
                print row
                raise StandardError

        if len(dist) == 1:
            dist = dist.values[0]
            return dist
        elif len(dist) > 1:
            print 'more res'
            print row
            raise StandardError

    distances = pd.Series()
    distances2 = pd.Series()
    for protein in links_site.Protein1.unique():
        if protein == 'P01024':
            continue
        pdb_single_prot = pdb_distances[protein].copy()
        pdb_single_prot = pdb_single_prot[pdb_single_prot['residue_1'] != pdb_single_prot['residue_2']]
        links_tt.loc[links_site.Protein1 == protein, 'dist'] = links_site[links_site.Protein1 ==
                                                                                            protein][
            [from_col, to_col]].apply(return_distance,
                                          axis=1)

        distances = pd.concat([distances, links_tt[links_site.Protein1 == protein]['dist']])
        joined_prot = [x[:-4] for x in diff_mark_tt if protein in x]

        distances2 = pd.concat([distances2, links_tt[(links_site.Protein1 == protein) & (links_tt['joined'].isin(joined_prot))]['dist']])
        plot_distogram(distances=links_tt[links_site.Protein1 == protein]['dist'].dropna(),
                       rnd=pdb_distances[protein][pdb_distances[protein].distance != 0].distance,
                       outpath=os.path.join(outpath, name + '_' + protein + '.svg'),
                       max_distance=max_distance,
                       distances_2=links_tt[(links_site.Protein1 == protein) & (links_tt['joined'].isin(joined_prot))]['dist'].dropna())

    links_tt.to_csv(os.path.join(outpath, name + '_fdrdf_dist.csv'))
    links_tt['dist'].to_csv(os.path.join(outpath, name + '_dist.csv'))


    links_prop = {'not_ld': len(distances[distances <= max_distance]),
                  'ld': len(distances[distances > max_distance]),
                  'NA': len(distances[distances.isnull()])}

    distances = distances.dropna()
    rnd = [pdb_distances[x][pdb_distances[x].distance != 0].distance for x in pdb_distances if not x == 'P01024']
    rnd = [item for sublist in rnd for item in sublist]

    plot_distogram(distances=distances, rnd=rnd, outpath=os.path.join(outpath, name + '.svg'), max_distance=max_distance, distances_2=distances2)
    pd.DataFrame(rnd).to_csv(os.path.join(outpath, 'random.csv'), index=None)
    return links_prop


def score_distribution(df, outpath, name, decoys=True):
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    try:
        scores_tt = df[df.isTT == True].Score
        if decoys:
            scores_dec = df[df.isTT == False].Score
    except AttributeError:
        scores_tt = df[df.decoy == '0']['pvalue'].astype(float)
        if decoys:
            scores_dec = df[df.decoy == '1']['pvalue'].astype(float)

    bins = np.arange(float(scores_tt.min()) - 2, float(scores_tt.max()) + 1, 0.25)
    fig, ax = plt.subplots()

    ax.hist(scores_tt, bins, alpha=0.9, normed=True, label='target')
    if decoys:
        ax.hist(scores_dec, bins, alpha=0.9, normed=True, label='decoy')

    ax.set_xlabel('Score')
    ax.set_ylabel('Frequency')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    sns.despine()
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, name))
    plt.close()
    return scores_tt.tolist()


def crystal_strc_comp(link_dir, settings, out_dir, diff_color='', psm=False):
    cl_distance = 30
    pdb_file = {
        'P02768-A': '1ao6.pdb',
        'P00563': '2crk.pdb',
        'P00432': '3J7U_5gkn.pdb',
        'P68082': '2frj.pdb',
        'P02789': '1ovt.pdb',
        'P00004': '3nbs.pdb'
    }
    chains = {
        'P01024': 'B', 'P02789': 'A', 'P02768-A': 'A', 'P00563': 'A', 'P00432': 'A', 'P00004': 'A',
        'P68082': 'X'
    }

    correction = {'P00432': 2, 'P02789': -19, 'P00563': 7,}

    # --------------------- main -----------------------
    distances = {}
    for protein in pdb_file:
        if protein == 'P01024': continue
        distances[protein] = distance_reader.get_distances_from_pdb(filename=pdb_file[protein], uniprot_id=protein,
                                                                    chain_id=chains[protein])
        nterm = distances[protein]['residue_1'].min()
        distances[protein] = distances[protein][((distances[protein]['residue_type_1'].isin(['LYS', 'SER', 'THR', 'TYR'])) | (distances[protein]['residue_1'] == nterm)) &
                                                ((distances[protein]['residue_type_2'].isin(['LYS', 'SER', 'THR', 'TYR']))| (distances[protein]['residue_2'] == nterm))
        ]
        if protein in correction:
            distances[protein]['residue_1'] -= correction[protein]
            distances[protein]['residue_2'] -= correction[protein]

    link_tables = []
    for single_setting in settings:
        if not psm:
            link_tables += [os.path.join(link_dir, single_setting, x) for x in
                            os.listdir(os.path.join(link_dir, single_setting)) if 'Links' in x]
        else:
            link_tables += [os.path.join(link_dir, single_setting, x) for x in
                            os.listdir(os.path.join(link_dir, single_setting)) if ('_PSM_' in x) and not ('Linear' in x)]
    region_names = [os.path.split(os.path.split(x)[0])[1] + '_' + os.path.split(x)[1].split('_')[0] for x in link_tables]

    link_df = {region_names[i]: pd.DataFrame.from_csv(link_tables[i]) for i in range(len(link_tables))}

    for region in region_names:
        score_distribution(df=link_df[region], outpath=os.path.join(out_dir, 'scores'), decoys=False,
                           name=region + '.png')

    overlength = {}
    for region in link_df:
        overlength[region] = get_link_distances(pdb_distances=distances, links=link_df[region],
                                                outpath=os.path.join(out_dir, 'distances'),
                                                name=region, max_distance=cl_distance, diff_color=diff_color,
                                                psm_level=psm)
    overlength_df = pd.DataFrame.from_dict(overlength, orient='index')
    overlength_df['ld %'] = overlength_df.apply(lambda x: x['ld'] * 100.0 / (x['ld'] + x['not_ld']), axis=1)
    overlength_df.to_csv(os.path.join(out_dir, 'distances', 'overlengths.csv'))

