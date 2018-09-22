import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib_venn as venn


def mip_kojak(row):
    if np.isnan(row['mass_raw']) or np.isnan(row['calc_neutral_pep_mass']):
        return None
    else:
        return int(np.round(row['calc_neutral_pep_mass'] - (row['npp_mz'] * row['kojak_z'] - row['kojak_z'] * 1.00727646658)))


def create_columns_for_xifdr(row):
    peptide1 = '.'.join(row['peptide'].split('+')[0].split('.')[:-1])[2:]
    peptide2 = '.'.join(row['peptide'].split('+')[1].split('.')[:-1])[2:]

    peptide_pos1 = row['peptide'].split('+')[0].split('(')[1][:-1]
    peptide_pos2 = row['peptide'].split('+')[1].split('(')[1][:-1]

    z = row['spectrum'].split('.')[-1]

    protein1 = row['protein'].split('+')[0].split('|')[1]
    protein2 = row['protein'].split('+')[1].split('|')[1]

    is_decoy1 = True if 'DECOY_' in row['protein'].split('+')[0] else False
    is_decoy2 = True if 'DECOY_' in row['protein'].split('+')[1] else False

    kojak_score = row['kojak_score'].split('(')[0]
    scan = row['spectrum'].split('.')[-2]
    run = row['spectrum'].split('.')[0]
    return pd.Series({'pep1': peptide1, 'pep2': peptide2, 'pep_pos1': peptide_pos1, 'pep_pos2': peptide_pos2,
                      'charge': z, 'protein1': protein1, 'protein2': protein2,
                      'is_decoy1': is_decoy1, 'is_decoy2': is_decoy2, 'probability': row['probability'],
                      'score': kojak_score, 'scan': scan, 'run': run})


base_dir = ''
kojak_file_dir = base_dir + ''
prec_df_path = ''

xi_precursor_df = pd.read_csv(prec_df_path, index_col=0)

file_ids = ['B160803_02', 'B160803_07', 'B160803_12']

xi_precursor_df['file'] = xi_precursor_df.apply(lambda x: '_'.join(x.name.split('_')[:2]), axis=1)
xi_precursor_df = xi_precursor_df[xi_precursor_df['file'].isin(file_ids)]
xi_precursor_df['scan'] = xi_precursor_df.index


res = []
for kojak_file in ['\B160803_02_interact.pep.csv',
                   '\B160803_07_interact.pep.csv',
                   '\B160803_12_interact.pep.csv']:

    df = pd.read_csv(kojak_file_dir + kojak_file)
    df['file'] = kojak_file[1:11]
    for_fdr = df[df.xlink_type == 'xl'].apply(create_columns_for_xifdr, axis=1)
    for_fdr.to_csv(kojak_file_dir + kojak_file.replace('.csv', '_for_fdr.csv'))
    res.append(df)

kojak_df = pd.concat(res)
kojak_df = kojak_df[kojak_df.xlink_type != 'na']

kojak_df['passed'] = kojak_df.probability.apply(lambda x: True if x >= 0.95 else False)
psm_num = []
for single_file in kojak_df.file.unique():
    a = len(kojak_df[(kojak_df.file == single_file) & (kojak_df.passed)])
    print('%s has %s PSMs passing the threshold.' % (single_file, a))
    psm_num.append(a)

kojak_df['scan_num'] = kojak_df.spectrum.apply(lambda x: (x.split('.')[-2]))
kojak_df['kojak_z'] = kojak_df.spectrum.apply(lambda x: int(x.split('.')[-1]))

kojak_df['scan'] = kojak_df.apply(lambda x: x['file'] + '_' + x['scan_num'], axis=1)

merged = xi_precursor_df.merge(kojak_df, on='scan', how='outer')
merged['kojak_mip'] = merged.apply(mip_kojak, axis=1)

# do FDR before manually

fdr_passed = []
for x in ['\B160803_12_PSM_xiFDR1.0.14.34.csv', '\B160803_02_PSM_xiFDR1.0.14.34.csv',
          '\B160803_07_PSM_xiFDR1.0.14.34.csv']:
    fdr_df = pd.read_csv(kojak_file_dir + x)
    fdr_df = fdr_df[fdr_df.isTT == True]
    fdr_df = fdr_df[fdr_df.fdrGroup.str.contains('Within')]
    fdr_df['scan'] = x[1:11] + '_' + fdr_df['scan'].astype(str)
    fdr_passed += fdr_df['scan'].tolist()

print sum(merged['passed'])
merged['passed'] = merged['scan'].apply(lambda x: True if x in fdr_passed else False)

merged.to_csv(kojak_file_dir + '\precursor_df.csv')
