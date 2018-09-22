import os
import plot_fnc
import pandas as pd

directory = ''
hsa_dir = os.path.join(directory, 'adam_HSA_5frags_SDA_Fusion')
pc_dir = os.path.join(directory, 'lars_PC_4frag_BS3_Lumos')
ribo_dir = os.path.join(directory, 'marta_ribo_HCD_BS3_Lumos')
chaet_dir = os.path.join(directory, 'chaetomium')
outpath = os.path.join(directory, 'eval_manuscript')

frag_dict = {'HCD': ['E151023_05', 'E151023_06', 'E151023_07', 'B160803_02', 'B160803_07', 'B160803_12',
                         'B160805_06', 'B160805_07', 'B160805_08', 'B160805_11', 'B160805_12', 'B160805_13',
                         'B160805_16', 'B160805_17', 'B160805_18', 'B160805_21', 'B160805_22', 'B160805_23'],
                 'EThcD': ['E151023_14', 'E151023_15', 'E151023_16', 'B160803_03', 'B160803_08', 'B160803_13'],
                 'ETciD': ['E151023_11', 'E151023_12', 'E151023_13', 'B160803_04', 'B160803_09', 'B160803_14'],
                 'CID': ['E151023_02', 'E151023_03', 'E151023_04', 'B160803_05', 'B160803_10', 'B160803_15'],
                 'ETD': ['E151023_08', 'E151023_09', 'E151023_10']}

hsa_counts = plot_fnc.xifdr_counts(xiFDR_folder=os.path.join(hsa_dir, 'Xi', 'xi_xiFDR'),
                                   out_path=os.path.join(outpath, 'hsa'))
pc_counts = plot_fnc.xifdr_counts(xiFDR_folder=os.path.join(pc_dir, 'Xi', 'xi_xiFDR'),
                                  out_path=os.path.join(outpath, 'pc'))
chaet_counts_1 = plot_fnc.xifdr_counts(xiFDR_folder=os.path.join(chaet_dir, 'fr3_6', 'Xi', 'xi_xiFDR'),
                                       out_path=os.path.join(outpath, 'chaet'), prefix='3-6_')
chaet_counts_2 = plot_fnc.xifdr_counts(xiFDR_folder=os.path.join(chaet_dir, 'fr30', 'Xi', 'xi_xiFDR'),
                                       out_path=os.path.join(outpath, 'chaet'), prefix='30_')
chaet_counts_2 = chaet_counts_2.rename(index={'mscon_PF_20': 'mscon_PF_20_100_0',
                                              'rel_1_mscon_PF_20': 'rel_1_mscon_PF_20_100_0',
                                              'rel_2_mscon_PF_20': 'rel_2_mscon_PF_20_100_0',
                                              'rel_3_mscon_PF_20': 'rel_3_mscon_PF_20_100_0',
                                              'rel_4_mscon_PF_20': 'rel_4_mscon_PF_20_100_0'})

chaet_all = pd.concat([chaet_counts_1, chaet_counts_2], axis=1)
chaet_all.to_csv(os.path.join(outpath, 'chaet', 'test_psm_merged.csv'))

plot_fnc.plt_fold_change_datasets(file_ls=[['E151023_05', 'E151023_06', 'E151023_07'],
                                           ['B160803_02', 'B160803_07', 'B160803_12'],
                                           chaet_all.columns.tolist()],
                                  counts_table=pd.concat([hsa_counts, pc_counts, chaet_all], axis=1),
                                  out_path=os.path.join(outpath, 'plots'), name='fig1A',
                                  setting_selection=['MQ_Top20_Win100_deiso0', 'OpenMS_PP0_FF51_MC0_SF3'],
                                  labels=['MaxQuant', 'OpenMS'])

df = plot_fnc.prec_info_table(setting_list=[('MQ_Top20_Win100_deiso0', 'MQ'), ('OpenMS_PP0_FF51_MC0_SF3', 'OMS'),
                                            ('rel_4_mscon_PF_20_100_0', 'Xi')],
                              xidir=os.path.join(pc_dir, 'Xi'), ids=frag_dict['HCD'],
                              outdir=os.path.join(outpath, 'pc'),
                              rawmgf_dir='', withins=True)

plot_fnc.correction_matches(prec_df=df, setting=['MQ'], outdir=os.path.join(outpath, 'plots'))

plot_fnc.overlap_scans2(xifdr_dir=os.path.join(pc_dir, 'Xi', 'xi_xiFDR'), file_ls=frag_dict['HCD'],
                        setting_ls=['other_NoPrePro', 'MQ_Top20_Win100_deiso0'],
                        outpath=os.path.join(outpath, 'plots'), withins=True)

plot_fnc.mip_dependency(prec_mz_table=df, setting='Xi', outdir=os.path.join(outpath, 'plots'))

plot_fnc.plt_fold_change(counts_table=pc_counts,
                         setting_selection=['other_NoPrePro', 'MQ_Top20_Win100_deiso0', '', 'rel_1_mscon_PF_20_100_0', 'rel_2_mscon_PF_20_100_0',
                                            'rel_3_mscon_PF_20_100_0', 'rel_4_mscon_PF_20_100_0', 'rel_5_mscon_PF_20_100_0'],
                         labels=['unproc.', 'MaxQuant-Xi', '', '<= -1', '<= -2', '<= -3', '<= -4', '<= -5'],
                         frag_dict=frag_dict, name='fig2A', group_frag=False,
                         out_path=os.path.join(outpath, 'plots'))

plot_fnc.plt_fold_change(counts_table=hsa_counts,
                         setting_selection=['other_NoPrePro', 'MQ_Top20_Win100_deiso0', '', 'rel_1_mscon_PF_20_100_0',
                                            'rel_2_mscon_PF_20_100_0',
                                            'rel_3_mscon_PF_20_100_0', 'rel_4_mscon_PF_20_100_0'],
                         labels=['unprocessed', 'MaxQuant', '', '-1 Da', '-2 Da', '-3 Da', '-4 Da'],
                         frag_dict=frag_dict, name='fig2A_hsa', group_frag=False,
                         out_path=os.path.join(outpath, 'plots'))

plot_fnc.plt_fold_change(counts_table=chaet_all,
                         setting_selection=['other_NoPrePro', 'MQ_Top20_Win100_deiso0', '',
                                            'rel_1_mscon_PF_20_100_0', 'rel_2_mscon_PF_20_100_0',
                                            'rel_3_mscon_PF_20_100_0', 'rel_4_mscon_PF_20_100_0'],
                         labels=['unprocessed', 'MaxQuant', '', '-1 Da', '-2 Da', '-3 Da', '-4 Da'],
                         frag_dict={'HCD': chaet_all.columns.tolist()}, name='fig2A_chaet', group_frag=False,
                         out_path=os.path.join(outpath, 'plots'))

plot_fnc.overlap_links(setting_ls=['4Da', 'OMS', 'MQ'],
                       outpath=os.path.join(outpath, 'plots'), fig_name='links_overlap',
                       xifdr_dir=os.path.join(pc_dir, 'links_fdr'), withins=True)

plot_fnc.crystal_strc_comp(link_dir='',
                           settings=['4Da'], out_dir=os.path.join(outpath, 'plots'),
                           diff_color='4Da.csv')
