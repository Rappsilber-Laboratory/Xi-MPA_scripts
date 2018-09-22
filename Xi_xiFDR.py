import os
import subprocess
import pandas as pd
from multiprocessing import Pool


def create_search_list(source_path, res_folder, rel=None):
    all_folders = [x for x in os.listdir(source_path) if '.conf' not in x]
    if rel is not None:
        rel_setting = 'rel_' + '_'.join([str(x) for x in rel[1:]]) + '_' + rel[0]
        all_folders.append(rel_setting)
    else: rel_setting = ''
    settings_done = os.listdir(res_folder)

    file_list = []
    for setting_folder in all_folders:
        if setting_folder not in settings_done:
            if setting_folder == rel_setting:
                file_dict = [{'full_path': os.path.join(source_path, rel[0], x),
                              'setting': setting_folder}
                             for x in os.listdir(os.path.join(source_path, rel[0]))]
            else:
                file_dict = [{'full_path': os.path.join(source_path, setting_folder, x),
                              'setting': setting_folder}
                             for x in os.listdir(os.path.join(source_path, setting_folder))]
            file_list += file_dict
        else:
            if setting_folder == rel_setting:
                setting_source = rel[0]
                orig_files = [x.split('.')[0] for x in
                              os.listdir(os.path.join(source_path, setting_source))]
                str_shift = len(rel_setting) - len(setting_source)
                files_done = [x.split('.')[0][str_shift:] for x in os.listdir(os.path.join(res_folder, setting_folder))]
            else:
                setting_source = setting_folder
                orig_files = [x.split('.')[0] for x in
                              os.listdir(os.path.join(source_path, setting_source))]
                files_done = [x.split('.')[0] for x in os.listdir(os.path.join(res_folder, setting_folder))]

            if not len(files_done) == len(orig_files):
                file_type = '.' + os.listdir(os.path.join(source_path, setting_source))[0].split('.')[1]
                file_dict = [{'full_path': os.path.join(source_path, setting_source, x + file_type),
                              'setting': setting_folder}
                             for x in orig_files if x not in files_done]
                file_list += file_dict

    return file_list


def create_search_list_xifdr(source_path, res_folder):
    all_folders = [x for x in os.listdir(source_path) if '.conf' not in x]
    settings_done = os.listdir(res_folder)

    file_list = []
    for setting_folder in all_folders:
        if setting_folder not in settings_done:
            file_dict = [{'full_path': os.path.join(source_path, setting_folder, x),
                          'setting': setting_folder}
                         for x in os.listdir(os.path.join(source_path, setting_folder))]
            file_list += file_dict

    return file_list


def xi_wrapper(arguments):
    xi = subprocess.Popen(arguments)
    xi.communicate()
    return


def run_Xi(out_path, Xi_par_dict, nthreads, peak_files_dir):
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    searchres_path = os.path.join(out_path, Xi_par_dict['searchname'])
    if not os.path.exists(searchres_path):
        os.makedirs(searchres_path)

    # get filenames from folder with results
    file_paths_list = create_search_list(source_path=peak_files_dir, res_folder=searchres_path, rel=Xi_par_dict['missing_peaks'])
    n_parallel = 3
    xi_cmds = ['java',
               Xi_par_dict['memory'], '-cp', Xi_par_dict['Xi_path'], 'rappsilber.applications.Xi',
               '--fasta=' + Xi_par_dict['fasta'],
               '--xiconf=UseCPUs:' + str(nthreads / n_parallel)]

    cmd_list = []
    for single_file in file_paths_list:
        file_name = os.path.split(single_file['full_path'])[1]
        file_identifier = '_'.join(file_name.split('.')[0].split('_')[-2:])
        setting_folder = os.path.join(searchres_path, single_file['setting'])

        if not os.path.exists(setting_folder):
            os.makedirs(setting_folder)
        try:
            cnf = Xi_par_dict['config'][file_identifier]
        except KeyError:
            cnf = Xi_par_dict['config']['all']
        # create command specific for file
        if Xi_par_dict['missing_peaks'] is not None:
            if 'rel_' + '_'.join([str(x) for x in Xi_par_dict['missing_peaks'][1:]]) + '_' + Xi_par_dict['missing_peaks'][0] == single_file['setting']:
                arg_missing_peaks = ['--xiconf=missing_isotope_peaks:' + str(Xi_par_dict['missing_peaks'][1])]
                if len(Xi_par_dict['missing_peaks']) == 3:
                    arg_missing_peaks += ['--xiconf=missing_isotope_peaks_unknown_charge:' + str(Xi_par_dict['missing_peaks'][2])]
                cmd_list.append(xi_cmds + ['--peaks=' + single_file['full_path'],
                                           '--config=' + cnf,
                                           '--output=' + os.path.join(setting_folder,
                                                                      'rel_' + str(Xi_par_dict['missing_peaks'][1]) + '_' +
                                                                      os.path.split(single_file['full_path'])[1].split('.')[0] + '.csv')] +
                                arg_missing_peaks)
        else:
            cmd_list.append(xi_cmds + ['--peaks=' + single_file['full_path'],
                                       '--config=' + cnf,
                                       '--output=' + os.path.join(setting_folder,
                                                                  os.path.split(single_file['full_path'])[1].split('.')[
                                                                      0] + '.csv')])

    pool = Pool(processes=n_parallel)
    pool.map(xi_wrapper, cmd_list)
    pool.close()
    pool.join()


def run_Xi_multiple_settings(out_path, Xi_par_dict, nthreads, prepro_settings, peak_files_dir, xi_cmds, prefix):
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    searchres_path = os.path.join(out_path, Xi_par_dict['searchname'])
    if not os.path.exists(searchres_path):
        os.makedirs(searchres_path)

    # get filenames from folder with results
    file_paths_list = []
    for folder in prepro_settings:
        file_dict = [{'full_path': os.path.join(peak_files_dir, folder, x),
                      'setting': folder}
                      for x in os.listdir(os.path.join(peak_files_dir, folder))]
        file_paths_list += file_dict
    n_parallel = 4
    xi_cmds = ['java',
               Xi_par_dict['memory'], '-cp', Xi_par_dict['Xi_path'], 'rappsilber.applications.Xi',
               '--fasta=' + Xi_par_dict['fasta'],
               '--xiconf=UseCPUs:' + str(nthreads / n_parallel)] + xi_cmds

    cmd_list = []
    for single_file in file_paths_list:
        file_name = os.path.split(single_file['full_path'])[1]
        file_identifier = '_'.join(file_name.split('.')[0].split('_')[-2:])
        setting_folder = os.path.join(searchres_path, prefix + '_' + single_file['setting'])

        if not os.path.exists(setting_folder):
            os.makedirs(setting_folder)
        try:
            cnf = Xi_par_dict['config'][file_identifier]
        except KeyError:
            cnf = Xi_par_dict['config']['all']
        # create command specific for file

        cmd_list.append(xi_cmds + ['--peaks=' + single_file['full_path'],
                                   '--config=' + cnf,
                                   '--output=' + os.path.join(setting_folder, prefix + '_' +
                                                              os.path.split(single_file['full_path'])[1].split('.')[
                                                                  0] + '.csv')])

    pool = Pool(processes=n_parallel)
    pool.map(xi_wrapper, cmd_list)
    pool.close()
    pool.join()


def filter_single_file(out_path, file_path, Xi_par_dict_filter):
    # read in xi out csv
    # dtype required because else some types not converted correctly and xiFDR does not work
    xi_out = pd.read_csv(file_path, index_col=None, dtype=object)

    # delete rows where pep1 or pep2 are not there
    xi_out_filtered = xi_out[xi_out['peptide1 unique matched conservative'].astype(float) >= Xi_par_dict_filter['pep1_count>=']]
    xi_out_filtered = xi_out_filtered[xi_out_filtered['peptide2 unique matched conservative'].astype(float) >=
                                      Xi_par_dict_filter['pep2_count>=']]
    join_run_scan = lambda x: '_'.join([str(x['Run']), str(x['Scan'])])
    xi_out_filtered['psmid'] = xi_out_filtered.apply(join_run_scan, axis=1)
    xi_out_filtered['match score'] = xi_out_filtered['match score'].astype(float)
    xi_out_filtered = xi_out_filtered.sort_values('match score', ascending=False)
    xi_out_filtered = xi_out_filtered.drop_duplicates('Scan')

    # save as csv
    xi_out_filtered.to_csv(out_path)


def filter_results(xi_path, xi_par_dict):
    xires_path = os.path.join(xi_path, xi_par_dict['searchname'])
    filterres_path = os.path.join(xi_path, xi_par_dict['searchname'] + '_filtered')
    if not os.path.exists(filterres_path):
        os.makedirs(filterres_path)

    # function to return filepaths and names of xi output to do
    file_list = create_search_list(source_path=xires_path, res_folder=filterres_path)

    # run in loop filter functions

    for single_file in file_list:
        setting_folder = os.path.join(filterres_path, single_file['setting'])
        if not os.path.exists(setting_folder):
            os.makedirs(setting_folder)

        filter_single_file(file_path=single_file['full_path'],
                           out_path=os.path.join(setting_folder, os.path.split(single_file['full_path'])[1]),
                           Xi_par_dict_filter=xi_par_dict['filter'])


def call_XiFDR(infile, XiFDR_par, out_path):
    standard_cmds = ['java', '-Xmx1g', '-cp', XiFDR_par['xiFDR_path'], 'org.rappsilber.fdr.CSVinFDR']
    name = os.path.split(infile)[1].split('.')[0]
    # --filter="..." leads to problems when arguments should be given as list --> give single string
    map_settings = '"' + ','.join(['run:Run',
                                   'scan:Scan',
                                   'peptide1:Peptide1',
                                   'peptide2:Peptide2',
                                   'peptide length 1:LengthPeptide1',
                                   'peptide length 2:LengthPeptide2',
                                   'peptide link 1:Link1',
                                   'peptide link 2:Link2',
                                   'is decoy 1:Protein1decoy',
                                   'is decoy 2:Protein2decoy',
                                   'precursor charge:PrecoursorCharge',
                                   'score:match score',
                                   'accession1:Protein1',
                                   'accession2:Protein2',
                                   'description1:Fasta1',
                                   'description2:Fasta2',
                                   'peptide position 1:Start1',
                                   'peptide position 2:Start2',
                                   'crosslinker:Crosslinker'
                                   ]) + '"'

    settings = ' '.join([subprocess.list2cmdline(XiFDR_par['arguments']),
                         '--csvOutDir="' + out_path + '"',
                         '--delimiter=","',
                         '--map=' + map_settings])

    # call XiFDR with different FDRs each
    for fdr_settings in XiFDR_par['fdrs']:
        cmd_string = ' '.join([subprocess.list2cmdline(standard_cmds),
                               settings,
                               subprocess.list2cmdline(XiFDR_par['fdrs'][fdr_settings]),
                               '--csvBaseName=' + name,
                               infile])
        # print(cmd_string)
        # raw_input()
        XiFDR = subprocess.Popen(cmd_string)
        XiFDR.communicate()


def run_xiFDR(xi_path, xi_par_dict, filtered=True):
    # define paths and folders
    if filtered:
        source_path = os.path.join(xi_path, xi_par_dict['searchname'] + '_filtered')
    else:
        source_path = os.path.join(xi_path, xi_par_dict['searchname'])
    xifdr_path = os.path.join(xi_path, xi_par_dict['searchname'] + '_xiFDR')
    if not os.path.exists(xifdr_path):
        os.makedirs(xifdr_path)

    # function to return filepaths and names of xi output to do
    file_list = create_search_list_xifdr(source_path=source_path, res_folder=xifdr_path)

    # call xiFDR for single file
    for single_file in file_list:
        setting_folder = os.path.join(xifdr_path, single_file['setting'])
        if not os.path.exists(setting_folder):
            os.makedirs(setting_folder)

        call_XiFDR(infile=single_file['full_path'],
                   XiFDR_par=xi_par_dict['FDR'],
                   out_path=setting_folder)
