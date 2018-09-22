from __future__ import print_function
import os
import subprocess
import fileinput
import itertools
import pandas


# create dict with entry for each node
# in each entry dataframe as saved in setting folder
def determine_settings_combos(algorithms, OpenMS_par_dict, out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    param_dict_combined = {}
    for algorithm in algorithms:
        # extract parameters of node to do from config variable
        parameters = OpenMS_par_dict[algorithm]
        # return dict of combinations of parameters for one node
        settings_dict = combine_settings(parameters)
        settings_filepath = os.path.join(out_path, algorithm + '_settings.csv')

        # if settings file already exists, combine settings so they won't be done again
        if os.path.exists(settings_filepath):
            # read in old settings file
            df_settings_old = pandas.DataFrame.from_csv(settings_filepath)
            # create dataframe from dict and concat old and new
            df_settings_new = pandas.DataFrame.from_dict(settings_dict, orient='index')
            df_settings_merged = pandas.concat([df_settings_old, df_settings_new])
            # read out identifiers not duplicated and subset to those
            duplicated = df_settings_merged.loc[:, df_settings_merged.columns != 'Identifier'].duplicated()
            non_duplicated = [not x for x in duplicated]
            df_settings_merged = df_settings_merged[non_duplicated]
            # read out identifiers in df, excluding empty entries
            identifier = df_settings_merged.Identifier[[not x for x in df_settings_merged.Identifier.isnull()]].tolist()
            # append new identifiers to df column, ex
            new_ids = [x for x in range(len(df_settings_merged)) if not x in identifier]
            # append new ids to entries where identifier is empty
            if len(new_ids) != 0:
                new_settings = df_settings_merged.Identifier.isnull()
                df_settings_merged.loc[new_settings, 'Identifier'] = new_ids
            # turn type of identifiers back to integer
            df_settings_merged['Identifier'] = df_settings_merged['Identifier'].astype(int)
            df_settings_merged.to_csv(settings_filepath)
            # merge new settings df with df of all settings
            df_current_settings = pandas.concat([df_settings_new, df_settings_merged])
            # select only duplicates of new settings
            df_settings = df_current_settings[
                df_current_settings.loc[:, df_current_settings.columns != 'Identifier'].duplicated()]
            df_settings['Identifier'] = df_settings['Identifier'].astype(int)

        # save settings as df to csv with identifiers assigned
        else:
            df_settings = pandas.DataFrame.from_dict(settings_dict, orient='index')
            df_settings['Identifier'] = range(len(df_settings))
            df_settings.to_csv(settings_filepath)

        # add df for single node to dict
        param_dict_combined[algorithm] = df_settings
    return param_dict_combined


# return dict of combinations of parameters for one node
def combine_settings(parameter_dict):
    # keep track of order
    parameter_list = []
    keys_sorted = sorted([x for x in parameter_dict])
    for key in keys_sorted:
        parameters = parameter_dict[key] if isinstance(parameter_dict[key], list) else [parameter_dict[key]]
        parameter_list.append(parameters)

    # set combinations via cartesian product
    combinations = list(itertools.product(*parameter_list))

    # save name and settings in separate dict
    counter = 0
    single_settings_dict = {}
    for i in combinations:
        setting = {}
        for j in range(len(keys_sorted)):
            setting[keys_sorted[j]] = i[j]

        single_settings_dict[str(counter)] = setting
        counter += 1
    return single_settings_dict


# function for calling MSConvert on files in raw_dir
def convert_to_mzML(raw_dir, out_path, filetype='raw'):
    # create list with all files in directory
    files_list = os.listdir(raw_dir)

    # subselect all raw files
    raw_files_list = [x for x in files_list if ('.' + filetype) in x]

    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # if output folder already there, check if same files
    # if all same end function, else give these to msconvert
    else:
        files_done = [x.split('.mzML')[0] for x in os.listdir(out_path)]
        raw_files_split = [x.split('.' + filetype)[0] for x in raw_files_list]
        not_done = [x for x in raw_files_split if not x in files_done]
        if len(not_done) == 0:
            return
        else:
            raw_files_list = [x + '.' + filetype for x in not_done]

    raw_files_list = [os.path.join(raw_dir, x) for x in raw_files_list]
    # write msconvert input file
    input_file_path = os.path.join(out_path, 'file_list.txt')
    input_file = open(input_file_path, 'w')
    for i in raw_files_list:
        input_file.write("%s\n" % i)
    input_file.close()

    # start process and print output
    convert = subprocess.Popen(["msconvert", '-f', input_file_path, '-o', out_path],
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    while True:
        line = convert.stdout.readline()
        if not line: break
        print(line, end='\r')


# wrapper for TOPP command line functions
# functions which are built 'standard' should work
#
# in_path: string - path to folder of step before (with subfolders of diff parameters)
# parameters: dict - name is name of option, after that value for option (numbers are converted to strings)
# out_path: string - path to folder where results should be stored
# algorithm: string - full name of algorithm which should be used
# features_path: string - needed for HighResPrecursorMassCorrector, path where features are stored
# out_type: string - needed for file converter, type for resulting file
def call_TOPP_tool(in_path, out_path, algorithm, num_threads, features_path='', out_type='.mgf', parameters=''):
    if os.path.exists(out_path):
        return
    else:
        os.makedirs(out_path)
        # generate input and output lists - return: mzml_files_list, out_files_list
        mzml_files_list, out_files_list = generate_in_out_list(mzml_in_path=in_path,
                                                               out_path=out_path,
                                                               algorithm=algorithm,
                                                               out_type=out_type)
        # return feature files for single folders (work on this to select specific feature files corresponding to mzml)
        if algorithm == 'HighResPrecursorMassCorrector':
            featxml_files_list = generate_feature_files_list(features_in_path=features_path,
                                                             mzml_files_list=mzml_files_list)
        else:
            featxml_files_list = ''

        # generate list of commands based on given parameternames and corresponding values
        commands_parameters = []
        if len(parameters) != 0:
            del parameters['Identifier']
            for i in range(len(parameters)):
                command = [parameters.index[i], str(parameters[i])]
                commands_parameters = commands_parameters + command

        for i in range(len(mzml_files_list)):
            # generate commands specific for single file
            all_commands = generate_commands(algorithm=algorithm, num_threads=num_threads,
                                             mzml_files_list=mzml_files_list,
                                             out_files_list=out_files_list, commands_parameters=commands_parameters,
                                             i=i, featxml_files_list=featxml_files_list)
            # call TOPP tool
            subprocess_TOPP(all_commands)


# generate lists with full paths to single files (in and out respectivly)
def generate_in_out_list(mzml_in_path, out_path, algorithm, out_type):
    files_list = os.listdir(mzml_in_path)
    mzml_names_list = [x for x in files_list if '.mzML' in x]

    # generate output list, use different endings for some tools
    if algorithm == 'FeatureFinderCentroided':
        out_files_list = [os.path.join(out_path, x.replace('.mzML', '.featureXML')) for x in mzml_names_list]
    elif algorithm == 'FileConverter':
        out_files_list = [os.path.join(out_path, x.replace('.mzML', out_type)) for x in mzml_names_list]
    else:
        out_files_list = [os.path.join(out_path, x) for x in mzml_names_list]
    mzml_files_list = [os.path.join(mzml_in_path, x) for x in mzml_names_list]

    return mzml_files_list, out_files_list


# return list of feature files (needed for HighResPrecursorMassCorrector)
def generate_feature_files_list(features_in_path, mzml_files_list):
    feature_list = os.listdir(features_in_path)
    featxml_names_list = [x for x in feature_list if '.featureXML' in x]
    # match order to mzml
    # extract names from mzml_files_list
    mzml_names_list = list()
    for i in range(len(mzml_files_list)):
        name = os.path.split(mzml_files_list[i])[1]
        name = name.split('.mzML')[0] + '.featureXML'
        mzml_names_list = mzml_names_list + [name]

    # match order of feature files to order of mzml files
    featxml_names_list_sorted = list()
    for i in mzml_names_list:
        match = [x for x in featxml_names_list if i in x]
        featxml_names_list_sorted = featxml_names_list_sorted + match

    featxml_files_list = [os.path.join(features_in_path, x) for x in featxml_names_list_sorted]
    return featxml_files_list


# return list of commands to give to calling of TOPP tool
def generate_commands(algorithm, num_threads, mzml_files_list, out_files_list, commands_parameters, i,
                      featxml_files_list):
    if algorithm == 'FileConverter':
        algorithm = 'C:\\Program Files\\OpenMS-2.1.0\\bin\\' + algorithm
        commands_general = [algorithm,
                            '-threads', str(num_threads),
                            '-in', mzml_files_list[i],
                            '-out', out_files_list[i]]
    else:
        commands_general = [algorithm,
                            '-threads', str(num_threads),
                            '-in', mzml_files_list[i],
                            '-out', out_files_list[i]]

        if algorithm == 'HighResPrecursorMassCorrector':
            commands_general = commands_general + ['-feature:in', featxml_files_list[i]]
    all_commands = commands_general + commands_parameters

    return all_commands


# start subprocess for TOPP tool
def subprocess_TOPP(all_commands):
    TOPP_tool = subprocess.Popen(all_commands, stdout=subprocess.PIPE, universal_newlines=True)
    # print output of console
    while True:
        line = TOPP_tool.stdout.readline()
        if not line: break
        print(line, end='\r')


def modify_scans(path, identifier, prefix='scan=', OpenMS_prefix=True):
    files = os.listdir(path)
    for filename in files:
        if filename[:2] == 'PP': continue
        if OpenMS_prefix:
            newfilename = '_'.join(['OpenMS', identifier, '_'.join(filename.split('_')[:2])])
        else:
            newfilename = '_'.join(['other', identifier, '_'.join(filename.split('_')[:2])])
        for line in fileinput.FileInput(os.path.join(path, filename), inplace=True):
            if prefix in line:
                Scan = line[line.find(prefix) + len(prefix):].split('_')
                Scan[1] = newfilename
                newline = line[:line.find(prefix)] + prefix + '_'.join(Scan)
                print(newline)
            else:
                print(line, end="")
        os.rename(os.path.join(path, filename),
                  os.path.join(path, newfilename + '.mgf'))


def write_to_ID_table(file_dir, setting, ID_table):
    result_dict = {}
    # create list with all files in directory
    files_list = os.listdir(file_dir)

    files = [x for x in files_list if '.mgf' in x]

    for i in files:
        # for j in settings:
        result_dict[i + setting] = {'File': i, 'Algorithm': 'OpenMS',
                                    'Setting': setting, 'ID': ''}

    df_settings = pandas.DataFrame.from_dict(result_dict, orient='index')
    df_ID_table = pandas.read_csv(ID_table,
                                  index_col=False, header=0)

    df_merged = pandas.concat([df_settings, df_ID_table])
    df_merged = df_merged[['File', 'Algorithm', 'Setting', 'ID']]
    df_merged.to_csv(ID_table, index=False)


def move_files(fileconv_path, out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    files_list = os.listdir(fileconv_path)

    for single_file in files_list:
        os.rename(os.path.join(fileconv_path, single_file),
                  os.path.join(out_path, single_file))


def peakpickingonly(result_path, raw_dir, mzml_dir, nthreads, Top, Win):
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    param_dict_combined = determine_settings_combos(
        algorithms=['SpectraFilter'],
        OpenMS_par_dict={'SpectraFilter': {'-algorithm:windowsize': Win,
                                           '-algorithm:peakcount': Top,
                                           '-algorithm:movetype': ['jump']}},
        out_path=os.path.join(result_path, 'Settings'))

    convert_to_mzML(raw_dir=raw_dir, out_path=mzml_dir)
    for SF_setting in param_dict_combined['SpectraFilter'].iterrows():
        call_TOPP_tool(algorithm='SpectraFilterWindowMower',
                       in_path=mzml_dir,
                       num_threads=nthreads,
                       parameters=SF_setting[1],
                       out_path=os.path.join(result_path, 'SpectraFilter'))

        # FileConverter
        call_TOPP_tool(algorithm='FileConverter',
                       in_path=os.path.join(result_path, 'SpectraFilter'),
                       num_threads=nthreads,
                       out_type='.mgf',
                       out_path=os.path.join(result_path, 'FileConverter'))

        modify_scans(path=os.path.join(result_path, 'FileConverter'),
                     identifier='PPOnlyWin' + str(Win) + 'Top' + str(Top), OpenMS_prefix=False)

        move_files(fileconv_path=os.path.join(result_path, 'FileConverter'),
                   out_path=os.path.join(os.path.split(result_path)[0], 'All_prepro_peakfiles',
                                         ('other_' + 'PPOnlyWin' + str(Win) + 'Top' + str(Top))))


def main(result_path, OpenMS_par_dict, raw_dir, mzml_dir, nthreads, keep_res):
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    # save paths to all node folders
    result_folders = {}
    for folder in ['PeakPicker', 'FeatureFinder', 'MassCorrector', 'SpectraFilter', 'FileConverter']:
        result_folders[folder] = os.path.join(result_path, folder)
        if not os.path.exists(os.path.join(result_path, folder)):
            os.makedirs(os.path.join(result_path, folder))

    param_dict_combined = determine_settings_combos(
        algorithms=['PeakPicker', 'FeatureFinder', 'MassCorrector', 'SpectraFilter'],
        OpenMS_par_dict=OpenMS_par_dict,
        out_path=os.path.join(result_path, 'Settings'))

    # convert files
    convert_to_mzML(raw_dir=raw_dir, out_path=mzml_dir)

    # peakpickerhires
    for PP_setting in param_dict_combined['PeakPicker'].iterrows():
        PP_folder_ident = 'PP' + str(PP_setting[1].Identifier)
        call_TOPP_tool(algorithm='PeakPickerHiRes',
                       in_path=mzml_dir, num_threads=nthreads,
                       parameters=PP_setting[1],
                       out_path=os.path.join(result_folders['PeakPicker'], PP_folder_ident))

        # featurefinderCentroided
        for FF_setting in param_dict_combined['FeatureFinder'].iterrows():
            FF_folder_ident = PP_folder_ident + '_FF' + str(FF_setting[1].Identifier)
            call_TOPP_tool(algorithm='FeatureFinderCentroided',
                           in_path=os.path.join(result_folders['PeakPicker'], PP_folder_ident),
                           num_threads=nthreads,
                           parameters=FF_setting[1],
                           out_path=os.path.join(result_folders['FeatureFinder'], FF_folder_ident))

            # HighResPrecursorMassCorrector
            for MC_setting in param_dict_combined['MassCorrector'].iterrows():
                MC_folder_ident = FF_folder_ident + '_MC' + str(MC_setting[1].Identifier)
                call_TOPP_tool(algorithm='HighResPrecursorMassCorrector',
                               in_path=os.path.join(result_folders['PeakPicker'], PP_folder_ident),
                               num_threads=nthreads,
                               features_path=os.path.join(result_folders['FeatureFinder'], FF_folder_ident),
                               parameters=MC_setting[1],
                               out_path=os.path.join(result_folders['MassCorrector'], MC_folder_ident))

                # SpectralFilterWindowMower
                for SF_setting in param_dict_combined['SpectraFilter'].iterrows():
                    SF_folder_ident = MC_folder_ident + '_SF' + str(SF_setting[1].Identifier)
                    call_TOPP_tool(algorithm='SpectraFilterWindowMower',
                                   in_path=os.path.join(result_folders['MassCorrector'], MC_folder_ident),
                                   num_threads=nthreads,
                                   parameters=SF_setting[1],
                                   out_path=os.path.join(result_folders['SpectraFilter'], SF_folder_ident))

                    # FileConverter
                    call_TOPP_tool(algorithm='FileConverter',
                                   in_path=os.path.join(result_folders['SpectraFilter'], SF_folder_ident),
                                   num_threads=nthreads,
                                   out_type='.mgf',
                                   out_path=os.path.join(result_folders['FileConverter'], SF_folder_ident))
                    if os.path.exists(os.path.join(result_folders['FileConverter'], SF_folder_ident)):
                        modify_scans(path=os.path.join(result_folders['FileConverter'], SF_folder_ident),
                                     identifier=SF_folder_ident)
                        move_files(fileconv_path=os.path.join(result_folders['FileConverter'], SF_folder_ident),
                                   out_path=os.path.join(os.path.split(result_path)[0], 'All_prepro_peakfiles',
                                                         ('OpenMS_' + SF_folder_ident)))
