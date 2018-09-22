import os
from lxml import etree
import subprocess
import time
import itertools
import pandas
import zipfile
import shutil


# -------- functions ----------


def write_MQ_settings(parameter_dict):
    # sort keys (settings to vary) of parameter dict and create list of variations in same order
    parameter_list = []
    keys_sorted = sorted([x for x in parameter_dict])
    for key in keys_sorted:
        parameter_list.append(parameter_dict[key])

    # combinations of settings via cartesian product
    combinations = list(itertools.product(*parameter_list))

    # save name and settings in separate dict
    # loop for each combination through each variable
    settings_combined_dict = {}
    for combination in combinations:
        setting = {}
        name = ''
        for j in range(len(keys_sorted)):
            setting[keys_sorted[j]] = combination[j]
            if keys_sorted[j] == 'deiso':
                name = name + keys_sorted[j] + str(int(combination[j])) + '_'
            else:
                name = name + keys_sorted[j] + str(combination[j]) + '_'

        # delete last _ from name
        # save settings in dict under name (/identifier)
        name = name[:-1]
        settings_combined_dict[name] = setting

    return settings_combined_dict


# wrapper for calling MQ from writing new xml to moving results
# MaxQuantCmd_file, mqpar_template_file, MQ_FASTA_file,result_path, raw_dir: string, paths to respective files / folders
# identifier: int or string, identifier for folder / settings
def call_MaxQuant(maxquant_param, MaxQuantCmd_file, mqpar_template_file, MQ_FASTA_file, identifier,
                  result_path, raw_dir, num_threads):
    # create folder, if alreadz there, skip MQ step for this setting
    # not checked how many files in folder, can also be empty (bec results can optionally be removed)
    current_dir = result_path + '\\MQ_' + str(identifier)
    if not os.path.exists(current_dir):
        os.makedirs(current_dir)
    else:
        print(str(identifier) + ' already done, not done again.')
        return False

    print(str(identifier) + ' started.')

    # rename raw files so setting appears in apl file, keep original names in varaible
    original_names = rename_raw_files(raw_dir=raw_dir, identifier=identifier)

    # write new xml file from template
    write_param_xml(template_par_file=mqpar_template_file,
                    raw_dir=raw_dir,
                    maxquant_param=maxquant_param,
                    out_path=current_dir,
                    num_threads=num_threads,
                    MQ_FASTA_file=MQ_FASTA_file)

    # call maxquant with new xml
    MaxQuant = subprocess.Popen([MaxQuantCmd_file, (current_dir + '\\mqpar_mod.xml')],
                                shell=True)  # shell=True wegen killing process?

    # terminate MQ if files '' are in directory (function)
    if terminate_MQ(current_dir):
        subprocess.call(['taskkill', '/F', '/T', '/PID', str(MaxQuant.pid)])

    # move files from raw-directory
    move_MQ_results(raw_dir=raw_dir, current_dir=current_dir)

    # rename raw files back
    files_list = os.listdir(raw_dir)
    raw_files_list_settings = [x for x in files_list if '.raw' in x]
    for i in range(len(raw_files_list_settings)):
        os.rename(raw_dir + '\\' + raw_files_list_settings[i],
                  raw_dir + '\\' + original_names[i])

    return True


def rename_raw_files(raw_dir, identifier):
    # create list with all files in directory
    files_list = os.listdir(raw_dir)
    # subselect all raw files
    raw_files_list_old = [x for x in files_list if '.raw' in x]
    raw_old_splitted = [x.split('_') for x in raw_files_list_old]

    raw_files_list_new = ['_'.join(['MQ', identifier] + x[:2]) + '.raw' for x in raw_old_splitted]

    # raw_files_list = [(raw_dir + '\\' + x) for x in raw_files_list_old]
    for i in range(len(raw_files_list_old)):
        os.rename(raw_dir + '\\' + raw_files_list_old[i],
                  raw_dir + '\\' + raw_files_list_new[i])
    pandas.Series(raw_files_list_old).to_csv(os.path.join(raw_dir, 'original_names.csv'), index=False)

    return raw_files_list_old


# function to edit xml file with specific parameters
# template_par_file: string, complete path to xml template with MQ parameters
# raw_dir: string, directory with raw files
# maxquant_param: dict, name for attribute and value - hard coded, all values must be given, also not more
# out_path: string, comnplete path to where result should be saved
def write_param_xml(template_par_file, raw_dir, maxquant_param, out_path, num_threads, MQ_FASTA_file):
    # create list with all files in directory
    files_list = os.listdir(raw_dir)
    # subselect all raw files, replace / by \ --> will be pasted to xml
    raw_files_list = [x for x in files_list if '.raw' in x]
    raw_files_list = [(raw_dir + '\\' + x) for x in raw_files_list]

    # parse xml with option (needed to enable pretty_print option, else elements in one row)
    parser = etree.XMLParser(remove_blank_text=True)
    template_xml = etree.parse(template_par_file, parser)

    # modify xml
    new_xml = modify_xml_parameters(template_xml=template_xml, maxquant_param=maxquant_param,
                                    raw_files_list=raw_files_list, out_path=out_path, num_threads=num_threads,
                                    MQ_FASTA_file=MQ_FASTA_file)

    # write new xml, parse (with parser), write --> else some elements in one row, bug?
    mqpar_mod_file = out_path + '\\mqpar_mod.xml'
    new_xml.write(mqpar_mod_file)
    new_xml_edit = etree.parse(mqpar_mod_file, parser)
    new_xml_edit.write(mqpar_mod_file, pretty_print=True)


# function to modify template xml
# template_xml: xml parsed with etree
# maxquant_param: dict with parameters to change, elements hard coded, not changeable
# raw_files_list: list with complete paths to raw files - will be pasted to xml
def modify_xml_parameters(template_xml, maxquant_param, raw_files_list, out_path, num_threads, MQ_FASTA_file):
    template_xml_new = template_xml
    # change parameters depending on number of raw files
    for i in range(0, len(raw_files_list)):
        # parse paths to raw files
        etree.SubElement(template_xml_new.findall('filePaths')[0], 'string')
        template_xml_new.findall('filePaths')[0][i].text = raw_files_list[i]

        # parse fixed elements
        etree.SubElement(template_xml_new.findall('experiments')[0], 'string')

        etree.SubElement(template_xml_new.findall('fractions')[0], 'short')
        template_xml_new.findall('fractions')[0][i].text = '32767'

        etree.SubElement(template_xml_new.findall('paramGroupIndices')[0], 'int')
        template_xml_new.findall('paramGroupIndices')[0][i].text = '0'

    # parse FASTA file with changed seperator
    etree.SubElement(template_xml_new.findall('fastaFiles')[0], 'string')
    template_xml_new.findall('fastaFiles')[0][0].text = MQ_FASTA_file

    # change editable parameters
    template_xml_new.findall('numThreads')[0].text = str(num_threads)
    template_xml_new.findall('fixedCombinedFolder')[0].text = out_path

    template_xml_new.findall('topxWindow')[0].text = str(maxquant_param['Win'])

    # only change for FTMS
    for i in range(0, len(template_xml_new.findall('msmsParamsArray')[0])):
        if template_xml_new.findall('msmsParamsArray')[0][i].get('Name') == 'FTMS':
            template_xml_new.findall('msmsParamsArray')[0][i].set('Topx', str(maxquant_param['Top']))
            template_xml_new.findall('msmsParamsArray')[0][i].set('Deisotope', str(maxquant_param['deiso']).lower())

    return template_xml_new


# function to determine if file 'Preparing_searches...finished.txt' is in proc (to terminate MQ)
def terminate_MQ(MQ_results_path):
    # add time to avoid 'folder not found'
    time.sleep(10)
    MaxQ_proc_path = MQ_results_path + '/combined/proc'
    finished = False
    while finished == False:
        files_MaxQ = os.listdir(MaxQ_proc_path)
        ## which file
        file_finished = [x for x in files_MaxQ if
                         'Combining_apl_files_' in x]
        file_finished = [x for x in file_finished if 'finished.txt' in x]
        if len(file_finished) > 0:
            finished = True
            # time.sleep(10)
        else:
            time.sleep(60)

    return finished


# moves files from raw file directory to other directory (all except .raw)
def move_MQ_results(raw_dir, current_dir):
    # create list with all files in directory
    files_list = os.listdir(raw_dir)

    # subselect all NOT raw files
    not_raw_files_list = [x for x in files_list if not '.raw' in x]
    formats_excl = [x for x in not_raw_files_list if not x in ['mgf', 'mzML']]
    for singlefile in formats_excl:
        print(raw_dir, singlefile, current_dir)
        os.rename((raw_dir + '\\' + singlefile), (current_dir + '\\' + singlefile))


# return zip file with correct name
def return_zipfile(path, single_group):
    group_zip = os.path.join(path, single_group)
    i = 1
    while os.path.exists(group_zip + '.zip'):
        group_zip = group_zip + str(i)
        i = i + 1
    zip_file = zipfile.ZipFile(file=group_zip + '.zip', mode='w')
    return zip_file


# write apl files to zip
def zip_apl_files(folder, zip_file):
    # select folder 'p0'
    result_files = os.listdir(os.path.join(folder, 'p0'))

    # select apl files with FTMS in them which are >0
    result_files = [x for x in result_files if '.apl' == x[-4:]]
    result_files = [x for x in result_files if 'FTMS' in x]
    result_files = [x for x in result_files if os.path.getsize(os.path.join(folder, 'p0', x)) > 0]

    if len(result_files) > 2: print('More than 2 apl files resulting!')

    for j in result_files:
        zip_file.write(os.path.join(folder, 'p0', j), arcname=j)


def zip_move_files(path_mq_res, out_path, keep_res):
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    files_list = os.listdir(path_mq_res)
    folder_list = [x for x in files_list if not '.index' in x]
    folder_list = [x for x in folder_list if not x in ['combined', 'mqpar_mod.xml', 'original_names.csv']]

    for folder_single_file in folder_list:
        # create zip_file
        zip_file = zipfile.ZipFile(file=os.path.join(out_path, folder_single_file + '.zip'), mode='w')

        # move files to zip
        zip_apl_files(folder=os.path.join(path_mq_res, folder_single_file), zip_file=zip_file)

    if not keep_res:
        for i in files_list:
            try:
                os.remove(os.path.join(path_mq_res, i))
            except OSError:
                shutil.rmtree(os.path.join(path_mq_res, i), ignore_errors=False)


def main(parameter_dict, MaxQuantCmd_file, mqpar_template_file, MQ_FASTA_file,
         result_path, raw_dir, num_threads, keep_res=True):
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    # set settings to loop over
    settings_dict = write_MQ_settings(parameter_dict=parameter_dict)  # config_run.MQ_par['MQ_peakpicking_var']

    # call MQ
    for settings in settings_dict:
        # correct terminate searches? takes forever.
        new = call_MaxQuant(maxquant_param=settings_dict[settings],
                            MaxQuantCmd_file=MaxQuantCmd_file,
                            mqpar_template_file=mqpar_template_file,
                            MQ_FASTA_file=MQ_FASTA_file,
                            identifier=settings,
                            result_path=result_path,
                            raw_dir=raw_dir,
                            num_threads=num_threads
                            )

        if new:
            # move correct apl files to folder and zip them
            # option whether to keep results
            zip_move_files(path_mq_res=(result_path + '\\MQ_' + str(settings)),
                           out_path=os.path.join(os.path.split(result_path)[0], 'All_prepro_peakfiles',
                                                 'MQ_' + settings),
                           keep_res=keep_res)


if __name__ == '__main__':
    main(parameter_dict, MaxQuantCmd_file, mqpar_template_file, MQ_FASTA_file, result_path, raw_dir, num_threads)
