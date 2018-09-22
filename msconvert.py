from __future__ import print_function
import os
import subprocess
import shutil


# function for calling MSConvert on files in raw_dir
def call_msconvert(raw_dir, out_path, filter_cmds, filetype='raw'):
    # create list with all files in directory
    files_list = os.listdir(raw_dir)

    # subselect all raw files
    raw_files_list = [x for x in files_list if ('.' + filetype) in x]

    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # if output folder already there, check if same files
    # if all same end function, else give these to msconvert
    else:
        files_done = [x.split('.mgf')[0] for x in os.listdir(out_path)]
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

    filter_formatted = []
    for i in range(len(filter_cmds)):
        filter_formatted.append('--filter')
        filter_formatted.append(filter_cmds[i])

    # start process and print output
    convert = subprocess.Popen(["msconvert.exe", '-f',
                                input_file_path, '-o', out_path, '--mgf'] + filter_formatted,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    while True:
        line = convert.stdout.readline()
        if not line: break
        print(line, end='\r')

    os.remove(input_file_path)


def main(result_path, raw_dir, msconvert_par_dict):

    for filter_setting in msconvert_par_dict:
        identifier = 'mscon_' + filter_setting
        # call msconvert for single setting
        call_msconvert(raw_dir=raw_dir, out_path=os.path.join(result_path, identifier),
                       filter_cmds=msconvert_par_dict[filter_setting], filetype='mgf')

        # move to folder all_prepro_files and rename
        target_dir = os.path.join(os.path.join(os.path.split(result_path)[0], 'All_prepro_peakfiles', identifier))

        if not os.path.exists(target_dir):
            shutil.copytree(os.path.join(result_path, identifier), target_dir)
            target_filelist = os.listdir(target_dir)
            for i in range(len(target_filelist)):
                old_name = target_filelist[i]
                new_name = identifier + '_' + '_'.join(old_name.split('_')[:2]) + '.mgf'
                os.rename(os.path.join(target_dir, old_name), os.path.join(target_dir, new_name))
