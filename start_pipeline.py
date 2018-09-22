import MQ
import OpenMS
import Xi_xiFDR
import sys
import getopt
import os
import logging
import time
import subprocess
import shutil
import bullseye
import msconvert
from multiprocessing import Pool, freeze_support

# ---------- functions -----------


def read_cmd_line_arg():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['outdir=', 'nthreads=', 'filedir=', 'log='])
    except getopt.GetoptError:
        print('start_pipeline.py --outdir <directory with config and for output> ' +
              '--nthreads <number of threads to use> --filedir <directory with files in different formats> ' +
              '--log <write to logfile, default False>')
        sys.exit()

    log = False
    for opt, arg in opts:
        if opt == '--outdir':
            dir = arg

        elif opt == '--nthreads':
            nthr = int(arg)

        elif opt == '--filedir':
            rawdir = arg

        elif opt == '--log':
            log = arg

    return dir, nthr, rawdir, log


def return_filelist(namefile):
    filenames_file = open(namefile)
    files = []
    for line in filenames_file:
        files.append(line.replace('\n', ''))

    return files


def start_logger(path, name):
    logger1 = logging.getLogger('basic_logger')
    handler1 = logging.FileHandler(os.path.join(path, name))
    handler1.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
    logger1.setLevel(logging.DEBUG)
    logger1.addHandler(handler1)
    return logger1


def no_preprocessing(mgf_folder, raw_folder, directory):
    # convert to mgf
    if not os.path.exists(mgf_folder):
        os.makedirs(mgf_folder)

        # create list with all files in directory
        files_list = os.listdir(raw_folder)
        files_list = [os.path.join(raw_folder, x) for x in files_list]

        # write msconvert input file
        input_file_path = os.path.join(mgf_folder, 'file_list.txt')
        input_file = open(input_file_path, 'w')
        for i in files_list:
            input_file.write("%s\n" % i)
        input_file.close()

        # start process and print output
        convert = subprocess.Popen(["msconvert", '-f',
                                    input_file_path, '-o', mgf_folder, '--mgf'],
                                   stdout=subprocess.PIPE,
                                   universal_newlines=True)
        convert.communicate()

        os.remove(input_file_path)

    # move to folder all_prepro_files
    target_dir = os.path.join(directory, 'All_prepro_peakfiles', 'other_NoPrePro')
    if not os.path.exists(target_dir):
        shutil.copytree(mgf_folder, target_dir)
        target_filelist = os.listdir(target_dir)
        for i in range(len(target_filelist)):
            old_name = target_filelist[i]
            new_name = 'other_NoPrePro_'+'_'.join(old_name.split('_')[:2]) + '.mgf'
            os.rename(os.path.join(target_dir, old_name), os.path.join(target_dir, new_name))


# ------------ main --------------
# cmd line argument of folder to do it in
if __name__ == '__main__':
    directory, nthreads, file_dir, log = read_cmd_line_arg()
    freeze_support()
    # import config specific for run
    sys.path.append(os.path.join(directory, 'config'))
    import config_run

    # define paths
    mq_res_dir = os.path.join(directory, 'MQ')
    OpenMS_res_dir = os.path.join(directory, 'OpenMS')
    Xi_res_dir = os.path.join(directory, 'Xi')
    mscon_dir = os.path.join(directory, 'mscon')
    raw_dir = os.path.join(file_dir, 'raw')
    mzml_dir = os.path.join(file_dir, 'mzML')
    mgf_dir = os.path.join(file_dir, 'mgf')
    mzxml_dir = os.path.join(file_dir, 'mzXML')

    if log: logger = start_logger(directory, 'log_all.log')

    # ------------ MaxQuant --------------
    if config_run.MQ_par is not None:
        start = time.clock()
        MQ.main(parameter_dict=config_run.MQ_par['MQ_prepro_var'],
                MaxQuantCmd_file=config_run.MQ_par['MQCmd_exe'],
                mqpar_template_file=config_run.MQ_par['mqpar_template'],
                MQ_FASTA_file=config_run.MQ_par['FASTA'],
                result_path=mq_res_dir,
                raw_dir=raw_dir,
                num_threads=nthreads,
                keep_res=config_run.MQ_par['keep_results'])
        if log: logger.info(str(config_run.MQ_par))
        if log: logger.info('MaxQuant took ' + str((time.clock() - start) / 60) + 'min')

    # ----------- OpenMS --------------
    if config_run.OpenMS_par is not None:
        start = time.clock()
        OpenMS.main(result_path=OpenMS_res_dir,
                    OpenMS_par_dict=config_run.OpenMS_par['settings_nodes'],
                    raw_dir=raw_dir,
                    mzml_dir=mzml_dir,
                    nthreads=nthreads,
                    keep_res=config_run.OpenMS_par['keep_fileconv_res'])
        if log: logger.info(str(config_run.OpenMS_par))
        if log: logger.info('OpenMS took ' + str((time.clock() - start) / 60) + ' min')

    if config_run.mscon_par is not None:
        msconvert.main(result_path=mscon_dir, raw_dir=mgf_dir, msconvert_par_dict=config_run.mscon_par)

    # ------- other prepro ------------
    if config_run.do_other_prepro:
        no_preprocessing(mgf_folder=mgf_dir, raw_folder=raw_dir, directory=directory)

    # ------- Search and FDR ------
    # Xi search
    start = time.clock()
    Xi_xiFDR.run_Xi(out_path=Xi_res_dir, peak_files_dir=os.path.join(directory, 'All_prepro_peakfiles'),
                    Xi_par_dict=config_run.Xi_par, nthreads=nthreads)
    if log: logger.info('Xi took ' + str((time.clock() - start) / 60) + ' min')

    # filter Xi results
    Xi_xiFDR.filter_results(xi_path=Xi_res_dir, xi_par_dict=config_run.Xi_par)

    # xiFDR
    Xi_xiFDR.run_xiFDR(xi_path=Xi_res_dir, xi_par_dict=config_run.Xi_par)

