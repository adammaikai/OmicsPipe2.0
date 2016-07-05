#!/usr/bin/env python

import os
import sys
import stat
import urllib
from optparse import OptionParser
import webbrowser as browser
import argparse
from omics_pipe.parameters.default_parameters import default_parameters 
import yaml
from omics_pipe.utils import *
import runpy
import subprocess
import os.path
import csv
from raven import Client

def main():
    '''make command line interface, read in analysis type, execute correct analysis pipeline script'''  
    client = Client('http://44fd4bee8b9b4d6fa33e29d297c70cec:6e1f13b4911d4a5f915c25edc328c381@sentry.sulab.org/2')
    
    parser = argparse.ArgumentParser(prog = 'omics_pipe', description = 'Run omics_pipe')
    parser.add_argument('analysis_type', action = "store",
                        choices = ['RNAseq_Tuxedo', 'RNAseq_count_based', 'RNAseq_cancer_report', 'TCGA_download', 'Variant_annotation', 'RNAseq_TCGA', 'RNAseq_TCGA_counts', 
                                   'Tumorseq_MUTECT', 'miRNAseq_count_based', 'miRNAseq_tuxedo', 'WES_GATK_report', 'WES_GATK', 'WES_GATK_group_calling', 'WGS_GATK_optimized', 'WGS_GATK', 'WGS_GATK_group_calling',
                                    'SomaticInDels', 'ChIPseq_MACS', 'ChIPseq_HOMER',  'test', 'custom'], 
                        help = 'type of analysis to run: RNAseq_Tuxedo, RNAseq_count_based, RNAseq_cancer_report, TCGA_download, Variant_annotation, RNAseq_TCGA, WGS_GATK_optimized, RNAseq_TCGA_counts, Tumorseq_MUTECT, WES_GATK_report, miRNAseq_count_based, miRNAseq_tuxedo, WES_GATK, WGS_GATK, WES_GATK_group_calling, SomaticInDels, ChIPseq_MACS, ChIPseq_HOMER,  custom')
    parser.add_argument('parameter_file', action = "store", help = 'specify parameter file to use for analysis')
    parser.add_argument('--custom_script_path', action = "store", help = 'specify custom script file with full path (/example/script.py) to use for analysis if you specify analysis type as custom')
    parser.add_argument('--custom_script_name', action = "store", help = 'specify custom script file with full path (/example/script.py) to use for analysis if you specify analysis type as custom')
    parser.add_argument('--compression', action = "store", help = 'select bzip or gzip if your fastq files are compressed. Leave this option off if your files are uncompressed', choices = ['gzip', 'bzip'])

    args = parser.parse_args()
    
    print args
    print args.analysis_type
    print args.custom_script_path
    print args.custom_script_name
    print args.parameter_file
    print args.compression
    
    analysis = args.analysis_type     
    parameters = os.path.abspath(args.parameter_file)
    path = args.custom_script_path
    script = args.custom_script_name
    compression = args.compression
    
    stream = file(parameters, 'r')
    params = yaml.load(stream)
         
    default_parameters.update(params)   #Update default parameters to user defined parameter file
    p = Bunch(default_parameters)

    if type(p.SAMPLE_LIST) == str:                          #Handles reading a list of files from a text file
            sample_file = open(p.SAMPLE_LIST, 'r')        
            reader = csv.reader(sample_file)
            sample_list = [row for row in reader] 
            sample_list2 = [item for sublist in sample_list for item in sublist]
            default_parameters.update(SAMPLE_LIST = sample_list2)   #Update default parameters to user defined parameter file
            p = Bunch(default_parameters)
     
    check_create_dir(p.LOG_PATH)
    check_create_dir(p.FLAG_PATH)
    
    current_cwd = os.getcwd()
    os.chdir(p.WORKING_DIR)
    for x in os.listdir(p.WORKING_DIR):
        os.chmod(x,0755)
    
    os.environ["DRMAA_LIBRARY_PATH"] = p.DRMAA_PATH 
    
    start_time = time.time()
    print start_time

    decompress(p.RAW_DATA_DIR, args.compression)     #Check if files zipped, if so, unzip them   

    record, project = sumatra_start(p.REPOSITORY, p.SUMATRA_DB_PATH, p.RESULTS_PATH, p.WORKING_DIR, p.HG_USERNAME, p.SUMATRA_RUN_NAME, parameters) #Create repo and sumatra project, start recording
 
    os.chdir(p.WORKING_DIR)
    
    if args.analysis_type == 'RNAseq_Tuxedo':
        runpy.run_module('omics_pipe.RNAseq_Tuxedo', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'RNAseq_count_based':
        runpy.run_module('omics_pipe.RNAseq_count_based', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'RNAseq_cancer_report':
        runpy.run_module('omics_pipe.RNAseq_cancer_report', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'RNAseq_TCGA_counts':
        runpy.run_module('omics_pipe.RNAseq_TCGA_counts', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'RNAseq_TCGA':
        runpy.run_module('omics_pipe.RNAseq_TCGA', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'Tumorseq_MUTECT':
        runpy.run_module('omics_pipe.Tumorseq_MUTECT', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'miRNAseq_count_based':
        runpy.run_module('omics_pipe.miRNAseq_count_based', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'miRNAseq_tuxedo':
        runpy.run_module('omics_pipe.miRNAseq_tuxedo', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'WES_GATK':
        runpy.run_module('omics_pipe.WES_GATK', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'WES_GATK_group_calling':
        runpy.run_module('omics_pipe.WES_GATK_group_calling', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'WGS_GATK':
        runpy.run_module('omics_pipe.WGS_GATK', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'WES_GATK_report':
        runpy.run_module('omics_pipe.WES_GATK_report', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'WGS_GATK_group_calling':
        runpy.run_module('omics_pipe.WGS_GATK_group_calling', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'SomaticInDels':
        runpy.run_module('omics_pipe.SomaticInDels', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'ChIPseq_MACS':
        runpy.run_module('omics_pipe.ChIPseq_MACS', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'ChIPseq_HOMER':
        runpy.run_module('omics_pipe.ChIPseq_HOMER', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'TCGA_download':
        runpy.run_module('omics_pipe.TCGA_download', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'Variant_annotation':
        runpy.run_module('omics_pipe.Variant_annotation', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'WGS_GATK_optimized':
        runpy.run_module('omics_pipe.WGS_GATK_optimized', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'test':
        runpy.run_module('omics_pipe.test', run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)
    elif args.analysis_type == 'custom':
        os.chdir(path)
        sys.path.append(path)
        runpy.run_module(script, run_name="__main__", alter_sys = True)
        compress(p.RAW_DATA_DIR, args.compression)
        sumatra_end(start_time, record, project)
        sys.exit(0)   
    else:
        print 'Error: unsupported analysis type. Please try again.'       
    return

    
if __name__ == '__main__':
    main()

