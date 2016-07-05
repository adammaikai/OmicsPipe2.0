#!/usr/bin/env python

#from sumatra.projects import load_project
#from sumatra.parameters import build_parameters
#from sumatra.decorators import capture
from ruffus import *
import sys 
import os
import time
import datetime 
import drmaa
from omics_pipe.utils import *
from omics_pipe.parameters.default_parameters import default_parameters 
from omics_pipe.modules.fastqc import fastqc
from omics_pipe.modules.bwa import bwa_mem
from omics_pipe.modules.picard_mark_duplicates import picard_mark_duplicates
from omics_pipe.modules.GATK_preprocessing_WES import GATK_preprocessing_WES
from omics_pipe.modules.GATK_variant_discovery import GATK_variant_discovery
from omics_pipe.modules.GATK_variant_filtering import GATK_variant_filtering


p = Bunch(default_parameters)
os.chdir(p.WORKING_DIR)
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    

print p

for step in p.STEPS:
    vars()['inputList_' + step] = []
    for sample in p.SAMPLE_LIST:
        vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
    print vars()['inputList_' + step]

#FASTQC
@parallel(inputList_fastqc)
@check_if_uptodate(check_file_exists)
def run_fastqc(sample, fastqc_flag):
    fastqc(sample, fastqc_flag)
    return

#BWA
@parallel(inputList_bwa_mem)
@check_if_uptodate(check_file_exists)
def run_bwa_mem(sample, bwa_mem_flag):
    bwa_mem(sample, bwa_mem_flag)
    return

#picard_mark_duplicates
@parallel(inputList_picard_mark_duplicates)
@check_if_uptodate(check_file_exists)
@follows(run_bwa_mem)
def run_picard_mark_duplicates(sample, picard_mark_duplicates_flag):
    picard_mark_duplicates(sample, picard_mark_duplicates_flag)
    return

#GATK_preprocessing
@parallel(inputList_GATK_preprocessing_WES)
@check_if_uptodate(check_file_exists)
@follows(run_picard_mark_duplicates)
def run_GATK_preprocessing_WES(sample, GATK_preprocessing_WES_flag):
    GATK_preprocessing_WES(sample, GATK_preprocessing_WES_flag)
    return


#GATK_variant_discovery
@parallel(inputList_GATK_variant_discovery)
@check_if_uptodate(check_file_exists)
@follows(run_GATK_preprocessing_WES)
def run_GATK_variant_discovery(sample, GATK_variant_discovery_flag):
    GATK_variant_discovery(sample, GATK_variant_discovery_flag)
    return

#GATK_filter_variants
@parallel(inputList_GATK_variant_filtering)
@check_if_uptodate(check_file_exists)
@follows(run_GATK_variant_discovery)
def run_GATK_variant_filtering(sample, GATK_variant_filtering_flag):
    GATK_variant_filtering(sample, GATK_variant_filtering_flag)
    return

    
@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_GATK_variant_filtering, run_fastqc)
def last_function(sample, last_function_flag):
    print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
    pipeline_graph_output = p.FLAG_PATH + "/pipeline_" + sample + "_" + str(date) + ".pdf"
    #pipeline_printout_graph (pipeline_graph_output,'pdf', step, no_key_legend=False)
    stage = "last_function"
    flag_file = "%s/%s_%s_completed.flag" % (p.FLAG_PATH, stage, sample)
    open(flag_file, 'w').close()
    return   


if __name__ == '__main__':

    pipeline_run(p.STEP, multiprocess = p.PIPE_MULTIPROCESS, verbose = p.PIPE_VERBOSE, gnu_make_maximal_rebuild_mode = p.PIPE_REBUILD)

