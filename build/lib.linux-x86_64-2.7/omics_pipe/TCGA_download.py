#!/usr/bin/env python

#from sumatra.projects import load_project
from ruffus import *
import sys 
import os
import time
import datetime 
import drmaa
from omics_pipe.utils import *
from omics_pipe.modules.TCGA_download import TCGA_download
from omics_pipe.parameters.default_parameters import default_parameters 
p = Bunch(default_parameters)


os.chdir(p.WORKING_DIR)
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    

print p

if not p.SAMPLE_LIST:
    sample_list = get_TCGA_ID(p.TCGA_XML_FILE)
    for step in p.STEPS:
         vars()['inputList_' + step] = []
         for sample in sample_list:
             vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
         print vars()['inputList_' + step]
else:    
    for step in p.STEPS:
        vars()['inputList_' + step] = []
        for sample in p.SAMPLE_LIST:
            vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
        print vars()['inputList_' + step]

@parallel(inputList_TCGA_download)
@check_if_uptodate(check_file_exists)
def run_TCGA_download(sample, TCGA_download_flag):
    TCGA_download(sample, TCGA_download_flag)
    return


@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_TCGA_download)
def last_function(sample, last_function_flag):
    print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
    pipeline_graph_output = p.FLAG_PATH + "/pipeline_" + sample + "_" + str(date) + ".pdf"
    pipeline_printout_graph (pipeline_graph_output,'pdf', step, no_key_legend=False)
    stage = "last_function"
    flag_file = "%s/%s_%s_completed.flag" % (p.FLAG_PATH, stage, sample)
    open(flag_file, 'w').close()
    return   


if __name__ == '__main__':

    pipeline_run(p.STEP, multiprocess = p.PIPE_MULTIPROCESS, verbose = p.PIPE_VERBOSE, gnu_make_maximal_rebuild_mode = p.PIPE_REBUILD)


