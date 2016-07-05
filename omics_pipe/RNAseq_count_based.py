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
import csv
from omics_pipe.utils import *
from omics_pipe.modules.fastqc import fastqc
from omics_pipe.modules.star import star
from omics_pipe.modules.htseq import htseq
from omics_pipe.modules.RNAseq_report_counts import RNAseq_report_counts

from omics_pipe.parameters.default_parameters import default_parameters 
p = Bunch(default_parameters)

os.chdir(p.WORKING_DIR)
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    

with open(p.DESEQ_META,"rb") as infile:
    next(infile, None)
    reader=csv.reader(infile) 
    samples = [x[0] for x in reader]

for step in p.STEPS:
    vars()['inputList_' + step] = []
    for sample in samples:
        vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
    
@parallel(inputList_fastqc)
@check_if_uptodate(check_file_exists)
def run_fastqc(sample, fastqc_flag):
    fastqc(sample, fastqc_flag)
    return

@parallel(inputList_star)
@check_if_uptodate(check_file_exists)
def run_star(sample, star_flag):
    star(sample, star_flag)
    return

@follows(run_star)
@parallel(inputList_htseq)
@check_if_uptodate(check_file_exists)
def run_htseq(sample, htseq_flag):
    htseq(sample, htseq_flag)
    return

@follows(run_fastqc, run_htseq)
@parallel([["report", "%s/RNAseq_report_report_completed.flag" % (p.FLAG_PATH)]])
@check_if_uptodate(check_file_exists)
def run_RNAseq_report_counts(sample, RNAseq_report_counts_flag):
    RNAseq_report_counts(sample, RNAseq_report_counts_flag)
    return

@follows(run_RNAseq_report_counts)
@parallel([["combined", "%s/last_function_combined_completed.flag" % (p.FLAG_PATH)]])
@check_if_uptodate(check_file_exists)
def last_function(sample, last_function_flag):
    print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
    pipeline_graph_output = p.FLAG_PATH + "/pipeline_" + sample + "_" + str(date) + ".pdf"
    pipeline_printout_graph (pipeline_graph_output,'pdf', step, no_key_legend=False)
    stage = "last_function"
    flag_file = "%s/%s_completed.flag" % (p.FLAG_PATH, stage)
    open(flag_file, 'w').close()
    return   


if __name__ == '__main__':

    pipeline_run(p.STEP, multiprocess = p.PIPE_MULTIPROCESS, verbose = p.PIPE_VERBOSE, gnu_make_maximal_rebuild_mode = p.PIPE_REBUILD)


