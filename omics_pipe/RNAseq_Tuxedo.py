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
from omics_pipe.modules.fastqc import fastqc
from omics_pipe.modules.tophat import tophat
from omics_pipe.modules.cufflinks import cufflinks
from omics_pipe.modules.rseqc import rseqc
from omics_pipe.modules.cuffmerge import cuffmerge
from omics_pipe.modules.cuffmergetocompare import cuffmergetocompare
from omics_pipe.modules.cuffdiff import cuffdiff
from omics_pipe.modules.RNAseq_report_tuxedo import RNAseq_report_tuxedo


from omics_pipe.parameters.default_parameters import default_parameters 
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

for step in p.STEPS_DE:
    vars()['inputList_' + step] = []
    vars()['inputList_' + step].append([step, "%s/%s_completed.flag" % (p.FLAG_PATH, step)])
    print vars()['inputList_' + step]
     

@parallel(inputList_fastqc)
@check_if_uptodate(check_file_exists)
def run_fastqc(sample, fastqc_flag):
    fastqc(sample, fastqc_flag)
    return

@parallel(inputList_tophat)
@check_if_uptodate(check_file_exists)
@follows(run_fastqc)   
def run_tophat(sample, tophat_flag):
    tophat(sample, tophat_flag)
    return

@parallel(inputList_cufflinks)
@check_if_uptodate(check_file_exists)
@follows(run_tophat)
def run_cufflinks(sample, cufflinks_flag):
    cufflinks(sample, cufflinks_flag)
    return


@parallel([["cuffmerge", "%s/cuffmerge_cuffmerge_completed.flag" % (p.FLAG_PATH)]])
@check_if_uptodate(check_file_exists)
@follows(run_cufflinks)
def run_cuffmerge(step, cuffmerge_flag):
    cuffmerge(step, cuffmerge_flag)
    return

@parallel([["cuffmergetocompare", "%s/cuffmergetocompare_cuffmergetocompare_completed.flag" % (p.FLAG_PATH)]])
@check_if_uptodate(check_file_exists)
@follows(run_cuffmerge)
def run_cuffmergetocompare(step, cuffmergetocompare_flag):
    cuffmergetocompare(step, cuffmergetocompare_flag)
    return

@check_if_uptodate(check_file_exists)
@parallel([["cuffdiff", "%s/cuffdiff_cuffdiff_completed.flag" % (p.FLAG_PATH)]])
@follows(run_cuffmergetocompare)
def run_cuffdiff(step, cuffdiff_flag):
    cuffdiff(step, cuffdiff_flag)
    return

@check_if_uptodate(check_file_exists)
@parallel([["report", "%s/RNAseq_report_tuxedo_report_completed.flag" % (p.FLAG_PATH)]])
@follows(run_cuffdiff)
def run_RNAseq_report_tuxedo(step, RNAseq_report_tuxedo_flag):
    RNAseq_report_tuxedo(step, RNAseq_report_tuxedo_flag)
    return

@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_RNAseq_report_tuxedo)
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

