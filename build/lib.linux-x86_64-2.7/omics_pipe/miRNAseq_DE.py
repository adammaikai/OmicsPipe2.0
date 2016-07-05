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
from omics_pipe.modules.cuffmerge_miRNA import cuffmerge_miRNA
from omics_pipe.modules.cuffmergetocompare_miRNA import cuffmergetocompare_miRNA
from omics_pipe.modules.cuffdiff_miRNA import cuffdiff_miRNA


from omics_pipe.parameters.default_parameters import default_parameters 
p = Bunch(default_parameters)

os.chdir(p.WORKING_DIR)
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    

print p

for step in p.STEPS_DE:
    vars()['inputList_' + step] = []
    vars()['inputList_' + step].append([step, "%s/%s_completed.flag" % (p.FLAG_PATH, step)])
    print vars()['inputList_' + step]
    

@parallel(inputList_cuffmerge_miRNA)
@check_if_uptodate(check_file_exists)
def run_cuffmerge_miRNA(step, cuffmerge_miRNA_flag):
    cuffmerge_miRNA(step, cuffmerge_miRNA_flag)
    return

@parallel(inputList_cuffmergetocompare_miRNA)
@check_if_uptodate(check_file_exists)
@follows(run_cuffmerge_miRNA)
def run_cuffmergetocompare_miRNA(step, cuffmergetocompare_miRNA_flag):
    cuffmergetocompare_miRNA(step, cuffmergetocompare_miRNA_flag)
    return

@check_if_uptodate(check_file_exists)
@parallel(inputList_cuffdiff_miRNA)
@follows(run_cuffmergetocompare_miRNA)
def run_cuffdiff_miRNA(step, cuffdiff_miRNA_flag):
    cuffdiff_miRNA(step, cuffdiff_miRNA_flag)
    return


@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_cuffdiff_miRNA)
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

