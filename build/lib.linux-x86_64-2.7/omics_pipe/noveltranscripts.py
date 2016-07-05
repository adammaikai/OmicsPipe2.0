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
from omics_pipe.modules.tophat_ncRNA import tophat_ncRNA
from omics_pipe.modules.cufflinks_ncRNA import cufflinks_ncRNA

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
    



@parallel(inputList_tophat_ncRNA)
@check_if_uptodate(check_file_exists)
def run_tophat_ncRNA(sample, tophat_ncRNA_flag):
    tophat_ncRNA(sample, tophat_ncRNA_flag)
    return

@parallel(inputList_cufflinks_ncRNA)
@check_if_uptodate(check_file_exists)
@follows(run_tophat_ncRNA)
def run_cufflinks_ncRNA(sample, cufflinks_ncRNA_flag):
    cufflinks_ncRNA(sample, cufflinks_ncRNA_flag)
    return



@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_cufflinks_ncRNA)
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
