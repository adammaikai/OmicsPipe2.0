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
    
def cleanup():
    for sample in p.SAMPLE_LIST:
        spawn_job(jobname = 'cleanup', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "2gb", script = "/cleanup.sh", args_list = [sample, p.STAR_RESULTS, p.BWA_RESULTS, p.SNPIR_RESULTS, p.FUSIONCATCHER_RESULTS. p.TEMP_DIR, p.RAW_DATA_DIR, p.TCGA_OUTPUT_PATH])
    return


if __name__ == '__main__':
    cleanup()
