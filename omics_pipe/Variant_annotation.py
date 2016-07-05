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
from omics_pipe.modules.intogen import intogen
from omics_pipe.modules.annovar import annovar
from omics_pipe.modules.filter_variants_DNA import filter_variants_DNA
from omics_pipe.modules.DNAseq_variant_report import DNAseq_variant_report

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


#Filter Variants
@parallel(inputList_filter_variants_DNA)
@check_if_uptodate(check_file_exists)
def run_filter_variants_DNA(sample, filter_variants_DNA_flag):
    filter_variants_DNA(sample, filter_variants_DNA_flag)
    return

#Annovar
@parallel(inputList_annovar)
@check_if_uptodate(check_file_exists)
def run_annovar(sample, annovar_flag):
    annovar(sample, annovar_flag)
    return

#Intogen
@parallel(inputList_intogen)
@check_if_uptodate(check_file_exists)
@follows(run_filter_variants_DNA)
def run_intogen(sample, intogen_flag):
    intogen(sample, intogen_flag)
    return


#DNAseq Report
@parallel(inputList_DNAseq_variant_report)
@check_if_uptodate(check_file_exists)
@follows(run_intogen, run_annovar)
def run_DNAseq_variant_report(sample, DNAseq_variant_report_flag):
    DNAseq_variant_report(sample,DNAseq_variant_report_flag)
    return

    
@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_DNAseq_variant_report)
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

