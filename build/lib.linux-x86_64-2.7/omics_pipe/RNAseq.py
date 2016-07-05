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
from omics_pipe.modules.star import star
from omics_pipe.modules.cufflinks import cufflinks
from omics_pipe.modules.htseq import htseq
from omics_pipe.modules.fusion_catcher import fusion_catcher
from omics_pipe.modules.call_variants import call_variants
from omics_pipe.modules.annotate_variants import annotate_variants
from omics_pipe.modules.RNAseq_report import RNAseq_report
from omics_pipe.parameters.default_parameters import default_parameters 
from omics_pipe.modules.bwa import bwa1
from omics_pipe.modules.bwa import bwa2
from omics_pipe.modules.bwa import bwa
from omics_pipe.modules.snpir_variants import snpir_variants
from omics_pipe.modules.filter_variants import filter_variants
from omics_pipe.modules.rseqc import rseqc
from omics_pipe.modules.intogen import intogen
from omics_pipe.modules.htseq_gencode import htseq_gencode

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

@parallel(inputList_star)
@check_if_uptodate(check_file_exists)
@follows(run_fastqc)
def run_star(sample, star_flag):
    star(sample, star_flag)
    return

@parallel(inputList_cufflinks)
@check_if_uptodate(check_file_exists)
@follows(run_tophat)
def run_cufflinks(sample, cufflinks_flag):
    cufflinks(sample, cufflinks_flag)
    return

@parallel(inputList_htseq)
@check_if_uptodate(check_file_exists)
@follows(run_star)
def run_htseq(sample, htseq_flag):
    htseq(sample, htseq_flag)
    return

@parallel(inputList_fusion_catcher)
@check_if_uptodate(check_file_exists)
def run_fusion_catcher(sample, fusion_catcher_flag):
    fusion_catcher(sample, fusion_catcher_flag)
    return


@parallel(inputList_bwa)
@check_if_uptodate(check_file_exists)
def run_bwa(sample, bwa_flag):
    bwa(sample, bwa_flag)
    return

@parallel(inputList_htseq_gencode)
@check_if_uptodate(check_file_exists)
@follows(run_star)
def run_htseq_gencode(sample, htseq_gencode_flag):
    htseq_gencode(sample, htseq_gencode_flag)
    return

@parallel(inputList_snpir_variants)
@check_if_uptodate(check_file_exists)
@follows(run_bwa)
def run_snpir_variants(sample, snpir_variants_flag):
    snpir_variants(sample, snpir_variants_flag)
    return


@parallel(inputList_filter_variants)
@check_if_uptodate(check_file_exists)
@follows(run_snpir_variants)
def run_filter_variants(sample, filter_variants_flag):
    filter_variants(sample, filter_variants_flag)
    return

@parallel(inputList_intogen)
@check_if_uptodate(check_file_exists)
@follows(run_filter_variants)
def run_intogen(sample, intogen_flag):
    intogen(sample, intogen_flag)
    return


@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_intogen, run_htseq_gencode)
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

