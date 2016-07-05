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
from omics_pipe.modules.star import star
from omics_pipe.modules.htseq import htseq
from omics_pipe.modules.htseq_gencode import htseq_gencode
from omics_pipe.modules.fusion_catcher import fusion_catcher
from omics_pipe.modules.call_variants import call_variants
from omics_pipe.modules.annotate_variants import annotate_variants
from omics_pipe.modules.BreastCancer_RNA_report import BreastCancer_RNA_report
from omics_pipe.modules.bwa import bwa1
from omics_pipe.modules.bwa import bwa2
from omics_pipe.modules.snpir_variants import snpir_variants
from omics_pipe.modules.filter_variants import filter_variants
from omics_pipe.modules.RNAseq_QC import RNAseq_QC
from omics_pipe.modules.intogen import intogen
from omics_pipe.modules.TCGA_download import TCGA_download
from omics_pipe.modules.oncofuse import oncofuse

from omics_pipe.parameters.default_parameters import default_parameters 
p = Bunch(default_parameters)


os.chdir(p.WORKING_DIR)
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    

print p


# sample_list = get_TCGA_ID(p.TCGA_XML_FILE)
# 
# for step in p.STEPS:
#     vars()['inputList_' + step] = []
#     for sample in sample_list:
#         vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
#     print vars()['inputList_' + step]

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
    

@parallel(inputList_fastqc)
@check_if_uptodate(check_file_exists)
@follows(run_TCGA_download)
def run_fastqc(sample, fastqc_flag):
    fastqc(sample, fastqc_flag)
    return

@parallel(inputList_star)
@check_if_uptodate(check_file_exists)
@follows(run_TCGA_download)
def run_star(sample, star_flag):
    star(sample, star_flag)
    return

@parallel(inputList_RNAseq_QC)
@check_if_uptodate(check_file_exists)
@follows(run_star)
def run_RNAseq_QC(sample, RNAseq_QC_flag):
    RNAseq_QC(sample, RNAseq_QC_flag)
    return

@parallel(inputList_fusion_catcher)
@check_if_uptodate(check_file_exists)
@follows(run_TCGA_download)
def run_fusion_catcher(sample, fusion_catcher_flag):
    fusion_catcher(sample, fusion_catcher_flag)
    return

@parallel(inputList_oncofuse)
@check_if_uptodate(check_file_exists)
@follows(run_fusion_catcher)
def run_oncofuse(sample, oncofuse_flag):
    oncofuse(sample, oncofuse_flag)
    return


@parallel(inputList_bwa1)
@check_if_uptodate(check_file_exists)
@follows(run_TCGA_download)
def run_bwa1(sample, bwa1_flag):
    bwa1(sample, bwa1_flag)
    return

@parallel(inputList_bwa2)
@check_if_uptodate(check_file_exists)
@follows(run_TCGA_download)
def run_bwa2(sample, bwa2_flag):
    bwa2(sample, bwa2_flag)
    return

@parallel(inputList_htseq)
@check_if_uptodate(check_file_exists)
@follows(run_star)
def run_htseq(sample, htseq_flag):
    htseq(sample, htseq_flag)
    return

@parallel(inputList_htseq_gencode)
@check_if_uptodate(check_file_exists)
@follows(run_star)
def run_htseq_gencode(sample, htseq_gencode_flag):
    htseq_gencode(sample, htseq_gencode_flag)
    return

@parallel(inputList_snpir_variants)
@check_if_uptodate(check_file_exists)
@follows(run_bwa1, run_bwa2)
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

@parallel(inputList_BreastCancer_RNA_report)
@check_if_uptodate(check_file_exists)
@follows(run_fastqc, run_RNAseq_QC, run_oncofuse, run_htseq, run_intogen)
def run_BreastCancer_RNA_report(sample, BreastCancer_RNA_report_flag):
    BreastCancer_RNA_report(sample,BreastCancer_RNA_report_flag)
    return


@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_BreastCancer_RNA_report, run_htseq_gencode)
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


