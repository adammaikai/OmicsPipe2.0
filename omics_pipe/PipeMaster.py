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
from multiprocessing import Pool
import itertools
from omics_pipe.utils import *
from omics_pipe.parameters.default_parameters import default_parameters 
from omics_pipe.modules.fastqc import fastqc
from omics_pipe.modules.Star import Star
from omics_pipe.modules.mojo import mojo
from omics_pipe.modules.WGS_preprocess import WGS_preprocess
from omics_pipe.modules.Samtools_pileup import Samtools_pileup
from omics_pipe.modules.Varscan_somatic import Varscan_somatic
from omics_pipe.modules.GATK_BQSR import GATK_BQSR
from omics_pipe.modules.MuTect import MuTect
from omics_pipe.modules.bbduk import bbduk
from omics_pipe.modules.snpir_variants import snpir_variants_plus_bwa
from omics_pipe.modules.GATK_variant_calling import GATK_variant_calling
from omics_pipe.modules.Varscan import Varscan
from omics_pipe.modules.GATK_Varscan_Annotation import GATK_Varscan_Annotation
from omics_pipe.modules.Somatic_Annotation import Somatic_Annotation
from omics_pipe.modules.Qualimap import Qualimap
p = Bunch(default_parameters)
os.chdir(p.OMICSPIPE["WORKING_DIR"])
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d %H:%M")    


ext_list = []
dna_ext_list = []
if p.DNA["TUMOR_EXT"]:
    ext_list.extend(p.DNA["TUMOR_EXT"])
    dna_ext_list.extend(p.DNA["TUMOR_EXT"])
if p.DNA["NORMAL_EXT"]:
    ext_list.extend(p.DNA["NORMAL_EXT"])
    dna_ext_list.extend(p.DNA["NORMAL_EXT"])
if p.RNA["TUMOR_EXT"]:
    ext_list.extend(p.RNA["TUMOR_EXT"])
if len(ext_list) == 0:
    ext_list = [None]
print p

for step in p.STEPS:
    vars()['inputList_' + step] = []
    for sample in p.SAMPLE_LIST:
        vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.OMICSPIPE["FLAG_PATH"], step, sample)])
    print vars()['inputList_' + step]


#Bbduk
@parallel(inputList_bbduk)
@check_if_uptodate(check_file_exists)
def run_bbduk(sample, bbduk_flag):
    # for ext in ext_list:
    #     bbduk(sample, ext, bbduk_flag)
    return 

#FASTQC
@parallel(inputList_fastqc)
@check_if_uptodate(check_file_exists)
@follows(run_bbduk)
def run_fastqc(sample, fastqc_flag):
    # for ext in ext_list:
    #     fastqc(sample, ext, fastqc_flag)
    return

#Star
@parallel(inputList_Star)
@check_if_uptodate(check_file_exists)
@follows(run_bbduk)
def run_Star(sample, Star_flag):
    # for ext in p.RNA["TUMOR_EXT"]:
    #     Star(sample, ext, Star_flag)
    return

#mojo
@parallel(inputList_mojo)
@check_if_uptodate(check_file_exists)
@follows(run_Star)
def run_mojo(sample, mojo_flag):
    # for ext in p.RNA["TUMOR_EXT"]:
    #     mojo(sample, ext, mojo_flag)
    return

#Preprocess (bwa-mem, sort, markdup, realign, sort, index)
@parallel(inputList_WGS_preprocess)
@check_if_uptodate(check_file_exists)
@follows(run_bbduk)
def run_WGS_preprocess(sample, WGS_preprocess_flag):
    # for ext in dna_ext_list:
    #     WGS_preprocess(sample, ext, WGS_preprocess_flag)
    return

#Qualimap
@parallel(inputList_Qualimap)
@check_if_uptodate(check_file_exists)
@follows(run_WGS_preprocess, run_Star)
def run_Qualimap(sample, Qualimap_flag):
    # for ext in ext_list:
    #     Qualimap(sample, ext, Qualimap_flag)
    return

# #GATK_BQSR
# @parallel(inputList_GATK_BQSR)
# @check_if_uptodate(check_file_exists)
# @follows(run_WGS_preprocess)
# def run_GATK_BQSR(sample, GATK_BQSR_flag):
#     for ext in dna_ext_list:
#         GATK_BQSR(sample, ext, GATK_BQSR_flag)
#     return

#MuTect
@parallel(inputList_MuTect)
@check_if_uptodate(check_file_exists)
@follows(run_WGS_preprocess)
def run_MuTect(sample, MuTect_flag):
    # for ext in p.DNA["TUMOR_EXT"]:
    #     MuTect(sample, ext, p.DNA["NORMAL_EXT"][0], MuTect_flag)
    return

#Samtools_pileup
@parallel(inputList_Samtools_pileup)
@check_if_uptodate(check_file_exists)
@follows(run_WGS_preprocess)
def run_Samtools_pileup(sample, Samtools_pileup_flag):
    # for ext in dna_ext_list:
    #     Samtools_pileup(sample, ext, Samtools_pileup_flag)
    return

#Varscan_somatic
@parallel(inputList_Varscan_somatic)
@check_if_uptodate(check_file_exists)
@follows(run_Samtools_pileup)
def run_Varscan_somatic(sample, Varscan_somatic_flag):
    # for ext in p.DNA["TUMOR_EXT"]:
    #     Varscan_somatic(sample, ext, p.DNA["NORMAL_EXT"][0], Varscan_somatic_flag)
    return

#GATK_variant_calling
@parallel(inputList_GATK_variant_calling)
@check_if_uptodate(check_file_exists)
@follows(run_WGS_preprocess)
def run_GATK_variant_calling(sample, GATK_variant_calling_flag):
    # for ext in dna_ext_list:
    #     GATK_variant_calling(sample, ext, GATK_variant_calling_flag)
    return

#Varscan
@parallel(inputList_Varscan)
@check_if_uptodate(check_file_exists)
@follows(run_Samtools_pileup)
def run_Varscan(sample, Varscan_flag):
    # for ext in dna_ext_list:
    #     Varscan(sample, ext, Varscan_flag)
    return

#Annotation
@parallel(inputList_GATK_Varscan_Annotation)
@check_if_uptodate(check_file_exists)
@follows(run_GATK_variant_calling, run_Varscan)
def run_GATK_Varscan_Annotation(sample, GATK_Varscan_Annotation_flag):
    # for ext in dna_ext_list:
    #     GATK_Varscan_Annotation(sample, ext, GATK_Varscan_Annotation_flag)
    return

#Somatic Annotation
@parallel(inputList_Somatic_Annotation)
@check_if_uptodate(check_file_exists)
@follows(run_MuTect, run_Varscan_somatic)
def run_Somatic_Annotation(sample, Somatic_Annotation_flag):
    # for ext in p.DNA["TUMOR_EXT"]:
    #     Somatic_Annotation(sample, ext, p.DNA["NORMAL_EXT"][0], Somatic_Annotation_flag)
    return

#snpir rna-seq variant calling
@parallel(inputList_snpir_variants_plus_bwa)
@check_if_uptodate(check_file_exists)
@follows(run_Star, run_Somatic_Annotation)
def run_snpir_variants_plus_bwa(sample, snpir_variants_plus_bwa_flag):
    # for ext in p.RNA["TUMOR_EXT"]:
    #     snpir_variants_plus_bwa(sample, ext, snpir_variants_plus_bwa_flag)
    return
    
@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_bbduk, run_fastqc, run_Qualimap, run_mojo, run_snpir_variants_plus_bwa, run_GATK_Varscan_Annotation, run_Somatic_Annotation)
def last_function(sample, last_function_flag):
    print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
    pipeline_graph_output = p.OMICSPIPE["FLAG_PATH"] + "/pipeline_" + sample + "_" + str(date) + ".pdf"
    #pipeline_printout_graph (pipeline_graph_output,'pdf', step, no_key_legend=False)
    stage = "last_function"
    flag_file = "%s/%s_%s_completed.flag" % (p.OMICSPIPE["FLAG_PATH"], stage, sample)
    open(flag_file, 'w').close()
    return   


if __name__ == '__main__':

    pipeline_run(p.STEP, multiprocess = p.OMICSPIPE["PIPE_MULTIPROCESS"], verbose = p.OMICSPIPE["PIPE_VERBOSE"], gnu_make_maximal_rebuild_mode = p.OMICSPIPE["PIPE_REBUILD"])


