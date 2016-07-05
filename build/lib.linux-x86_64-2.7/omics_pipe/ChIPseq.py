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
from omics_pipe.modules.fastqc import fastqc
from omics_pipe.modules.ChIP_trim import ChIP_trim
from omics_pipe.modules.bowtie import bowtie
from omics_pipe.modules.read_density import read_density
from omics_pipe.modules.macs import macs
from omics_pipe.modules.homer_peaks import homer_peaks
from omics_pipe.modules.peak_track import peak_track
from omics_pipe.modules.annotate_peaks import annotate_peaks
from omics_pipe.modules.find_motifs import find_motifs

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


for steps in p.STEPS_PAIRS:
    vars()['inputList_' + steps] = []
    vars()['inputList_' + steps].append([steps, "%s/%s_completed.flag" % (p.FLAG_PATH, steps)])
    print vars()['inputList_' + steps]
    

@parallel(inputList_fastqc)
@check_if_uptodate(check_file_exists)
def run_fastqc(sample, fastqc_flag):
    fastqc(sample, fastqc_flag)
    return

@parallel(inputList_ChIP_trim)
@check_if_uptodate(check_file_exists)
def run_ChIP_trim(sample, ChIP_trim_flag):
    ChIP_trim(sample, ChIP_trim_flag)
    return


@parallel(inputList_bowtie)
@check_if_uptodate(check_file_exists)
@follows(run_ChIP_trim)
def run_bowtie(sample, bowtie_flag):
    bowtie(sample, bowtie_flag)
    return


@parallel(inputList_read_density)
@check_if_uptodate(check_file_exists)
@follows(run_bowtie)
def run_read_density(sample, read_density_flag):
    read_density(sample, read_density_flag)
    return

@parallel(inputList_homer_peaks)
@check_if_uptodate(check_file_exists)
@follows(run_read_density)
def run_homer_peaks(step, homer_peaks_flag):
    homer_peaks(step, homer_peaks_flag)
    return

@parallel(inputList_peak_track)
@check_if_uptodate(check_file_exists)
@follows(run_homer_peaks)
def run_peak_track(step, peak_track_flag):
    peak_track(step, peak_track_flag)
    return

@parallel(inputList_annotate_peaks)
@check_if_uptodate(check_file_exists)
@follows(run_homer_peaks)
def run_annotate_peaks(step, annotate_peaks_flag):
    annotate_peaks(step, annotate_peaks_flag)
    return

@parallel(inputList_find_motifs)
@check_if_uptodate(check_file_exists)
@follows(run_homer_peaks)
def run_find_motifs(step, find_motifs_flag):
    find_motifs(step, find_motifs_flag)
    return

@parallel(inputList_macs)
@check_if_uptodate(check_file_exists)
@follows(run_bowtie)
def run_macs(step, macs_flag):
    macs(step, macs_flag)
    return



@parallel(inputList_last_function)
@check_if_uptodate(check_file_exists)
@follows(run_fastqc, run_find_motifs, run_peak_track, run_annotate_peaks, run_macs)
def last_function(sample, last_function_flag):
    print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
    pipeline_graph_output = p.FLAG_PATH + "/pipeline_" + sample + "_" + str(date) + ".pdf"
    #pipeline_printout_graph(pipeline_graph_output,'pdf', [step,steps], no_key_legend=False)
    stage = "last_function"
    flag_file = "%s/%s_%s_completed.flag" % (p.FLAG_PATH, stage, sample)
    open(flag_file, 'w').close()
    return   


if __name__ == '__main__':

    pipeline_run(p.STEP, multiprocess = p.PIPE_MULTIPROCESS, verbose = p.PIPE_VERBOSE, gnu_make_maximal_rebuild_mode = p.PIPE_REBUILD)

