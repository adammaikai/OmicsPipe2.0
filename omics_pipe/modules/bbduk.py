#!/usr/bin/env python
# -*- coding: utf-8 -*-
from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
import os.path
p = Bunch(default_parameters)


def bbduk(sample, extension, bbduk_flag):
    ''' 
    
    input: 
        fastq
    output: 
        
    citation: 
        
    link: 
        
    parameters from parameters file: 

        '''
    if extension is not None:
        if "RNA" in extension:
            sample_type="rna"
        else:
            sample_type="dna"
        sample_full = sample + extension
        if os.path.isfile(p.BBDUK["RAW_DATA"] + "/" + sample +  "/" + sample_type + "/" + sample_full + "_1.fastq.gz"):
            RAW_DATA_DIR = p.BBDUK["RAW_DATA"] + "/" + sample +  "/" + sample_type
        else:
            RAW_DATA_DIR = p.BBDUK["RAW_DATA"]
    else:
        RAW_DATA_DIR = p.BBDUK["RAW_DATA"]
        sample_full = sample

    spawn_job(jobname = 'bbduk', SAMPLE = sample_full, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.BBDUK["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.BBDUK["NODES"], ppn = p.BBDUK["CPU"], memory = p.BBDUK["MEMORY"], script = "/bbduk.sh", args_list = [sample_full, RAW_DATA_DIR, p.BBDUK["VERSION"], p.BBDUK["RESULTS"], p.BBDUK["REF_ADAPTER"], p.BBDUK["OPTIONS"], p.BBDUK["MEMORY"], str(p.BBDUK["CPU"])])
    job_status(jobname = 'bbduk', resultspath = p.BBDUK["RESULTS"] + "/" + sample_full + "/stats/", SAMPLE = sample_full, outputfilename = sample_full + ".stats", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    bbduk(sample, extension, bbduk_flag)
    sys.exit(0)
