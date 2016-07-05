#!/usr/bin/env python
# -*- coding: utf-8 -*-
from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def fusion_catcher(sample, fusion_catcher_flag):
    '''Detects fusion genes in paired-end RNAseq data. 
    
    input: 
        paired end .fastq files
    output: 
        list of candidate fusion genes
    citation: 
        S. Kangaspeska, S. Hultsch, H. Edgren, D. Nicorici, A. Murumgi, O.P. Kallioniemi, Reanalysis of RNA-sequencing data reveals several additional fusion genes with multiple isoforms, PLOS One, Oct. 2012. http://dx.plos.org/10.1371/journal.pone.0048745
    link: 
        https://code.google.com/p/fusioncatcher
    parameters from parameters file: 
        ENDS:
        
        RAW_DATA_DIR:
        
        FUSION_RESULTS:
        
        FUSIONCATCHERBUILD_DIR:
        
        TEMP_DIR:
        
        SAMTOOLS_VERSION:
        
        FUSIONCATCHER_VERSION:
        
        FUSIONCATCHER_OPTIONS:
        
        TISSUE:
        
        PYTHON_VERSION:
        '''
    if p.ENDS == "PE":
        spawn_job(jobname = 'fusion_catcher', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/fusion_catcher_drmaa.sh", args_list = [p.RAW_DATA_DIR, p.FUSION_RESULTS, p.FUSIONCATCHERBUILD_DIR, p.TEMP_DIR, p.SAMTOOLS_VERSION, p.FUSIONCATCHER_VERSION, p.FUSIONCATCHER_OPTIONS, sample, p.TISSUE, p.PYTHON_VERSION])
        job_status(jobname = 'fusion_catcher', resultspath = p.FUSION_RESULTS, SAMPLE = sample, outputfilename = sample + "/final-list_candidate-fusion-genes.txt", FLAG_PATH = p.FLAG_PATH)
    else:
        print "Fusion Catcher can only be run on paired-end samples. Input samples are single-ended. Skipping step."
        flag_file = "%s/fusion_catcher_completed.flag" % p.FLAG_PATH
        open(flag_file, 'w').close() 
    return

if __name__ == '__main__':
    fusion_catcher(sample, fusion_catcher_flag)
    sys.exit(0)