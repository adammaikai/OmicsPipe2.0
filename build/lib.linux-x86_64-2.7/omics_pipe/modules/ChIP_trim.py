#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def ChIP_trim(sample, ChIP_trim_flag):
    '''Runs Homer Tools to trim adapters from .fastq files.
    
    input: 
        .fastq file
    output: 
        .fastq file
    citation: 
        Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
    link: 
        http://homer.salk.edu/homer/
    parameters from parameters file: 
        ENDS: 
        
        RAW_DATA_DIR:
        
        HOMER_TRIM_OPTIONS:
        
        TRIMMED_DATA_PATH:
        
        HOMER_VERSION:
        '''
    spawn_job(jobname = 'ChIP_trim', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/ChIP_trim_" + p.ENDS + ".sh", args_list = [sample,p.RAW_DATA_DIR,p.HOMER_TRIM_OPTIONS, p.TRIMMED_DATA_PATH, p.HOMER_VERSION])
    if p.ENDS == "PE":
        job_status(jobname = 'ChIP_trim', resultspath = p.TRIMMED_DATA_PATH, SAMPLE = sample, outputfilename = sample + "_1.fastq", FLAG_PATH = p.FLAG_PATH)
    else:
        job_status(jobname = 'ChIP_trim', resultspath = p.TRIMMED_DATA_PATH, SAMPLE = sample, outputfilename = sample + ".fastq", FLAG_PATH = p.FLAG_PATH)

    return

if __name__ == '__main__':
    ChIP_trim(sample, ChIP_trim_flag)
    sys.exit(0)