#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def cutadapt_miRNA(sample, cutadapt_miRNA_flag):
    '''Runs Cutadapt to trim adapters from reads.
    
    input: 
        .fastq
    output: 
        .fastq
    citation: 
        Martin 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17: 10-12.
    link: 
        https://code.google.com/p/cutadapt/
    parameters from parameters file:       
        RAW_DATA_DIR:
        
        ADAPTER:
        
        TRIMMED_DATA_PATH:
        
        PYTHON_VERSION'''
    spawn_job(jobname = 'cutadapt_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/cutadapt_drmaa.sh", args_list = [sample, p.RAW_DATA_DIR, p.ADAPTER, p.TRIMMED_DATA_PATH, p.PYTHON_VERSION])
    job_status(jobname = 'cutadapt_miRNA', resultspath = p.TRIMMED_DATA_PATH, SAMPLE = sample, outputfilename = sample + "_trimmed.fastq", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    cutadapt_miRNA(sample, cutadapt_miRNA_flag)
    sys.exit(0)