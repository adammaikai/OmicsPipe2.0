#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)



def fastq_length_filter_miRNA(sample, fastq_length_filter_miRNA_flag):
    '''Runs custom Python script to filter miRNA reads by length.
    
    input: 
        .fastq
    output: 
        .fastq
    parameters from parameter file: 
        TRIMMED_DATA_PATH:
    '''
    sample1 = p.TRIMMED_DATA_PATH + "/" + sample + '_trimmed.fastq'
    script = p.WORKING_DIR + "/" + "fastq_length_filter.py"
    spawn_job(jobname = 'fastq_length_filter_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/run_py_script.sh", args_list = [script, sample1,p.TRIMMED_DATA_PATH])
    job_status(jobname = 'fastq_length_filter_miRNA', resultspath = p.TRIMMED_DATA_PATH, SAMPLE = sample, outputfilename = sample + ".fastq", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    fastq_length_filter_miRNA(sample, fastq_length_filter_miRNA_flag)
    sys.exit(0)