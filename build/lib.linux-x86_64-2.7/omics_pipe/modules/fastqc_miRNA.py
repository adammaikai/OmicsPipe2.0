#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def fastqc_miRNA(sample, fastqc_miRNA_flag):
    '''QC check of raw .fastq files using FASTQC.
    
        input: 
            .fastq file
        output: 
            folder and zipped folder containing html, txt and image files
        citation: 
            Babraham Bioinformatics
        link: 
            http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
        parameters from parameters file: 
            RAW_DATA_DIR:
            
            QC_PATH:
            
            FASTQC_VERSION:
            ''' 
    print "sample name is: ", sample 
    if p.ENDS == "PE":
        SAMPLE1 = sample + "_1"
        SAMPLE2 = sample + "_2"        
        spawn_job(jobname = 'fastqc_miRNA', SAMPLE = SAMPLE1, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/fastqc_drmaa.sh", args_list = [SAMPLE1, p.RAW_DATA_DIR,p.QC_PATH, p.FASTQC_VERSION])
        spawn_job(jobname = 'fastqc_miRNA', SAMPLE = SAMPLE2, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/fastqc_drmaa.sh", args_list = [SAMPLE2, p.RAW_DATA_DIR,p.QC_PATH, p.FASTQC_VERSION])
        job_status(jobname = 'fastqc_miRNA', resultspath = p.QC_PATH, SAMPLE = SAMPLE1,  outputfilename = SAMPLE1 + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
        job_status(jobname = 'fastqc_miRNA', resultspath = p.QC_PATH, SAMPLE = SAMPLE2,  outputfilename = SAMPLE2 + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
    else:
        spawn_job(jobname = 'fastqc_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/fastqc_drmaa.sh", args_list = [sample, p.TRIMMED_DATA_PATH,p.QC_PATH, p.FASTQC_VERSION])
        job_status(jobname = 'fastqc_miRNA', resultspath = p.QC_PATH, SAMPLE = sample, outputfilename = sample + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
    
    return

if __name__ == '__main__':
    fastqc_miRNA(sample, fastqc_miRNA_flag)
    sys.exit(0)