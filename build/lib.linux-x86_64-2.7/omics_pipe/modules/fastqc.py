#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def fastqc(sample, extension, fastqc_flag):
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
            
            COMPRESSION:
            ''' 
    print "sample name is: ", sample

    sample = sample + extension    
    spawn_job(jobname = 'fastqc', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.FASTQC["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.FASTQC["NODES"], ppn = p.FASTQC["CPU"], memory = p.FASTQC["MEMORY"], script = "/fastqc_drmaa.sh", args_list = [sample, p.FASTQC["PATH"],p.FASTQC['RESULTS'], p.FASTQC["VERSION"], p.FASTQC["COMPRESSION"]])
    job_status(jobname = 'fastqc', resultspath = p.FASTQC["RESULTS"] + "/" + sample, SAMPLE = sample, outputfilename = sample + "_2" + "_fastqc.html", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    fastqc(sample, extension, fastqc_flag)
    sys.exit(0)
