#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def TCGA_download(sample, TCGA_download_flag):  
    '''Downloads and unzips TCGA data from Manifest.xml downloaded from CGHub.
    input:
        TGCA XML file
    output:
        downloaded files from TCGA
    citation: 
        The Cancer Genome Atlas
    link: 
        https://cghub.ucsc.edu/software/downloads.html
    parameters from parameters file: 
        TCGA_XML_FILE:
        
        TCGA_KEY:
        
        TCGA_OUTPUT_PATH:
        
        CGATOOLS_VERSION:
        '''
    spawn_job(jobname = 'TCGA_download', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "31gb", script = "/TCGA_download.sh", args_list = [str(sample),p.TCGA_KEY,p.TCGA_OUTPUT_PATH,p.CGATOOLS_VERSION])
    job_status(jobname = 'TCGA_download', resultspath = p.TCGA_OUTPUT_PATH, SAMPLE = sample, outputfilename = sample + "_1.fastq", FLAG_PATH = p.FLAG_PATH) 
    return

if __name__ == '__main__':
    TCGA_download(sample, TCGA_download_flag)
    sys.exit(0)