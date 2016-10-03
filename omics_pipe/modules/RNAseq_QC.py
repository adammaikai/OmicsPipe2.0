#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def RNAseq_QC(sample, RNAseq_QC_flag):   
    '''Runs picard rnaseqmetrics and insertsize estimation
    
    input: 
        .bam
    output: 
        pdf plot
    link: 
        
    parameters from parameters file: 
        STAR_RESULTS:
        
        QC_PATH:
        
        BAM_FILE_NAME:
        
        RSEQC_REF:
              
        TEMP_DIR:
        
        PICARD_VERSION:
        
        R_VERSION:
        '''

    spawn_job(jobname = 'RNAseq_QC', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.RNASEQQC["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.RNASEQQC["NODES"], ppn = p.RNASEQQC["CPU"], memory = p.RNASEQQC["MEMORY"], script = "/RNAseq_QC.sh", args_list = [p.RNASEQQC["ALIGNMENT_DIR"], p.RNASEQQC["RESULTS"], p.RNASEQQC["BAM_FILE_NAME"], p.RNASEQQC["REFFLAT"], p.RNASEQQC["TEMP_DIR"], sample, p.RNASEQQC["PICARD_VERSION"], p.RNASEQQC["R_VERSION"]])
    job_status(jobname = 'RNAseq_QC', resultspath = p.RNASEQQC["RESULTS"], SAMPLE = sample, outputfilename = sample + "/insertSizeHist.pdf", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return
if __name__ == '__main__':
    RNAseq_QC(sample, RNAseq_QC_flag)
    sys.exit(0)
