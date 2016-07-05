#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def RNAseq_QC(sample, RNAseq_QC_flag):   
    '''Runs rseqc to determine insert size as QC for alignment.
    
    input: 
        .bam
    output: 
        pdf plot
    link: 
        http://rseqc.sourceforge.net/
    parameters from parameters file: 
        STAR_RESULTS:
        
        QC_PATH:
        
        BAM_FILE_NAME:
        
        RSEQC_REF:
              
        TEMP_DIR:
        
        PICARD_VERSION:
        
        R_VERSION:
        '''

    spawn_job(jobname = 'RNAseq_QC', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 1, memory = "12gb", script = "/RNAseq_QC.sh", args_list = [p.STAR_RESULTS, p.QC_PATH, p.BAM_FILE_NAME, p.RSEQC_REF, p.TEMP_DIR, sample, p.PICARD_VERSION, p.R_VERSION])
    job_status(jobname = 'RNAseq_QC', resultspath = p.QC_PATH, SAMPLE = sample, outputfilename = sample + "/insertSizeHist.pdf", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    RNAseq_QC(sample, RNAseq_QC_flag)
    sys.exit(0)
