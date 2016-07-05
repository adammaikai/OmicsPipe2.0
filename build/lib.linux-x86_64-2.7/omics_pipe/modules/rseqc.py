#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def rseqc(sample, rseqc_flag):   
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
        
        RSEQC_VERSION:
        
        TEMP_DIR:
        '''
    spawn_job(jobname = 'rseqc', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/rseqc_drmaa.sh", args_list = [p.STAR_RESULTS, p.QC_PATH, p.BAM_FILE_NAME, p.RSEQC_REF, p.RSEQC_VERSION, p.TEMP_DIR, sample])
    job_status(jobname = 'rseqc', resultspath = p.QC_PATH, SAMPLE = sample, outputfilename = sample + "/insert_size.inner_distance_plot.pdf", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    rseqc(sample, rseqc_flag)
    sys.exit(0)