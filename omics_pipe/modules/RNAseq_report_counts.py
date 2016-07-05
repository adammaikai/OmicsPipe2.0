#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def RNAseq_report_counts(sample, RNAseq_report_counts_flag):   
    '''Runs R script with knitr to produce report from RNAseq pipeline.
     
    input: 
        results from other steps in RNAseq pipelines
    output: 
        html report
    citation: 
        T. Meissner
        Love MI, Huber W and Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, pp. 550. http://dx.doi.org/10.1186/s13059-014-0550-8
    parameters from parameter file:
        WORKING_DIR:
        
        R_VERSION:
        
        REPORT_RESULTS:
        
        PARAMS_FILE:
        '''
    spawn_job(jobname = 'RNAseq_report_counts', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "2gb", script = "/RNAseq_report_counts.sh", args_list = [sample,p.WORKING_DIR,p.R_VERSION, p.REPORT_RESULTS, p.DESEQ_RESULTS, p.PARAMS_FILE])
    job_status(jobname = 'RNAseq_report_counts', resultspath = p.REPORT_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + ".html", FLAG_PATH = p.FLAG_PATH)
    return
    
if __name__ == '__main__':
    RNAseq_report_counts(sample, RNAseq_report_counts_flag)
    sys.exit(0)