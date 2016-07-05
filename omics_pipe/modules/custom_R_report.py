#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def custom_R_report(sample, custom_R_report_flag):   
    '''Runs R script with knitr to produce report from omics pipeline.
     
    input: 
        results from other steps in RNAseq pipelines
    output: 
        html report
    citation: 
        T. Meissner
    parameters from parameter file:
        REPORT_SCRIPT:
        
        R_VERSION:
        
        REPORT_RESULTS:
        
        R_MARKUP_FILE:
        
        DPS_VERSION:
        
        PARAMS_FILE:
        '''
    spawn_job(jobname = 'custom_R_report', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "31gb", script = "/custom_R_report.sh", args_list = [sample,p.REPORT_SCRIPT,p.R_VERSION, p.REPORT_RESULTS, p.R_MARKUP_FILE,  p.DPS_VERSION, p.PARAMS_FILE])
    job_status(jobname = 'custom_R_report', resultspath = p.REPORT_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + ".html", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    custom_R_report(sample, call_variants_flag)
    sys.exit(0)