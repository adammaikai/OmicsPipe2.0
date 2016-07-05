#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def DNAseq_variant_report(sample, DNAseq_variant_report_flag):   
    '''Runs R script with knitr to produce report from DNAseq pipeline.
     
    input: 
        results from other steps in DNAseq pipelines
    output: 
        html report
    citation: 
        K. Fisch, T. Meissner
    parameters from parameter file:
        WORKING_DIR:
        
        R_VERSION:
        
        REPORT_RESULTS:
        
        PARAMS_FILE:
        
        TABIX_VERSION:
        
        TUMOR_TYPE: 
        
        GENELIST: 
        
        COSMIC: 
        
        CLINVAR: 
        
        PHARMGKB_rsID: 
        
        PHARMGKB_Allele: 
        
        DRUGBANK: 
        
        CADD: 
        '''
    spawn_job(jobname = 'DNAseq_variant_report', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/DNAseq_report.sh", args_list = [sample,p.WORKING_DIR,p.R_VERSION, p.REPORT_RESULTS, p.PARAMS_FILE, p.TABIX_VERSION])
    job_status(jobname = 'DNAseq_variant_report', resultspath = p.REPORT_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + ".html", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    DNAseq_variant_report(sample, DNAseq_variant_report_flag)
    sys.exit(0)
    
