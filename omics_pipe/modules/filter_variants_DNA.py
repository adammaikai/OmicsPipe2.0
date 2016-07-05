#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def filter_variants_DNA(sample, filter_variants_DNA_flag):   
    '''Filters variants to remove common variants.
    
    input: 
        .bam or .sam file
    output: 
        .vcf file 
    citation: 
        
    link:
             
    parameters from parameters file: 

    VARIANT_RESULTS:
    
    TEMP_DIR:
    
    BEDTOOLS_VERSION:
    
    SNPEFF_VERSION:
    
    dbNSFP:
    
    VCFTOOLS_VERSION:
    
    WORKING_DIR:
    
    SNP_FILTER_OUT_REF:
    
    R_VERSION:
        '''
    spawn_job(jobname = 'filter_variants_DNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/filter_variants.sh", args_list = [p.VARIANT_RESULTS,sample,p.TEMP_DIR,p.BEDTOOLS_VERSION,p.SNPEFF_VERSION,p.dbNSFP,p.VCFTOOLS_VERSION,p.WORKING_DIR,p.SNP_FILTER_OUT_REF,p.R_VERSION])
    job_status(jobname = 'filter_variants_DNA', resultspath = p.VARIANT_RESULTS, SAMPLE = sample, outputfilename = sample + "/intogen_input.vcf", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    filter_variants_DNA(sample, filter_variants_DNA_flag)
    sys.exit(0)
    
