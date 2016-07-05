#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def call_variants(sample, call_variants_flag):    
    '''Calls variants from alignment .bam files using Varcall.
    
    input: 
        Aligned.out.sort.bam or accepted_hits.bam
    output: 
        .vcf file 
    citation: 
        Erik Aronesty (2011). ea-utils : "Command-line tools for processing biological sequencing data";
    link: 
        https://code.google.com/p/ea-utils/wiki/Varcall
    parameters from parameters file: 
        STAR_RESULTS:
        
        GENOME:
        
        VARSCAN_PATH:
        
        VARSCAN_OPTIONS:
        
        VARIANT_RESULTS:
        
        TEMP_DIR:
        
        SAMTOOLS_VERSION:
        
        ANNOVAR_VERSION:
        
        VCFTOOLS_VERSION:
        
        VARSCAN_VERSION:
        
        SAMTOOLS_OPTIONS:
        '''
    spawn_job(jobname = 'call_variants', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "31gb", script = "/varcall_drmaa_bash.sh", args_list = [sample,p.STAR_RESULTS,p.GENOME, p.VARSCAN_PATH, p.VARSCAN_OPTIONS, p.VARIANT_RESULTS, p.TEMP_DIR, p.SAMTOOLS_VERSION, p.ANNOVAR_VERSION, p.VCFTOOLS_VERSION, p.VARSCAN_VERSION, p.SAMTOOLS_OPTIONS])
    job_status(jobname = 'call_variants', resultspath = p.VARIANT_RESULTS, SAMPLE = sample, outputfilename = sample + "/variants.vcf", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    call_variants(sample, call_variants_flag)
    sys.exit(0)