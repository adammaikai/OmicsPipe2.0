#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def intogen(sample, intogen_flag):   
    '''Runs Intogen to rank mutations and implication for cancer phenotype. Follows variant calling. 
    
    input: 
        .vcf
    output: 
        variant list
    citation: 
        Gonzalez-Perez et al. 2013. Intogen mutations identifies cancer drivers across tumor types. Nature Methods 10, 1081-1082. 
    link: 
        http://www.intogen.org/
    parameters from parameter file: 
        VCF_FILE:
        
        INTOGEN_OPTIONS:
        
        INTOGEN_RESULTS:
        
        INTOGEN_VERSION:
              
        USERNAME:
        
        WORKING_DIR:
        
        TEMP_DIR:
        
        SCHEDULER:
        
        VARIANT_RESULTS:
        '''
    #vcf_file = p.SNPIR_RESULTS + "/" +sample + "/final_variants.vcf"
    spawn_job(jobname = 'intogen', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/intogen_drmaa.sh", args_list = [p.VCF_FILE, p.INTOGEN_OPTIONS, p.INTOGEN_RESULTS, p.INTOGEN_VERSION,  p.USERNAME, p.WORKING_DIR, p.TEMP_DIR, sample, p.SCHEDULER, p.VARIANT_RESULTS])
    job_status(jobname = 'intogen', resultspath = p.INTOGEN_RESULTS, SAMPLE = sample, outputfilename = sample + "/variant_genes.tsv", FLAG_PATH = p.FLAG_PATH)
    return
if __name__ == '__main__':
    intogen(sample, intogen_flag)
    sys.exit(0)