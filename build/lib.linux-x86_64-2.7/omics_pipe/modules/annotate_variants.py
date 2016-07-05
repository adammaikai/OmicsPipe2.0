#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def annotate_variants(sample, annotate_variants_flag):
    '''Annotates variants with ANNOVAR variant annotator. Follows VarCall.
    input: 
        .vcf
    output: 
        .vcf
    citation: 
        Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010
    link: 
        http://www.openbioinformatics.org/annovar/
    parameters from parameters file:
        VARIANT_RESULTS:
        
        ANNOVARDB:
        
        ANNOVAR_OPTIONS:
        
        ANNOVAR_OPTIONS2:
        
        TEMP_DIR:
        
        ANNOVAR_VERSION:
        
        VCFTOOLS_VERSION:
        '''
    spawn_job(jobname = 'annotate_variants', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "31gb", script = "/varannotate_drmaa_bash.sh", args_list = [sample,p.VARIANT_RESULTS,p.ANNOVARDB,p.ANNOVAR_OPTIONS, p.ANNOVAR_OPTIONS2, p.TEMP_DIR, p.ANNOVAR_VERSION, p.VCFTOOLS_VERSION])
    job_status(jobname = 'annotate_variants', resultspath = p.VARIANT_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + ".vcf.gz", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    annotate_variants(sample, annotate_variants_flag)
    sys.exit(0)