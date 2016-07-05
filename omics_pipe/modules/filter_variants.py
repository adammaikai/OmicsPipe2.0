#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def filter_variants(sample, filter_variants_flag):   
    '''Filters variants to remove common variants.
    
    input: 
        .bam or .sam file
    output: 
        .vcf file 
    citation: 
        Piskol et al. 2013. Reliable identification of genomic variants from RNA-seq data. The American Journal of Human Genetics 93: 641-651.
    link:
        http://lilab.stanford.edu/SNPiR/        
    parameters from parameters file: 
        VARIANT_RESULTS:
        
        TEMP_DIR:
        
        SAMTOOLS_VERSION:
        
        BWA_VERSION:
        
        PICARD_VERSION:
        
        GATK_VERSION:
        
        BEDTOOLS_VERSION:
        
        UCSC_TOOLS_VERSION:
        
        GENOME:
        
        REPEAT_MASKER:
        
        SNPIR_ANNOTATION:
        
        RNA_EDIT:
        
        DBSNP:
        
        MILLS:
        
        G1000:
        
        WORKING_DIR:
        
        BWA_RESULTS:
        
        SNPIR_VERSION:
        
        SNPIR_CONFIG:
        
        SNPIR_DIR:
        
        SNPEFF_VERSION:
        
        dbNSFP:
        
        VCFTOOLS_VERSION:
        
        WORKING_DIR:
        
        SNP_FILTER_OUT_REF:
        '''
    spawn_job(jobname = 'filter_variants', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "28gb", script = "/filter_snpir_drmaa.sh", args_list = [p.VARIANT_RESULTS, p.TEMP_DIR, p.SAMTOOLS_VERSION, p.BWA_VERSION, p.PICARD_VERSION, p.GATK_VERSION, p.BEDTOOLS_VERSION, p.UCSC_TOOLS_VERSION, p.GENOME, p.REPEAT_MASKER, p.SNPIR_ANNOTATION, p.RNA_EDIT, p.DBSNP, p.MILLS, p.G1000, p.WORKING_DIR, sample, p.BWA_RESULTS, p.SNPIR_VERSION, p.SNPIR_CONFIG, p.SNPIR_DIR, p.SNPEFF_VERSION, p.dbNSFP, p.VCFTOOLS_VERSION, p.WORKING_DIR, p.SNP_FILTER_OUT_REF])
    job_status(jobname = 'filter_variants', resultspath = p.VARIANT_RESULTS, SAMPLE = sample, outputfilename = sample + "/intogen_input.vcf", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    filter_variants(sample, filter_variants_flag)
    sys.exit(0)