#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def snpir_variants(sample, snpir_variants_flag):    
    '''Calls variants using SNPIR pipeline.
    
    input: 
        Aligned.out.sort.bam or accepted_hits.bam
    output: 
        final_variants.vcf file
    citation: 
        Piskol, R., et al. (2013). "Reliable Identification of Genomic Variants from RNA-Seq Data." The American Journal of Human Genetics 93(4): 641-651. 
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
        
        ENCODING:

	    LOCIBED:

        CALLABLELOCI:
     
        '''
    spawn_job(jobname = 'snpir_variants', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.SNPIR["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.SNPIR["NODES"], ppn = p.SNPIR["CPU"], memory = p.SNPIR["MEMORY"], script = "/snpir_drmaa.sh", args_list = [p.SNPIR["RESULTS"], p.SNPIR["TEMP"], p.SNPIR["SAMTOOLS_VERSION"], p.SNPIR["BWA_VERSION"], p.SNPIR["PICARD_VERSION"], p.SNPIR["GATK_VERSION"], p.SNPIR["BEDTOOLS_VERSION"], p.SNPIR["UCSC_TOOLS_VERSION"], p.SNPIR['GENOME'], p.SNPIR["REPEAT_MASKER"], p.SNPIR["ANNOTATION"], p.SNPIR["RNA_EDIT"], p.SNPIR["DBSNP"], p.SNPIR["MILLS"], p.SNPIR["G1000"], p.OMICSPIPE["WORKING_DIR"], sample, p.SNPIR["BWA_RESULTS"], p.SNPIR["VERSION"], p.SNPIR["CONFIG"], p.SNPIR["DIR"], p.SNPIR["ENCODING"], p.SNPIR["LOCIBED"], p.SNPIR["CALLABLELOCI"]])
    job_status(jobname = 'snpir_variants', resultspath = p.SNPIR["RESULTS"], SAMPLE = sample, outputfilename = sample + "/final_variants.vcf", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

def snpir_variants_plus(sample, snpir_variants_plus_flag):
    '''Calls variants using improved SNPIR pipeline.
    
    input: 
        Aligned.out.sort.bam or accepted_hits.bam
    output: 
        final_variants.vcf file
    citation: 
        Piskol, R., et al. (2013). "Reliable Identification of Genomic Variants from RNA-Seq Data." The American Journal of Human Genetics 93(4): 641-651. 
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
        
        ENCODING:

        LOCIBED:

	    CALLABLELOCI:
        
        '''
    spawn_job(jobname = 'snpir_variants_plus', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.SNPIR["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.SNPIR["NODES"], ppn = p.SNPIR["CPU"], memory = p.SNPIR["MEMORY"], script = "/snpir_plus_drmaa.sh", args_list = [p.SNPIR["RESULTS"], p.SNPIR["TEMP"], p.SNPIR["SAMTOOLS_VERSION"], p.SNPIR["BWA_VERSION"], p.SNPIR["PICARD_VERSION"], p.SNPIR["GATK_VERSION"], p.SNPIR["BEDTOOLS_VERSION"], p.SNPIR["UCSC_TOOLS_VERSION"], p.SNPIR['GENOME'], p.SNPIR["REPEAT_MASKER"], p.SNPIR["ANNOTATION"], p.SNPIR["RNA_EDIT"], p.SNPIR["DBSNP"], p.SNPIR["MILLS"], p.SNPIR["G1000"], p.OMICSPIPE["WORKING_DIR"], sample, p.SNPIR["BWA_RESULTS"], p.SNPIR["VERSION"], p.SNPIR["CONFIG"], p.SNPIR["DIR"], p.SNPIR["ENCODING"], p.SNPIR["LOCIBED"], p.SNPIR["CALLABLELOCI"]])
    job_status(jobname = 'snpir_variants_plus', resultspath = p.SNPIR["RESULTS"], SAMPLE = sample, outputfilename = sample + "/final_variants.vcf", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

def snpir_variants_plus_bwa(sample, extension, snpir_variants_plus_bwa_flag):
    '''Calls variants using improved SNPIR pipeline.
    
    input: 
        fastq.gz
    output: 
        final_variants.vcf file
    citation: 
        Piskol, R., et al. (2013). "Reliable Identification of Genomic Variants from RNA-Seq Data." The American Journal of Human Genetics 93(4): 641-651. 
    link:
        http://lilab.stanford.edu/SNPiR/
    parameters from parameters file:
        RESULTS:
        
        TEMP:
        
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

        SAMPLE:
        
        RAW_DATA:
        
        SNPIR_VERSION:
        
        SNPIR_CONFIG:
        
        SNPIR_DIR:
        
        ENCODING:

        LOCIBED:

        CALLABLELOCI:

        RGLB:

        RGPL:

        CPU:

        MEMORY:

        SAMBLASTER_VERSION:

        SAMBAMBA_VERSION:

        R_VERSION:
        
        '''
  
    sample = sample + extension
    spawn_job(jobname = 'snpir_variants_plus_bwa', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.SNPIR["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.SNPIR["NODES"], ppn = p.SNPIR["CPU"], memory = p.SNPIR["MEMORY"], script = "/snpir_plus_bwa.sh", args_list = [p.SNPIR["RESULTS"], p.SNPIR["TEMP"], p.SNPIR["SAMTOOLS_VERSION"], p.SNPIR["BWA_VERSION"], p.SNPIR["PICARD_VERSION"], p.SNPIR["GATK_VERSION"], p.SNPIR["BEDTOOLS_VERSION"], p.SNPIR["UCSC_TOOLS_VERSION"], p.SNPIR["GENOME"], p.SNPIR["REPEAT_MASKER"], p.SNPIR["ANNOTATION"], p.SNPIR["RNA_EDIT"], p.SNPIR["DBSNP"], p.SNPIR["MILLS"], p.SNPIR["G1000"], p.OMICSPIPE["WORKING_DIR"], sample, p.SNPIR["RAW_DATA"], p.SNPIR["VERSION"], p.SNPIR["CONFIG"], p.SNPIR["DIR"], p.SNPIR["ENCODING"], p.SNPIR["LOCIBED"], p.SNPIR["CALLABLELOCI"], p.SNPIR["RGLB"], p.SNPIR["RGPL"], str(p.SNPIR["CPU"]), p.SNPIR["MEMORY"], p.SNPIR["SAMBLASTER_VERSION"], p.SNPIR["SAMBAMBA_VERSION"], p.SNPIR["R_VERSION"], p.SNPIR["BWA_INDEX"], p.STAR["RESULTS"]])
    job_status(jobname = 'snpir_variants_plus_bwa', resultspath = p.SNPIR["RESULTS"] + "/" + sample, SAMPLE = sample, outputfilename =  sample + "_final_variants_sort.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

def snpir_variants_pe(sample, snpir_variants_pe_flag):
    '''Calls variants using improved SNPIR pipeline.
    
    input: 
        .bam file
    output: 
        final_variants.vcf file
    citation: 
        Piskol, R., et al. (2013). "Reliable Identification of Genomic Variants from RNA-Seq Data." The American Journal of Human Genetics 93(4): 641-651. 
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
        
        ENCODING:

        LOCIBED:

        CALLABLELOCI:
        
        '''
    spawn_job(jobname = 'snpir_variants_pe', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.SNPIRPE["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.SNPIRPE["NODES"], ppn = p.SNPIRPE["CPU"], memory = p.SNPIRPE["MEMORY"], script = "/snpir_pe.sh", args_list = [p.SNPIRPE["RESULTS"], p.SNPIRPE["TEMP"], p.SNPIRPE["SAMTOOLS_VERSION"], p.SNPIRPE["BWA_VERSION"], p.SNPIRPE["PICARD_VERSION"], p.SNPIRPE["GATK_VERSION"], p.SNPIRPE["BEDTOOLS_VERSION"], p.SNPIRPE["UCSC_TOOLS_VERSION"], p.SNPIRPE['GENOME'], p.SNPIRPE["REPEAT_MASKER"], p.SNPIRPE["ANNOTATION"], p.SNPIRPE["RNA_EDIT"], p.SNPIRPE["DBSNP"], p.SNPIRPE["MILLS"], p.SNPIRPE["G1000"], p.OMICSPIPE["WORKING_DIR"], sample, p.SNPIRPE["INPUT"], p.SNPIRPE["VERSION"], p.SNPIRPE["CONFIG"], p.SNPIRPE["DIR"], p.SNPIRPE["ENCODING"], p.SNPIRPE["LOCIBED"], p.SNPIRPE["CALLABLELOCI"], p.SNPIRPE["BAM_FILE_NAME"]])
    job_status(jobname = 'snpir_variants_pe', resultspath = p.SNPIRPE["RESULTS"], SAMPLE = sample, outputfilename = sample + "/final_variants.vcf", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    snpir_variants(sample, snpir_variants_flag)
    snpir_variants_plus(sample, snpir_variants_plus_flag)
    snpir_variants_plus_bwa(sample, extension, snpir_variants_plus_bwa_flag)
    snpir_variants_pe(sample, snpir_variants_pe_flag)
    sys.exit(0)
