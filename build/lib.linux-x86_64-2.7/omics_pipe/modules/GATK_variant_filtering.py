#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_variant_filtering(sample, GATK_variant_filtering_flag):
    '''GATK_variant_filtering.
    
    input: 
        sorted.rg.md.bam
    output: 
        .ready.bam
    citation: 
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
    link: GATK_variant_filtering
        http://www.broadinstitute.org/gatk/
    parameters from parameters file:
        VARIANT_RESULTS:
            
        TEMP_DIR:
            
        GATK_VERSION:
        
        GENOME:
        
        DBSNP:
        
        MILLS:
        
        OMNI:
        
        HAPMAP:
        
        R_VERSION:
        
        G1000_SNPs: 
        
        G1000_Indels:
        
        
        ''' 
    
    spawn_job(jobname = 'GATK_variant_filtering', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/GATK_variant_filtering_" + p.GATK_VERSION + ".sh", args_list = [p.VARIANT_RESULTS, sample, p.TEMP_DIR,p.GATK_VERSION, p.GENOME, p.DBSNP,p.MILLS_G1000,p.OMNI,p.HAPMAP, p.R_VERSION, p.G1000, p.G1000_SNPs, p.G1000_Indels])
    job_status(jobname = 'GATK_variant_filtering', resultspath = p.BWA_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + ".vqsr.vcf", FLAG_PATH = p.FLAG_PATH)
    return


def GATK_variant_filtering_group(sample, GATK_variant_filtering_group_flag):
    '''GATK_variant_filtering.
    
    input: 
        sorted.rg.md.bam
    output: 
        .ready.bam
    citation: 
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
    link: GATK_variant_filtering
        http://www.broadinstitute.org/gatk/
    parameters from parameters file:
        
        VARIANT_RESULTS:            
        
        TEMP_DIR:
            
        GATK_VERSION:
        
        GENOME:
        
        DBSNP:
        
        MILLS_G1000:
        
        OMNI:
        
        HAPMAP:
        
        R_VERSION:
        
        G1000:

        ''' 
    
    spawn_job(jobname = 'GATK_variant_filtering_group', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "30gb", script = "/GATK_variant_filtering_group_" + p.GATK_VERSION + ".sh", args_list = [p.VARIANT_RESULTS, sample, p.TEMP_DIR,p.GATK_VERSION, p.GENOME, p.DBSNP,p.MILLS_G1000,p.OMNI,p.HAPMAP, p.R_VERSION,p.G1000_SNPs, p.G1000_Indels])
    job_status(jobname = 'GATK_variant_filtering_group', resultspath = p.VARIANT_RESULTS, SAMPLE = sample,  outputfilename = "GenotypeGVCF_output.vqsr.vcf", FLAG_PATH = p.FLAG_PATH)
    return


if __name__ == '__main__':
    GATK_variant_filtering(sample, GATK_variant_filtering_flag)
    GATK_variant_filtering_group(sample, GATK_variant_filtering_group_flag)
    sys.exit(0)
    
    
    