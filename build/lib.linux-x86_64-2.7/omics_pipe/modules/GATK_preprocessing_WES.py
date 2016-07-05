#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_preprocessing_WES(sample, GATK_preprocessing_WES_flag):
    '''GATK preprocessing steps for whole exome sequencing.
    
    input: 
        sorted.rg.md.bam
    output: 
        .ready.bam
    citation: 
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
    link: 
        http://www.broadinstitute.org/gatk/
    parameters from parameters file:
        BWA_RESULTS:
            
        TEMP_DIR:
            
        GATK_VERSION:
        
        GENOME:
        
        DBSNP:
        
        MILLS:
        
        G1000:
        
        CAPTURE_KIT_BED:
        
        SAMTOOLS_VERSION:
        ''' 
    
    spawn_job(jobname = 'GATK_preprocessing_WES', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/GATK_preprocessing_WES_" + p.GATK_VERSION + ".sh", args_list = [p.BWA_RESULTS, sample, p.TEMP_DIR,p.GATK_VERSION, p.GENOME, p.DBSNP,p.MILLS,p.G1000,p.CAPTURE_KIT_BED, p.SAMTOOLS_VERSION])
    job_status(jobname = 'GATK_preprocessing_WES', resultspath = p.BWA_RESULTS, SAMPLE = sample,  outputfilename = sample + "/" + sample + ".ready.bam", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    GATK_preprocessing_WES(sample, GATK_preprocessing_WES_flag)
    sys.exit(0)