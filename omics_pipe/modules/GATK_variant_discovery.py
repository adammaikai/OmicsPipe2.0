#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_variant_discovery(sample, GATK_variant_discovery_flag):
    '''GATK_variant_discovery.
    
    input: 
        sorted.rg.md.bam
    output: 
        .ready.bam
    citation: 
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
    link: GATK_variant_discovery
        http://www.broadinstitute.org/gatk/
    parameters from parameters file:
        BWA_RESULTS:
            
        TEMP_DIR:
            
        GATK_VERSION:
        
        GENOME:
        
        DBSNP:
        
        VARIANT_RESULTS:
        ''' 
    
    spawn_job(jobname = 'GATK_variant_discovery', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/GATK_variant_discovery_" + p.GATK_VERSION + ".sh", args_list = [p.BWA_RESULTS, sample, p.TEMP_DIR,p.GATK_VERSION, p.GENOME, p.DBSNP, p.VARIANT_RESULTS])
    job_status(jobname = 'GATK_variant_discovery', resultspath = p.BWA_RESULTS, SAMPLE = sample,  outputfilename = sample + "/" +  sample + ".raw.vcf", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    GATK_variant_discovery(sample, GATK_variant_discovery_flag)
    sys.exit(0)
