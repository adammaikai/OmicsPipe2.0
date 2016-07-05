#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_variant_joint_analysis(sample, GATK_variant_joint_analysis_flag):
    '''GATK Group Variant Discovery.
    
    input: 
        sorted.rg.md.bam
    output: 
        .ready.bam
    citation: 
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
    link: GATK_variant_discovery
        http://www.broadinstitute.org/gatk/
    parameters from parameters file:
 
        TEMP_DIR:
            
        GATK_VERSION:
        
        GENOME:
        
        SAMPLE_LIST:
        
        VARIANT_RESULTS:
        
        DBSNP:
        ''' 
    sample_list = ('\t'.join(map(str,p.SAMPLE_LIST))) 
    spawn_job(jobname = 'GATK_variant_joint_analysis', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/GATK_genotypeGVCFs_" + p.GATK_VERSION + ".sh", args_list = [p.VARIANT_RESULTS, sample_list, p.TEMP_DIR, p.GATK_VERSION, p.GENOME, p.DBSNP])
    job_status(jobname = 'GATK_variant_joint_analysis', resultspath = p.VARIANT_RESULTS, SAMPLE = sample,  outputfilename = "GenotypeGVCF_output.vcf", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    GATK_variant_joint_analysis(sample, GATK_variant_joint_analysis_flag)
    sys.exit(0)