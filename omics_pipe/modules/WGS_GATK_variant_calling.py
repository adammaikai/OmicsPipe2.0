#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def WGS_GATK_variant_calling(sample, WGS_GATK_variant_calling_flag):
    '''Calling variants with GATK for whole exome sequencing.
        
        input:
        _WES_realigned_sorted.bam
        output:
        _vqsr.vcf
        citation:
        
        link:
        
        parameters from parameters file:
        
        TEMP_DIR:
        
        GENOME:
        
        GATK_VERSION:
        
        R_VERSION:
        
        CAPTURE_KIT_BED:
        
        ALIGNMENT_DIR:
        
        RESULTS:
        
        OMNI:
        
        HAPMAP:
        
        DBSNP:
        
        MILLS:
        
        G1000:
        
        KG:
        '''
    
    spawn_job(jobname = 'WGS_GATK_variant_calling', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = "240:00:00", queue = p.OMICSPIPE["QUEUE"], nodes = 1, ppn = 32, memory = "58gb", script = "/WGS_GATK_variant_calling.sh", args_list = [sample, p.OMICSPIPE["TEMP_DIR"], p.PREPROCESS["GENOME"], p.GATK["VERSION"], p.GATK["R_VERSION"], p.PREPROCESS["ALIGNMENT_DIR"], p.GATK["RESULTS"], p.GATK["OMNI"], p.GATK["HAPMAP"], p.GATK["DBSNP"], p.GATK["MILLS"], p.GATK["G1000"], p.GATK["KG"]])
    job_status(jobname = 'WGS_GATK_variant_calling', resultspath = p.GATK["RESULTS"], SAMPLE = sample,  outputfilename = sample + "_vqsr.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    WGS_GATK_variant_calling(sample, WGS_GATK_variant_calling_flag)
    sys.exit(0)
