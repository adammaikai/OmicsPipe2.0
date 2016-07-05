#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_BQSR(sample, extension, GATK_BQSR_flag):
    '''Recalibrate base quality scores.
        
        input:
        _WES_realigned_sorted.bam
        output:
        _gatk_recal.bam
        citation:
        
        link:
        
        parameters from parameters file:
        
        TEMP_DIR:
        
        GENOME:
        
        GATK_VERSION:
        
        R_VERSION:
        
        CAPTURE_KIT_BED:
        
        ALIGNMENT_DIR:
        
        DBSNP:
        
        MILLS:
        
        G1000:
        '''
    sample = sample + extension
    spawn_job(jobname = 'GATK_BQSR', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.BQSR["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.BQSR["NODES"], ppn = p.BQSR["CPU"], memory = p.BQSR["MEMORY"], script = "/GATK_BQSR.sh", args_list = [sample, p.OMICSPIPE["TEMP_DIR"], p.SAMTOOLS["GENOME"], p.BQSR["VERSION"], p.VARSCAN["R_VERSION"], p.CAPTURE_KIT_BED, p.BQSR["ALIGNMENT_DIR"], p.BQSR["DBSNP"], p.BQSR["MILLS"], p.BQSR["G1000"]])
    job_status(jobname = 'GATK_BQSR', resultspath = p.BQSR["ALIGNMENT_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + "_gatk_recal.bam", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    GATK_BQSR(sample, extension, GATK_BQSR_flag)
    sys.exit(0)
