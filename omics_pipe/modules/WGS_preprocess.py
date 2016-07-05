#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def WGS_preprocess(sample, extension, WGS_preprocess_flag):
    '''Preprocessing steps for whole genome sequencing.
        
        input:
        fastq
        output:
        _gatk_recal.bam
        citation:
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
        link:
        http://www.broadinstitute.org/gatk/
        parameters from parameters file:
        
        TEMP_DIR:
        
        GENOME:
        
        BWA_VERSION:
        
        ABRA_VERSION:
        
        SAMTOOLS_VERSION:
        
        SAMBAMBA_VERSION:
        
        SAMBLASTER_VERSION:
    
        SAMTOOLS_VERSION:
        '''
    
    sample = sample + extension
    spawn_job(jobname = 'WGS_preprocess', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = "240:00:00", queue = p.OMICSPIPE["QUEUE"], nodes = 1, ppn = 32, memory = "58gb", script = "/WGS_preprocess.sh", args_list = [sample, p.OMICSPIPE["TEMP_DIR"], p.PREPROCESS["GENOME"], p.PREPROCESS["BWA_INDEX"], p.PREPROCESS["BWA_VERSION"], p.PREPROCESS["SAMTOOLS_VERSION"], p.PREPROCESS["SAMBAMBA_VERSION"], p.PREPROCESS["SAMBLASTER_VERSION"], p.GATK["VERSION"], p.PREPROCESS["ALIGNMENT_DIR"], p.PREPROCESS["FASTQ_PATH"], p.GATK["DBSNP"], p.GATK["MILLS"], p.GATK["G1000"], p.CAPTURE_KIT_BED])
    job_status(jobname = 'WGS_preprocess', resultspath = p.PREPROCESS["ALIGNMENT_DIR"], SAMPLE = sample,  outputfilename = sample + "_gatk_recal.bam", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    WGS_preprocess(sample, extension, WGS_preprocess_flag)
    sys.exit(0)
