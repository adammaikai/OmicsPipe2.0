#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Varscan(sample, extension, Varscan_flag):
    '''Calling variants with Varscan for whole genome sequencing.
        
        input:
        _gatk_recal.bam
        output:
        _varscan_somatic.4.1.vcf.gz
        citation:
        
        link:
        
        parameters from parameters file:
        
        TEMP_DIR:
        
        GENOME:
        
        SAMTOOLS_VERSION:
        
        VARSCAN_VERSION:
        
        CAPTURE_KIT_BED:
        
        ALIGNMENT_DIR:
        
        RESULTS:
        
        OPTIONS:
        '''
 
    dna_sample = sample + extension
    spawn_job(jobname = 'Varscan', SAMPLE = dna_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.VARSCAN["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.VARSCAN["NODES"], ppn = p.VARSCAN["CPU"], memory = p.VARSCAN["MEMORY"], script = "/Varscan.sh", args_list = [dna_sample, p.OMICSPIPE["TEMP_DIR"], p.VARSCAN["GENOME"], p.VARSCAN["SAMTOOLS_VERSION"], p.VARSCAN["VERSION"], p.VARSCAN["ALIGNMENT_DIR"], p.VARSCAN["RESULTS"], p.VARSCAN["OPTIONS"], p.VARSCAN["R_VERSION"], p.VARSCAN["VCFTOOLS_VERSION"]])
    job_status(jobname = 'Varscan', resultspath = p.VARSCAN["RESULTS"] + "/" + dna_sample, SAMPLE = dna_sample,  outputfilename = dna_sample + "_varscan.4.1.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])

    return

if __name__ == '__main__':
   Varscan(sample, extension, Varscan_flag)
   sys.exit(0)
