#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Varscan_somatic(sample, tumor, normal, Varscan_somatic_flag):
    '''Runs Varscan2 on paired tumor/normal samples to detect somatic point mutations in cancer genomes.
        
        input:
        .pileup
        output:
        _varscan_somatic.vcf.gz
        citation:
        parameters from parameters file:

        ALIGNMENT_DIR:
        
        TEMP_DIR:
        
        VARSCAN_VERSION:
        
        GENOME:
        
        CAPTURE_KIT_BED:
        '''
    
    spawn_job(jobname = 'Varscan_somatic', SAMPLE = sample + tumor  + normal, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.VARSCAN["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.VARSCAN["NODES"], ppn = p.VARSCAN["CPU"], memory = p.VARSCAN["MEMORY"], script = "/Varscan_somatic.sh", args_list = [sample, p.VARSCAN["VERSION"], p.VARSCAN["R_VERSION"], p.VARSCAN["VCFTOOLS_VERSION"], p.VARSCAN["SOMATIC_RESULTS"], p.VARSCAN["ALIGNMENT_DIR"], normal, tumor, p.VARSCAN["VCFLIB_VERSION"]])
    job_status(jobname = 'Varscan_somatic', resultspath = p.VARSCAN["RESULTS"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + tumor + normal + "_varscan_somatic.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    Varscan_somatic(sample, tumor, normal, Varscan_somatic_flag)
    sys.exit(0)
