#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def mojo(sample, extension, mojo_flag):
    '''Runs mojo on a processed .bam file.
        
        input:
        fastq.gz
        output:
    
        citation:
        parameters from parameters file:

        RAW_DATA:
        
        RESULTS:
        
        MOJO_VERSION:

        CONFIG:
        '''
    rna_tumor_sample = sample + extension
    spawn_job(jobname = 'mojo', SAMPLE = rna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.MOJO["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.MOJO["NODES"], ppn = p.MOJO["CPU"], memory = p.MOJO["MEMORY"], script = "/mojo.sh", args_list = [rna_tumor_sample, p.MOJO["RAW_DATA"], p.MOJO["RESULTS"], p.MOJO["VERSION"], p.MOJO["CONFIG"], p.MOJO["CPU"], p.MOJO["MEMORY"], p.MOJO["FUSIONANNOTATOR_VERSION"], p.MOJO["FUSIONANNOTATOR_LIB"]])
    job_status(jobname = 'mojo', resultspath = p.MOJO["RESULTS"] + "/" + rna_tumor_sample, SAMPLE = rna_tumor_sample,  outputfilename = "genome_results.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    mojo(sample, extension, mojo_flag)
    sys.exit(0)
