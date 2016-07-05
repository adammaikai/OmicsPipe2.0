#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def star(sample, star_flag):
    '''Runs STAR to align .fastq files.
    
    input: 
        .fastq file
    output: 
        Aligned.out.bam
    citation: 
        A. Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635 "STAR: ultrafast universal RNA-seq aligner"
    link: 
        https://code.google.com/p/rna-star/
    parameters from parameters file: 
        ENDS:
        
        RAW_DATA_DIR:
        
        STAR_INDEX:
        
        STAR_OPTIONS:
        
        STAR_RESULTS:
        
        SAMTOOLS_VERSION:
        
        STAR_VERSION:

	COMPRESSION:

	REF_GENES:
        ''' 
    spawn_job(jobname = 'star', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 32, memory = "40gb", script = "/star_drmaa_" + p.ENDS + ".sh", args_list = [sample, p.RAW_DATA_DIR, p.STAR_INDEX, p.STAR_OPTIONS, p.STAR_RESULTS, p.SAMTOOLS_VERSION, p.STAR_VERSION, p.COMPRESSION, p.REF_GENES])
    job_status(jobname = 'star', resultspath = p.STAR_RESULTS, SAMPLE = sample, outputfilename = sample + "/SJ.out.tab", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    star(sample, star_flag)
    sys.exit(0)
