#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def star_miRNA(sample, star_miRNA_flag):
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
        
        TRIMMED_DATA_PATH:
        
        STAR_INDEX:
        
        STAR_OPTIONS:
        
        STAR_RESULTS:
        
        SAMTOOLS_VERSION:
        
        STAR_VERSION:
        ''' 
    spawn_job(jobname = 'star_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/star_drmaa_" + p.ENDS + ".sh", args_list = [sample, p.TRIMMED_DATA_PATH, p.STAR_INDEX, p.STAR_OPTIONS, p.STAR_RESULTS, p.SAMTOOLS_VERSION, p.STAR_VERSION])
    job_status(jobname = 'star_miRNA', resultspath = p.STAR_RESULTS, SAMPLE = sample, outputfilename = sample + "/SJ.out.tab", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    star_miRNA(sample, star_miRNA_flag)
    sys.exit(0)