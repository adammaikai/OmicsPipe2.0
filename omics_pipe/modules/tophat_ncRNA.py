#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def tophat_ncRNA(sample, tophat_ncRNA_flag):
    '''Runs TopHat to align .fastq files.
    
    input: 
        .fastq file
    output: 
        accepted_hits.bam
    citation: 
        Trapnell C, Pachter L, Salzberg SL. TopHat: discovering splice junctions with RNA-Seq. Bioinformatics doi:10.1093/bioinformatics/btp120
    link: 
        http://tophat.cbcb.umd.edu/
    parameters from parameters file: 
        RAW_DATA_DIR:
        
        REF_GENES:
        
        TOPHAT_RESULTS:
        
        NONCODE_BOWTIE_INDEX:
        
        TOPHAT_VERSION:
        
        TOPHAT_OPTIONS:
        
        BOWTIE_VERSION:
        
        SAMTOOLS_VERSION:
        '''
    spawn_job(jobname = 'tophat_ncRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/tophat_drmaa_" + p.ENDS + "_ncrna.sh", args_list = [sample,p.RAW_DATA_DIR, p.TOPHAT_RESULTS, p.NONCODE_BOWTIE_INDEX, p.TOPHAT_VERSION, p.TOPHAT_OPTIONS, p.BOWTIE_VERSION, p.SAMTOOLS_VERSION])
    job_status(jobname = 'tophat_ncRNA', resultspath = p.TOPHAT_RESULTS, SAMPLE = sample,  outputfilename = sample + "/accepted_hits.bam", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    tophat_ncRNA(sample, tophat_ncRNA_flag)
    sys.exit(0)
    
    
    