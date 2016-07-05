#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def tophat_miRNA(sample, tophat_miRNA_flag):
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
        
        miRNA_GTF:
        
        TOPHAT_RESULTS:
        
        miRNA_BOWTIE_INDEX:
        
        TOPHAT_VERSION:
        
        TOPHAT_OPTIONS:
        
        BOWTIE_VERSION:
        
        SAMTOOLS_VERSION:
        '''
    spawn_job(jobname = 'tophat_miRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/tophat_drmaa_miRNA.sh", args_list = [sample,p.TRIMMED_DATA_PATH, p.TOPHAT_RESULTS, p.miRNA_BOWTIE_INDEX, p.TOPHAT_VERSION, p.TOPHAT_OPTIONS, p.miRNA_GTF, p.BOWTIE_VERSION, p.SAMTOOLS_VERSION])
    job_status(jobname = 'tophat_miRNA', resultspath = p.TOPHAT_RESULTS, SAMPLE = sample, outputfilename = sample + "/accepted_hits.bam", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    tophat_miRNA(sample, tophat_miRNA_flag)
    sys.exit(0)