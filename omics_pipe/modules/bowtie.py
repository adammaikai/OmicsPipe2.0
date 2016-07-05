#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def bowtie(sample, bowtie_flag):
    '''Runs Bowtie to align .fastq files.
    
    input:
        .fastq file
    output:
        sample.bam
    citation: 
        Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 10:R25
    link: 
        http://bowtie-bio.sourceforge.net/index.shtml
    parameters from parameters file: 
        ENDS: 
        
        TRIMMED_DATA_PATH:
        
        BOWTIE_OPTIONS:
        
        BOWTIE_INDEX:
        
        BOWTIE_RESULTS:
        
        BOWTIE_VERSION:
        
        SAMTOOLS_VERSION:
        
        BEDTOOLS_VERSION:
        
        TEMP_DIR:
        '''
    spawn_job(jobname = 'bowtie', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/bowtie_drmaa_" + p.ENDS + ".sh", args_list = [sample,p.TRIMMED_DATA_PATH,p.BOWTIE_OPTIONS, p.BOWTIE_INDEX, p.BOWTIE_RESULTS, p.BOWTIE_VERSION, p.SAMTOOLS_VERSION, p.BEDTOOLS_VERSION, p.TEMP_DIR])
    job_status(jobname = 'bowtie', resultspath = p.BOWTIE_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + ".bam", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    bowtie(sample, bowtie_flag)
    sys.exit(0)