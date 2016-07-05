#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def htseq_ncRNA(sample, htseq_ncRNA_flag):
    '''Runs htseq-count to get raw count data from alignments.
    
    input: 
        Aligned.out.sort.bam
    output: 
        counts.txt
    citation: 
        Simon Anders, EMBL
    link: 
        http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
    parameters from parameters file: 
        ALIGNMENT_RESULTS:
        
        HTSEQ_OPTIONS:
        
        ncRNA_GFF:
        
        HTSEQ_ncRNA_RESULTS: 
        
        TEMP_DIR:
        
        SAMTOOLS_VERSION:
        
        BAM_FILE_NAME:
        '''
    spawn_job(jobname = 'htseq_ncRNA', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/htseq_drmaa.sh", args_list = [sample,p.ALIGNMENT_RESULTS,p.HTSEQ_OPTIONS,p.ncRNA_GFF,p.HTSEQ_ncRNA_RESULTS,p.TEMP_DIR,p.SAMTOOLS_VERSION, p.BAM_FILE_NAME, p.PYTHON_VERSION])
    job_status(jobname = 'htseq_ncRNA', resultspath = p.HTSEQ_RESULTS, SAMPLE = sample, outputfilename = sample + "_counts.txt", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    htseq_ncRNA(sample, htseq_ncRNA_flag)
    sys.exit(0)