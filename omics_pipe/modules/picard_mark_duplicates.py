#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def picard_mark_duplicates(sample, picard_mark_duplicates_flag):
    '''Picard tools Mark Duplicates.
    
    input: 
        sorted.bam
    output: 
        _sorted.rg.md.bam
    citation: 
        http://picard.sourceforge.net/
    link: 
        http://picard.sourceforge.net/
    parameters from parameters file:
        BWA_RESULTS:
            
        TEMP_DIR:
            
        PICARD_VERSION:
        
        SAMTOOLS_VERSION:
    ''' 
    spawn_job(jobname = 'picard_mark_duplicates', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "29gb", script = "/picard_mark_duplicates.sh", args_list = [p.BWA_RESULTS, sample, p.TEMP_DIR, p.PICARD_VERSION, p.SAMTOOLS_VERSION])
    job_status(jobname = 'picard_mark_duplicates', resultspath = p.BWA_RESULTS, SAMPLE = sample,  outputfilename = sample + "/" + sample + "_sorted.rg.md.bam", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    picard_mark_duplicates(sample, picard_mark_duplicates_flag)
    sys.exit(0)