#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def macs(step, macs_flag):
    '''Runs MACS to call peaks from ChIPseq data.
    input: 
        .fastq file
    output: 
        peaks and .bed file
    citation: 
        Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137
    link: 
        http://liulab.dfci.harvard.edu/MACS/
    parameters from parameters file: 
        PAIR_LIST:
        
        BOWTIE_RESULTS:
        
        CHROM_SIZES:
        
        MACS_RESULTS:
        
        MACS_VERSION:
        
        TEMP_DIR:
        
        BEDTOOLS_VERSION:
        
        PYTHON_VERSION:
        '''
    spawn_job(jobname = 'macs', SAMPLE = step, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/macs_drmaa.sh", args_list = [p.PAIR_LIST, p.BOWTIE_RESULTS, p.CHROM_SIZES, p.MACS_RESULTS, p.MACS_VERSION, p.TEMP_DIR, p.BEDTOOLS_VERSION, p.PYTHON_VERSION])
    job_status(jobname = 'macs', resultspath = p.MACS_RESULTS, SAMPLE = step, outputfilename = p.PAIR_LIST[1] + "_macs_enrichment.bed.gz", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    macs(step, macs_flag)
    sys.exit(0)