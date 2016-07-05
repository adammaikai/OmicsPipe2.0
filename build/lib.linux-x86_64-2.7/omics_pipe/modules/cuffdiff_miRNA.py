#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def cuffdiff_miRNA(step, cuffdiff_miRNA_flag):
    '''Runs Cuffdiff to perform differential expression. Runs after Cufflinks. Part of Tuxedo Suite.
    
    input: 
        .bam files
    output: 
        differential expression results
    citation: 
        Trapnell C, et al. Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation Nature Biotechnology doi:10.1038/nbt.1621
    link:
         http://cufflinks.cbcb.umd.edu/
    parameters from parameters file:
        CUFFDIFF_RESULTS:
        
        GENOME:
        
        CUFFDIFF_OPTIONS:
        
        CUFFMERGE_RESULTS: 
        
        CUFFDIFF_INPUT_LIST_COND1:
        
        CUFFDIFF_INPUT_LIST_COND2: 
        
        CUFFLINKS_VERSION:
        '''
    spawn_job(jobname = 'cuffdiff_miRNA', SAMPLE = step, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/cuffdiff_drmaa.sh", args_list = [p.CUFFDIFF_RESULTS,p.GENOME,p.CUFFDIFF_OPTIONS,p.CUFFMERGE_RESULTS, p.CUFFDIFF_INPUT_LIST_COND1, p.CUFFDIFF_INPUT_LIST_COND2, p.CUFFLINKS_VERSION])
    job_status(jobname = 'cuffdiff_miRNA', resultspath = p.CUFFDIFF_RESULTS, SAMPLE = step, outputfilename = "gene_exp.diff", FLAG_PATH = p.FLAG_PATH)
    return


if __name__ == '__main__':
    cuffdiff_miRNA(step, cuffdiff_miRNA_flag)
    sys.exit(0)