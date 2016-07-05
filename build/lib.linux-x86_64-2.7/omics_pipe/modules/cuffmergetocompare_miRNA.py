#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def cuffmergetocompare_miRNA(step, cuffmergetocompare_miRNA_flag):
    '''Runs cuffcompare to annotate merged .gtf files from Cufflinks. 
    
    input: 
        assembly_GTF_list.txt
    output: 
        merged.gtf
    citation: 
        Trapnell C, et al. Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation  Nature Biotechnology doi:10.1038/nbt.1621
    link: 
        http://cufflinks.cbcb.umd.edu/
    parameters from parameters file: 
        CUFFMERGE_RESULTS:
        
        miRNA_GTF:
        
        GENOME:
        
        CUFFMERGETOCOMPARE_OPTIONS:
        
        CUFFLINKS_VERSION:
        '''
    spawn_job(jobname = 'cuffmergetocompare_miRNA', SAMPLE = step, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/cuffmergetocompare_drmaa.sh", args_list = [p.CUFFMERGE_RESULTS,p.miRNA_GTF,p.GENOME,p.CUFFMERGETOCOMPARE_OPTIONS,p.CUFFLINKS_VERSION])
    job_status(jobname = 'cuffmergetocompare_miRNA', resultspath = p.CUFFMERGE_RESULTS, SAMPLE = step, outputfilename = "merged.gtf", FLAG_PATH = p.FLAG_PATH)
    return


if __name__ == '__main__':
    cuffmergetocompare_miRNA(step, cuffmergetocompare_miRNA_flag)
    sys.exit(0)