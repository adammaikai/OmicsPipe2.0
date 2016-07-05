#!/usr/bin/env python
import os
from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def cuffmerge_miRNA(step, cuffmerge_miRNA_flag):
    '''Runs cuffmerge to merge .gtf files from Cufflinks. 
    
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
        
        CUFFMERGE_OPTIONS:
        
        CUFFLINKS_VERSION:
        '''
    check_create_dir(p.CUFFMERGE_RESULTS)
    f = open(p.CUFFMERGE_RESULTS + '/assembly_GTF_list.txt', 'w')
    for sample in p.SAMPLE_LIST:
        f.write(p.CUFFLINKS_RESULTS + "/" + sample + "/transcripts.gtf\n")
    f.close()
    gtf_list = p.CUFFMERGE_RESULTS + '/assembly_GTF_list.txt'
    spawn_job(jobname = 'cuffmerge_miRNA', SAMPLE = step, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "120:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/cuffmerge_drmaa.sh", args_list = [gtf_list,p.CUFFMERGE_RESULTS,p.miRNA_GTF,p.GENOME,p.CUFFMERGE_OPTIONS,p.CUFFLINKS_VERSION])
    job_status(jobname = 'cuffmerge_miRNA', resultspath = p.CUFFMERGE_RESULTS, SAMPLE = step,  outputfilename = "merged.gtf", FLAG_PATH = p.FLAG_PATH)
    return



if __name__ == '__main__':
    cuffmerge_miRNA(step, cuffmerge_miRNA_flag)
    sys.exit(0)