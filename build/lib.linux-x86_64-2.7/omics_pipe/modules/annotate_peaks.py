#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)



def annotate_peaks(step, annotate_peaks_flag):
    '''Runs HOMER to annotate peak track from ChIPseq data.
    
    input: 
        .tag input file
    output: 
        .txt file
    citation:
         Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
    link: 
        http://homer.salk.edu/homer/
    parameters from parameters file:
        PAIR_LIST:
        
        HOMER_RESULTS:
        
        HOMER_VERSION:
        
        TEMP_DIR:
        
        HOMER_GENOME:
        
        HOMER_ANNOTATE_OPTIONS:
        '''
    spawn_job(jobname = 'annotate_peaks', SAMPLE = step, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/annotate_peaks.sh", args_list = [p.PAIR_LIST, p.HOMER_RESULTS, p.HOMER_VERSION, p.TEMP_DIR, p.HOMER_GENOME, p.HOMER_ANNOTATE_OPTIONS])
    job_status(jobname = 'annotate_peaks', resultspath = p.HOMER_RESULTS, SAMPLE = step, outputfilename =  "*_regions_annot.txt", FLAG_PATH = p.FLAG_PATH)
    return


if __name__ == '__main__':
    annotate_peaks(step, annotate_peaks_flag)
    sys.exit(0)