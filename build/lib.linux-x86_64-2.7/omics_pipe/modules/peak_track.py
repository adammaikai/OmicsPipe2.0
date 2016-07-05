#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def peak_track(step, peak_track_flag):
    '''Runs HOMER to create peak track from ChIPseq data.
    
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
        '''
    spawn_job(jobname = 'peak_track', SAMPLE = step, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/peak_track.sh", args_list = [p.PAIR_LIST, p.HOMER_RESULTS, p.HOMER_VERSION, p.TEMP_DIR])
    job_status(jobname = 'peak_track', resultspath = p.HOMER_RESULTS, SAMPLE = step, outputfilename =  "*_peaks.bed", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    peak_track(step, peak_track_flag)
    sys.exit(0)