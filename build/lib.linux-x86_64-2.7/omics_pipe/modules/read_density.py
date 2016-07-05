#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def read_density(sample, read_density_flag):
    '''Runs HOMER to visualize read density from ChIPseq data.
    
    input: 
        .bam file
    output: 
        .txt file
    citation: 
        Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
    link: 
        http://homer.salk.edu/homer/
    parameters from parameters file: 
        BOWTIE_RESULTS:
        
        CHROM_SIZES:
        
        HOMER_RESULTS:
        
        HOMER_VERSION:
        
        TEMP_DIR:
        '''
    spawn_job(jobname = 'read_density', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "15gb", script = "/read_density.sh", args_list = [sample, p.BOWTIE_RESULTS, p.CHROM_SIZES, p.HOMER_RESULTS, p.HOMER_VERSION, p.TEMP_DIR])
    job_status(jobname = 'read_density', resultspath = p.HOMER_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + "_trackInfo.txt", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    read_density(sample, read_density_flag)
    sys.exit(0)