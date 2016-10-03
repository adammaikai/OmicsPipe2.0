#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Star(sample, extension, Star_flag):
    '''Runs Star aligner for RNA-seq.
        
        input:
        fastq.gz
        output:
        _Aligned.sortedByCoord.out.bam
        citation:
        parameters from parameters file:

        TEMP_DIR:
        
        RAW_DATA_DIR:
        
        ALIGNMENT_DIR:
        
        STAR_VERSION:
        
        GENOME:
        '''
    
    if extension is None:
        sample = sample
    else:
        sample = sample + extension
    spawn_job(jobname = 'Star', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.STAR["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.STAR["NODES"], ppn = p.STAR["CPU"], memory = p.STAR["MEMORY"], script = "/Star.sh", args_list = [sample, p.OMICSPIPE["TEMP_DIR"], p.STAR["RAW_DATA"], p.STAR["RESULTS"], p.STAR["VERSION"], p.STAR["GENOME"], p.STAR["REFGTF"], p.STAR["STARFUSION_LIB"], p.STAR["STARFUSION_VERSION"], p.STAR["SAMBAMBA_VERSION"]])
    job_status(jobname = 'Star', resultspath = p.STAR["RESULTS"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + "_Aligned.sortedByCoord.out.bam", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    Star(sample, extension, Star_flag)
    sys.exit(0)
