#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Qualimap(sample, extension, Qualimap_flag):
    '''Runs Qualimap on a processed .bam file.
        
        input:
        _gatk_recal.bam
        output:
    
        citation:
        parameters from parameters file:

        ALIGNMENT_DIR:
        
        RESULTS:
        
        QUALIMAP_VERSION:

        QUALIMAP_OPTIONS:
        '''
    if "DNA" in extension:
        options = p.QUALIMAP["DNA_OPTIONS"]
        sample_type = "DNA"
    else:
        options = p.QUALIMAP["RNA_OPTIONS"]
        sample_type = "RNA"

    sample = sample + extension
    spawn_job(jobname = 'Qualimap', SAMPLE =  sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.QUALIMAP["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.QUALIMAP["NODES"], ppn = p.QUALIMAP["CPU"], memory = p.QUALIMAP["MEMORY"], script = "/Qualimap.sh", args_list = [sample, p.QUALIMAP["ALIGNMENT_DIR"], p.QUALIMAP["RESULTS"], p.QUALIMAP["VERSION"], sample_type, options, p.CAPTURE_KIT_BED])
    job_status(jobname = 'Qualimap', resultspath = p.QUALIMAP["RESULTS"] + "/" + sample, SAMPLE = sample,  outputfilename = "genome_results.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    Qualimap(sample, extension, Qualimap_flag)
    sys.exit(0)
