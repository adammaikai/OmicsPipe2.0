#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Samtools_pileup(sample, extension, Samtools_pileup_flag):
    '''Runs Samtools mpileup on tumor/normal .bam files.
        
        input:
        _realigned_sorted.bam
        output:
        .pileup
        citation:
        parameters from parameters file:

        ALIGNMENT_DIR:
        
        SAMTOOLS VERSION:
        
        GENOME:
        
        CAPTURE_KIT_BED:
        '''
    sample = sample + extension
    spawn_job(jobname = 'Samtools_pileup', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.SAMTOOLS["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.SAMTOOLS["NODES"], ppn = p.SAMTOOLS["CPU"], memory = p.SAMTOOLS["MEMORY"], script = "/Samtools_pileup.sh", args_list = [sample, p.SAMTOOLS["ALIGNMENT_DIR"], p.SAMTOOLS["VERSION"], p.SAMTOOLS["GENOME"], p.CAPTURE_KIT_BED, p.SAMTOOLS["SAMBAMBA_VERSION"]])
    job_status(jobname = 'Samtools_pileup', resultspath = p.SAMTOOLS["ALIGNMENT_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + ".flt.pileup", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    Samtools_pileup(sample, extension, Samtools_pileup_flag)
    sys.exit(0)
