#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def MuTect(sample, tumor, normal, MuTect_flag):
    '''Runs MuTect on paired tumor/normal samples to detect somatic point mutations in cancer genomes.
        
        input:
        _gatk_recal.bam
        output:
        _mutect.filt.vcf.gz
        citation:
        Cibulskis, K. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnology (2013).doi:10.1038/nbt.2514
        link:
        http://www.broadinstitute.org/cancer/cga/mutect
        parameters from parameters file:
        
        TEMP_DIR:
        
        GATK_VERSION:
        
        GENOME:
        
        DBSNP:
        
        COSMIC:
        
        CAPTURE_KIT_BED:
        '''
    
    spawn_job(jobname = 'MuTect', SAMPLE = sample + tumor + normal, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"],  walltime = p.MUTECT["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.MUTECT["NODES"], ppn = p.MUTECT["CPU"], memory = p.MUTECT["MEMORY"], script = "/MuTect.sh", args_list = [sample, p.MUTECT["RESULTS"], p.MUTECT["COSMIC"], p.BQSR["DBSNP"], p.CAPTURE_KIT_BED, p.BQSR["ALIGNMENT_DIR"], p.SAMTOOLS["GENOME"], p.MUTECT["VERSION"], normal, tumor, p.MUTECT["OPTIONS"]])
    job_status(jobname = 'MuTect', resultspath = p.MUTECT["RESULTS"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + tumor + normal + "mutect.filt.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    MuTect(sample, tumor, normal, MuTect_flag)
    sys.exit(0)
