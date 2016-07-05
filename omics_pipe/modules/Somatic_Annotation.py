#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Somatic_Annotation(sample, tumor, normal, Somatic_Annotation_flag):
    '''Intersect MuTect and Varscan Somatic VCF files and annotating with MyVariant.info.
        
        input:
        _mutect.filt.vcf.gz, _varscan_somatic.4.1.vcf.gz
        output:
        _somatic_annotations.txt
        citation:
        
        link:
        
        parameters from parameters file:
        
        sample:
        
        VARIANT_DIR:
        '''
    
    spawn_job(jobname = 'Somatic_Annotation', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.ANNOTATION["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.ANNOTATION["NODES"], ppn = p.ANNOTATION["CPU"], memory = p.ANNOTATION["MEMORY"], script = "/Somatic_Annotation.sh", args_list = [sample, p.ANNOTATION["VARIANT_DIR"], p.VARSCAN["R_VERSION"], tumor, normal])
    job_status(jobname = 'Somatic_Annotation', resultspath = p.ANNOTATION["VARIANT_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + "_somatic_annotations.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
   Somatic_Annotation(sample, tumor, normal, Somatic_Annotation_flag)
   sys.exit(0)

