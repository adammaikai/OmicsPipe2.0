#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_Varscan_Annotation(sample, extension, GATK_Varscan_Annotation_flag):
    '''Intersect GATK and Varscan VCF files and annotating with MyVariant.info.
        
        input:
        _vqsr.vcf.gz, _varscan.4.1.vcf.gz
        output:
        _merged.annnotated.vcf.gz
        citation:
        
        link:
        
        parameters from parameters file:
        
        sample:
        
        VARIANT_DIR:
        '''
    sample = sample + extension
    spawn_job(jobname = 'GATK_Varscan_Annotation', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.ANNOTATION["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.ANNOTATION["NODES"], ppn = p.ANNOTATION["CPU"], memory = p.ANNOTATION["MEMORY"], script = "/GATK_Varscan_Annotation.sh", args_list = [sample, p.ANNOTATION["VARIANT_DIR"], p.ANNOTATION["GENOME"], p.VARSCAN["R_VERSION"], p.VARSCAN["VCFLIB_VERSION"], p.VARSCAN["VCFTOOLS_VERSION"], p.ANNOTATION["ANNOTATE_SNPS"], p.ANNOTATION["ANNOTATE_INDELS"]])
    job_status(jobname = 'GATK_Varscan_Annotation', resultspath = p.ANNOTATION["VARIANT_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + "_merged.annotated.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
   GATK_Varscan_Annotation(sample, extension, GATK_Varscan_Annotation_flag)
   sys.exit(0)

