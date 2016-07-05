#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def snpeff(sample, snpeff_flag):
    '''Annotates variants with SNPeff variant annotator. 
    input: 
        .vcf
    output: 
        .vcf
    citation: 
        "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672
    link: 
        http://snpeff.sourceforge.net
    parameters from parameters file:

    VARIANT_RESULTS:
    
    TEMP_DIR:
    
    WORKING_DIR:
    
    SNPEFF_VERSION:
        '''
    spawn_job(jobname = 'snpeff', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "28gb", script = "/snpeff.sh", args_list = [p.VARIANT_RESULTS, sample, p.TEMP_DIR, p.WORKING_DIR, p.SNPEFF_VERSION])
    job_status(jobname = 'snpeff', resultspath = p.VARIANT_RESULTS, SAMPLE = sample, outputfilename = sample + "/" + sample + "_final_variants_filt_snpeff.vcf", FLAG_PATH = p.FLAG_PATH)
    return

if __name__ == '__main__':
    snpeff(sample, snpeff_flag)
    sys.exit(0)