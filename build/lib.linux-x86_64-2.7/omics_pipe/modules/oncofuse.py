#!/usr/bin/env python
# -*- coding: utf-8 -*-
from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def oncofuse(sample, oncofuse_flag):
    '''Predicts oncogenic potential of fusion genes. 
    
    input: 
        list of candidate fusion genes
    output: 
        list of candidate fusion genes with oncogenic potential
    citation: 
        Mikhail Shugay, Inigo Ortiz de Mend�bil, Jose L. Vizmanos and Francisco J. Novo. Oncofuse: a computational framework for the prediction of the oncogenic potential of gene fusions. Bioinformatics. 16 Aug 2013. doi:10.1093/bioinformatics/btt445.     
    link: 
        http://www.unav.es/genetica/oncofuse.html
    parameters from parameters file: 
        
        FUSION_RESULTS:
              
        TEMP_DIR:
        
        ONCOFUSE_VERSION:

        '''
    if p.ENDS == "PE":
        spawn_job(jobname = 'oncofuse', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, SCHEDULER = p.SCHEDULER, walltime = "240:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "2gb", script = "/oncofuse.sh", args_list = [sample, p.FUSION_RESULTS, p.TEMP_DIR, p.ONCOFUSE_VERSION])
        job_status(jobname = 'oncofuse', resultspath = p.FUSION_RESULTS, SAMPLE = sample, outputfilename = sample + "/oncofuse_res.txt", FLAG_PATH = p.FLAG_PATH)
    else:
        print "Oncofuse can only be run on paired-end samples. Input samples are single-ended. Skipping step."
        flag_file = "%s/oncofuse_completed.flag" % p.FLAG_PATH
        open(flag_file, 'w').close() 
    return

if __name__ == '__main__':
    oncofuse(sample, oncofuse_flag)
    sys.exit(0)