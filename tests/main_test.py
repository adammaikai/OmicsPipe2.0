#!/usr/bin/env python

#first register command
#read analysis type
#from sumatra.decorators import capture
import os
import sys
import urllib
from optparse import OptionParser
import webbrowser as browser
import argparse
from omics_pipe.parameters.default_parameters import default_parameters 
import yaml
from omics_pipe.utils import *
import runpy

#@capture  #sumatra capture. Might not need if running from command line within script
def main():
    '''make command line interface, read in analysis type, execute correct analysis pipeline script'''
    parser = argparse.ArgumentParser(prog = 'omics_pipe', description = 'Run omics_pipe')
    parser.add_argument('analysis_type', action = "store",  choices = ['RNAseq', 'RNAseq_sanford', 'miRNAseq', 'ChIPseq', 'GWAS', 'modular', 'custom'], help = 'type of analysis to run: RNAseq, ChIPseq, GWAS, modular, custom')
    parser.add_argument('parameter_file', action = "store", help = 'specify parameter file to use for analysis')
    parser.add_argument('--custom_script_path', action = "store", help = 'specify custom script file with full path (/example/script.py) to use for analysis if you specify analysis type as custom')
    parser.add_argument('--custom_script_name', action = "store", help = 'specify custom script file with full path (/example/script.py) to use for analysis if you specify analysis type as custom')

    args = parser.parse_args()
    
    print args
    print args.analysis_type
    print args.custom_script_path
    print args.custom_script_name
    print args.parameter_file
    
    analysis = args.analysis_type     
    parameters = args.parameter_file
    path = args.custom_script_path
    script = args.custom_script_name
    
    stream = file(parameters, 'r')
    params = yaml.load(stream)
        
    default_parameters.update(params)   #Update default parameters to user defined parameter file
    p = Bunch(default_parameters)
       
    check_create_dir(p.LOG_PATH)
    check_create_dir(p.FLAG_PATH)
#   check if sumatra project is present, if not, create one
  
    if args.analysis_type == 'RNAseq':
        runpy.run_module('omics_pipe.RNAseq', run_name="__main__", alter_sys = True)
        sys.exit(0)
    elif args.analysis_type == 'RNAseq_sanford':
        runpy.run_module('omics_pipe.RNAseq_sanford', run_name="__main__", alter_sys = True)
        sys.exit(0)
    elif args.analysis_type == 'miRNAseq':
        runpy.run_module('omics_pipe.miRNAseq', run_name="__main__", alter_sys = True)
        sys.exit(0)
    elif args.analysis_type == 'ChIPseq':
        runpy.run_module('omics_pipe.ChIPseq', run_name="__main__", alter_sys = True)
        sys.exit(0)
    elif args.analysis_type == 'GWAS':
        runpy.run_module('omics_pipe.GWAS', run_name="__main__", alter_sys = True)
        sys.exit(0)
    elif args.analysis_type == 'custom':
        os.chdir(path)
        runpy.run_module(script, run_name="__main__", alter_sys = True)
        sys.exit(0)   
    else:
        print 'Error: unsupported analysis type'       
    return




    
if __name__ == '__main__':
    main()

