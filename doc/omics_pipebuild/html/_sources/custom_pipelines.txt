.. index:: omics_pipe; Custom Pipeline

=======================
Omics Pipe Tutorial -- Creating a Custom Pipeline Script
=======================

A pipeline script is a .py file that has the steps in the pipeline that you want to run in your analysis. 
To create a custom pipeline, you will create a new Python script (*.py) file and place it in your working 
directory (or wherever you want).  The available analysis steps built-in to omics_pipe are available in the
(List of currently available omics_pipe analysis steps).  

You may add new modules directly to the module directory (see :doc:`Adding Custom Modules <adding_modules>`), 
and they will become available steps that you can use in your custom pipeline.  

There are three steps to creating a custom pipeline:
	1.	Designing the structure of your pipeline
	2.	Creating the script
	3.	Updating your parameters file

The section below details each of these steps. 

Designing the structure of the pipeline
=========================
Omics_pipe depends upon the pipelining module Ruffus to handle the automation. Please read the documentation
at the Ruffus website http://www.ruffus.org.uk/ for more information.  To design your pipeline, you need to decide
 - What analysis modules you want to run, 
 - What order you want the analysis modules to run in, 
 - Which, if any, of the analysis modules depend upon the results from another analysis module. 
 
For example, we will create a custom pipeline that runs fastqc, star and htseq (depends on output from star).

Creating the script
==========================
To create the script, create a text file name custom_script.py (or whatever name you choose).  
At the top of the file, cut and copy this text: ::

	#!/usr/bin/env python
	from ruffus import *
	import sys 
	import os
	import time
	import datetime 
	import drmaa
	from omics_pipe.utils import *
	from omics_pipe.parameters.default_parameters import default_parameters 
	p = Bunch(default_parameters)
	os.chdir(p.WORKING_DIR)
	now = datetime.datetime.now()
	date = now.strftime("%Y-%m-%d %H:%M")    
	print p
	for step in p.STEPS:
		vars()['inputList_' + step] = []
		for sample in p.SAMPLE_LIST:
			vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (p.FLAG_PATH, step, sample)])
		print vars()['inputList_' + step]

After this has been completed, you will need to import each of the analysis modules that you will use in your pipeline. 
For each analysis module, write the line below (fill in analysis_name with the name of the analysis module): ::

	from omics_pipe.modules.analysis_name import analysis_name

In our example, you will have three lines (see below): ::

	from omics_pipe.modules.fastqc import fastqc
	from omics_pipe.modules.star import star
	from omics_pipe.modules.htseq import htseq

Now you are ready to write the functions to run each of these steps in the analysis.  For each step in our analysis pipeline, 
we will need to write a function.  You can cut and copy these from the pre-packaged analysis pipeline scripts (or eventually
a function reference) or create them.  Each function needs to have two decorators from Ruffus: @parallel(inputList_analysis_name) 
to specify that the pipeline should run in parallel for more than one sample and @check_if_uptodate(check_file_exists) to call a 
function to check if that step in the pipeline is up to date.  Name each function with the name of the analysis prefixed by “run_.”  
The function input should always be (sample, analysis_name_flag).  Within the function, you will call the analysis module that you 
loaded above.  If you want an analysis module to run only after a module it depends upon finishes, you must add the @follows() 
Ruffus decorator before the function, with the name of the step that it depends upon. For example, if htseq needs to run after star, 
you would put @follows(run_star) above the run_htseq function. If you have steps that do not have functions that are dependent 
upon them, you can create a more complex pipeline structure by creating a “Last Function” that ties together all steps of your pipeline. 
The last function below is an example of such a function, and it also produces a PDF diagram of your pipeline when it completes. 
The functions for our example are below. ::


	@parallel(inputList_fastqc)
	@check_if_uptodate(check_file_exists)
	def run_fastqc(sample, fastqc_flag):
		fastqc(sample, fastqc_flag)
		return

	@parallel(inputList_star)
	@check_if_uptodate(check_file_exists)
	def run_star(sample, star_flag):
		star(sample, star_flag)
		return

	@parallel(inputList_htseq)
	@check_if_uptodate(check_file_exists)
	@follows(run_star)
	def run_htseq(sample, htseq_flag):
		htseq(sample, htseq_flag)
		return


	@parallel(inputList_last_function)
	@check_if_uptodate(check_file_exists)
	@follows(run_fastqc, run_htseq)
	def last_function(sample, last_function_flag):
		print "PIPELINE HAS FINISHED SUCCESSFULLY!!! YAY!"
		pipeline_graph_output = p.FLAG_PATH + "/pipeline_" + sample + "_" + str(date) + ".pdf"
		pipeline_printout_graph (pipeline_graph_output,'pdf', step, no_key_legend=False)
		stage = "last_function"
		flag_file = "%s/%s_%s_completed.flag" % (p.FLAG_PATH, stage, sample)
		open(flag_file, 'w').close()
		return  
	
Once you have created all of the functions for each step of your pipeline, cut and copy the code below to the bottom of your script: ::

	if __name__ == '__main__':

		pipeline_run(p.STEP, multiprocess = p.PIPE_MULTIPROCESS, verbose = p.PIPE_VERBOSE, gnu_make_maximal_rebuild_mode = p.PIPE_REBUILD)

At this point, please save your script and move on to step 3.  

Updating your parameters file
==========================
In order for your script to run successfully, you need to configure your parameter file so that each analysis module has the necessary parameters
to execute successfully. The full list of parameters for all modules in the current version of omics_pipe are located in the 
omics_pipe/parameters/default_parameters.py file (and eventually organized somewhere).  You can view the list of necessary parameters for each
analysis module by importing the analysis module into an interactive python session (from omics_pipe.modules.analysis_module import analysis_module) 
and typing analysis_module.__doc__.  The parameters necessary for that analysis module will be listed under “parameters from parameters file.” These
parameters must be put into your parameters.yaml file and spelled exactly as shown (including all caps).   Below is the list of parameters that are 
necessary to run omics_pipe in addition to the module specific parameters. ::

	SAMPLE_LIST: [test, test1]
	STEP: run_last_function
	STEPS: [fastqc, star, htseq, last_function]
	RAW_DATA_DIR: /gpfs/group/sanford/patient/SSKT/test_patient/RNA/RNA_seq/data
	FLAG_PATH: /gpfs/group/sanford/patient/SSKT/test_patient/RNA/RNA_seq/logs/flags
	LOG_PATH: /gpfs/group/sanford/patient/SSKT/test_patient/RNA/RNA_seq/logs
	WORKING_DIR: /gpfs/home/kfisch/virt_env/virt2/lib/python2.6/site-packages/omics_pipe-1.0.7-py2.6.egg/omics_pipe/scripts
	ENDS: PE
	PIPE_MULTIPROCESS: 100
	PIPE_REBUILD: 'True'
	PIPE_VERBOSE: 5
	RESULTS_EMAIL: kfisch@scripps.edu
	TEMP_DIR: /scratch/kfisch
	DPS_VERSION: '1.3.1111'
	QUEUE: bigmem
	PARAMS_FILE: /gpfs/home/kfisch/omics_pipe_docs/test_params.yaml 
	USERNAME: kfisch 
	GENOME: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
	CHROM: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
	REF_GENES: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
	STAR_INDEX: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/star_genome
	BOWTIE_INDEX: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

Once you have all of the necessary parameters in your parameters.yaml file, for your custom script you will need to change the STEP and STEPS parameters.
In the STEP parameter, you will write the name of the last function in your pipeline that you want to run, which should be configured so that it captures all steps
in the pipeline (as in the example above). Make sure to put run_ in front of this, since you are calling the function, not the analysis module.  In order for 
omics_pipe to know what steps you have in your pipeline, you need to list each analysis module name in the STEPS parameter separated with commas 
(without run_ in the prefix).  You are now ready to run your custom script.

Running omics_pipe with a custom pipeline script
When you call the omics_pipe function, you will specify the path to your custom script using the command ::

	omics_pipe custom  --custom_script_path ~/path/to/the/script --custom_script_name name_of_customscript /path/to/parameters.yaml  

This will automatically load your custom script and run through the steps in your pipeline using the default modules available in omics_pipe.
