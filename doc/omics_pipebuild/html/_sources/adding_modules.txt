.. index:: omics_pipe; Adding Modules

==========================================
Omics Pipe Tutorial -- Adding a New Module (Tool)
==========================================

Users can easily create new analysis modules for use within omics_pipe.  
The user has two options for creating new analysis modules: 
- Adding analysis modules directly within the omics_pipe/scripts installation directory 
- Creating a new working directory where all analysis modules scripts are located 
(this can be changed in the parameters file by changing the WORKING_DIR parameter 
to the desired location).  If you want to use option 2, in order to use pre-installed analysis 
modules, for the time being you must copy these analysis modules to your new working 
directory.  If you choose option 1, you can simply add additional analysis modules and they
will be accessible along with the pre-installed analysis modules.  

To create a new analysis module, you need to perform four steps: 
1.	Create a Bash script with the command to be sent to the cluster
2.	Create a Python module that calls the Bash script 
3.	Add  your module to your custom pipeline
4.	Add new module parameters to parameters file

The section below details each of these steps. 

1.	Create a Bash script
===================

The first step in creating your custom module is to create the Bash script with the command you 
would like to run. If you are unsure how to write a Bash script, you can look at the examples in 
omics_pipe/scripts or work through this tutorial (http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html).   
In many cases, this will be a simple script with a one line command to call the analysis program.  You
should name your script something that will be easily identifiable and it should have the 
suffix .sh (e.g. analysis_script.sh).  At the beginning of your analysis script, you should put the following lines: ::

	#!/bin/bash
	set -x
	#Source modules for current shell
	source $MODULESHOME/init/bash
	#Make output directory if it doesn't exist
	mkdir -p ${variable} #RESULTS_DIR
	#Move tmp dir to scratch 
	export TMPDIR=${variable} #TEMP_DIR
	#Load specified software version 
	module load fastqc/${variable} #VERSION

The ${variable} will be changed to ${number} (e.g. $1) based on the location of the variable in the input
script (more on this below).   These settings are assuming you are working on a cluster with a modular 
structure. If not, “module load” may not be appropriate to load the software, so please ask your system 
administrator to provide assistance with this if your cluster has a different system.  After you specify the 
software and other configuration variables, you can write the commands for the software you would like
to use. When you are finished with the commands, exit the script with ‘exit 0.’  An example script for running
the software program FASTQC is below. ::

	#Runs fastqc with $1=SAMPLE, $2=RAW_DATA_DIR, $3=QC_PATH
	fastqc -o $3 $2/$1.fastq

	exit 0

Substitute all variables that you would like to change from the parameter file with a variable notation, in the 
form of $1, $2, $3, etc for the first, second, third, etc input parameter that will be passed to the script.  Once 
you have appropriately parameterized the script, save the script either in your working directory (along will all 
the other scripts you will need, possibly copied from omics_pipe/scripts) or in the omics_pipe/scripts directory.  

2.	Create a Python module
===================

Now that you have created your custom script, you can create the Python module that will handle that script and 
schedule a job on the compute cluster using DRMAA (https://code.google.com/p/drmaa-python/wiki/Tutorial).  
You should name the Python module the same name as your custom analysis module, but with the extension .py.  
In this example, your Python module would be named analysis_script.py and the function within it would also be 
called analysis_script. Save your custom Python module within the same directory as your custom pipeline script.  
At the top of your Python module, cut and copy the text below.  ::

	#!/usr/bin/env python

	import drmaa
	from omics_pipe.parameters.default_parameters import default_parameters
	from omics_pipe.utils import *
	p = Bunch(default_parameters)

	You will then write a simple Python function that take the form of the function below. You can directly cut and copy 
	this function and then change the necessary names/parameters to fit your custom analysis.  ::

	def fastqc(sample, fastqc_flag):
		'''QC check of raw .fastq files using FASTQC
			input: .fastq file
			output: folder and zipped folder containing html, txt and image files
			citation: Babraham Bioinformatics
			link: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
			parameters from parameters file: RAW_DATA_DIR,QC_PATH, FASTQC_VERSION''' 
			  spawn_job(jobname = 'fastqc', SAMPLE = sample, LOG_PATH = p.LOG_PATH, RESULTS_EMAIL = p.RESULTS_EMAIL, walltime = "12:00:00", queue = p.QUEUE, nodes = 1, ppn = 8, memory = "16gb", script = "/fastqc_drmaa.sh", args_list = [sample, p.RAW_DATA_DIR,p.QC_PATH, p.FASTQC_VERSION])
			job_status(jobname = 'fastqc', resultspath = p.QC_PATH, SAMPLE = sample, outputfilename = sample + "_fastq/" + "fastqc_data.txt", FLAG_PATH = p.FLAG_PATH)
		return


Name your function the same as the names of both the Bash and Python scripts you just created for consistency. 
In our example, the first line would look like: “def analysis_script(sample, analysis_script_flag):”.  As you can see,
I changed the name of the function as well as the name of the flag input file.  The document string should be change
to describe what your analysis module does, what type of input file it takes, a citation and link to the tool that you are
calling, as well as the parameters that are needed in the parameters file that will be passed to the Bash script that you 
created.  After you are done documenting your function, you will change a few items within the spawn_job and job_status 
functions that are called from the omics_pipe.utils module.  In the spawn_job function, you should change the job name to 
match the name of your function, you can customize the resources your job will request from the compute cluster, you will
need to change the name of the script to match that of the Bash script that you just created, and then you will change the 
parameters listed in the variable “args_list.”  The variable “sample” is lower case because it is passed to this function from
omics_pipe, but input parameters coming from the parameters file must be prefixed with “p.” List the parameters that you 
need to feed into your custom analysis script in the order that you numbered them in the Bash script. In the example above,
$1 corresponds to ‘sample’ $2 corresponds to p.RAW_DATA_DIR, etc.  Once you have the spawn_job function updated, you
will update the job_status function with the job name, results path and a name of an output file that will be produced from your 
Bash script. This can be any file that is created.  This function will check that this file exists in the specified results directory, 
check that its size is greater than zero, and then it will create a flag file if it exists.  Once you complete this, you are finished
creating your custom Python module.  


3.	Add custom Python module to your custom pipeline
======================================

In order to use your custom analysis module, you will need to create a custom pipeline with your custom analysis module included
as a step in the pipeline. For a tutorial on how to create a custom pipeline, see Section “Creating a Custom Pipeline Script.”  Once 
you have a custom pipeline script, please make sure your custom analysis module and custom pipeline script are in the same directory. 

4.	Add new parameters to parameters file
==============================
Add the parameters necessary for your custom analysis module to run into the 
parameters file.  Simply add the parameters to your parameters script, save it, and then run your custom pipeline.

5.	Create Custom Module Directory
==============================
In order for Omics Pipe to find your custom modules and scripts, create a local directory and put your new scripts, module files and custom pipeline within this folder. In your parameter file,
specify this directory as the WORKING_DIR: parameter. In your custom pipeline script, you can import custom modules at the beginning of the script that are within the same folder as the custom pipeline. 
For Omics Pipe modules, import those modules as they appear in the Omics Pipe pipeline scripts. 

6. Run your pipeline
======================
Run your custom modules/pipeline as you would from the command line for running a custom pipeline. 