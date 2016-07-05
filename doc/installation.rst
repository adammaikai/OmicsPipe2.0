.. index:: omics_pipe; installation

======================================   
Using Omics Pipe
======================================

Omics Pipe is a Python framework for automating 'best practice' next generation sequencing pipelines. 
Omics Pipe can be run from the command-line by providing it with a YAML parameter file specifying your 
directory structure and software specific parameters. This executes a parallel automated pipeline 
on a Distributed Resource Management system (local cluster or Amazon Web Services (AWS)) that efficiently handles job resource allocation, monitoring and 
restarting. The goals of Omics Pipe are to provide researchers with an open-source computational solution to 
implement 'best practice' pipelines with minimal development overhead and providing visual outputs to aid the
researcher in biological interpretation.

To install Omics Pipe, first determine if you are going to be using it on a local compute cluster or on AWS. If you are 
going to be installing it on your local cluster, follow the directions below (or have your system administrator install it globally). 
If you are going to create a local installation in your home directory on your cluster but you do not have administrative permissions, 
you can create a `Python Virtual Environment`_ and then follow the instructions below within the virtual environment.  

Requirements
==========
* HPC Cluster or AWS Star Cluster (:doc:`Resource Requirements <resources>`)

* `Python`_ >=2.6

* `Modules`_

* `R`_ >= 2.15
	*  R Packages (:doc:`Third Party Software Dependencies <dependencies>`)

* Third Party Software Dependencies (:doc:`Third Party Software Dependencies <dependencies>`)

* Reference Databases (:doc:`databases`)

Installation
=============

* **Option 1**: Install from pypi using pip::
    
	pip install omics_pipe

* **Option 2**: Install from pypi using easy_install::
	
	easy_install omics_pipe
	
* **Option 3**: Install from source: Download/extract the source code and run::

	python setup.py install

* **Option 4**: Install the latest code directly from the repository::

	pip install -e hg+https://bitbucket.org/sulab/omics_pipe#egg=omics_pipe

* **Option 5**: If you do not have administrator privileges on your system::

    Step 1: Set up a `Python Virtual Environment`_ 
    Step 2: Use one of the Options (1-4) above to install Omics Pipe within your virtual environment.


Usage
==========
Once you have successfully installed Omics Pipe, you can run a pipeline by typing the command::

	omics_pipe [-h] [--custom_script_path CUSTOM_SCRIPT_PATH]
                  [--custom_script_name CUSTOM_SCRIPT_NAME]
				  [--compression {gzip, bzip}]
                  {RNAseq_Tuxedo, RNAseq_count_based, RNAseq_cancer_report, RNAseq_TCGA, RNAseq_TCGA_counts, Tumorseq_MUTECT, miRNAseq_count_based, miRNAseq_tuxedo, WES_GATK, WGS_GATK, SomaticInDels, ChIPseq_MACS, ChIPseq_HOMER,  custom} 
                  parameter_file


Running Omics Pipe on Amazon Web Services (AWS)
==================================
:doc:`AWS Installation Instructions <AWS_installation>`
	Installation instructions for setting up the AWS Omics Pipe AMI

				  
Tutorial
============
:doc:`Tutorial <tutorial>`
	Step-by-step tutorial for running Omics Pipe

:doc:`Creating a custom pipeline <custom_pipelines>`
	Tutorial for creating and running a custom pipeline in Omics Pipe using existing modules

:doc:`Adding new modules/tools <adding_modules>`	
	Tutorial for adding new modules to Omics Pipe be used in a custom pipeline 
	  
				  
Version history
=============

:doc:`history`	

Documentation
============
The latest copy of this documentation should always be available at:
		`<http://packages.python.org/omics_pipe>`_
	
Questions
===============
Email: kfisch@scripps.edu

Twitter: @kathleenfisch


.. _R: http://cran.r-project.org/
.. _Modules: 	http://modules.sourceforge.net/
.. _Python:  https://www.python.org/
.. _Python Virtual Environment: http://docs.python-guide.org/en/latest/dev/virtualenvs/