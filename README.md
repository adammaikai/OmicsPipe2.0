Omics Pipe - A computational framework for reproducible next generation sequencing analysis
=================================================================================

Summary
---------
Welcome to the documentation for Omics Pipe! Omics pipe is an open-source, modular 
computational platform that automates best practice multi-omics data analysis pipelines 
published in Nature Protocols and other commonly used pipelines, such as GATK.  It currently
automates and provides summary reports for two RNA-seq pipelines, variant calling from whole
exome sequencing (WES), variant calling and copy  number variation analysis from whole genome 
sequencing (WGS), two ChIP-seq pipelines and a custom RNA-seq pipeline for personalized genomic medicine 
reporting.  It also provides automated support for interacting with the The Cancer Genome 
Atlas (TCGA) datasets, including automatic download and processing of the samples in this database.  

![omics_pipe_overview.png](https://bitbucket.org/repo/nKR7eL/images/2606129465-omics_pipe_overview.png)

![omics_pipe_pipelines_20140402.png](https://bitbucket.org/repo/nKR7eL/images/2365251253-omics_pipe_pipelines_20140402.png)

Requirements
------
HPC Cluster or AWS Star Cluster

Python 2.6 or newer

Software Dependencies installed as Modules

Reference Databases

Installation
--------------

Option 1
    
    pip install omics_pipe
    
Option 2

    easy_install omics_pipe

Option 3 -- download/extract the source code

    python setup.py install

Option 4 -- install the latest code directly from the repository

    pip install -e hg+https://bitbucket.org/sulab/omics_pipe#egg=omics_pipe

Option 5 -- if you do not have administrator privileges 

    *Step 1: Set up a [Python virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/)
    *Step 2: Use one of the Options (1-4) above to install Omics Pipe within your virtual environment.

Usage
------
command-line usage::
	
    omics_pipe [-h] [--custom_script_path CUSTOM_SCRIPT_PATH]
                  [--custom_script_name CUSTOM_SCRIPT_NAME]
                  [--compression {gzip,bzip}]
                  {RNAseq_Tuxedo, RNAseq_count_based, RNAseq_cancer_report, RNAseq_TCGA, RNAseq_TCGA_counts, Tumorseq_MUTECT, miRNAseq_count_based, miRNAseq_tuxedo, WES_GATK, WGS_GATK, SomaticInDels, ChIPseq_MACS, ChIPseq_HOMER,  custom} 
                  parameter_file


Documentation
------

[Tutorial for omics_pipe](http://pythonhosted.org/omics_pipe "Tutorial for Omics Pipe")

Website
---------
[Omics Pipe home page](http://sulab.scripps.edu/omicspipe "Omics Pipe Home Page")


Repository
-------------
[Bitbucket Repo](https://bitbucket.org/sulab/omics_pipe "Bitbucket Repo")


Distribution
-------------
[PyPI](https://pypi.python.org/pypi/omics_pipe/1.1.3 "PyPi")

bioRxiv Preprint
-------------
[bioRxiv](http://biorxiv.org/content/early/2014/08/23/008383 "Preprint")

AWS
-------------

Current DB Snapshot: snap-2db684e3 (27 Aug 2014)

Current AMI: ami-2b676b6e


Contact
--------
Feedback welcome at: omics_pipe@googlegroups.com

Twitter:
Kathlessn Fisch: @kathleenfisch
Adam Mark: @AdamMaikai


Google Group: https://groups.google.com/forum/#!forum/omics_pipe
