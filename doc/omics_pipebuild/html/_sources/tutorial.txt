.. index:: omics_pipe; tutorial

====================
Omics Pipe Tutorial
====================

Installation
=======

:doc:`Installation instructions <installation>`

Test your installation by typing::

	omics_pipe
	
on the command line. If you get the omics_pipe help readout, it has been installed correctly and you can continue. 

Before Running Omics Pipe: Configuring Parameters File
=============================

.. note::
	Before running omics_pipe, you must configure the parameters file, 
	which is a YAML document.   Follow the instructions here: :doc:`Configuring the parameters file <parameter_file>`

Running Omics Pipe
===================
When you are ready to run omics pipe, simply type the command::

	omics_pipe RNAseq_count_based /path/to/parameter_file.yaml  

To run the basic RNAseq_count_based pipeline with your parameter file. Additional usage instructions below and are available by typing omics_pipe –h.::  

	omics_pipe [-h] [--custom_script_path CUSTOM_SCRIPT_PATH]
                  [--custom_script_name CUSTOM_SCRIPT_NAME]
				  [--compression {gzip, bzip}]
                  {RNAseq_Tuxedo, RNAseq_count_based, RNAseq_cancer_report, RNAseq_TCGA, RNAseq_TCGA_counts, 
				  Tumorseq_MUTECT, miRNAseq_count_based, miRNAseq_tuxedo, WES_GATK, WGS_GATK, SomaticInDels, ChIPseq_MACS, ChIPseq_HOMER,  custom} 
                  parameter_file

If your .fastq files are compressed, please use the compression option and indicate the type of compression used for your files. Currently supported compression types are gzip and bzip.  


Running Omics Pipe with the Test Data and Parameters
====================================
To run Omics Pipe with the test parameter files and data, type the commands below to run each pipeline. 

.. note::
	Replace the ~ with the path to your Omics Pipe installation. 


**RNA-seq (Tuxedo)**::

	omics_pipe RNAseq_Tuxedo ~/tests/test_params_RNAseq_Tuxedo.yaml

**RNA-seq(Anders 2013)**::

	omics_pipe RNAseq_count_based ~/tests/test_params_RNAseq_counts.yaml

**Whole Exome Sequencing (GATK)**::

	omics_pipe WES_GATK ~/tests/test_params_WES_GATK.yaml

**Whole Genome Sequencing (GATK)**::

	omics_pipe WGS_GATK ~/tests/test_params_WGS_GATK.yaml

**Whole Genome Sequencing (MUTECT)**::

	omics_pipe Tumorseq_MUTECT ~/tests/test_params_MUTECT.yaml


**ChIP-seq (MACS)**::

	omics_pipe ChIPseq_MACS ~/tests/test_params_MACS.yaml

**ChIP-seq (HOMER)**::

	omics_pipe ChIPseq_HOMER ~/tests/test_params_HOMER.yaml

**Breast Cancer Personalized Genomics Report- RNAseq**::

	omics_pipe RNAseq_cancer_report ~/tests/test_params_RNAseq_cancer.yaml

**TCGA Reanalysis Pipeline - RNAseq**::

	omics_pipe RNAseq_TCGA ~/tests/test_params_RNAseq_TCGA.yaml

**TCGA Reanalysis Pipeline - RNAseq Counts**::

	omics_pipe RNAseq_TCGA_counts ~/tests/test_params_RNAseq_TCGA_counts.yaml

**miRNAseq Counts (Anders 2013)**::

	omics_pipe miRNAseq_count_based ~/tests/test_params_miRNAseq_counts.yaml

**miRNAseq (Tuxedo)**::

	omics_pipe miRNAseq_tuxedo ~/tests/test_params_miRNAseq_Tuxedo.yaml
	

Running Omics Pipe with your own data
=====================================

1. Copy the test parameter file for the pipeline that you want to run into your home directory::
	
	cp ~/tests/test_params_RNAseq_counts.yaml ~/my_params.yaml


2. Configure the parameter file to point to the path to your data (fastq files), results directories, correct software versions, third party software tool parameters and the correct genome/annotations as described here: :doc:`Configuring the parameters file <parameter_file>`.


3. Ensure that your fastq files follow the naming convention Sample1_1.fastq Sample1_2.fastq for paired end samples. 


4. Type the Omics Pipe command corresponding to your parameter file/pipeline of interest to run the pipeline::

	omics_pipe RNAseq_count_based ~/my_params.yaml
	

5. Omics Pipe will log out to the screen as it is running through the steps in the pipeline. 

	*The pipeline will log out to the screen details regarding the progress of the analysis, including the analysis and status of each step in the pipeline. 
	
	*Individual log files for each job will be located in /data/results/logs (LOG_PATH parameter in parameter file)
	
	*If there are flag files present in the /data/results/flags (FLAG_PATH parameter) folder, these steps in the pipeline will be skipped as they have completed successfully. To redo these steps on the next run of the pipeline, simply delete the flag files and rerun the pipeline. 
	
	*Monitor the progress of the pipeline through the standard output from the Omics Pipe command along with the individual log files for each job to ensure completion. 
 	
	*Each job (step in the pipeline for each sample) will be completed on one of the slave nodes, and Omics Pipe (controller script) will be run on the master node. 
	
	*To check the status of your jobs, type the qstat command. 


6. Wait for the pipeline to finish completely and check the results folders you specified for result files. 

