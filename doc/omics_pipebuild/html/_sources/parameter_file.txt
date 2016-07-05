.. index:: omics_pipe; parameters

====================
Omics Pipe Tutorial -- Configuring the Parameter File
====================

Before running Omics Pipe, you must configure the parameters file, 
which is a `YAML`_ document.  Example parameters files are located within 
the omics_pipe/test folder for each pipeline.  Copy one of these 
parameters files into your working directory, and edit the parameters to 
work with your sample names, directory structure, software options and software versions. 
Make sure to keep the formatting and parameter names exactly the same
as in the example parameters files. 

.. note::
	Make sure to follow the YAML format exactly. Ensure that there is only one space after each colon.
	
.. note::
	For parameters in quotes in the test parameters file, please make sure to keep them in quotes in your custom parameter file. 

The STEP parameter should be the function name of the last step in the pipeline 
that you want to run (e.g. run_tophat). To run the pre-installed pipelines all the 
way through, this should be “last_function.”  

.. warning::
	Do not change the STEPS or STEPS_DE parameters for a pre-installed pipeline.   


.. note::
	Fastq files: paired end: 2 files, “Name_1.fastq” and “Name_2.fastq” representing read 1 and read 2.  
	Have all fastq files in same raw data folder

.. warning::
	Default parameters have been included for each third party software tool included in each of the pipelines. Before running, please view the documentation for each software tool to determine if these parameters are appropriate for your analysis. We do not advise using the default parameters included in Omics Pipe without a full understanding of the tools/parameters. 

Example Omics Pipe Parameter File
=======================

test_params.yaml in omics_pipe/tests::

	SAMPLE_LIST: [test1, test2, test3]
	
	STEP: last_function
	
	STEPS: [fastqc, star, htseq, last_function]

	RAW_DATA_DIR: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests

	FLAG_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/flags

	HTSEQ_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/counts

	LOG_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/logs

	QC_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run

	RESULTS_PATH: /gpfs/home/kfisch/test

	STAR_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/star

	WORKING_DIR: /gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts

	REPORT_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run

	ENDS: SE

	FASTQC_VERSION: '0.10.1'

	GENOME: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

	HTSEQ_OPTIONS: -m intersection-nonempty -s no -t exon

	PIPE_MULTIPROCESS: 100

	PIPE_REBUILD: 'True'

	PIPE_VERBOSE: 5

	REF_GENES: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf

	RESULTS_EMAIL: kfisch@scripps.edu

	STAR_INDEX: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/star_genome

	STAR_OPTIONS: --readFilesCommand cat --runThreadN 8 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical

	STAR_VERSION: '2.3.0'

	TEMP_DIR: /scratch/kfisch

	QUEUE: workq

	USERNAME: kfisch 

	DRMAA_PATH: /opt/applications/pbs-drmaa/current/gnu/lib/libdrmaa.so

	DPS_VERSION: '1.3.1111'

	BAM_FILE_NAME: Aligned.out.bam

	PARAMS_FILE: '/gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_params_RNAseq_counts.yaml'

	DESEQ_META: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/counts_meta.csv

	DESIGN: '~ condition'

	PVAL: '0.05'

	DESEQ_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/DESEQ

	SUMATRA_DB_PATH: /gpfs/home/kfisch/sumatra

	SUMATRA_RUN_NAME: test_counts_sumatra_project

	REPOSITORY: https://kfisch@bitbucket.org/sulab/omics_pipe

	HG_USERNAME: Kathleen Fisch <kfisch@scripps.edu>

	
Explanation of Variables in Omics Pipe Parameter File
=======================================================

Parameters vary by pipeline and the correct parameter file for each pipeline must be used. See examples in the /tests/ folder. 

RNAseq Count Based Pipeline 
----------------------------

test_params_RNAseq_counts.yaml in omics_pipe/tests::

	#sample names ie “Name” for paired and single end reads. So, “Name” for paired-end would expect two fastq files named “Name_1.fastq” and Name_2.fastq”
	SAMPLE_LIST: [test1, test2, test3]	
	
	#Function to be run within pipeline. If you want to run the whole pipeline, leave this as last_function
	STEP: last_function	
	
	#All steps within the pipeline. DO NOT CHANGE this parameter for pre-installed pipelines. If you create your own pipeline, you will need to modify this by listing all of the steps in your pipeline. 
	STEPS: [fastqc, star, htseq, last_function]	

	#Directory where your raw .fastq files are located.
	RAW_DATA_DIR: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests	

	#Directory where you would like to have the flag files created. Flag files are empty files that indicate if a step in the pipeline has completed successfully. 
	FLAG_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/flags	

	#Directory for HTSEQ results
	HTSEQ_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/counts	

	#Directory where log files will be written
	LOG_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/logs	

	#Directory for QC results
	QC_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run	

	#Upper level results directory. Sumatra will check all subfolders of this directory for new files to add to the run tracking database. 
	RESULTS_PATH: /gpfs/home/kfisch/test	

	#Directory where STAR results will be written
	STAR_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/star	
	
	#Where omics_pipe is installed, this path will be pointing to ~/omics_pipe/scripts. 
	WORKING_DIR: /gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts	

	#Directory for the summary report
	REPORT_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run	

	#SE is single end, PE is paired-end sequencing reads
	ENDS: SE	

	#Version number of FASTQC
	FASTQC_VERSION: '0.10.1'	

	#Full path to Genome fasta file
	GENOME: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa	

	#Options for HTSEQ. Please see HTSEQ-count documentation for parameter options. http://www-huber.embl.de/users/anders/HTSeq/doc/count.html#options
	HTSEQ_OPTIONS: -m intersection-nonempty -s no -t exon	

	#Number of multiple processes you want Ruffus to spawn at once
	PIPE_MULTIPROCESS: 100

	#Ruffus parameter. No need to change.
	PIPE_REBUILD: 'True'		

	#Ruffus parameter. No need to change. 
	PIPE_VERBOSE: 5		

	#Full path to reference gene annotations
	REF_GENES: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf	

	#Your email. 
	RESULTS_EMAIL: kfisch@scripps.edu	

	#Directory pointing to STAR_INDEX (you may have to create this)
	STAR_INDEX: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/star_genome	

	#Options for STAR. Please read parameter options https://code.google.com/p/rna-star/
	STAR_OPTIONS: --readFilesCommand cat --runThreadN 8 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical	

	#Version number of STAR
	STAR_VERSION: '2.3.0'	

	#Path to temporary directory
	TEMP_DIR: /scratch/kfisch	

	#Name of the queue on your local cluster you wish to use
	QUEUE: workq	

	#Username for local cluster
	USERNAME: kfisch		

	#Path to your local cluster installation of DRMAA (ask your sys admin for this)
	DRMAA_PATH: /opt/applications/pbs-drmaa/current/gnu/lib/libdrmaa.so		

	#Name of created Bam file. Will be Aligned.out.bam if you are using STAR and accepted_hits.bam if you are using TopHat
	BAM_FILE_NAME: Aligned.out.bam		

	#Full path to your parameter file. Make sure to include the single quotes. 
	PARAMS_FILE: '/gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_params_RNAseq_counts.yaml' 	

	#Location of the meta data csv file for DESEQ. See tests/counts_meta.csv for an example. This file contains the Design file for your study. http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html
	DESEQ_META: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/counts_meta.csv		

	#Design for DESEQ differential expression. Leave as is if you use the exact design as in the counts_meta.csv file. 
	DESIGN: '~ condition'		

	#P-value threshold
	PVAL: '0.05'		

	#Directory for DESEQ results
	DESEQ_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/DESEQ		

	#Directory where you want to store your Sumatra database. Once you run this once, you do not have to change this. 
	SUMATRA_DB_PATH: /gpfs/home/kfisch/sumatra		

	#Name of your project. You do not need to change this for subsequent runs of the pipeline, but you can if you wish. 
	SUMATRA_RUN_NAME: test_counts_sumatra_project		

	#Location of omics pipe repository (you can leave this)
	REPOSITORY: https://kfisch@bitbucket.org/sulab/omics_pipe		

	#Your Mercurial username
	HG_USERNAME: Kathleen Fisch <kfisch@scripps.edu>		
	
	#Version of Python installed on system
	PYTHON_VERSION: 2.6.5
	
	#Type of cluster scheduler (options: PBS, SGE)
	SCHEDULER: PBS
	
	#Full path to WORKING_DIR/reporting on your system
	R_SOURCE_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting
	
	#Are your raw fastq files compressed? If not, leave this parameter as none. If so, please type 'GZIP' and it will automatically process your gzip files. 
	COMPRESSION: none
	
RNAseq Tuxedo Pipeline 
----------------------------

test_params_RNAseq_Tuxedo.yaml in omics_pipe/tests::

	#sample names ie “Name” for paired and single end reads. So, “Name” for paired-end would expect two fastq files named “Name_1.fastq” and Name_2.fastq”
	SAMPLE_LIST: [test1, test2]
	
	#Function to be run within pipeline. If you want to run the whole pipeline, leave this as last_function
	STEP: last_function
	
	#All steps within the pipeline. DO NOT CHANGE this parameter for pre-installed pipelines. If you create your own pipeline, you will need to modify this by listing all of the steps in your pipeline. 
	STEPS: [fastqc, tophat, rseqc, cufflinks]
	
	#All steps within the pipeline. DO NOT CHANGE this parameter for pre-installed pipelines. If you create your own pipeline, you will need to modify this by listing all of the steps in your pipeline. 
	STEPS_DE: [cuffmerge, cuffmergetocompare, cuffdiff, RNAseq_report_tuxedo, last_function]
	
	#SE is single end, PE is paired-end sequencing reads
	ENDS: SE
	
	#Your email address
	RESULTS_EMAIL: kfisch@scripps.edu
	
	#Path to temporary directory (make sure this is a large, writable directory)
	TEMP_DIR: /scratch/kfisch
	
	#Name of the queue on your cluster
	QUEUE: workq
	
	#Your username on the cluster
	USERNAME: kfisch 
	
	#Full path to your raw data files
	RAW_DATA_DIR: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests
	
	#Full path to where you want the Flag files written
	FLAG_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run
	
	#Full path to where you want the log files for each step written
	LOG_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run
	
	#Full path to where you want QC results written
	QC_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run
	
	#Full path to upper level results path
	RESULTS_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run
	
	#Where omics_pipe is installed, this path will be pointing to ~/omics_pipe/scripts
	WORKING_DIR: /gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts
	
	#Full path to where you want Tophat results written
	TOPHAT_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/alignments
	
	#Full path to where you want Cufflinks results written
	CUFFLINKS_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/assemblies
	
	#Full path to where you want Cuffmerge results written
	CUFFMERGE_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/cuffmerge
	
	#Full path to where you want Cuffdiff results written
	CUFFDIFF_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/cuffdiff
	
	#List of full paths to alignment files divided by condition. Each sample will have the path TOPHAT_RESULTS/SAMPLE_NAME/accepted_hits.bam. Divide these up into your two conditions for differential expression analysis. See http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/#cuffdiff-arguments for more details. 
	CUFFDIFF_INPUT_LIST_COND1: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/alignments/test1/accepted_hits.bam
	
	CUFFDIFF_INPUT_LIST_COND2: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/alignments/test2/accepted_hits.bam
	
	#Options for Tophat. Please read the TopHat documentation to set these options for your analysis. http://ccb.jhu.edu/software/tophat/manual.shtml#toph
	TOPHAT_OPTIONS: -p 8 -a 5 --microexon-search --library-type fr-secondstrand
	
	CUFFLINKS_OPTIONS: -u -N  
	
	CUFFMERGE_OPTIONS: -p 8
	
	CUFFMERGETOCOMPARE_OPTIONS: -CG
	
	CUFFDIFF_OPTIONS: -p 8 -FDR 0.01 -L Group1,Group2 -N --compatible-hits-norm
	
	#Software versions installed on your system
	FASTQC_VERSION: '0.10.1'
	TOPHAT_VERSION: '2.0.9'
	CUFFLINKS_VERSION: '2.1.1'
	R_VERSION: '3.0.1'
	BOWTIE_VERSION: 2.2.3
	SAMTOOLS_VERSION: 0.1.19
	PYTHON_VERSION: 2.6.5
	BOWTIE_VERSION: 2.2.3
	SAMTOOLS_VERSION: 0.1.19
	
	#Ruffus specific parameters. See above or documentation for details. http://www.ruffus.org.uk/pipeline_functions.html#pipeline-functions-pipeline-run
	PIPE_MULTIPROCESS: 100
	PIPE_REBUILD: 'True'
	PIPE_VERBOSE: 5
	
	#Full path to gene annotation gtf file
	REF_GENES: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
	
	#Full path to genome file
	GENOME: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
	
	#Full path to BOWTIE index
	BOWTIE_INDEX: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
	
	#Location of chromosomes folder
	CHROM: /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
	
	#Path to your local cluster installation of DRMAA (ask your sys admin for this)
	DRMAA_PATH: /opt/applications/pbs-drmaa/current/gnu/lib/libdrmaa.so
	
	#Full path to directory where you want report results written
	REPORT_RESULTS: /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests
	
	#Full path to your parameter file. Make sure to include the single quotes. 
	PARAMS_FILE: '/gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_params_RNAseq_Tuxedo.yaml'
	
	#Gene IDs of interest for visualization with CummeRbund
	GENEIDS: [GAPDH, COL2A1, BRCA2]
	
	#Name of samples in Condition 1
	COND1: Group1
	
	#Name of samples in Condition 2
	COND2: Group2
	
	#Name of your project. You do not need to change this for subsequent runs of the pipeline, but you can if you wish.
	NAME: Test_run_date
	
	#Directory where you want to store your Sumatra database. Once you run this once, you do not have to change this. 
	SUMATRA_DB_PATH: /gpfs/home/kfisch/sumatra
	
	#Name of your project. You do not need to change this for subsequent runs of the pipeline, but you can if you wish. 
	SUMATRA_RUN_NAME: default_sumatra_project
	
	#Location of omics pipe repository (you can leave this)
	REPOSITORY: https://kfisch@bitbucket.org/sulab/omics_pipe
	
	#Your Mercurial username
	HG_USERNAME: Kathleen Fisch <kfisch@scripps.edu>
	
	#Type of cluster scheduler (options: PBS, SGE)
	SCHEDULER: PBS
	
	#Full path to WORKING_DIR/reporting on your system
	R_SOURCE_PATH: /gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting
	
	#Are your raw fastq files compressed? If not, leave this parameter as none. If so, please type 'GZIP' and it will automatically process your gzip files. 
	COMPRESSION: none

.. _`YAML`: http://www.yaml.org/