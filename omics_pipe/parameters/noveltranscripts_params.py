#!/usr/bin/env python
#PATHS
WORKING_DIR = "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts"   ###PATH TO ../omics-pipeline/omics_pipe/scripts
RAW_DATA_DIR = "/gpfs/group/su/kfisch/OA/data/rawdata/NormalvsOA20130808"
RESULTS_PATH = "/gpfs/group/su/kfisch/OA/results"
QC_PATH = "/gpfs/group/su/kfisch/OA/data/QC"
LOG_PATH = "/gpfs/group/su/kfisch/OA/logs/ncrna"
TOPHAT_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/ncRNA/alignments"
STAR_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/STAR_alignments"
HTSEQ_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/counts"
CUFFLINKS_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/ncRNA/assemblies"
CUFFMERGE_RESULTS = "/gpfs/group/su/kfisch/Bodymap/results/assemblies/cuffmerge"
CUFFDIFF_RESULTS = "/gpfs/group/su/kfisch/Bodymap/results/cuffdiff"
TEMP_DIR = "/scratch/kfisch"


RESULTS_EMAIL = "kfisch@scripps.edu"

#SAMPLE INFO
ENDS = "SE"   #PE=paired-ends, SE=single ends  If paired end, sample file names must be samplename_1.fastq and samplename_2.fastq
#SAMPLE = "test"

# SAMPLE_LIST = [
# 				'OA_Cart_DZ_1_9',
# 				'OA_Cart_DZ_2_10',
# 				'OA_Cart_DZ_3_11',
# 				'OA_Cart_DZ_4_12',
# 				'OA_Cart_MZ_1_5',
# 				'OA_Cart_MZ_2_6',
# 				'OA_Cart_MZ_3_7',
# 				'OA_Cart_MZ_4_8',
# 				'OA_Cart_SZ_1_1',
# 				'OA_Cart_SZ_2_2',
# 				'OA_Cart_SZ_3_3',
# 				'OA_Cart_SZ_4_4',
# 				
# 				]

SAMPLE_LIST =[
  'Normal_Cart_10_8',
 	'Normal_Cart_7_3',
 	'Normal_Cart_9_7',
 	'OA_Cart_10_9',
 	'OA_Cart_6_1',
 	'OA_Cart_7_2',
 	'OA_Cart_8_5',
 	'OA_Cart_9_6',
 	'Normal_Cart_2_2',
	'Normal_Cart_3_3',
	'Normal_Cart_4_4',
	'Normal_Cart_5_5',
	'Normal_Cart_6_6',
	'OA_Cart_1_7',
	'OA_Cart_2_8',
	'OA_Cart_3_9',
	'OA_Cart_4_10',
	'OA_Cart_5_5'
 	
 	]


#GENOME INFO
GENOME = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
REF_GENES = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
CHROM = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes"
BOWTIE_INDEX = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
STAR_INDEX = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/star_genome"


LNCRNA_GTF = "/gpfs/group/su/kfisch/references/broadinstitute/lincRNAs_transcripts.gtf"
NONCODE_BOWTIE_INDEX = "/gpfs/group/su/kfisch/references/noncode/human/noncode"
NONCODE_FASTA = "/gpfs/group/su/kfisch/references/noncode/human/noncode.fa"


#SOFTWARE VERSIONS
FASTQC_VERSION = "0.10.1"
TOPHAT_VERSION = "2.0.9"
STAR_VERSION = "2.3.0"
CUFFLINKS_VERSION = "2.1.1"
R_VERSION = "3.0.1"
SAMTOOLS_VERSION = "0.1.19"

#SOFTWARE OPTIONS
TOPHAT_OPTIONS = "-p 8 -a 5 --microexon-search --library-type fr-secondstrand"  #-p8 -G options hardcoded in script. These options go after.
STAR_OPTIONS = "--readFilesCommand cat --runThreadN 8 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"
HTSEQ_OPTIONS = "-m intersection-nonempty -s no -t exon"
CUFFLINKS_OPTIONS = "-u -N" #'-p8 -g REF_GENES -b GENOME' options hardcoded in script
CUFFMERGE_OPTIONS = ""
CUFFDIFF_OPTIONS = ""


#STEP TO RUN IN PIPELINE
#			fastqc
#tophat				STAR
#cufflinks			htseqcount
#		last_function
STEP = "last_function" #enclose name in brackets (eg [last_function]) if you only want to rerun that step in the pipeline. If you want to run between certain steps in the pipeline, write steps enclosed in brackets separated by a comma (eg [STAR], [htseqcount])

PIPE_VERBOSE = 5
PIPE_MULTIPROCESS = 16  #sample number * number parallel tasks in pipeline
PIPE_REBUILD = "True" #Add gnu_make_maximal_rebuild_mode = False if you want to run only one task


