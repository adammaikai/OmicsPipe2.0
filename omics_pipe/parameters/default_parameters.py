#!/usr/bin/env python

default_parameters = dict(                                                                                   #MAKE SURE TO KEEP dict() and , after each variable
#PATHS
WORKING_DIR = "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts"   ,   ###PATH TO ../omics-pipeline/omics_pipe/scripts
RAW_DATA_DIR = "/gpfs/group/su/kfisch/OA/data/rawdata/NormalvsOA20130808"           ,   
RESULTS_PATH = "/gpfs/group/su/kfisch/OA/results"                                   ,
QC_PATH = "/gpfs/group/su/kfisch/OA/data/QC"                                        ,
LOG_PATH = "/gpfs/group/su/kfisch/OA/logs"                                          ,
FLAG_PATH = "/gpfs/group/su/kfisch/OA/logs"                                         ,
TOPHAT_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/alignments"                 ,
STAR_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/STAR_alignments"              ,
HTSEQ_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/counts"                      ,
HTSEQ_GENCODE_RESULTS = "",
CUFFLINKS_RESULTS = "/gpfs/group/su/kfisch/OA/results/mRNA/assemblies"              ,
CUFFMERGE_RESULTS = "/gpfs/group/su/kfisch/Bodymap/results/assemblies/cuffmerge"    ,
CUFFDIFF_RESULTS = "/gpfs/group/su/kfisch/Bodymap/results/cuffdiff"                 ,
TEMP_DIR = "/scratch/kfisch"                                                        ,
VARIANT_RESULTS = "/gpfs/group/su/kfisch/OA/results/variants"                       ,
FUSIONCATCHER_RESULTS = "/gpfs/group/su/kfisch/OA/results/fusions"                  ,
TOPHAT_FUSION_RESULTS = "/gpfs/group/su/kfisch/OA/results/fusions"                  ,


RESULTS_EMAIL = "kfisch@scripps.edu"                                                ,
QUEUE = "workq"                                                                   ,
SCHEDULER = "PBS",

#SAMPLE INFO
ENDS = "SE"                                                                         ,     #PE=paired-ends, SE=single ends  If paired end, sample file names must be samplename_1.fastq and samplename_2.fastq
COMPRESSION = "None",
SAMPLE_LIST =[
    'Normal_Cart_2_2',
    'Normal_Cart_3_3',
    'Normal_Cart_4_4',
    'Normal_Cart_5_5',
    'Normal_Cart_6_6',
    'OA_Cart_1_7',
    'OA_Cart_2_8',
    'OA_Cart_3_9',
    'OA_Cart_4_10',
    'OA_Cart_5_5',
    'Normal_Cart_10_8',
     'Normal_Cart_7_3',
     'Normal_Cart_9_7',
     'OA_Cart_10_9',
     'OA_Cart_6_1',
     'OA_Cart_7_2',
     'OA_Cart_8_5',
     'OA_Cart_9_6']                                                                 ,



#GENOME INFO
GENOME = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"     ,
REF_GENES = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"           ,
CHROM = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes"                     ,
BOWTIE_INDEX = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"      ,
STAR_INDEX = "/gpfs/group/databases/Homo_sapiens/UCSC/hg19/star_genome"                         ,
REF_GENES_GENCODE = "",    
BWA_INDEX = "",
#SOFTWARE VERSIONS
FASTQC_VERSION = "0.10.1"                                                           ,
TOPHAT_VERSION = "2.0.9"                                                            ,
STAR_VERSION = "2.3.0"                                                              ,
CUFFLINKS_VERSION = "2.1.1"                                                         ,            
R_VERSION = "3.0.1"                                                                 ,
SAMTOOLS_VERSION = "0.1.19"                                                         ,    
ANNOVAR_VERSION = "07292013"                                                        ,
VCFTOOLS_VERSION = "0.1.10"                                                         ,
VARSCAN_VERSION = "2.3.6"                                                           ,
FUSIONCATCHER_VERSION = "0.98"                                                      ,   
TRIMGALORE_VERSION = "0.3.3",

#SOFTWARE OPTIONS
TOPHAT_OPTIONS = "-a 5 --microexon-search --library-type fr-secondstrand --no-discordant"                                               ,    #-p8 -G options hardcoded in script. These options go after.
STAR_OPTIONS = "--readFilesCommand cat --runThreadN 8 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"       ,  
HTSEQ_OPTIONS = "-m intersection-nonempty -s no -t exon"                                                                                ,
CUFFLINKS_OPTIONS = "-u -N"                                                                                                             ,   #'-p8 -g REF_GENES -b GENOME' options hardcoded in script

#STEP TO RUN IN PIPELINE
#            fastqc
#tophat                STAR
#cufflinks            htseqcount
#        last_function
STEP = "last_function"                                                                                                                  , #enclose name in brackets (eg [last_function]) if you only want to rerun that step in the pipeline. If you want to run between certain steps in the pipeline, write steps enclosed in brackets separated by a comma (eg [STAR], [htseqcount])
STEPS = [
        'fastqc',
        'tophat', 
        'star', 
        'htseq',
        'cufflinks',
        'fusion_catcher',
        'call_variants',
        'annotate_variants',
        'last_function'
        ]                                                                                                                               ,    



PIPE_VERBOSE = 5                                                                                                                        ,
PIPE_MULTIPROCESS = 100                                                                                                                 ,   #sample number * number parallel tasks in pipeline
PIPE_REBUILD = "True"                                                                                                                   ,   #Add gnu_make_maximal_rebuild_mode = False if you want to run only one task

#Variant calling options
ANNOVARDB = "/gpfs/group/databases/annovar/humandb"                                                                                                                                     ,
VARSCAN_PATH = "/opt/applications/varscan/2.3.6/VarScan.jar"                                                                                                                            ,
VARSCAN_OPTIONS = "--min-var-freq 0.5 --min-avg-qual 30 --p-value 0.995 --output-vcf 1"                                                                                                 ,
ANNOVAR_OPTIONS = "-buildver hg19 -protocol nonsyn_splicing,1000g2012apr_all,esp6500_ea,esp6500_aa,snp135NonFlagged,cg46,ljb_sift,ljb_pp2,dominant -operation g,f,f,f,f,f,f,f,m"        ,
ANNOVAR_OPTIONS2 = "-genetype knowngene --remove"                                                                                                                                       ,
SAMTOOLS_OPTIONS = "-C 500",
#Fusion Catcher options
FUSIONCATCHERBUILD_DIR = "/gpfs/group/databases/ensembl_v72"                                                                                                                            ,
FUSIONCATCHER_OPTIONS = ""                                                                                                                                                              ,

REPORT_SCRIPT = "/gpfs/group/sanford/src/RNA/knitMeR.R",
REPORT_RESULTS = "/gpfs/group/sanford/patient/SSKT/test_patient/RNA/RNA_seq",
R_MARKUP_FILE = "/gpfs/group/sanford/src/RNA/sanfordRNASeqReport.Rmd",
BAM_FILE_NAME = "accepted_hits.bam",
R_SOURCE_PATH = "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting",

#Trim galore options
ADAPTER = "TGGAATTCTCGGGTGCCAAGG",
TRIM_LENGTH_MIN = 10,
TRIMMED_DATA_PATH = "/gpfs/home/kfisch",

miRNA_BOWTIE_INDEX = "/gpfs/group/su/kfisch/references/bowtieindexes/miRNA",
miRNA_GTF=  "/gpfs/group/su/kfisch/references/miRBase/hsa.gff3",
CUFFMERGE_OPTIONS=  "-p 8",
CUFFMERGETOCOMPARE_OPTIONS = "-CG",
CUFFDIFF_OPTIONS=  "-p 8 -FDR 0.01 -L Normal,OA -N --compatible-hits-norm",
CUFFDIFF_INPUT_LIST_COND1= "/gpfs/group/su/kfisch/OA/results/miRNA/bamfiles/sample1.bam,/gpfs/group/su/kfisch/OA/results/miRNA/bamfiles/sample2.bam" ,
CUFFDIFF_INPUT_LIST_COND2 = "/gpfs/group/su/kfisch/OA/results/miRNA/bamfiles/sample3.bam",
#TCGA download
TCGA_XML_FILE = "/gpfs/group/sanford/TCGA/tcga_test.xml",
TCGA_KEY = "/gpfs/group/sanford/TCGA/cghub.key",
TCGA_OUTPUT_PATH = "/gpfs/group/sanford/TCGA",
SSH_USER_NAME = "kfisch@hpcdata.scripps.edu",

#BWA and SNPIR
SNPIR_RESULTS = "",
BWA_VERSION ="0.7.4",
PICARD_VERSION ="1.92",
GATK_VERSION ="3.1-1",
BEDTOOLS_VERSION ="2.17.0",
UCSC_TOOLS_VERSION ="273",
REPEAT_MASKER = "",
SNPIR_ANNOTATION ="",
RNA_EDIT  ="",
DBSNP  ="",
MILLS  ="",
G1000  ="",
BWA_RESULTS ="",
GATK_READ_GROUP_INFO ="",
SNPIR_VERSION = "1.0",
PARAMS_FILE = "",
VCF_FILE= "/gpfs/group/sanford/patient/SSKT/SSKT_2/RNA/RNA_seq/results/SNPIR/final_variants.vcf",
INTOGEN_OPTIONS= "--single-tumor ",
INTOGEN_RESULTS= "/gpfs/group/sanford/patient/SSKT/SSKT_2/RNA/RNA_seq/results/intogen",
INTOGEN_VERSION= '2.3.0',
INTOGEN_CONFIG= "/gpfs/group/sanford/patient/SSKT/SSKT_2/RNA/RNA_seq/results/intogen.conf",
USERNAME= "kfisch",
SNPIR_CONFIG = "/gpfs/group/su/meissto/cancer_report/src/RNA",
SNPIR_DIR = "/opt/applications/snpir/1.0/bin",
RSEQC_VERSION = "2.3.7",
RSEQC_REF = "/gpfs/group/su/meissto/data/knowngene.bed",
SNPEFF_VERSION = "3.5",
dbNSFP = "/gpfs/group/su/meissto/data/dbNSFP/dbNSFP2.0.txt",
SNP_FILTER_OUT_REF = "/gpfs/group/sanford/src/RNA/vcf/common_no_known_medical_impact_00-latest.vcf",
TISSUE= "EPI",
#ChIP-seq parameters
BOWTIE_VERSION = '1.0.0',
BOWTIE_OPTIONS = '',
HOMER_VERSION = '4.5',
HOMER_TRIM_OPTIONS = '-3 GATCGGAAGAGCACACGTCT -mis 1 -minMatchLength 6 -min 45',
MACS_VERSION = '1.4.2',
STEPS_PAIRS = ['macs'], 
PAIR_LIST = 'test1_chip-test1_input test2_chip-test2_input',
CHROM_SIZES = '/gpfs/group/su/kfisch/references/hg19.chrom.sizes',
MACS_RESULTS = '',
HOMER_PEAKS_OPTIONS = '-style factor -o auto' ,
HOMER_ANNOTATE_OPTIONS = '',
HOMER_MOTIFS_OPTIONS = '-size 200 -mask',
HOMER_RESULTS = "/gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run/homer",

#OMICS_PIPE
DRMAA_PATH = "/opt/applications/pbs-drmaa/current/gnu/lib/libdrmaa.so",
DPS_VERSION= '1.3.1111' ,
PYTHON_VERSION="2.6.5",

#WGS parameters
BWA_OPTIONS = "-t 8 -M",
RESOURCES = "-l nodes=1:ppn=8 -l mem=31gb",

#SNIPR
TABIX_VERSION = '0.2.6',
ENCODING = "phred64",
TUMOR_TYPE= 'BRCA',
GENELIST= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/brca.txt",
COSMIC= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/cosmic.tsv",
CLINVAR="/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/clinvar.txt",
PHARMGKB_rsID= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/pharmgkbRSID.csv",
PHARMGKB_Allele= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/pharmgkbAllele.tsv",
DRUGBANK= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/drugbank.tsv",
CADD= "/gpfs/home/kfisch/scripts/omics_pipeline-devel/omics_pipe/scripts/reporting/ref/cadd.tsv.gz"

) #Please do not remove ). It is needed to create a dictionary from the parameters