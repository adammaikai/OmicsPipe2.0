#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 

export TMPDIR=$6  #TEMP_DIR

#Make output directory if it doesn't exist
mkdir -p $2/$1/reduce

#Load specific module versions
module load annovar/$7
module load vcftools/$8


#Variant annotation with $1=SAMPLE, $2=VARIANT_RESULTS, $3=ANNOVARDB, $4=ANNOVAR_OPTIONS, $5=ANNOVAR_OPTIONS2, $6=TEMP_DIR, $7=ANNOVAR_VERSION, $8=VCFTOOLS_VERSION

convert2annovar.pl $2/$1/$1.vcf -format vcf4 -outfile $2/$1/$1.avinput --includeinfo


#Since the input file is in hg19 coordinate, we added "-buildver hg19" in every command above. Similarly, if you generated variant calls from human GRCh38 coordinate, add -buildver hg38' in every command

#gene based annotation
annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/ 
annotate_variation.pl -downdb -buildver hg19 -webfrom annovar knownGene humandb/ 
annotate_variation.pl -downdb -buildver hg19 -webfrom annovar ensGene humandb/ 
annotate_variation.pl -geneanno -buildver hg19 -hgvs example/ex1.avinput humandb/ -dbtype knownGene

#region based annotation
annotate_variation.pl -downdb -buildver hg19 cytoBand humandb/
annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 example/ex1.avinput humandb/ 

#filter-based annotation
annotate_variation.pl -downdb -buildver hg19 1000g2014oct humandb/
annotate_variation.pl -filter -dbtype 1000g2014oct -buildver hg19 example/ex1.avinput humandb/

#Download hg19 annotation databases
annotate_variation.pl -buildver hg19 -downdb  cytoBand humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct humandb/
annotate_variation.pl -buildver hg19 -downdb phastConsElements46way humandb/ #most conserved element annotation
annotate_variation.pl -buildver hg19 -downdb tfbsConsSites humandb/ #transcription factor binding site annotation
annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/ #variants located in segmental duplications (likely sequence alignment errors)
annotate_variation.pl -buildver hg19 -downdb wgRna humandb/ #variants disrupting miRNAs and snoRNAs
annotate_variation.pl -buildver hg19 -downdb targetScanS humandb/ #variants disrupting predicted miRNA  binding sites from TargetScanHuman
annotate_variation.pl -buildver hg19 -downdb gwasCatalog humandb/ #Published GWAS variants
annotate_variation.pl -buildver hg19 -downdb dgvMerged humandb/ #previously reported structural variants in DGV (Database of genomic variants)
annotate_variation.pl -buildver hg19 -downdb wgEncodeRegTfbsClustered humandb/ #variants in ENCODE transcription factor ChIP=seq data
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/ #dbSNP variants
#annotate_variation.pl -buildver hg19 -downdb -webfrom annovar wgEncodeBroadHmmGm12878HMM humandb/ #non-coding variants that disrupt enhancers, repressors, promoters using chromHMM predictions
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/ #dbNSFP annotations (sift, polyphen, CADD etc)
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500si_all humandb/ #exome sequencing project (
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cg46 humandb/ #CG (complete genomics) frequency annotations
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cg69 humandb/ #CG frequency annotation
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar popfreq_all humandb/ #Annovar compiled pop freq from 1000G2012, ESP6500si and CG46
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20140929 humandb/ #Clinvar
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cosmic70 humandb/ #Cosmic
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar nci60 humandb/ #NCI 60 cell lines
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac01 humandb/ #ExAC 65000 exome alelle freq data

#download hg18 annotation databases
annotate_variation.pl -buildver hg18 -downdb wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 humandb/ #variants within transcribed regions in RNAseq data for GM12878 cell lines
annotate_variation.pl -buildver hg18 -downdb wgEncodeBroadChipSeqPeaksGm12878H3k4me1 humandb/ #variants in enhancer regions based on H3K4Me1
annotate_variation.pl -buildver hg18 -downdb wgEncodeRegDnaseClustered humandb/ #variants in DNase I hypersensitivity regions
annotate_variation.pl -buildver hg18 -downdb wgEncodeBroadChipSeqPeaksGm12878Ctcf humandb/ #variants in CTCF binding sites
annotate_variation.pl -buildver hg18 -downdb wgEncodeRegTfbsClustered humandb/ #variants in ENCODE transcription factor ChIP=seq data
annotate_variation.pl -buildver hg18 -downdb  cytoBand humandb/
annotate_variation.pl -buildver hg18 -downdb -webfrom annovar 1000g2012apr humandb/
#annotate_variation.pl -buildver hg18 -downdb phastConsElements46way humandb/ #most conserved element annotation
annotate_variation.pl -buildver hg18 -downdb tfbsConsSites humandb/ #transcription factor binding site annotation
annotate_variation.pl -buildver hg18 -downdb genomicSuperDups humandb/ #variants located in segmental duplications (likely sequence alignment errors)
annotate_variation.pl -buildver hg18 -downdb wgRna humandb/ #variants disrupting miRNAs and snoRNAs
annotate_variation.pl -buildver hg18 -downdb targetScanS humandb/ #variants disrupting predicted miRNA  binding sites from TargetScanHuman
annotate_variation.pl -buildver hg18 -downdb gwasCatalog humandb/ #Published GWAS variants
annotate_variation.pl -buildver hg18 -downdb dgvMerged humandb/ #previously reported structural variants in DGV (Database of genomic variants)
annotate_variation.pl -buildver hg18 -downdb wgEncodeRegTfbsClustered humandb/ #variants in ENCODE transcription factor ChIP=seq data
annotate_variation.pl -buildver hg18 -downdb -webfrom annovar snp138 humandb/ #dbSNP variants
#annotate_variation.pl -buildver hg18 -downdb -webfrom annovar wgEncodeBroadHmmGm12878HMM humandb/ #non-coding variants that disrupt enhancers, repressors, promoters using chromHMM predictions
annotate_variation.pl -buildver hg18 -downdb -webfrom annovar ljb26_all humandb/ #dbNSFP annotations (sift, polyphen, CADD etc)
annotate_variation.pl -buildver hg18 -downdb -webfrom annovar esp6500si_all humandb/ #exome sequencing project (
annotate_variation.pl -buildver hg18 -downdb -webfrom annovar cg46 humandb/ #CG (complete genomics) frequency annotations
annotate_variation.pl -buildver hg18 -downdb -webfrom annovar cg69 humandb/ #CG frequency annotation
#annotate_variation.pl -buildver hg18 -downdb -webfrom annovar popfreq_all humandb/ #Annovar compiled pop freq from 1000G2012, ESP6500si and CG46
#annotate_variation.pl -buildver hg18 -downdb -webfrom annovar clinvar_20140929 humandb/ #Clinvar
#annotate_variation.pl -buildver hg18 -downdb -webfrom annovar cosmic70 humandb/ #Cosmic
#annotate_variation.pl -buildver hg18 -downdb -webfrom annovar nci60 humandb/ #NCI 60 cell lines
#annotate_variation.pl -buildver hg18 -downdb -webfrom annovar exac01 humandb/ #ExAC 65000 exome alelle freq data


#output to csv
table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all -operation g,r,r,f,f,f,f -nastring . -csvout

#output to vcf
table_annovar.pl example/ex2.vcf humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all -operation g,r,r,f,f,f,f -nastring . -vcfinput


#all hg 19
table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove \
-protocol knownGene,phastConsElements46way,tfbsConsSites,cytoBand,wgRna,targetScanS,\
genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeRegTfbsClustered,\
esp6500si_all,ALL.sites.2014_10,snp138,ljb26_all,cg46,cg69,popfreq_all,clinvar_20140929,cosmic70,nci60 \
-operation g,r,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -csvout

#all hg 18
table_annovar.pl example/ex1.avinput humandb/ -buildver hg18 -out myanno -remove \
-protocol knownGene,tfbsConsSites,cytoBand,wgRna,targetScanS,\
genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75,\
wgEncodeBroadChipSeqPeaksGm12878H3k4me1,wgEncodeRegDnaseClustered,wgEncodeBroadChipSeqPeaksGm12878Ctcf,\
wgEncodeRegTfbsClustered,esp6500si_all,ALL.sites.2012_04,snp138,ljb26_all,cg46,cg69 \
-operation g,r,r,r,r,r,r,r,r,r,r,r,r,f,f,f,f,f,f -nastring . -csvout

perl convert2annovar.pl /Users/kfisch/Dropbox/Research/UCSD/research/autism_Sebat/vcf/74-0075/vqsr/combined-74-0075-vqsr-recalibrated-snps-only.vcf -format vcf4 -allsample -withfreq -outfile /Users/kfisch/Dropbox/Research/UCSD/research/autism_Sebat/vcf/74-0075/vqsr/combined-74-0075-vqsr-recalibrated-snps-only.avinput --includeinfo

perl table_annovar.pl /Users/kfisch/Dropbox/Research/UCSD/research/autism_Sebat/vcf/74-0075/vqsr/combined-74-0075-vqsr-recalibrated-snps-only.avinput humandb/ -buildver hg18 -out sebat/74-0075_snps -remove -protocol knownGene,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeRegTfbsClustered,esp6500si_all,ALL.sites.2012_04,snp138,ljb26_all,cg46,cg69 -operation g,r,r,r,r,r,r,r,r,f,f,f,f,f,f -nastring . -csvout