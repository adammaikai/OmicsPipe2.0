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
# convert vcf to annovar specific input

convert2annovar.pl $2/$1/$1.vcf -format vcf4 -outfile $2/$1/$1.avinput --includeinfo

# annotate using anovar: identify a small subset of variants/genes that are likely to be related to diseases
# Identifying nonsynonymous and splicing variants
# Removing variants in segmental duplication regions
# Removing variants not observed in 1000 Genomes Project 2012 April release or ESP5400 European Americans or ESP5400 African Americans, removing variants observed in dbSNP135 Non Flagged set
# Remove any variants observed in the CG46 database, 
# MAF threhsold of 0.01 will be applied to all the 1000G, ESP6500 and CG46 databases. 
# Variants believed to be likely benign by SIFT or PolyPhen are removed
# Apply a dominant disease model.
# UCSC Known Gene will be used for gene-based annotation.

variants_reduction.pl $2/$1/$1.avinput $3 $4 -outfile $2/$1/reduce $5

# sort, zip & index
#bgzip $2/$1/$1_rename.vcf 
#bgzip $2/$1/$1.vcf 
#tabix -p vcf $2/$1/$1.vcf.gz
        
#Clean up
#rm $2/$1/*.avinput


exit 0

