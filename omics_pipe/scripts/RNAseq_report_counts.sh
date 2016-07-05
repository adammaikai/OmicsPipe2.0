#!/bin/bash
set -x

#source modules for current shell
source $MODULESHOME/init/bash

#Make report directory if it doesn't exist
mkdir -p $4/$1

#Make DEseq directory if it doesn't exist
mkdir -p $5

#Load specific modules 
module load R/$3

#Runs Knitr R script to generate RNAseq report $1=SAMPLE, $2=WORKING_DIR, $3=R_version, $4=REPORT_RESULTS, $5=DESEQ_RESULTS, $6 PARAMS_FILE

Rscript $2/reporting/src/knitMeR.R $2/reporting/src/RNAseq_counts.Rmd $6 $1 $4/$1

exit 0
