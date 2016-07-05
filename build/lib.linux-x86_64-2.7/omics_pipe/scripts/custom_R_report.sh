#!/bin/bash
set -x

#source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $4/$1

#Load specific modules 
module load R/$3
module load dps/$6

#Runs Knitr R script to generate RNAseq report $1=SAMPLE, $2=REPORT_SCRIPT, $3=R_version, $4=REPORT_RESULTS, $5=R_MARKUP_FILE, $6 DPS_VERSION $7 PARAMS_FILE

Rscript $2 $5 $7 $1 $4/$1


exit 0
