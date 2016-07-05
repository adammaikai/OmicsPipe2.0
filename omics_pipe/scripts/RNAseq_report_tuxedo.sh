#!/bin/bash
set -x

#source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $4/$1

#Load specific modules 
module load R/$3
module load dps/$5

#Runs Knitr R script to generate RNAseq report $1=SAMPLE, $2=WORKING_DIR, $3=R_version, $4=REPORT_RESULTS, $5 DPS_VERSION $6 PARAMS_FILE

Rscript $2/reporting/src/knitMeR.R $2/reporting/src/RNAseq_tuxedo.Rmd $6 $1 $4/$1


exit 0
