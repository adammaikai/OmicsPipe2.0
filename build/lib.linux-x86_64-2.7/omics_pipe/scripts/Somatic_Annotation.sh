#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load R/$3

####INPUTS: $1=sample $2=VARIANT_DIR $3=R_VERSION 

normal=$(echo $5 | sed -e 's/-DNA//' | cut -c 2-)
tumor=$(echo $4 | sed -e 's/-DNA//' | cut -c 2-)

Rscript ~/.virtualenvs/pm/omics_pipe/omics_pipe/scripts/somaticAnnotation.R $2/$1/$1\_$tumor\_$normal $2/$1/$1\_$tumor\_$normal\_varscan_somatic.vcf.gz $2/$1/$1\_$tumor\_$normal\_mutect.filt.vcf.gz TRUE

exit 0
