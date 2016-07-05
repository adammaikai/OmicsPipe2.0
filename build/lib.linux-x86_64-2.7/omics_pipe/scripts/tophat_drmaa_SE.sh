#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make alignment output directory if it doesn't exist
mkdir -p $4

#Loads specified tophat version
module load tophat/$6
module load bowtie/$8
module load samtools/$9
module load bowtie2/$8

#Runs Tophat with $1 = SAMPLE, $2=RAW_DATA_DIR, $3=REF_GENES, $4=RESULTS_PATH, $5=BOWTIE_INDEX, $6=TOPHAT_VERSION, $7=TOPHAT_OPTIONS $8 BOWTIE_VERSION, $9 SAMTOOLS_VERSION

tophat $7 -G $3 -o $4/$1 $5 $2/$1.fastq

exit 0
