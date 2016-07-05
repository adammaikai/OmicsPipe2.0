#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make alignment output directory if it doesn't exist
mkdir -p $3

#Loads specified tophat version
module load tophat/$5
module load bowtie/$7
module load samtools/$8
module load bowtie2/$7

#Runs Tophat with $1 = SAMPLE, $2=RAW_DATA_DIR, $3=RESULTS_PATH, $4=BOWTIE_INDEX, $5=TOPHAT_VERSION, $6=TOPHAT_OPTIONS $7BOWTIE_VERSION, $8 SAMTOOLS_VERSION

tophat $6 -o $3/$1 $4 $2/$1.fastq

exit 0
