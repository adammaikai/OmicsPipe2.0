#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make alignment output directory if it doesn't exist
mkdir -p $4

#Loads specified module versions
module load homer/$5

#Trim fastq files using homerTools $1 Sample $2 RAW_DATA_DIR $3 HOMER_TRIM_OPTIONS $4 TRIMMED_DATA_PATH $5HOMER_VERSION

homerTools trim $3 $2/$1.fastq 

mv $2/$1.fastq.trimmed $4/$1.fastq

exit 0
