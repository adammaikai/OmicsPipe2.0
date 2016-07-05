#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Loads specified python addons module for cutadapt
module load python/$5

#Make assemblies output director if it doesn't exist
mkdir -p $4

#Runs cutadapt with $1=SAMPLE, $2=RAW_DATA_DIR, $3=ADAPTER, $4=TRIMMED_DATA_PATH $5PYTHON_VERSION
cutadapt -a $3 $2/$1.fastq > $4/$1_trimmed.fastq

exit 0
