#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $4/$1

#Move tmp dir to scratch 
export TMPDIR=$6  #TEMP_DIR
 
#Loads specified module versions
module load homer/$5
module load ucsc-tools

#Calls scripts from HOMER to make TagDirectory and UCSC file
#1sample $2BOWTIE_RESULTS $3CHROM_SIZES $4HOMER_RESULTS $5HOMER_VERSION $6 TEMP_DIR

makeTagDirectory $4/$1_tag $2/$1/$1.bam

makeUCSCfile $4/$1_tag -o auto -bigWig $3 -fsize 1e20 > $4/$1/$1_trackInfo.txt

exit 0