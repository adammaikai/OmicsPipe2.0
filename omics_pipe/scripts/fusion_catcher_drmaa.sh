#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make output directory if it doesn't exist
mkdir -p $2

#Move tmp dir to scratch 
export TMPDIR=$4/$BASHPID  #TEMP_DIR
scratch=$TMPDIR

#Load specified samtools version and python
module load samtools/$5
module load python/$10
module load fusioncatcher/$6

#Updated command to be compatible with recent releases 7/11/14
#Runs fusion catcher with $1= RAW_DATA_DIR $2=FUSIONCATCHER_RESULTS, $3=FUSIONCATCHERBUILD_DIR, $4 TEMP_DIR, $5 SAMTOOLS_VERSION, $6 FUSIONCATCHER_VERSION. $7 FUSIONCATCHER_OPTIONS, $8 Sample $9 TISSUE $10 PYTHON_VERSION
fusioncatcher -d $3 -i $1/$8_1.fastq,$1/$8_2.fastq -o $2/$8 $7

cat $2/$8/final-list_candidate-fusion-genes.txt | cut -f 9,10 | tail -n +2 | sed 's/:+//g' | sed 's/:-//g' | awk -v tissue=$9 '{print "chr"$1, "\t", "chr"$2, "\t", tissue}' | sed -e 's+:+\t+g' > $2/$8/final_list_coord.txt

exit 0
