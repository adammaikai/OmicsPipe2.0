#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash


#Make HTSEQ count output directory if it doesn't exist
mkdir -p $5

#Move tmp dir to scratch 
mkdir -p $6/$BASHPID
export TMPDIR=$6/$BASHPID  

#Load specified samtools version and python
module load samtools/$7
module load python/$9

#Runs HTSEQ count with $1=SAMPLE, $2=STAR_RESULTS, $3=HTSEQ_OPTIONS, $4REF_GENES, $5=HTSEQ_RESULTS, $6=TEMP_DIR, $7=SAMTOOLS_VERSION, $8 BAM_FILE_NAME $9 PYTHON_VERSION
samtools view $2/$1/$8 | sort -s -k 1,1 - |  htseq-count $3 - $4 > $5/$1_counts.txt

exit 0
