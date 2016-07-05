#!/bin/bash
set -x

####INPUTS: $1: BWA_RESULTS $2: TEMP_DIR $3 SAMTOOLS_VERSION $4 BWA_VERSION $5 GENOME $6 SAMPLE $7 RAW_DATA_DIR $8 BWA_OPTIONS $9 COMPRESSION

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/$6

#Move tmp dir to scratch 
export TMPDIR=$2  #TEMP_DIR
 
#Load specified module versions
module load samtools/$3
module load bwa/$4


######################## MAPPING WITH BWA #######################
comp=$9
#create config file
if [ "$comp" = "GZIP" ]; then
	bwa mem $8 $5 $7/$6.fastq.gz > $1/$6/$6.sam
else
	bwa mem $8 $5 $7/$6.fastq > $1/$6/$6.sam
fi



######################## SAMTOOLS SORT AND INDEX #######################
samtools view $1/$6/$6.sam -bS -o $1/$6/$6.bam
samtools sort $1/$6/$6.bam $1/$6/$6_sorted
samtools index $1/$6/$6_sorted.bam
exit 0