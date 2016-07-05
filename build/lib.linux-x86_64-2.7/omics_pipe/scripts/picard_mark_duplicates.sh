#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR

#Loads specified module versions
module load picard/$4
module load samtools/$5


####INPUTS: $1: BWA_RESULTS $2 SAMPLE $3: TEMP_DIR $4 PICARD_VERSION $5 SAMTOOLS_VERSION


######################## PICARD AddOrReplaceReadGroups #################
java -Xms454m -Xmx3181m -jar `which AddOrReplaceReadGroups.jar` \
	I=$1/$2/$2_sorted.bam \
	O=$1/$2/$2_sorted.rg.bam \
	SM=$2 \
	VALIDATION_STRINGENCY=LENIENT \
	ID=1 \
	LB=unknown \
	PL=illumina \
	PU=1234

######################## PICARD MarkDuplicates #########################
java -Xms454m -Xmx3181m -jar `which MarkDuplicates.jar` \
	I=$1/$2/$2_sorted.rg.bam \
	O=$1/$2/$2_sorted.rg.md.bam \
	M=$1/$2/$2_sorted.rg.md.metrics \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=LENIENT

samtools index $1/$2/$2_sorted.rg.md.bam

exit 0