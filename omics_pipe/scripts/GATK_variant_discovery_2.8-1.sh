#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR

#Loads specified module versions
module load gatk/$4

#Make output directory if it doesn't exist
mkdir -p $7/$2

####INPUTS: $1: BWA_RESULTS $2 SAMPLE $3: TEMP_DIR $4 GATK_VERSION $5 GENOME $6 DBSNP $7 VARIANT_RESULTS

dbsnp=$6
temp=$3

######################## GATK HaplotypeCaller ###############################
java -Xms454m -Xmx3181m -Djava.io.tmpdir=$temp \
	-jar `which GenomeAnalysisTK.jar` \
	-T HaplotypeCaller \
	-nct 8 \
	-R $5 \
	--dbsnp $dbsnp \
	-I $1/$2/$2.ready.bam \
	-o $7/$2/$2.raw.vcf

exit 0