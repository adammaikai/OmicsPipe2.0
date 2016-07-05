#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR

#Loads specified module versions
module load mutect/$4

#Make output directory if it doesn't exist
mkdir -p $8

####INPUTS: $1: BWA_RESULTS $2 SAMPLE $3: TEMP_DIR $4 MUTECT_VERSION $5 GENOME $6 DBSNP $7 COSMIC  $8 VARIANT_RESULTS $9 TUMOR_SAMPLE_NAME $10 NORMAL_SAMPLE_NAME


java -Xms454m -Xmx3181m -jar `which muTect.jar` \
--analysis_type MuTect \
--reference_sequence $5 \
--cosmic $7 \
--dbsnp $6 \
--intervals $1/${10}/${10}_sorted.rg.md.intervals \
--input_file:normal $1/${10}/${10}.ready.bam \
--input_file:tumor $1/$9/$9.ready.bam \
--out $8/$2/call_stats.txt \
--coverage_file $8/$2/coverage.wig.txt \
--vcf $8/$2/$2.raw.vcf \

exit 0

