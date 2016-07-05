#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR

#Loads specified module versions
module load gatk/$4
module load samtools/${10}


####INPUTS: $1: BWA_RESULTS $2 SAMPLE $3: TEMP_DIR $4 GATK_VERSION $5 GENOME $6 DBSNP $7 MILLS $8 G1000 $9 CAPTURE_KIT_BED regions.bed $10 SAMTOOLS_VERSION 

dbsnp=$6
mills=$7
g1000=$8


######################## GATK RealignerTargetCreator ###################
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
    -T RealignerTargetCreator -nt 8 \
    -R $5 \
    -I $1/$2/$2_sorted.rg.md.bam \
    -o $1/$2/$2_sorted.rg.md.intervals \
    --fix_misencoded_quality_scores \
    --filter_mismatching_base_and_quals \
    -L $9 \
    -known $mills \
    -known $dbsnp \

    
    
######################## GATK IndelRealigner ###########################
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
    -T IndelRealigner \
    -R $5 \
    -I $1/$2/$2_sorted.rg.md.bam \
    -targetIntervals $1/$2/$2_sorted.rg.md.intervals \
    -o $1/$2/$2_sorted.rg.md.ir.bam \
    --fix_misencoded_quality_scores \
    --filter_bases_not_stored \
    --filter_mismatching_base_and_quals \
    -known $mills \
    -known $dbsnp \

######################## ANOTHER SAMTOOLS INDEX ########################
samtools index $1/$2/$2_sorted.rg.md.ir.bam

######################## GATK BaseRecalibrator #########################
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
    -T BaseRecalibrator -nct 8 \
    -R $5 \
    -I $1/$2/$2_sorted.rg.md.ir.bam \
    -o $1/$2/$2_sorted.rg.md.ir.grp \
    -knownSites $mills \
    -knownSites $dbsnp \
    -knownSites $g1000

######################## GATK PrintReads ###############################
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
    -T PrintReads -nct 8 \
    -R $5 \
    -I $1/$2/$2_sorted.rg.md.ir.bam \
    -BQSR $1/$2/$2_sorted.rg.md.ir.grp \
    -o $1/$2/$2.ready.bam
    
exit 0
