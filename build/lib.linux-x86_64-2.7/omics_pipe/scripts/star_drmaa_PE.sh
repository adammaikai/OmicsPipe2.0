#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make star alignment output directory if it doesn't exist
mkdir -p $5/$1

#Load specified versions of modules 
module load samtools/$6
module load star/$7

#Runs STAR with $1=SAMPLE, $2=RAW_DATA_DIR, $3 STAR_INDEX, $4=STAR_OPTIONS, $5 = STAR_RESULTS, #6 SAMTOOLS_VERSION, $7 STAR_VERSION $8 COMPRESSION $9 REF_GENES
comp=$8

if [ "$comp" = "GZIP" ]; then
	STAR \
	--genomeDir $3 \
	--sjdbGTFfile $9 \
        --readFilesIn $2/$1_1.fastq $2/$1_2.fastq \
        --readFilesCommand zcat \
        $4 \
        --outFileNamePrefix $5/$1/
else
	STAR \
	--genomeDir $3 \
	--sjdbGTFfile $9 \
	--readFilesIn $2/$1_1.fastq $2/$1_2.fastq \
	$4 \
	--outFileNamePrefix $5/$1/
fi

samtools view -Sb $5/$1/Aligned.out.sam > $5/$1/Aligned.out.bam
samtools sort -m 2G -@ 8 $5/$1/Aligned.out.bam $5/$1/Aligned.out.sorted
samtools index $5/$1/Aligned.out.sorted.bam

#rm $5/$1/Aligned.out.sam
#rm $5/$1/Aligned.out.bam

# remove the temp folder
rm -rf $5/$1/_STARtmp

exit 0
