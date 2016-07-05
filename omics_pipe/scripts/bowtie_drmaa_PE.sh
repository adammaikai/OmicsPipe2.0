#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make alignment output directory if it doesn't exist
mkdir -p $5/$1

#Move tmp dir to scratch 
mkdir -p $9/$BASHPID
export TMPDIR=$9/$BASHPID  

#Loads specified module versions
module load bowtie/$6
module load bowtie2/$6
module load samtools/$7
module load bedtools/$8

#Runs Bowtie with $1 = SAMPLE, $2 = RAW_DATA_DIR, $3 = BOWTIE_OPTIONS, $4 BOWTIE_INDEX, $5=BOWTIE_RESULTS, $6 = BOWTIE_VERSION, $7 SAMTOOLS_VERSION, $8 BEDTOOLS_VERSION, $9 TEMP_DIR

bowtie $3 $4 -1 $2/$1_1.fastq -2 $2/$1_2.fastq > $5/$1/$1.sam

#Convert file from SAM to BAM format
samtools view -Sb $5/$1/$1.sam > $5/$1/$1_nonSorted.bam

#Sort BAM file
samtools sort $5/$1/$1_nonSorted.bam $5/$1/$1

#Create index file (BAI)
samtools index $5/$1/$1.bam

#Remove intermediate files
rm $5/$1/$1.sam $5/$1/$1_nonSorted.bam

#QC for read mapping
raw=$(samtools view $5/$1/$1.bam | wc -l)
bamToBed -i $5/$1/$1.bam | awk -vRAW=$raw '{coordinates=$1":"$2"-"$3;total++;count[coordinates]++}END{for(coordinates in count){if(!max||count[coordinates]>max){max=count[coordinates];maxCoor= coordinates};if(count[coordinates]==1){unique++}};print RAW,total,total*100/RAW,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}' > $5/$1/$1_num_reads.txt
samtools view -f 0x0004 $5/$1/$1.bam | awk '{read=$10;total++;count[read]++}END{print "Total_non-mapped_reads",total;for(read in count){print read,count[read]+0}}' | sort -k2,2nr | head -11 > $5/$1/$1_nonmappedreads.txt

exit 0
