#!/bin/bash
set -x

####INPUTS: $1: BWA_RESULTS $2: TEMP_DIR $3 SAMTOOLS_VERSION $4 BWA_VERSION $5 BWA_INDEX $6 SAMPLE $7 RAW_DATA_DIR $8 BWA_OPTIONS $9 COMPRESSION $10 SAMBAMBA_VERSION $11 SAMBLASTER_VERSION $12 SAMBAMBA_OPTIONS

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/$6

#Move tmp dir to scratch 
export TMPDIR=$2  #TEMP_DIR
 
#Load specified module versions
module load samtools/$3
module load bwa/$4
module load sambamba/${10}
module load samblaster/${11}


######################## MAPPING WITH BWA #######################
comp=$9
#create config file
if [ "$comp" = "GZIP" ]; then
	bwa mem $8 $5 $7/$6_1.fastq.gz $7/$6_2.fastq.gz \
	| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
	| sambamba sort ${12} --tmpdir $2 -o $1/$6/$6-sr.bam /dev/stdin) \
	--discordantFile >(samtools view -S -u /dev/stdin \
	| sambamba sort ${12} --tmpdir $2 -o $1/$6/$6-disc.bam /dev/stdin) \
	| samtools view -S -u /dev/stdin \
	| sambamba sort ${12} --tmpdir $2 -o $1/$6/$6_sorted.bam /dev/stdin
else
	bwa mem $8 $5 $7/$6_1.fastq $7/$6_2.fastq \
	| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
	| sambamba sort ${12} --tmpdir $2 -o $1/$6/$6-sr.bam /dev/stdin) \
	--discordantFile >(samtools view -S -u /dev/stdin \
	| sambamba sort ${12} --tmpdir $2 -o $1/$6/$6-disc.bam /dev/stdin) \
	| samtools view -S -u /dev/stdin \
	| sambamba sort ${12} --tmpdir $2 -o $1/$6/$6_sorted.bam /dev/stdin
fi


exit 0


#SAMBAMBA_OPTIONS: -t 28 -m 682M