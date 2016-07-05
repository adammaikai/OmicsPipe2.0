#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load star/$5
module load star-fusion/$9

####INPUTS $1=sample $2=TEMP_DIR $3=RAW_DATA_DIR $4=ALIGNMENT_DIR $5=STAR_VERSION $6=GENOME $7=REF_GTF $8=STARFUSION_LIB $9=STARFUSION_VERSION

mkdir -p "$4/$1"

cpu=$(grep -c "processor" /proc/cpuinfo)

rgid="1234"
rgpl="ILLUMINA"
rglb="TrueSeq"

if [ -f $3/$1_1.fastq.gz ]
then
	rawdata_dir=$3
else
	rawdata_dir=$3/$1
fi

# make sure temp is empty
if [ -d $2/$1 ]
then
    rm -rf $2/$1
fi

mkdir -p $2/$1 && cd $2/$1

STAR \
--genomeDir $6 \
--readFilesIn $rawdata_dir/$1_1.fastq.gz $rawdata_dir/$1_2.fastq.gz \
--readFilesCommand zcat \
--twopassMode Basic \
--quantMode GeneCounts \
--sjdbGTFfile $7 \
--alignIntronMax 200000 \
--alignMatesGapMax 200000 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattrRGline ID:$rgid SM:$1 PL:$rgpl LB:$rglb \
--outSAMstrandField intronMotif \
--runThreadN $cpu \
--chimSegmentMin 25 \
--chimJunctionOverhangMin 25 \
--outFileNamePrefix $1_

cp $2/$1/* $4/$1/

#clean up
rm -rf $2/$1
rm -rf $2/$1

#STAR-FUSION
STAR-Fusion \
--genome_lib_dir $8 \
-J $4/$1/$1\_Chimeric.out.junction \
--output_dir $4/$1


exit 0
