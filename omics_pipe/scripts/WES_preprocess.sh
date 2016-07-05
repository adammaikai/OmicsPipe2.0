#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load bwa/$5
module load samtools/$6
module load sambamba/$7
module load samblaster/$8
module load abra/$9
cpu=$(grep -c "processor" /proc/cpuinfo)

####INPUTS: $1=sample $2=TEMP_DIR $3=GENOME $4=BWA_INDEX $5=BWA_VERSION $6=SAMTOOLS_VERSION $7=SAMBAMBA_VERSION $8=SAMBLASTER_VERSION $9=ABRA_VERSION $10=CAPTURE_KIT_BED $11=ALIGNMENT_DIR $12=FASTQ_PATH $13=LOG_PATH 

mkdir -p "${11}"
mkdir -p "${11}/$1"

RGR="@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:$1"

if [ -f ${12}/$1_1.fastq.gz ]
then
	rawdata_dir=${12}
else
	rawdata_dir=${12}/$1
fi

## align to hg19
bwa mem -M -t $cpu -R $RGR $4 $rawdata_dir/$1_1.fastq.gz $rawdata_dir/$1_2.fastq.gz \
| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 10G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 10G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin) \
| samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 10G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin

## index
samtools index $2/$1\_sorted.bam

## mark duplicates
if [ -f $2/$1\_sorted.bam ]; then
	sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 $2/$1\_sorted.bam $2/$1\_dedup.bam
fi

## index the bam file
if [ -f $2/$1\_dedup.bam ]; then
	samtools index  $2/$1\_dedup.bam
	rm $2/$1\_sorted*
fi

## realign with abra
java -Xmx16G -jar `which abra.jar` \
--in $2/$1\_dedup.bam \
--out $2/$1\_realigned.bam \
--ref $3 \
--targets ${10} \
--threads $cpu \
--working $2/$1\_abra_temp_dir > $2/$1\_abra.log 2>&1

if [ -f $2/$1\_realigned.bam ]; then
    rm $2/$1\_dedup.ba*
fi

## sort and index again
if [ -f $2/$1\_realigned.bam ]; then
    sambamba sort -t $cpu -m 10G --tmpdir $2 -o ${11}/$1/$1\_realigned_sorted.bam $2/$1\_realigned.bam
    samtools index ${11}/$1/$1\_realigned_sorted.bam
fi

if [ -f ${11}/$1/$1\_realigned_sorted.bam ]; then
    rm $2/$1\_realigned.bam
fi
