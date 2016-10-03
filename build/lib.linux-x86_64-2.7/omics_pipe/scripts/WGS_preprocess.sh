#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load bwa/$5
module load samtools/$6
module load sambamba/$7
module load samblaster/$8
module load gatk/$9
cpu=$(grep -c "processor" /proc/cpuinfo)

####INPUTS 
####		$1=sample $2=TEMP_DIR $3=GENOME $4=BWA_INDEX $5=BWA_VERSION 
####        $6=SAMTOOLS_VERSION $7=SAMBAMBA_VERSION $8=SAMBLASTER_VERSION $9=GATK_VERSION $10=ALIGNMENT_DIR 
####		$11=FASTQ_PATH $12=DBSNP $13=MILLS $14=G1000 $15=CAPTURE_KIT_BED

mkdir -p "${10}/$1"

RGR="@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:$1"

# ## align to hg19
bwa mem -M -t $cpu -R $RGR $4 ${11}/$1/$1\_1.fastq.gz ${11}/$1/$1\_2.fastq.gz \
| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 16G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 16G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin) \
| samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 16G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin

## index
samtools index $2/$1\_sorted.bam

## mark duplicates
if [ -f $2/$1\_sorted.bam ]; then
	sambamba markdup \
	-t $cpu \
	--overflow-list-size 1000000 \
	--hash-table-size 1000000 \
	$2/$1\_sorted.bam \
	${10}/$1/$1\_dedup.bam
fi

## index the bam file
if [ -f ${10}/$1/$1\_dedup.bam ]; then
	samtools index  ${10}/$1/$1\_dedup.bam
    rm $2/$1\_sorted*
fi

if [ -z "$15" ]; then

	## GATK Realignment
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T RealignerTargetCreator -nt $cpu \
	-R $3 -I ${10}/$1/$1\_dedup.bam \
	-known ${13} \
	-known ${14} \
	-o $2/$1\.intervals

	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T IndelRealigner -R $3 -I ${10}/$1/$1\_dedup.bam \
	-targetIntervals $2/$1\.intervals \
	-known ${13} \
	-known ${14} \
	-o ${10}/$1/$1\_realigned.bam

	if [ -f ${10}/$1/$1\_realigned.bam ]; then
	    rm ${10}/$1/$1\_dedup.ba*
	fi

	## BQSR
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T BaseRecalibrator -nct $cpu \
	-R $3 \
	-I ${10}/$1/$1\_realigned.bam  \
	-knownSites ${12} \
	-knownSites ${13} \
	-knownSites ${14} \
	-o $2/$1\_recal_data.table

	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T PrintReads -nct $cpu \
	-R $3 \
	-I ${10}/$1/$1\_realigned.bam \
	-BQSR $2/$1\_recal_data.table \
	-o ${10}/$1/$1\_gatk_recal.bam
else
	## GATK Realignment
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T RealignerTargetCreator -nt $cpu \
	-L ${15} \
	-R $3 -I ${10}/$1/$1\_dedup.bam \
	-known ${13} \
	-known ${14} \
	-o $2/$1\.intervals

	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T IndelRealigner -R $3 -I ${10}/$1/$1\_dedup.bam \
	-targetIntervals $2/$1\.intervals \
	-known ${13} \
	-known ${14} \
	-o ${10}/$1/$1\_realigned.bam

	if [ -f ${10}/$1/$1\_realigned.bam ]; then
	    rm ${10}/$1/$1\_dedup.ba*
	fi

	## BQSR
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T BaseRecalibrator -nct $cpu \
	-R $3 \
	-L ${15} \
	-I ${10}/$1/$1\_realigned.bam  \
	-knownSites ${12} \
	-knownSites ${13} \
	-knownSites ${14} \
	-o $2/$1\_recal_data.table

	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T PrintReads -nct $cpu \
	-R $3 \
	-I ${10}/$1/$1\_realigned.bam \
	-BQSR $2/$1\_recal_data.table \
	-o ${10}/$1/$1\_gatk_recal.bam
fi

## clean up
if [ -f ${10}/$1/$1\_gatk_recal.bam ]; then
    rm ${10}/$1/$1\_realigned.ba*
    rm $2/$1\_recal_data.table
    rm $2/$1.intervals
fi


exit 0




