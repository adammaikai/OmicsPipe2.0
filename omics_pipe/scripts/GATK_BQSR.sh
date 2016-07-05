#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load gatk/$4
module load samtools
cpu=$(grep -c "processor" /proc/cpuinfo)

mkdir -p "$7"

####INPUTS: $1=sample $2=TEMP_DIR $3=GENOME $4=GATK_VERSION $5=R_VERSION $6=CAPTURE_KIT_BED $7=ALIGNMENT_DIR $8=DBSNP $9=MILLS $10=G1000


java -Xmx28g -Djava.io.tmpdir=$2 -jar `which GenomeAnalysisTK.jar` \
-T BaseRecalibrator -nct $cpu \
-R $3 \
-I $7/$1/$1\_realigned_sorted.bam \
-o $2/$1\_gatk_recal_data.table \
-L $6 \
-ip 1 \
-rf BadCigar \
-knownSites ${8} \
-knownSites ${9} \
-knownSites ${10}

java -Xmx28g -Djava.io.tmpdir=$2 -jar `which GenomeAnalysisTK.jar` \
-T PrintReads -nct $cpu \
-R $3 \
-I $7/$1/$1\_realigned_sorted.bam \
-BQSR $2/$1\_gatk_recal_data.table \
-o $7/$1/$1\_gatk_recal.bam

if [ -f $$7/$1/$1\_gatk_recal.bam ]; then
	rm $2/$1\_gatk.intervals
	rm $2/$1\_gatk_recal_data.table
	rm $7/$1/$1\_realigned_sorted.ba*
fi

exit 0


