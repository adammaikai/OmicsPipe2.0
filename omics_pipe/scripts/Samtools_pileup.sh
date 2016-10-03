#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load samtools/$3
module load sambamba/$6

cpu=$(grep -c "processor" /proc/cpuinfo)

####INPUTS $1=sample $2=ALIGNMENT_DIR $3=SAMTOOLS_VERSION $4=GENOME $5=CAPTURE_KIT_BED $6=SAMBAMBA_VERSION

## samtools mpileup
# if [ -z "$5" ]; then
#     samtools mpileup -f $4 $2/$1/$1\_gatk_recal.bam > $2/$1/$1\.pileup
#     #sambamba mpileup -t $cpu $2/$1/$1\_gatk_recal.bam -o $2/$1/$1\.pileup --samtools -f $4
# else
#     samtools mpileup -f $4 --positions $5 $2/$1/$1\_gatk_recal.bam > $2/$1/$1\.pileup
#     #sambamba mpileup -t $cpu -L $5 $2/$1/$1\_gatk_recal.bam -o $2/$1/$1\.pileup --samtools -f $4
# fi

if [ -z "$5" ]; then
    samtools view -b -q 20 $2/$1/$1\_gatk_recal.bam | samtools mpileup -f $4 - > $2/$1/$1\.pileup
else
    samtools view -b -q 20 $2/$1/$1\_gatk_recal.bam | samtools mpileup -f $4 --positions $5 - > $2/$1/$1\.pileup
fi

# filter for alleles with coverage > 0
awk -F"\t" '$4 > 0' $2/$1/$1\.pileup > $2/$1/$1\.flt.pileup
rm $2/$1/$1\.pileup

exit 0
