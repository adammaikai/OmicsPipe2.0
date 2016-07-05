#!/bin/bash
set -x

####INPUTS: $1: BWA_RESULTS $2: TEMP_DIR $3 SAMTOOLS_VERSION $4 BWA_VERSION $5 BWA_INDEX $6 SAMPLE $7 RAW_DATA_DIR $8 GATK_READ_GROUP_INFO $9 COMPRESSION

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/$6

#Move tmp dir to scratch 
export TMPDIR=$2/$BASHPID   #TEMP_DIR
 
#Load specified module versions
module load samtools/$3
module load bwa/$4


# SNPiR="/gpfs/home/meissto/bin/SNPiR" # this needs to go on garibaldi $X
# export PERL5LIB="/gpfs/home/meissto/bin" # point to one dir upstream of SNPiR config.pm $X1

ncpu=$(grep -c "processor" /proc/cpuinfo) 
#nthreads=$((ncpu/2))


bwa_index=$5 # needs to go on garibaldi

out_dir=$1/$6

#workaround ..
#cp $SNPiR/convertCoordinates.* $out_dir # seems to run only from the working dir, maybe works after its on garibaldi..

# input fastq
in1=$7/$6.fastq

# read group information for read1 & read2, is required for GATK
RGR1=$8

comp=$9
if [ "$comp" = "GZIP" ]; then
	bwa aln -t $ncpu $bwa_index $in1.gz > $out_dir/$6.sai 
	bwa samse -n4 -r $RGR1 $bwa_index $out_dir/$6.sai $in1.gz > $out_dir/$6.sam
else
	bwa aln -t $ncpu $bwa_index $in1 > $out_dir/$6.sai
	bwa samse -n4 -r $RGR1 $bwa_index $out_dir/$6.sai $in1 > $out_dir/$6.sam
fi

exit 0
