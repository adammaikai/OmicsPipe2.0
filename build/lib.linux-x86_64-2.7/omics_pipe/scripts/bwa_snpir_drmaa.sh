#!/bin/bash
set -x

####INPUTS: $1: BWA_RESULTS $2: TEMP_DIR $3 SAMTOOLS_VERSION $4 BWA_VERSION $5 BWA_INDEX_FOLDER $6 SAMPLE $7 RAW_DATA_DIR $8 GATK_READ_GROUP_INFO

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/$6

#Move tmp dir to scratch 
export TMPDIR=$2/$BASHPID   #TEMP_DIR
 
#Load specified module versions
module load samtools/$3
module load bwa/$4


SNPiR="/gpfs/home/meissto/bin/SNPiR" # this needs to go on garibaldi $X
export PERL5LIB="/gpfs/home/meissto/bin" # point to one dir upstream of SNPiR config.pm $X1

ncpu=$(grep -c "processor" /proc/cpuinfo) 
nthreads=$((ncpu/2))


bwa_index=$5 # needs to go on garibaldi

out_dir=$1/$6

#workaround ..
cp $SNPiR/convertCoordinates.* $out_dir # seems to run only from the working dir, maybe works after its on garibaldi..

# input fastq
in1=$7/$6.fastq

# read group information for read1 & read2, is required for GATK
RGR1=$8


# select bwa index based on sequence length of the input
zcat in1=$7/$6.fastq | awk '{if(NR%4==2) print length($1)}' > $out_dir/input.readslength.txt
length=$(sort -r $out_dir/input.readslength.txt | uniq -c | head -n 1 | cut -d ' ' -f 2)

if [ $length -lt 75 ]; then
	# align with bwa as single reads
	bwa aln -t $nthreads $bwa_index/hg19_genome_junctions_45.fa $in1 > $out_dir/$6.sai 	
	bwa samse -n4 -r $RGR1 $bwa_index/hg19_genome_junctions_45.fa $out_dir/$6.sai $in1 > $out_dir/$6.sam 
fi

if [ $length -eq 75 ] && [ $length -lt 100 ]; then
	# align with bwa as single reads
	bwa aln -t $nthreads $bwa_index/hg19_genome_junctions_75.fa $in1 > $out_dir/$6.sai 	
	bwa samse -n4 -r $RGR1 $bwa_index/hg19_genome_junctions_75.fa $out_dir/$6.sai $in1 > $out_dir/$6.sam 
fi

if [ $length -eq 100 ] && [ $length -lt 150 ]; then
	# align with bwa as single reads
	bwa aln -t $nthreads $bwa_index/hg19_genome_junctions_95.fa $in1 > $out_dir/$6.sai 	
	bwa samse -n4 -r $RGR1 $bwa_index/hg19_genome_junctions_95.fa $out_dir/$6.sai $in1 > $out_dir/$6.sam 
fi

if [ $length -eq 150 ]]; then
	# align with bwa as single reads
	bwa aln -t $nthreads $bwa_index/hg19_genome_junctions_140.fa $in1 > $out_dir/$6.sai 	
	bwa samse -n4 -r $RGR1 $bwa_index/hg19_genome_junctions_140.fa $out_dir/$6.sai $in1 > $out_dir/$6.sam 
fi

else 
	echo "Sequence length not supported!"
fi

exit 0
