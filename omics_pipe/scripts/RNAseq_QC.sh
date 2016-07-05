#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $2/$6

#Move tmp dir to scratch 
export TMPDIR=$5  #TEMP_DIR
 
#Load specified module versions
module load python
module load picard/$7
module load R/$8

####INPUTS: $1: STAR_RESULTS $2 QC_PATH $3 BAM_FILE_NAME $4 RSEQC_REF  $5 TMP_DIR $6 SAMPLE $7 PICARD_VERSION $8 R_VERSION


java -Xmx2g -jar `which CollectRnaSeqMetrics.jar` REF_FLAT=$4 STRAND_SPECIFICITY=NONE INPUT=$1/$6/$3 OUTPUT=$2/$6/rnaseqmetrics &

java -Xmx2g -jar `which CollectInsertSizeMetrics.jar` HISTOGRAM_FILE=$2/$6/insertSizeHist.pdf INPUT= $1/$6/$3 OUTPUT=$2/$6/insertSize.txt &

wait

exit 0