#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $2

#Move tmp dir to scratch 
export TMPDIR=$6  #TEMP_DIR
 
#Load specified module versions
module load python
module load rseqc/$5


####INPUTS: $1: STAR_RESULTS $2 QC_PATH $3 BAM_FILE_NAME $4 RSEQC_REF $5 RSEQC_VERSION $6 TMP_DIR $7 SAMPLE 


inner_distance.py -i $1/$7/$3 -o $2/$7/insert_size -r $4 

wait

exit 0