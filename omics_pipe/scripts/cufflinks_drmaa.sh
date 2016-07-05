#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $3
 
#Load specified cufflinks version
module load cufflinks/$7


#Run Cufflinks with $1=SAMPLE, $2=TOPHAT_RESULTS, $3=CUFFLINKS_RESULTS, $4 =REF_GENES, $5=GENOME, $6 = CUFFLINKS_OPTIONS, $7=CUFFLINKS_VERSION 
cufflinks $6 -g $4 -b $5 -o $3/$1 $2/$1/accepted_hits.bam

exit 0
