#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1
 
#Load specified cufflinks version
module load cufflinks/$5


#Run cuffmerge with  $1=CUFFMERGE_RESULTS, $2 =REF_GENES, $3=GENOME, $4 = CUFFLINKS_OPTIONS, $5=CUFFLINKS_VERSION 


cuffcompare $4 -o $1 -s $3 -r $2 $1/merged.gtf

exit 0
