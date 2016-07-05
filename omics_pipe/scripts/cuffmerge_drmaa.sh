#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $2
 
#Load specified cufflinks version
module load cufflinks/$6


#Run cuffmerge with  $1=GTF_LIST, $2=CUFFMERGE_RESULTS, $3 =REF_GENES, $4=GENOME, $5 = CUFFLINKS_OPTIONS, $6=CUFFLINKS_VERSION 


cuffmerge $5 -g $3 -o $2 -s $4 $1  


exit 0
