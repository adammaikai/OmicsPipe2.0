#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1
 
#Load specified cufflinks version
module load cufflinks/$7


#Run cuffdiff $1=CUFFDIFF_RESULTS $2=GENOME, $3 =CUFFDIFF_OPTIONS (-p 8 -FDR 0.01 -L Normal, OA -N --compatible-hits-norm); $4=CUFFMERGE_RESULTS $5=CUFFDIFF_INPUT_LIST_COND1 $6CUFFDIFF_INPUT_LIST_COND2 $7 CUFFLINKS_VERSION
cuffdiff $3 -o $1 -b $2 -u $4/merged.gtf \
$5 \
$6


exit 0
