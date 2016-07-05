#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$5  #TEMP_DIR
 
#Loads specified module versions
module load homer/$4

#Calls scripts from HOMER to call peaks
#1PAIR_LIST  $2HOMER_RESULTS $3 HOMER_PEAKS_OPTIONS $4HOMER_VERSION $5 TEMP_DIR

#PAIR LIST NEEDS TO BE IN THIS FORMAT 'pair1_chip-pair1_input pair2_chip-pair2_input'

array=($1)  

for pair in ${array[*]}; do
    echo -en $pair"\t"
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
	findPeaks $2/${chip}_tag $3 -i $2/${input}_tag
done
exit 0