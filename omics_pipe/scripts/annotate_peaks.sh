#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$4  #TEMP_DIR
 
#Loads specified module versions
module load homer/$3

#Calls scripts from HOMER to call peaks
#1PAIR_LIST  $2HOMER_RESULTS $3HOMER_VERSION $4 TEMP_DIR $5 GENOME $6 HOMER_ANNOTATE_OPTIONS
#PAIR LIST NEEDS TO BE IN THIS FORMAT 'pair1_chip-pair1_input pair2_chip-pair2_input'
array=($1)  

for pair in ${array[*]}; do
    echo -en $pair"\t"
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
# peak annotation
	annotatePeaks.pl $2/${chip}_peaks.txt $5 $6 -go $2/${chip}_regions_GO -genomeOntology $2/${chip}_regions_GenomeO > $2/${chip}_regions_annot.txt
done




exit 0