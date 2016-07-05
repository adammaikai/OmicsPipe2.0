#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $3
 
#Load specified 
module load python

#Run script  $1: script to run $2 input parameters
python $1 $2

exit 0
