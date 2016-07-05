#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load R/$2

####INPUTS: $1=sample 

Rscript ~/.virtualenvs/pm/omics_pipe/omics_pipe/scripts/gep.R  

exit 0
