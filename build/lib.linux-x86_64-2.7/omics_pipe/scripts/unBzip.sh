#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$2  #TEMP_DIR

for FILE in $1;
do
	bunzip2 ${FILE}.bz2
	echo "${FILE} has been unzipped"
done

exit 0
