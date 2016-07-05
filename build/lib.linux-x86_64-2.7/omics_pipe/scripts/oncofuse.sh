#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make output directory if it doesn't exist
mkdir -p $2/$1

#Move tmp dir to scratch 
export TMPDIR=$3/$BASHPID  #TEMP_DIR
scratch=$TMPDIR

#Load specified oncofuse version
module load oncofuse/$4


#Runs oncofuse with $1= sample $2=FUSION_RESULTS, $3 TEMP_DIR, $4 ONCOFUSE_VERSION

cat $2/$1/final-list_candidate-fusion-genes.txt | cut -f 9,10 | tail -n +2 | sed 's/:+//g' | sed 's/:-//g' | awk -v tissue=EPI '{print "chr"$1, "\t", "chr"$2, "\t", tissue}' | sed -e 's+:+\t+g' > $2/$1/final_list_coord.txt


java -Xmx1G -jar `which Oncofuse.jar` $2/$1/final_list_coord.txt coord - $2/$1/oncofuse_res.txt

exit 0
