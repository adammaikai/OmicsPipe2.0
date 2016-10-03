#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load qualimap/$4

mkdir -p "$3/$1"

####INPUTS $1=sample $2=ALIGNMENT_DIR $3=RESULTS $4=QUALIMAP_VERSION $5=ANALYSIS_TYPE $6=OPTIONS $7=CAPTURE_KIT_BED_DNA

if [ "$5" == "DNA" ]; then
qualimap $6 \
-gff $7 \
-bam $2/$1/$1\_gatk_recal.bam \
-outdir $3/$1
else
qualimap $6 \
-bam $2/$1/$1\_Aligned.sortedByCoord.out.bam \
-outdir $3/$1
fi

exit 0
