#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load mojo/$4
module load fusionannotator/$9
cpu=$(grep -c "processor" /proc/cpuinfo)

mkdir -p "$3/$1"

####INPUTS $1=sample $2=RAW_DATA $3=RESULTS $4=MOJO_VERSION $5=CONFIG $6=CPU $7=MEMORY $8=ALIGNMENT_DIR $9=FUSIONANNOTOR_VERSION $10=FUSIONANNOTOR_LIB

MOJO \
--config $5 \
--sample_name $1 \
--output_dir $3/$1 \
--fq1 $2/$1/$1_1.fastq.gz \
--fq2 $2/$1/$1_2.fastq.gz \
--cores $6 \
--mem 50

## Extract intersection of fusions called from Star-Fusion and Mojo
cat <(awk '{print $1}' $8/$1/star-fusion.fusion_candidates.final.abridged | sort | uniq -c | awk '{print $2}') \
<(cat $3/$1/$1/$1\.fusions | awk -F '\t' '{print $1}' | sed -e 's/_/--/g' | sort | uniq -c | awk '{print $2}') \
| sort | uniq -c | grep -v '^ *1 ' | awk '{print $2}' > $3/$1/$1/$1\_fusions_merged.txt

#Fusion Annotator
/opt/software/fusionannotator/FusionAnnotator-0.0.2/FusionAnnotator \
--fusion_annot_lib ${10} \
--annotate $3/$1/$1/$1\_fusions_merged.txt \
>  $3/$1/$1/$1\_fusion_annotations.txt



exit 0