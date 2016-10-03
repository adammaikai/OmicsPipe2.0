#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load R/$3
module load vcflib/$7

####INPUTS: $1=sample $2=VARIANT_DIR $3=R_VERSION $4=tumor $5=normal $6=VCFLIB_VERSION $7=GENOME

normal=$(echo $5 | sed -e 's/-DNA//' | cut -c 2-)
tumor=$(echo $4 | sed -e 's/-DNA//' | cut -c 2-)

cat $2/$1/$1\_$tumor\_$normal\_varscan_somatic.snp.Somatic.hc.vcf | cut -f 11 | awk -F ":" '{ print $7}' > $2/$1/$1\_$tumor\_$normal\_varscan_somatic.snp.strand.dp.csv
cat $2/$1/$1\_$tumor\_$normal\_varscan_somatic.indel.Somatic.hc.vcf | cut -f 11 | awk -F ":" '{ print $7}' > $2/$1/$1\_$tumor\_$normal\_varscan_somatic.indel.strand.dp.csv

vcfintersect -r $6 -u $2/$1/$1\_$tumor\_$normal\_varscan_somatic.vcf.gz $2/$1/$1\_$tumor\_$normal\_mutect.filt.vcf.gz > $2/$1/$1\_$tumor\_$normal\_somatic_merged.vcf
bgzip -f $2/$1/$1\_$tumor\_$normal\_somatic_merged.vcf 
tabix -p vcf $2/$1/$1\_$tumor\_$normal\_somatic_merged.vcf.gz

#Rscript ~/.virtualenvs/op2/OmicsPipe2.0/omics_pipe/scripts/somaticAnnotation.R $2/$1/$1\_$tumor\_$normal\_varscan_somatic.vcf.gz $2/$1/$1\_$tumor\_$normal\_mutect.filt.vcf.gz TRUE

exit 0
