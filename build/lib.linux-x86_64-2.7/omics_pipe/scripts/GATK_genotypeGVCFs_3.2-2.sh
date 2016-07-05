#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR

#Loads specified module versions
module load gatk/$4


####INPUTS: $1: VARIANT_RESULTS $2 SAMPLE_LIST $3: TEMP_DIR $4 GATK_VERSION $5 GENOME $6 DBSNP
temp=$3
echo $2
echo type $2
string=$2
OIFS=$IFS;
IFS=$'\t';
SAMPLE_LIST=($string)
for ((i=0; i<${#SAMPLE_LIST[@]}; ++i));
do
    echo "sample $i: ${SAMPLE_LIST[$i]}";
done

IFS=$OIFS;

echo $SAMPLE_LIST
###Create variant argument for all samples in SAMPLE_LIST
for i in "${SAMPLE_LIST[@]}"; do
	sample="--variant $1/$i/$i.raw.g.vcf "
	variable=$variable$sample
done
echo $variable

###########GATK CombineGVCFs########################################
java -Xms454m -Xmx3181m -Djava.io.tmpdir=$temp \
	-jar `which GenomeAnalysisTK.jar` \
	-T CombineGVCFs \
	-R $5 \
	$variable \
	-o $1/mergeGvcf.vcf


######################## GATK GenotypeGVCFs ###############################
java -Xms454m -Xmx3181m -Djava.io.tmpdir=$temp \
	-jar `which GenomeAnalysisTK.jar` \
	-T GenotypeGVCFs \
	-R $5 \
	--variant $1/mergeGvcf.vcf \
    -o $1/GenotypeGVCF_output.vcf \
    --dbsnp $6
exit 0

