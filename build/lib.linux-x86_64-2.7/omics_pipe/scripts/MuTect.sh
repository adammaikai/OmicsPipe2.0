#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load mutect/$8

####INPUTS: $1=sample $2=RESULTS $3=COSMIC $4=DBSNP $5=CAPTURE_KIT_BED $6=ALIGNMENT_DIR $7=GENOME $8=MUTECT_VERSION $9=DNA_NORMAL_EXT $10=DNA_TUMOR_EXT

mkdir -p "$2"
mkdir -p "$2/$1"

normal=$(echo $9 | sed -e 's/-DNA//' | cut -c 2-)
tumor=$(echo ${10} | sed -e 's/-DNA//' | cut -c 2-)

if [ -z "$5" ]; then
    java -Xmx8g -jar `which mutect-1.1.7.jar` \
    --analysis_type MuTect \
    --reference_sequence $7 \
    --cosmic $3 \
    --dbsnp $4 \
    --input_file:normal $6/$1$9/$1$9_gatk_recal.bam \
    --input_file:tumor $6/$1${10}/$1${10}_gatk_recal.bam \
    --vcf $2/$1/$1\_$tumor\_$normal\_mutect.raw.vcf \
    -rf BadCigar \
    --out $2/$1/$1\_call_stats.txt \
    --coverage_file $2/$1/$1\_coverage.wig.txt
else
    java -Xmx8g -jar `which mutect-1.1.7.jar` \
    --analysis_type MuTect \
    --reference_sequence $7 \
    --cosmic $3 \
    --dbsnp $4 \
    --intervals $5 \
    --input_file:normal $6/$1$9/$1$9_gatk_recal.bam \
    --input_file:tumor $6/$1${10}/$1${10}_gatk_recal.bam \
    --vcf $2/$1/$1\_$tumor\_$normal\_mutect.raw.vcf \
    -rf BadCigar \
    --out $2/$1/$1\_call_stats.txt \
    --coverage_file $2/$1/$1\_coverage.wig.txt
fi

grep -v REJECT $2/$1/$1\_$tumor\_$normal\_mutect.raw.vcf > $2/$1/$1\_$tumor\_$normal\_mutect.filt.vcf

bgzip -f $2/$1/$1\_$tumor\_$normal\_mutect.raw.vcf
tabix -f -p vcf $2/$1/$1\_$tumor\_$normal\_mutect.raw.vcf.gz

bgzip -f $2/$1/$1\_$tumor\_$normal\_mutect.filt.vcf
tabix -f -p vcf $2/$1/$1\_$tumor\_$normal\_mutect.filt.vcf.gz

rm $2/$1/$1\_call_stats.txt

exit 0
