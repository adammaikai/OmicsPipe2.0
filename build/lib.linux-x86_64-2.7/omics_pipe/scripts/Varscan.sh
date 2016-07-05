#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load samtools/$4
module load varscan/$5
module load R/$9
module load vcftools/${10}
module load vcflib/1.0

####INPUTS: $1=sample $2=TEMP_DIR $3=GENOME $4=SAMTOOLS_VERSION $5=VARSCAN_VERSION $6=ALIGNMENT_DIR $7=RESULTS $8=OPTIONS $9=R_VERSION $10=VCFTOOLS_VERSION

mkdir -p "$7/$1"

## varscan2 call variants
# VarScan.v$5.jar \
# mpileup2snp \
# $6/$1/$1\.flt.pileup \
# $8 \
# > $2/$1\_varscan_snp.vcf

# VarScan.v$5.jar \
# mpileup2indel \
# $6/$1/$1\.flt.pileup \
# $8 \
# > $2/$1\_varscan_indel.vcf

# vcfcombine $2/$1\_varscan_snp.vcf $2/$1\_varscan_indel.vcf > $7/$1/$1\_varscan.4.1.vcf

# rm $2/$1\_varscan_snp.vcf
# rm $2/$1\_varscan_indel.vcf

## varscan2 call variants
VarScan.v$5.jar \
mpileup2cns \
$6/$1/$1\.flt.pileup \
$8 \
> $7/$1/$1\_varscan.vcf


bgzip -f $7/$1/$1\_varscan.vcf
tabix -f -p vcf $7/$1/$1\_varscan.vcf.gz

exit 0
