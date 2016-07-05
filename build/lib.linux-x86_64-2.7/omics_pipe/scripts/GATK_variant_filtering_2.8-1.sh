#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR

#Loads specified module versions
module load gatk/$4
module load R/${10}


####INPUTS: $1: VARIANT_RESULTS $2 SAMPLE $3: TEMP_DIR $4 GATK_VERSION $5 GENOME $6 DBSNP $7 MILLS_G1000 $8 OMNI $9 hapmap $10 R_VERSION $11 G1000_SNPs $12 G1000_Indels 

dbsnp=$6
mills=$7
omni=$8
hapmap=$9
temp=$3
g1000snps=${11}
g1000indels=${12}
######################## GATK VQSR FILTRATION ###############################
# Create a Gaussian mixture model for SNPs
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
  -nt 7 \
  -R $5 \
  -T VariantRecalibrator \
  --maxGaussians 4 \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000snps \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
  -mode SNP \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -input  $1/$2/$2.raw.vcf \
  -recalFile $1/$2/$2.snp.recal.vcf \
  -tranchesFile $1/$2/$2.snp.tranches \
  -rscriptFile $1/$2/$2.snp.r

# Create a Gaussian mixture model for INDELs
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
  -nt 7 \
  -R $5 \
  -T VariantRecalibrator \
  --maxGaussians 4 \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000indels \
  -an DP -an FS -an ReadPosRankSum -an MQRankSum \
  -mode INDEL \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -input $1/$2/$2.raw.vcf \
  -recalFile $1/$2/$2.indel.recal.vcf \
  -tranchesFile $1/$2/$2.indel.tranches \
  -rscriptFile $1/$2/$2.indel.r

# Apply the model for SNPs
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
  -nt 7 \
  -R $5 \
  -T ApplyRecalibration \
  -mode SNP \
  --ts_filter_level 99.0 \
  -input $1/$2/$2.raw.vcf \
  -recalFile $1/$2/$2.snp.recal.vcf \
  -tranchesFile $1/$2/$2.snp.tranches \
  -o $1/$2/$2.snpAr.vcf

# Apply the model for INDELS
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
  -nt 7 \
  -R $5 \
  -T ApplyRecalibration \
  -mode indel \
  --ts_filter_level 99.0 \
  -input $1/$2/$2.snpAr.vcf \
  -recalFile $1/$2/$2.indel.recal.vcf \
  -tranchesFile $1/$2/$2.indel.tranches \
  -o $1/$2/$2.snpAr.indelAr.vcf

# Select the variants that passed the models
java -Xms454m -Xmx3181m -jar `which GenomeAnalysisTK.jar` \
  -nt 7 \
  -R $5 \
  -T SelectVariants \
  --excludeNonVariants \
  --excludeFiltered \
  --variant $1/$2/$2.snpAr.indelAr.vcf \
  --out $1/$2/$2.vqsr.vcf
  
exit 0