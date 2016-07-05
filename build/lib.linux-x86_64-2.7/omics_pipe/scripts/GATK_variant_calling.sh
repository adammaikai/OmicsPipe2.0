#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load gatk/$4
module load R/$5
module load samtools
cpu=$(grep -c "processor" /proc/cpuinfo)

mkdir -p "$7/$1"

####INPUTS: $1=sample $2=TEMP_DIR $3=GENOME $4=GATKVERSION $5=R_VERSION $6=ALIGNMENT_DIR $7=RESULTS $8=OMNI $9=HAPMAP $10=DBSNP $11=MILLS $12=G1000 $13=KG $14=CAPTURE_KIT_BED 

## merge input bam with 1000 Genomes cohort for group variant calling and recalibration (need 30 exomes to perform correctly).
# samtools merge $7/$1/$1_G1000.merged.bam $7/$1/$1\_gatk_recal.bam /data/s3/averaprojects/G1000/alignments/G1000.merged.mapped.ILLUMINA.bwa.GBR.exome.bam
# samtools index $7/$1/$1_G1000.merged.bam

## call variants
if [ -z "${14}" ]; then
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T HaplotypeCaller \
	-nct $cpu \
	-R $3 \
	--dbsnp ${10} \
	-I $6/$1/$1_gatk_recal.bam \
	-rf BadCigar \
	--genotyping_mode DISCOVERY \
	-stand_emit_conf 10 \
	-stand_call_conf 30 \
	-o $2/$1\_gatk.vcf

	## recalibrate variant quality scores
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T VariantRecalibrator \
	-R $3 \
	-input $2/$1\_gatk.vcf \
	-nt 1 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $9 \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $8 \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${12} \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${10} \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-rf BadCigar \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	-recalFile $2/$1\_gatk_recalibrate_SNP.recal \
	-tranchesFile $2/$1\_gatk_recalibrate_SNP.tranches \
	-rscriptFile $2/$1\_gatk_recalibrate_SNP_plots.R

	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T ApplyRecalibration \
	-R $3 \
	-input $2/$1\_gatk.vcf \
	-tranchesFile $2/$1\_gatk_recalibrate_SNP.tranches \
	-recalFile $2/$1\_gatk_recalibrate_SNP.recal \
	-o $2/$1\_gatk_recalibrate_SNP.vcf \
	--ts_filter_level 99.5 \
	-mode SNP

	if [ -f $2/$1\_gatk_recalibrate_SNP.vcf ]; then
		rm $2/$1\_gatk.vcf
		rm $2/$1\_gatk_recalibrate_SNP.recal
		rm $2/$1\_gatk_recalibrate_SNP.tranches
		rm $2/$1\_gatk_recalibrate_SNP_plots.R
	else
		bgzip -f -c $2/$1\_gatk.vcf > $7/$1/$1\_gatk.vcf.gz
		tabix -f -p vcf $7/$1/$1\_gatk.vcf.gz
		rm $2/$1\_gatk.vcf
	fi

	java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
	-T VariantRecalibrator \
	-R $3 \
	-input $2/$1\_gatk_recalibrate_SNP.vcf \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 ${11} \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${10} \
	--maxGaussians 4 \
	--minNumBadVariants 5000 \
	-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	-recalFile $2/$1\_gatk_recalibrate_INDEL.recal \
	-tranchesFile $2/$1\_gatk_recalibrate_INDEL.tranches \
	-rscriptFile $2/$1\_gatk_recalibrate_INDEL_plots.R

	java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
	-T ApplyRecalibration \
	-R $3 \
	-input $2/$1\_gatk_recalibrate_SNP.vcf \
	-tranchesFile $2/$1\_gatk_recalibrate_INDEL.tranches \
	-recalFile $2/$1\_gatk_recalibrate_INDEL.recal \
	-o $2/$1\_gatk_vqsr.vcf \
	--ts_filter_level 99.0 \
	-mode INDEL

	if [ -f $2/$1\_gatk.vcf ]; then 
		java -Xmx4g -jar `which GenomeAnalysisTK.jar` \
		-T SelectVariants \
		-nt $cpu \
		-R $3 \
		--excludeNonVariants \
		--excludeFiltered \
		--variant $2/$1\_gatk.vcf \
		--out $7/$1/$1\_gatk.filt.vcf

		bgzip -f $7/$1/$1\_gatk.filt.vcf
		tabix -f -p vcf $7/$1/$1\_gatk.filt.vcf.gz
	fi

	if [ -f $2/$1\_gatk_vqsr.vcf ]; then
		rm $2/$1\_gatk_recalibrate_indel.vcf
		rm $2/$1\_gatk_recalibrate_indel.recal
		rm $2/$1\_gatk_recalibrate_indel.tranches
		rm $2/$1\_gatk_recalibrate_indel_plots.R

		java -Xmx4g -jar `which GenomeAnalysisTK.jar` \
		-T SelectVariants \
		-nt $cpu \
		-R $3 \
		--excludeNonVariants \
		--excludeFiltered \
		--variant $2/$1\_gatk_vqsr.vcf \
		--out $7/$1/$1\_gatk_vqsr.filt.vcf

		bgzip -f $7/$1/$1\_gatk_vqsr.filt.vcf
		tabix -f -p vcf $7/$1/$1\_gatk_vqsr.filt.vcf.gz
	fi

else
	java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
	-T HaplotypeCaller \
	-nct $cpu \
	-R $3 \
	--dbsnp ${10} \
	-I $6/$1/$1_gatk_recal.bam \
	-L ${14} \
	-ip 1 \
	-rf BadCigar \
	--genotyping_mode DISCOVERY \
	-stand_emit_conf 10 \
	-stand_call_conf 30 \
	-o $2/$1\_gatk.vcf

	java -Xmx8g -jar `which GenomeAnalysisTK.jar` \
	-T SelectVariants \
	-nt $cpu \
	-R $3 \
	--excludeNonVariants \
	--excludeFiltered \
	--variant $2/$1\_gatk.vcf \
	--out $7/$1/$1\_gatk.filt.vcf

	bgzip -f $7/$1/$1\_gatk.filt.vcf
	tabix -f -p vcf  $7/$1/$1\_gatk.filt.vcf.gz
	rm $2/$1\_gatk.vcf
fi


exit 0


