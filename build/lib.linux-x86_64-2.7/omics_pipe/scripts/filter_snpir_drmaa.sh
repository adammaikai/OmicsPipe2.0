#!/bin/bash
set -x

####INPUTS: $1: SNIPR_RESULTS $2: TEMP_DIR $3 SAMTOOLS_VERSION $4 BWA_VERSION $5 PICARD_VERSION $6 GATK_VERSION $7 BEDTOOLS_VERSION $8 UCSC_TOOLS_VERSION
#$9 GENOME $10 REPEAT_MASKER $11 REF_GENES $12 RNA EDIT $13 DBSNP $14 MILLS $15 G1000 $16 WORKING_DIR $17 SAMPLE $18 BWA_RESULTS $19 SNIPR_VERSION $20 SNIPR_CONFIG, $21 SNIPR_DIR


#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/${17}

#Move tmp dir to scratch 
export TMPDIR=$2  #TEMP_DIR
 
#Load specified module versions
module load samtools/$3
module load bwa/$4
module load picard/$5
module load gatk/$6
module load bedtools/$7
module load ucsc_tools/$8 
module load snpir/${19}
module load snpeff/${22} #3.3e ####################
module load vcftools/${24}
export PERL5LIB=${20}

ncpu=$(grep -c "processor" /proc/cpuinfo) 
nthreads=$((ncpu/2))

hg19_reference=$9
RepeatMasker=${10} 
gene_annotation=${11} 
rnaedit=${12}
dbsnp=${13}
mills=${14}
g1000=${15}
extractvcf=${16}/extractvcf.R
out_dir=$1/${17}
SNPiR=${21}
dbNSFP=${23}
script_path=${25}
no_impact=${26}
cp $SNPiR/convertCoordinates.* $out_dir


####INPUTS: $1: SNIPR_RESULTS $2: TEMP_DIR $3 SAMTOOLS_VERSION $4 BWA_VERSION $5 PICARD_VERSION $6 GATK_VERSION $7 BEDTOOLS_VERSION $8 UCSC_TOOLS_VERSION
#$9 GENOME $10 REPEAT_MASKER $11 REF_GENES $12 RNA EDIT $13 DBSNP $14 MILLS $15 G1000 $16 WORKING_DIR $17 SAMPLE $18 BWA_RESULTS $19 SNIPR_VERSION $20 SNIPR_CONFIG, $21 SNIPR_DIR
#$22 SNPEFF_VERSION $23 dbNSFP $24 VCFTOOLS_VERSION $25 WORKING_DIR $26 SNP_FILTER_OUT_REF 

#temporary fix for header
cat $out_dir/raw_variants.vcf | grep '#' > $out_dir/header.vcf

################################################################################
# filter variants
################################################################################
# filter by ncbi common variants
sed "s/chr//g" $out_dir/final_variants.vcf > $out_dir/final_variants_nochr.vcf
bedtools intersect -v -a $out_dir/final_variants_nochr.vcf -b $no_impact > $out_dir/final_variants_filt.temp.vcf
cat $out_dir/header.vcf $out_dir/final_variants_filt.temp.vcf > $out_dir/final_variants_filt.vcf

# use snpeff to annotate the variants
java -jar `which snpEff.jar` \
	-c $script_path/snpEff.config \
	-v hg19 \
	-s $out_dir/snpEffSummary.html \
	$out_dir/final_variants_filt.vcf > $out_dir/final_variants_filt_snpeff.vcf
	
# annotate using dbNSFP
java -jar `which SnpSift.jar` dbnsfp \
	-v $dbNSFP \
	$out_dir/final_variants_filt_snpeff.vcf > $out_dir/final_variants_filt_snpeff_dbnsfp.vcf
	
# filter based on snpeff high & moderate impact classes
#java -jar `which SnpSift.jar` filter \
	-f $out_dir/final_variants_filt_snpeff_dbnsfp.vcf \
	" ( EFF[*].IMPACT = 'HIGH' | EFF[*].IMPACT = 'MODERATE' ) " \
	> $out_dir/final_variants_filt_snpeff_dbnsfp_filtered.vcf
	
# remove 1000 genomes maf > 0.01
java -jar `which SnpSift.jar` filter \
	-f $out_dir/final_variants_filt_snpeff_dbnsfp.vcf \
	" ( exists dbNSFP_1000Gp1_AF ) " \
	> $out_dir/maf_1kg_temp.vcf

vcftools --vcf $out_dir/maf_1kg_temp.vcf --get-INFO dbNSFP_1000Gp1_AF --out $out_dir/maf
cut $out_dir/maf.INFO -f 5 > $out_dir/maf.txt
skipn=$(cat $out_dir/maf_1kg_temp.vcf | grep '#' | wc -l)
cat $out_dir/final_variants_filt_snpeff_dbnsfp_filtered.vcf | grep '#' > $out_dir/header.vcf 
Rscript $script_path/removeLines.R $out_dir $out_dir/maf.txt $out_dir/maf_1kg_temp.vcf $skipn

bedtools intersect -v -a $out_dir/final_variants_filt_snpeff_dbnsfp_filtered.vcf -b $out_dir/maf_filt.vcf > $out_dir/intogen_body.vcf
cat $out_dir/header.vcf $out_dir/intogen_body.vcf > $out_dir/intogen_input.vcf

#sort the file
vcf-sort $out_dir/intogen_input.vcf > $out_dir/intogen_input.vcf.temp && mv $out_dir/intogen_input.vcf.temp $out_dir/intogen_input.vcf
exit 0
