#!/bin/bash
set -x


#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/${2}

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR
 
#Load specified module versions
module load bedtools/$4
module load snpeff/$5 #3.3e ####################
module load vcftools/$7
module load R/${10}


ncpu=$(grep -c "processor" /proc/cpuinfo) 
nthreads=$((ncpu/2))


extractvcf=$8/extractvcf.R
out_dir=$1/$2
dbNSFP=$6
script_path=$8
no_impact=$9
sample=$2

####INPUTS: $1: VARIANT_RESULTS $2 SAMPLE $3: TEMP_DIR  $4 BEDTOOLS_VERSION  #5 SNPEFF_VERSION $6 dbNSFP $7 VCFTOOLS_VERSION $8 WORKING_DIR $9 SNP_FILTER_OUT_REF $10 R_VERSION

#temporary fix for header
cat $out_dir/$sample.vcf | grep '#' > $out_dir/header.vcf

################################################################################
# filter variants
################################################################################
# filter by ncbi common variants
sed "s/chr//g" $out_dir/$sample.vcf > $out_dir/final_variants_nochr.vcf
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
	$out_dir/final_variants_filt_snpeff.vcf > $out_dir/final_variants_filt_snpeff_dbnsfp.vcf  #this works

###TO DO: POSSIBLY ADD ANNOVAR ANNOTATION HERE
	
# remove 1000 genomes maf > 0.01
java -jar `which SnpSift.jar` filter \
	-f $out_dir/final_variants_filt_snpeff_dbnsfp.vcf \
	" ( exists dbNSFP_1000Gp1_AF ) " \
	> $out_dir/maf_1kg_temp.vcf

vcftools --vcf $out_dir/maf_1kg_temp.vcf --get-INFO dbNSFP_1000Gp1_AF --out $out_dir/maf
cut $out_dir/maf.INFO -f 5 > $out_dir/maf.txt
skipn=$(cat $out_dir/maf_1kg_temp.vcf | grep '#' | wc -l)
cat $out_dir/final_variants_filt_snpeff_dbnsfp.vcf | grep '#' > $out_dir/header.vcf 
Rscript $script_path/removeLines.R $out_dir $out_dir/maf.txt $out_dir/maf_1kg_temp.vcf $skipn

bedtools intersect -v -a $out_dir/final_variants_filt_snpeff_dbnsfp.vcf -b $out_dir/maf_filt.vcf > $out_dir/intogen_body.vcf
cat $out_dir/header.vcf $out_dir/intogen_body.vcf > $out_dir/intogen_input.vcf

#sort the file
vcf-sort $out_dir/intogen_input.vcf > $out_dir/intogen_input.vcf.temp && mv $out_dir/intogen_input.vcf.temp $out_dir/intogen_input.vcf
exit 0
