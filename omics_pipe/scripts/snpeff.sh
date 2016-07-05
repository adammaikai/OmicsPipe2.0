#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $1/$2

#Move tmp dir to scratch 
export TMPDIR=$3  #TEMP_DIR
 
#Load specified module versions

module load snpeff/$5 #3.3e ####################


out_dir=$1/$2
script_path=$4



####INPUTS: $1: VARIANT_RESULTS   $2 SAMPLE  $3: TEMP_DIR $4 WORKING_DIR #$5 SNPEFF_VERSION 

#temporary fix for header
#cat $out_dir/raw_variants.vcf | grep '#' > $out_dir/header.vcf


# use snpeff to annotate the variants
java -jar `which snpEff.jar` \
	-c $script_path/snpEff.config \
	-v hg19 \
	-s $out_dir/snpEffSummary.html \
	$out_dir/$2.vcf > $out_dir/$2_final_variants_filt_snpeff.vcf
	

exit 0
