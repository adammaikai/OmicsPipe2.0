#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Move tmp dir to scratch 
export TMPDIR=$7 #TEMP_DIR

#Make output directory if it doesn't exist
mkdir -p $6/$1

#Get number of threads from node
ncpu=$(grep -c "processor" /proc/cpuinfo) 

#Load specific module versions
module load samtools/$8
module load annovar/$9
module load vcftools/${10}
module load varscan/${11}


#Variant calling with $1=SAMPLE, $2=STAR_RESULTS (location of bam file), $3 GENOME, $4 VARSCAN_PATH, $5 VARSCAN_OPTIONS, $6=VARIANT_RESULTS, $7=TEMP_DIR, $8=SAMTOOLS_VERSION, $9=ANNOVAR_VERSION, $10=VCFTOOLS_VERSION, $11=VARSCAN_VERSION, $12=SAMTOOLS_OPTIONS
samtools view -H $2/$1/Aligned.out.sorted.bam | \
grep "\@SQ" | \
sed 's/^.*SN://g' | \
cut -f 1 | \
xargs -I {} -n 1 -P $ncpu sh -c \
"samtools mpileup ${12} -f $3 -r {} $2/$1/Aligned.out.sorted.bam | java -jar $4 mpileup2snp - $5 > $6/$1/tmpvar.{}.vcf" &

#Indel calling
samtools view -H $2/$1/Aligned.out.sorted.bam | \
grep "\@SQ" | \
sed 's/^.*SN://g' | \
cut -f 1 | \
xargs -I {} -n 1 -P $ncpu sh -c \
"samtools mpileup ${12} -f $3 -r {} $2/$1/Aligned.out.sorted.bam | java -jar $4 mpileup2indel - $5 > $6/$1/tmpind.{}.vcf" &



wait #Waits for variant calling & indel calling process to finish up, then joins the temporary files...

#Join the split chr together to one final variant .vcf file
list=$(samtools view -H $2/$1/Aligned.out.sorted.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1)
gawk '{if(FILENAME==ARGV[1]){print}else{if($1 !~ "^#"){print}}}' `for j in $list; do echo $6/$1/tmpvar.${j}.vcf; done ` > $6/$1/variants.vcf &      
        
#Join the split chr together to one final indel .vcf
list=$(samtools view -H $2/$1/Aligned.out.sorted.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1)
gawk '{if(FILENAME==ARGV[1]){print}else{if($1 !~ "^#"){print}}}' `for j in $list; do echo $6/$1/tmpind.${j}.vcf; done ` > $6/$1/indel.vcf &
        
#Wait for the processes to finish up
wait 

#Check if variants.vcf and indel.vcf have the header
if  ! head -n 1 $6/$1/variants.vcf | grep -ci '##fileformat=VCFv4.1'; then  
	header='##fileformat=VCFv4.1\n##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1'
	echo -e $header | cat - $6/$1/variants.vcf > $6/$1/variants.new.vcf
	mv $6/$1/variants.new.vcf $6/$1/variants.vcf
fi

if  ! head -n 1 $6/$1/indel.vcf | grep -ci '##fileformat=VCFv4.1'; then  
	header='##fileformat=VCFv4.1\n##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1'
	echo -e $header | cat - $6/$1/indel.vcf > $6/$1/indel.new.vcf
	mv $6/$1/indel.new.vcf $6/$1/indel.vcf
fi



#join the vcf files
vcf-concat $6/$1/variants.vcf $6/$1/indel.vcf > $6/$1/$1_conc.vcf
        
#Set sample name in vcf file to patient ID (VarScan default is Sample1)
cat $6/$1/$1_conc.vcf | sed -e "s/Sample1/$1/g" > $6/$1/$1_rename.vcf


vcf-sort -c $6/$1/$1_rename.vcf > $6/$1/$1.vcf


# Clean up 
rm $6/$1/tmp*.*.vcf



exit 0
