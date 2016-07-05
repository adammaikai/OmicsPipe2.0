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
#module load samtools/dnanexus-1.0
module load bwa/$4
module load picard/$5
module load gatk/$6
module load bedtools/$7
module load ucsc_tools/$8 
module load snpir/${19}
module load R
 
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
encoding=${22} #phred33 or #phred64
#encoding=phred64


cp $SNPiR/convertCoordinates.* $out_dir


################################################################################
# MAPPING
################################################################################


# merge alignments
cat ${18}/${17}_1/${17}_1.sam <(grep -v '^@' ${18}/${17}_2/${17}_2.sam) > $out_dir/merged.sam

# convert the position of reads that map across splicing junctions onto the genome
cd $out_dir
java -Xmx4g -cp . convertCoordinates < $out_dir/merged.sam > $out_dir/merged.conv.sam
cd ~

rm $out_dir/merged.sam # clean up

# rewrite the .sam header to drop the pseudochromosomes
samtools view -HS $out_dir/merged.conv.sam > $out_dir/header.sam
sed -n '1,25 p' $out_dir/header.sam > $out_dir/t1.sam
tail -n -2 $out_dir/header.sam > $out_dir/t2.sam
cat $out_dir/t1.sam $out_dir/t2.sam > $out_dir/new_header.sam

java -jar `which ReplaceSamHeader.jar` \
	INPUT=$out_dir/merged.conv.sam \
	HEADER=$out_dir/new_header.sam \
	OUTPUT=$out_dir/merged.conv.nh.sam 
	
rm $out_dir/merged.conv.sam $out_dir/new_header.sam # clean up

# sort file, filter out unmappes reads and reads with mapping quality < 20 and convert to .bam
# samtools view -bS -q 20 -F 4 $out_dir/merged.conv.nh.sam | samtools rocksort -@ $ncpu -m 1500M - $out_dir/merged.conv.nh.sort
samtools view -bS -q 20 -F 4 $out_dir/merged.conv.nh.sam | samtools sort - $out_dir/merged.conv.nh.sort

# remove duplicate reads ###THIS PART FINISHES SUCCESSFULLY
java -Xmx4g -jar `which MarkDuplicates.jar` \
	INPUT=$out_dir/merged.conv.nh.sort.bam \
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=LENIENT \
	AS=true \
	METRICS_FILE=$out_dir/SM1.dups \
	OUTPUT=$out_dir/merged.conv.nh.sort.rd.bam \
	TMP_DIR=$TMPDIR
	
rm $out_dir/merged.conv.nh.sort.bam $out_dir/merged.conv.nh.sam # clean up	

# index ###THIS PART WORKS
samtools index $out_dir/merged.conv.nh.sort.rd.bam

# indel realignment & base quality score recalibration ###THIS PART FINISHES SUCCESSFULLY
if [ $encoding == "phred64" ]; then
	java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
		-T RealignerTargetCreator \
		-R $hg19_reference \
		-I $out_dir/merged.conv.nh.sort.rd.bam \
		-o $out_dir/output.intervals \
		-known $mills \
		-known $g1000 \
		-nt 8 \
		--fix_misencoded_quality_scores \
		--filter_reads_with_N_cigar
	
	java -Xmx16g -Djava.io.tmpdir=$TMPDIR -jar `which GenomeAnalysisTK.jar` \
		-I $out_dir/merged.conv.nh.sort.rd.bam \
		-R $hg19_reference \
		-T IndelRealigner \
		-targetIntervals $out_dir/output.intervals \
		-o $out_dir/merged.conv.sort.rd.realigned.bam \
		-known $mills \
		-known $g1000 \
		--consensusDeterminationModel KNOWNS_ONLY \
		-LOD 0.4 \
		--fix_misencoded_quality_scores \
		--filter_reads_with_N_cigar
else
	java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
		-T RealignerTargetCreator \
		-R $hg19_reference \
		-I $out_dir/merged.conv.nh.sort.rd.bam \
		-o $out_dir/output.intervals \
		-known $mills \
		-known $g1000 \
		-nt 8 \
		--filter_reads_with_N_cigar
	
	java -Xmx16g -Djava.io.tmpdir=$TMPDIR -jar `which GenomeAnalysisTK.jar` \
		-I $out_dir/merged.conv.nh.sort.rd.bam \
		-R $hg19_reference \
		-T IndelRealigner \
		-targetIntervals $out_dir/output.intervals \
		-o $out_dir/merged.conv.sort.rd.realigned.bam \
		-known $mills \
		-known $g1000 \
		--consensusDeterminationModel KNOWNS_ONLY \
		-LOD 0.4 \
		--filter_reads_with_N_cigar
fi

###THIS PART FINISHES SUCCESSFULLY	
# Base Quality Score Recalibration
	#1 Analyze patterns of covariation in the sequence dataset
java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
    -T BaseRecalibrator \
    -R $hg19_reference \
    -I $out_dir/merged.conv.sort.rd.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $g1000 \
    -o $out_dir/recal_data.table \
    -nct $ncpu
    
	#2 Do a second pass to analyze covariation remaining after recalibration
###THIS PART FINISHES SUCCESSFULLY
java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
    -T BaseRecalibrator \
    -R $hg19_reference \
    -I $out_dir/merged.conv.sort.rd.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $g1000 \
    -BQSR $out_dir/recal_data.table \
    -o $out_dir/post_recal_data.table \
    -nct $ncpu

	#3 Generate before / after plots This part throws error
java -jar `which GenomeAnalysisTK.jar` \
    -T AnalyzeCovariates \
    -R $hg19_reference \
    -before $out_dir/recal_data.table \
    -after $out_dir/post_recal_data.table \
    -plots $out_dir/recalibration_plots.pdf
    
	#4 Apply the recalibration to your sequence data
java -jar `which GenomeAnalysisTK.jar` \
    -T PrintReads \
    -R $hg19_reference \
    -I $out_dir/merged.conv.sort.rd.realigned.bam \
    -BQSR $out_dir/recal_data.table \
    -o $out_dir/merged.conv.sort.rd.realigned.recal.bam
    
################################################################################
# Call/filter variants
################################################################################     
# call variants using UnifiedGenotyper
# as suggested by Robert Piskol by mail
java -Xmx16g -jar `which GenomeAnalysisTK.jar` -T UnifiedGenotyper \
	-R $hg19_reference \
	-I $out_dir/merged.conv.sort.rd.realigned.recal.bam \
	-stand_call_conf 0 \
	-stand_emit_conf 0 \
	--dbsnp $dbsnp \
	-out_mode EMIT_VARIANTS_ONLY \
	-rf BadCigar \
	-o $out_dir/raw_variants.vcf \
	-nt 1 \
	-nct 32

# do the filtering
# convert vcf format into custom SNPiR format and filter variants with quality <20
$SNPiR/convertVCF.sh $out_dir/raw_variants.vcf $out_dir/raw_variants.txt 20

# filter mismatches at read ends
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format
	$SNPiR/filter_mismatch_first6bp.pl \
	-infile $out_dir/raw_variants.txt \
	-outfile $out_dir/raw_variants.rmhex.txt \
	-bamfile $out_dir/merged.conv.sort.rd.realigned.bam

# filter variants in repetitive regions
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $out_dir/raw_variants.rmhex.txt | \
	intersectBed -a stdin -b $RepeatMasker -v | \
	cut -f1,3-7 > $out_dir/raw_variants.rmhex.rmsk.txt

# filter intronic sites that are within 4bp of splicing junctions
# make sure your gene annotation file is in UCSC text format and sorted by chromosome and 
# transcript start position
$SNPiR/filter_intron_near_splicejuncts.pl \
	-infile $out_dir/raw_variants.rmhex.rmsk.txt \
	-outfile $out_dir/raw_variants.rmhex.rmsk.rmintron.txt \
	-genefile $gene_annotation

# filter variants in homopolymers
$SNPiR/filter_homopolymer_nucleotides.pl \
	-infile $out_dir/raw_variants.rmhex.rmsk.rmintron.txt \
	-outfile $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.txt \
	-refgenome $hg19_reference

# filter variants that were caused by mismapped reads
# this may take a while depending on the number of variants to screen and the size of the reference genome
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format
	$SNPiR/BLAT_candidates.pl \
	-infile $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.txt \
	-outfile $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt \
	-bamfile $out_dir/merged.conv.sort.rd.realigned.bam \
	-refgenome $hg19_reference

# remove known RNA editing sites
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt | \
	intersectBed -a stdin -b $rnaedit -v  > \
	$out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed
	
# extract variants from the raw vcf
skipn=$(cat $out_dir/raw_variants.vcf | grep '#' | wc -l)
cat $out_dir/raw_variants.vcf | grep '#' > $out_dir/header.vcf 
Rscript $extractvcf $out_dir $skipn
cat $out_dir/header.vcf $out_dir/red_variants.vcf > $out_dir/final_variants.vcf

################################################################################
# clean up stuff..
################################################################################
#shopt -s extglob
#cd $out_dir
#rm !(final_variants.vcf)
#rm $out_dir/*.bam

exit 0
