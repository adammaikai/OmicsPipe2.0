#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

# set date
date1=$(date +%s.%N)

####INPUTS: 
#$1: SNIPR_RESULTS 
resDir=${1}
#$2: TEMP_DIR 
tempDir=${2}
#$3 SAMTOOLS_VERSION 
samtoolsVersion=${3}
#$4 BWA_VERSION
bwaVersion=${4}
#$5 PICARD_VERSION 
picardVersion=${5}
#$6 GATK_VERSION 
gatkVersion=${6}
#$7 BEDTOOLS_VERSION
bedtoolsVersion=${7} 
#$8 UCSC_TOOLS_VERSION
ucscVersion=${8}
#$9 GENOME 
hg19_reference=$9
#$10 REPEAT_MASKER
RepeatMasker=${10} 
#$11 REF_GENES 
gene_annotation=${11} 
#$12 RNA EDIT 
rnaedit=${12}
#$13 DBSNP 
dbsnp=${13}
#$14 MILLS 
mills=${14}
#$15 G1000 
g1000=${15}
#$16 WORKING_DIR 
workingDir=${16}
#$17 SAMPLE 
sample=${17}
#$18 RAW_DATA
rawDir=${18}
#$19 SNIPR_VERSION 
snpirVersion=${19}
#$20 SNIPR_CONFIG
snpirConfig=${20}
#$21 SNIPR_DIR
SNPiR=${21}
#$22 ENCODING #phred33 or #phred64
encoding=${22} 
# callable loci 
locibed=${23} # bed file for CallableLoci tool
callableloci=${24} # run CallableLoci, yes no
# read group info
rgid=$sample
rglb=${25}
rgpl=${26}
rgsm=$sample
# ncpu
ncpu=${27}
nthreads=${27}
# mem
mem=${28}
# samblaster
samblasterVersion=${29}
# sambamba
sambambaVersion=${30}
# R
rVersion=${31}
# BWA INDEX
bwa_index=${32}
#alignment directory
alignment_dir=${33}
#Make assemblies output director if it doesn't exist
out_dir=$resDir/$sample
mkdir -p $out_dir
cpu=$(grep -c "processor" /proc/cpuinfo)

#Move tmp dir to scratch 
export TMPDIR=$tempDir  #TEMP_DIR

# SET PERL lib
export PERL5LIB=$snpirConfig
 
#Load specified module versions
module load samtools/$samtoolsVersion
#module load samtools/dnanexus-1.0
module load bwa/$bwaVersion
module load picard/$picardVersion
module load gatk/$gatkVersion
module load bedtools/$bedtoolsVersion
module load ucsc-tools/$ucscVersion
module load snpir/$snpirVersion
module load samblaster/$samblasterVersion
module load sambamba/$sambambaVersion
module load R/$rVersion

extractvcf=$workingDir/extractvcf.R

cp $SNPiR/convertCoordinates.* $out_dir

# input read1 & read2
in1=$rawDir/${sample}/${sample}_1.fastq.gz
in2=$rawDir/${sample}/${sample}_2.fastq.gz

# read group information for read1 & read2, is required for GATK
RGR1="@RG\tID:$rgid\tLB:$rglb\tPL:$rgpl\tSM:$rgsm"
RGR2="@RG\tID:$rgid\tLB:$rglb\tPL:$rgpl\tSM:$rgsm"


################################################################################
# MAPPING
################################################################################
# align with bwa as single reads
# bwa mem -M -t $ncpu -R $RGR1 $bwa_index $in1 > $out_dir/out1.sam
# bwa mem -M -t $ncpu -R $RGR2 $bwa_index $in2 > $out_dir/out2.sam

# # merge alignments
# cat $out_dir/out1.sam <(grep -v '^@' $out_dir/out2.sam) > $out_dir/merged.sam
# rm $out_dir/out1.sam $out_dir/out2.sam # clean up

# # convert the position of reads that map across splicing junctions onto the genome
# cd $out_dir
# java -Xmx4g -cp . convertCoordinates < $out_dir/merged.sam > $out_dir/merged.conv.sam
# cd ~

# rm $out_dir/merged.sam $out_dir/convertCoordinates.* # clean up

# # rewrite the .sam header to drop the pseudochromosomes
# samtools view -HS $out_dir/merged.conv.sam > $out_dir/header.sam
# sed -n '1,25 p' $out_dir/header.sam > $out_dir/t1.sam
# tail -n -2 $out_dir/header.sam > $out_dir/t2.sam
# cat $out_dir/t1.sam $out_dir/t2.sam > $out_dir/new_header.sam

# java -jar `which ReplaceSamHeader.jar` \
# 	INPUT=$out_dir/merged.conv.sam \
# 	HEADER=$out_dir/new_header.sam \
# 	OUTPUT=$out_dir/merged.conv.nh.sam 
	
# rm $out_dir/merged.conv.sam $out_dir/new_header.sam $out_dir/header.sam $out_dir/t1.sam $out_dir/t2.sam # clean up

# # sort file, filter out unmappes reads and reads with mapping quality < 20 and convert to .bam
# # samtools view -bS -q 20 -F 4 $out_dir/merged.conv.nh.sam | samtools rocksort -@ $ncpu -m 1500M - $out_dir/merged.conv.nh.sort
# #cat $out_dir/merged.conv.nh.sam | samblaster -r | samtools view -bS -q 20 -F 4 - | sambamba sort -m $mem -t $ncpu -p -o $out_dir/merged.conv.nh.sort.rd
# samtools view -bS -q 20 -F 4 $out_dir/merged.conv.nh.sam | samtools sort -@ 8 -m 2G - $out_dir/merged.conv.nh.sort

# remove duplicate reads ###THIS PART FINISHES SUCCESSFULLY
# java -Xmx4g -jar `which MarkDuplicates.jar` \
# 	INPUT=$out_dir/merged.conv.nh.sort.bam \
# 	REMOVE_DUPLICATES=false \
# 	VALIDATION_STRINGENCY=LENIENT \
# 	AS=true \
# 	METRICS_FILE=$out_dir/$sample.dups \
# 	OUTPUT=$out_dir/merged.conv.nh.sort.rd.bam \
# 	TMP_DIR=$TMPDIR

# rm $out_dir/merged.conv.nh.sort.bam $out_dir/merged.conv.nh.sam # clean up

# # index ###THIS PART WORKS
# samtools index $out_dir/merged.conv.nh.sort.rd.bam

samtools index $alignment_dir/$sample/$sample\_Aligned.sortedByCoord.out.bam
sambamba markdup \
-t $cpu \
--overflow-list-size 1000000 \
--hash-table-size 1000000 \
$alignment_dir/$sample/$sample\_Aligned.sortedByCoord.out.bam \
$alignment_dir/$sample/$sample\_dedup.bam

samtools index $alignment_dir/$sample/$sample\_dedup.bam

#SplitNCigarReads  (from Broad RNASeq Variant best practice)
java -jar `which GenomeAnalysisTK.jar` \
	-T SplitNCigarReads \
	-R $hg19_reference \
	-I $alignment_dir/$sample/$sample\_dedup.bam \
	-o $out_dir/merged.conv.nh.sort.rd.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS \
    -fixNDN

# indel realignment & base quality score recalibration ###THIS PART FINISHES SUCCESSFULLY
if [ $encoding == "phred64" ]; then
	java -Xmx50g -jar `which GenomeAnalysisTK.jar` \
		-T RealignerTargetCreator \
		-R $hg19_reference \
		-I $out_dir/merged.conv.nh.sort.rd.split.bam \
		-o $out_dir/output.intervals \
		-known $mills \
		-known $g1000 \
		-nt $cpu \
		--fix_misencoded_quality_scores #\
		#--filter_reads_with_N_cigar
	
	java -Xmx50g -Djava.io.tmpdir=$TMPDIR -jar `which GenomeAnalysisTK.jar` \
		-I $out_dir/merged.conv.nh.sort.rd.split.bam \
		-R $hg19_reference \
		-T IndelRealigner \
		-targetIntervals $out_dir/output.intervals \
		-o $out_dir/merged.conv.sort.rd.split.realigned.bam \
		-known $mills \
		-known $g1000 \
		--consensusDeterminationModel KNOWNS_ONLY \
		-LOD 0.4 \
		--fix_misencoded_quality_scores #\
		#--filter_reads_with_N_cigar
else
	java -Xmx50g -jar `which GenomeAnalysisTK.jar` \
		-T RealignerTargetCreator \
		-R $hg19_reference \
		-I $out_dir/merged.conv.nh.sort.rd.split.bam \
		-o $out_dir/output.intervals \
		-known $mills \
		-known $g1000 \
		-nt $cpu #\
		#--filter_reads_with_N_cigar
	
	java -Xmx50g -Djava.io.tmpdir=$TMPDIR -jar `which GenomeAnalysisTK.jar` \
		-I $out_dir/merged.conv.nh.sort.rd.split.bam \
		-R $hg19_reference \
		-T IndelRealigner \
		-targetIntervals $out_dir/output.intervals \
		-o $out_dir/merged.conv.sort.rd.split.realigned.bam \
		-known $mills \
		-known $g1000 \
		--consensusDeterminationModel KNOWNS_ONLY \
		-LOD 0.4 #\
		#--filter_reads_with_N_cigar
fi

###THIS PART FINISHES SUCCESSFULLY	
# Base Quality Score Recalibration
	#1 Analyze patterns of covariation in the sequence dataset
java -Xmx50g -jar `which GenomeAnalysisTK.jar` \
    -T BaseRecalibrator \
    -R $hg19_reference \
    -I $out_dir/merged.conv.sort.rd.split.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $g1000 \
    -o $out_dir/recal_data.table \
    -nct $cpu
    
	#2 Do a second pass to analyze covariation remaining after recalibration
###THIS PART FINISHES SUCCESSFULLY
java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
    -T BaseRecalibrator \
    -R $hg19_reference \
    -I $out_dir/merged.conv.sort.rd.split.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $g1000 \
    -BQSR $out_dir/recal_data.table \
    -o $out_dir/post_recal_data.table \
    -nct $cpu

	#3 Generate before / after plots This part throws error
#java -jar `which GenomeAnalysisTK.jar` \
#    -T AnalyzeCovariates \
#    -R $hg19_reference \
#    -before $out_dir/recal_data.table \
#    -after $out_dir/post_recal_data.table \
#    -plots $out_dir/recalibration_plots.pdf
    
	#4 Apply the recalibration to your sequence data
java -jar `which GenomeAnalysisTK.jar` \
    -T PrintReads \
    -R $hg19_reference \
    -I $out_dir/merged.conv.sort.rd.split.realigned.bam \
    -BQSR $out_dir/recal_data.table \
    -o $out_dir/merged.conv.sort.rd.split.realigned.recal.bam

################################################################################
# CallableLoci
###############################################################################
if  [ $callableloci=="yes" ]; then
	java -jar `which GenomeAnalysisTK.jar` \
    		-T CallableLoci \
		-R $hg19_reference \
		-I $out_dir/merged.conv.sort.rd.split.realigned.recal.bam \
		-L $locibed \
		-summary $out_dir/${sample}_callable_loci.txt \
		-o $out_dir/${sample}_callable_loci.bed \
		--filter_reads_with_N_cigar
fi

################################################################################
# Call/filter variants
################################################################################     
# call variants using UnifiedGenotyper
# as suggested by Robert Piskol by mail
# modified using Broad best practices
java -Xmx28g -jar `which GenomeAnalysisTK.jar` -T HaplotypeCaller \
	-R $hg19_reference \
	-I $out_dir/merged.conv.sort.rd.split.realigned.recal.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20 \
	-stand_emit_conf 20 \
	--dbsnp $dbsnp \
	-out_mode EMIT_VARIANTS_ONLY \
	-rf BadCigar \
	-o $out_dir/raw_variants.vcf \
	-nct $cpu

# do the filtering
# GATK based variant filtering
java -jar `which GenomeAnalysisTK.jar` \
	-T VariantFiltration \
	-R $hg19_reference \
	-V $out_dir/raw_variants.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS \
	-filter "FS > 30.0" \
	-filterName QD \
	-filter "QD < 2.0" \
	-o $out_dir/raw_variants_gatk_filt.vcf 


# # convert vcf format into custom SNPiR format and filter variants with quality <20
# $SNPiR/convertVCF.sh $out_dir/raw_variants_gatk_filt.vcf $out_dir/raw_variants.txt 20

# # filter mismatches at read ends
# # note: add the -illumina option if your reads are in Illumina 1.3+ quality format
# 	$SNPiR/filter_mismatch_first6bp.pl \
# 	-infile $out_dir/raw_variants.txt \
# 	-outfile $out_dir/raw_variants.rmhex.txt \
# 	-bamfile $out_dir/merged.conv.sort.rd.split.realigned.bam

# # filter variants in repetitive regions
# awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $out_dir/raw_variants.rmhex.txt | \
# 	intersectBed -a stdin -b $RepeatMasker -v | \
# 	cut -f1,3-7 > $out_dir/raw_variants.rmhex.rmsk.txt

# # filter intronic sites that are within 4bp of splicing junctions
# # make sure your gene annotation file is in UCSC text format and sorted by chromosome and 
# # transcript start position
# $SNPiR/filter_intron_near_splicejuncts.pl \
# 	-infile $out_dir/raw_variants.rmhex.rmsk.txt \
# 	-outfile $out_dir/raw_variants.rmhex.rmsk.rmintron.txt \
# 	-genefile $gene_annotation

# # filter variants in homopolymers
# $SNPiR/filter_homopolymer_nucleotides.pl \
# 	-infile $out_dir/raw_variants.rmhex.rmsk.rmintron.txt \
# 	-outfile $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.txt \
# 	-refgenome $hg19_reference

# # filter variants that were caused by mismapped reads
# # this may take a while depending on the number of variants to screen and the size of the reference genome
# # note: add the -illumina option if your reads are in Illumina 1.3+ quality format
# 	$SNPiR/BLAT_candidates.pl \
# 	-infile $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.txt \
# 	-outfile $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt \
# 	-bamfile $out_dir/merged.conv.sort.rd.split.realigned.bam \
# 	-refgenome $hg19_reference

# # remove known RNA editing sites
# awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.txt | \
# 	intersectBed -a stdin -b $rnaedit -v  > \
# 	$out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed
	
# # extract variants from the raw vcf
# skipn=$(cat $out_dir/raw_variants.vcf | grep '#' | wc -l)
# cat $out_dir/raw_variants.vcf | grep '#' > $out_dir/header.vcf 
# Rscript $extractvcf $out_dir $skipn
# cat $out_dir/header.vcf $out_dir/red_variants.vcf > $out_dir/final_variants.vcf

# # sort vcf file
# java -jar `which SortVcf.jar` \
# 	I=$out_dir/final_variants.vcf \
# 	O=$out_dir/final_variants_sort.vcf

# # zip and index
# bgzip -f $out_dir/final_variants_sort.vcf 
# mv $out_dir/final_variants_sort.vcf.gz $out_dir/${sample}_final_variants_sort.vcf.gz
# tabix -f -p vcf $out_dir/${sample}_final_variants_sort.vcf.gz

# ################################################################################
# # clean up stuff..
# ################################################################################
# rm $out_dir/*.bam
# rm $out_dir/*.bai
# rm $out_dir/*.vcf
# rm $out_dir/*txt*
# rm $out_dir/*.idx
# rm $out_dir/output.intervals $out_dir/post_recal_data.table $out_dir/recal_data.table
# rm $out_dir/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed
# rm $out_dir/convertCoordinates*

# # annotate variants
# Rscript ~/.virtualenvs/pm/omics_pipe/omics_pipe/scripts/annotateVariantsFromVcf.R $out_dir/${sample}_final_variants_sort.vcf.gz FALSE

# record runtime
date2=$(date +%s.%N)
dt=$(echo "$date2 - $date1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds

exit 0
