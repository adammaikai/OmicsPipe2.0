#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make alignment output directory if it doesn't exist
mkdir -p $4

#Move tmp dir to scratch 
export TMPDIR=$6  #TEMP_DIR

#Loads specified module versions
module load python/$8
module load macs/$5
module load bedtools/$7


###Runs MACS for ChIPseq peak calling  $1PAIR_LIST $2 BOWTIE_RESULTS $3CHROM_SIZES $4 MACS_RESULTS $5MACS_VERSION $6 TEMP_DIR $7 BEDTOOLS_VERSION $8 PYTHON_VERSION

## Peak calling using MACS (http://liulab.dfci.harvard.edu/MACS)
# chip-input above is filename of the chip bam file plus filename of the input bam file, connected by a dash
#PAIR LIST NEEDS TO BE IN THIS FORMAT 'pair1_chip-pair1_input pair2_chip-pair2_input'

cd $4

array=($1)  

for pair in ${array[*]}; do
    echo -en $pair"\t"
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
    # Run MACS
    GEN_SIZE=$(awk '{size+=$2}END{print size}' $3)
    PVALUE=1e-5
    macs14 -t $2/${chip}/${chip}.bam -c $2/${input}/${input}.bam --name=${pair}_macs_p05 --format=BAM --gsize=$GEN_SIZE --pvalue=$PVALUE 2> $4/${pair}_macs_p05.log
    # Print genomic fragment length
    grep "# d = " $4/${pair}_macs_p05_peaks.xls | awk '{print $4}'
    # Check warnings
    grep "WARNING" $4/${pair}_macs_p05.log
    # Remove intermediate files
    rm $4/${pair}_macs_p05{_model.r,_negative_peaks.xls,_peaks.bed}
done
# Re-run MACS with a lower mfold parameter if WARNING "Fewer paired peaks (X) than 1000!"
# Number of peaks at different FDR thresholds
(echo -e "FDR\tAll\t5\t1\t0"
for pair in ${array[*]}; do
    echo -en $pair
    for fdr in 100 5 1 0; do
        echo -en "\t"$(grep -v "#" $4/${pair}_macs_p05_peaks.xls | awk -vFDR=$fdr '(NR>1&&$9<=FDR)' | wc -l)
    done
    echo
done)


# Define confident peaks (FDR), enriched regions (p-value<=10e-5) and control peaks
FDR=1
for pair in ${array[*]};do 
    # Confident peaks
    grep -v "#" $4/${pair}_macs_p05_peaks.xls | awk -vOFS='\t' -vFDR=$FDR '(NR>1&&$9<=FDR){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > $4/${pair}_macs_confident.txt
    # Regions with significant enrichment
    grep -v "#" $4/${pair}_macs_p05_peaks.xls | sed '/^$/d' | awk -vOFS='\t' '(NR>1){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > $4/${pair}_macs_enrichment.txt    #remove blank line instead of # grep returns a blank line before the header; (NR>1) alone would remove the blank line and return the header; (NR>1&&$9<=100) is a workaround. 
    # Control peaks
    shuffleBed -i $4/${pair}_macs_enrichment.txt -g $3 -chrom | sort -k1,1 -k2,2n > $4/${pair}_macs_control.txt
done

## Step 9: Peak visualization
for pair in ${array[*]}; do
    # Create BED files
    (echo -e "track name=\"${pair}_confident_peaks\" description=\"${pair}_confident_peaks\" visibility=2"
    sort -k5,5gr $4/${pair}_macs_confident.txt | awk -vOFS='\t' '{print $1,$2,$3,"PEAK_"NR,$5,"."}' | sort -k1,1 -k2,2n) | gzip > $4/${pair}_macs_confident.bed.gz
    (echo -e "track name=\"${pair}_enriched_regions\" description=\"${pair}_enriched_regions\" visibility=2"
    sort -k5,5gr $4/${pair}_macs_enrichment.txt | awk -vOFS='\t' '{print $1,$2,$3,"PEAK_"NR,$5,"."}' | sort -k1,1 -k2,2n) | gzip > $4/${pair}_macs_enrichment.bed.gz
done