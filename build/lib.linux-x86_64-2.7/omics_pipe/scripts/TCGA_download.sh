#!/bin/bash
set -x


#Source modules for current shell
source $MODULESHOME/init/bash

module load ucsc_tools/$4

#Make alignment output directory if it doesn't exist
mkdir -p $3

#Runs gtdownload with $0 TCGA_XML_FILE $1 = sample, $2 = TCGA_KEY, $3 = TCGA_OUTPUT_PATH, $4 =  CGATOOLS_VERSION
#Create XML file from https://browser.cghub.ucsc.edu/ 

sample=$1

echo "${sample}"

temp="${sample#\TCGA_}"

echo "${temp}"

gtdownload -d ${temp} -c $2 -p $3 

tar -zxvf $3/${temp}/*.tar.gz -C $3/${temp}/

mv $3/${temp}/*_1.fastq $3/${sample}_1.fastq
mv $3/${temp}/*_2.fastq $3/${sample}_2.fastq

rm $3/${temp}/*.tar.gz

exit 0

