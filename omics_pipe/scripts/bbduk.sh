#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Loads specified python addons module for cutadapt
module load bbmap/$3

#Make assemblies output director if it doesn't exist
mkdir -p $4/$1

#Runs bbduk with $1=PATIENT, $2=RAW_DATA_DIR, $3=BBDUK_VERSION, $4=RESULTS, $5=REF_ADAPTOR, $6=BBDUK_PARAMETERS, $7=MEMORY, $8=CPU

mem=$(expr $(echo $7 | sed 's/gb//'))

bbduk.sh \
	in=$2/$1_1.fastq.gz \
	in2=$2/$1_2.fastq.gz \
	out=$4/$1/$1_1.fastq.gz \
	out2=$4/$1/$1_2.fastq.gz \
	ref=$5 \
	stats=$4/$1/stats/$1.stats \
	bhist=$4/$1/stats/$1.bhist \
	qhist=$4/$1/stats/$1.qhist \
	aqhist=$4/$1/stats/$1.aqhist \
	lhist=$4/$1/stats/$1.lhist \
	gchist=$4/$1/stats/$1.gchist \
	gcbins=auto \
	threads=${8} \
	$6 
#	-Xmx${mem}g


exit 0
