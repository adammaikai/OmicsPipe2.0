#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Make assemblies output director if it doesn't exist
mkdir -p $3/$8

#Move tmp dir to scratch 
export TMPDIR=$7  #TEMP_DIR
 
#Load specific module versions
module load intogen/$4


#$1 p.VCF_FILE, 2 p.INTOGEN_OPTIONS, 3p.INTOGEN_RESULTS, 4p.INTOGEN_VERSION, 5p.USERNAME, 6p.WORKING_DIR, 7p.TEMP_DIR, 8 sample $9 SCHEDULER $10 VARIANT_RESULTS

sched=$9
#create config file
if [ "$sched" = "SGE" ]; then
	$6/intogen_config_AWS.sh $3/$8 > $3/$8/intogen.conf
else
	$6/intogen_config.sh $3/$8 > $3/$8/intogen.conf
fi

vcf=$1
if ["$vcf" = ""]; then
	run analysis -p $8 $2 -C $3/$8/intogen.conf ${10}/$8/$8.vcf
else
	run analysis -p $8 $2 -C $3/$8/intogen.conf ${10}/$8/$1
fi

unzip $3/$8/results/$5/default/projects/$8/results.zip -d $7

cp $7/$8/consequences.tsv $3/$8/consequences.tsv
cp $7/$8/variant_genes.tsv $3/$8/variant_genes.tsv

exit 0