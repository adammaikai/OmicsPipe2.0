===========================
Reference Databases Needed
===========================

To run the pipelines, you will need to have reference databases installed on your cluster. If you are using the
AWS installation, these databases are provided for you.  If you need to install your references, please install the ones below.
Omics Pipe is compatible with all species genome files.  Examples below are for hg19, but you can substitute them for the equivalent files from other species. 

All Pipelines
===========

Genome
------------------

* .fa file can be downloaded from: `<http://cufflinks.cbcb.umd.edu/igenomes.html>`_

Reference Annotation Files
------------------------------
You can use any reference annotations you would like, as long as they are GTF files. 

Examples include: 

* gencode.v18.annotation.gtf

* UCSC genes.gtf

Reference Data for Cancer Reporting Scripts (RNAseq cancer, TCGA pipelines)
===================================================================================================
For the cancer pipelines, please download the file from the link below, extract it and put the files in the respective directories.
`Reporting_data <http://sulab.scripps.edu/var/www/kfisch/omics_pipe_ref.tar.gz>`_

In your omics pipe installation directory under omics_pipe/scripts/reporting/ref place the files.
------------------------------------------------------------------------------------------------------

* `tcga_brca.R <http://sulab.scripps.edu/kfisch/omics_pipe/ref/tcga_brca.R>`_

* brca.txt <http://sulab.scripps.edu/kfisch/omics_pipe/ref/brca.txt>`_

In your omics pipe installation directory under omics_pipe/scripts/reporting/data place the remaining files.
-----------------------------------------------------------------------------------------------------------------


* brca_mol_class/*

* DoG/*

* geneLists/*

* SPIA

* deseq.tcga_brca.Rdata

* loggeoameansBRCA.Rdata



References for Variants (RNA-seq cancer, RNA-seq cancer TCGA, WES and WGS pipelines)
======================================================================================
For pipelines performing variant calling, please download the references below and put them in the specified directories.

You can put these files in any directory. You will point to their location in the parameters file. 
---------------------------------------------------------------------------------------------------

* `dbNSFP2.0.txt <https://www.firedrive.com/file/19D60333C6A3D3B8>`_

* `common_no_known_medical_impact_00-latest.vcf <ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/common_no_known_medical_impact-latest.vcf.gz>`_

Available within the `GATK recource bundle v.2.5 <http://www.broadinstitute.org/gatk/guide/article.php?id=1215>`_:

* dbsnp_137.hg19.vcf

* Mills_and_1000G_gold_standard.indels.hg19.vcf

* 1000G_phase1.indels.hg19.vcf

* hapmap_3.3.hg19.vcf

* 1000G_omni2.5.hg19.vcf


In your omics pipe installation directory under omics_pipe/scripts/reporting/ref place these files.
----------------------------------------------------------------------------------------------------
* cadd.tsv.gz from `http://cadd.gs.washington.edu/download <http://cadd.gs.washington.edu/download>`_

* `drugbank.tsv <http://sulab.scripps.edu/kfisch/omics_pipe/ref/drugbank.tsv>`_

* `cosmic.tsv <http://sulab.scripps.edu/kfisch/omics_pipe/ref/cosmic.txt>`_

* `clinvar.txt <http://sulab.scripps.edu/kfisch/omics_pipe/ref/clinvar.txt>`_

from `PharmGKB <https://www.pharmgkb.org/downloads/>`_:

* pharmgkbAllele.tsv

* pharmgkbRSID.csv


WES Pipeline 
===============

* `truseq_exome_targeted_regions.hg19.bed <http://supportres.illumina.com/documents/myillumina/5dfd7e70-c4a5-405a-8131-33f683414fb7/truseq_exome_targeted_regions.hg19.bed.chr.gz>`_


ChIP-seq Pipelines
========================
* `hg19.chrom.sizes <https://genome.ucsc.edu/goldenPath/help/hg19.chrom.sizes>`_


SNPiR Pipelines (RNA-seq cancer and RNA-seq cancer TCGA pipelines)
========================================================================

* BWA Index

* RNA editing sites (Human_AG_all_hg19.bed)

* RepeatMasker.bed

* anno_combined_sorted

* knowngene.bed