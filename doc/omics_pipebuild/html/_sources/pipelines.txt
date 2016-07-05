.. index:: omics_pipe; pipelines

============
Omics Pipe Available Pipelines
============

.. rst-class:: html-toggle
RNA-seq (Tuxedo)
===================
:doc:`RNAseq_Tuxedo_Modules`
	Modules included in the Tuxedo RNA-seq pipeline.
	
* FASTQC

* TopHat

* Cufflinks

* Cuffmerge

* Cuffmergetocompare

* Cuffdiff

* R Summary Report - CummeRbund


.. rst-class:: html-toggle
RNA-seq(Anders 2013)
=======================
:doc:`RNAseq_counts_modules`
	Modules included in the count-based RNA-seq pipeline.

* FASTQC

* STAR

* HTSEQ

* R Summary Report - DESEQ2


.. rst-class:: html-toggle
Whole Exome Sequencing (GATK)
============================
:doc:`WGS_WES_modules`
	Modules included in the whole exome sequencing pipeline. 
	
* FASTQC

* BWA-MEM

* PICARD Mark Duplicates

* GATK Preprocessing

* GATK Variant Discovery

* GATK Variant Filtering


.. rst-class:: html-toggle
Whole Genome Sequencing (GATK)
=================================
:doc:`WGS_WES_modules`
	Modules included in the whole genome sequencing pipeline. 
	
* FASTQC

* BWA-MEM

* PICARD Mark Duplicates

* GATK Preprocessing

* GATK Variant Discovery

* GATK Variant Filtering

.. rst-class:: html-toggle
Whole Genome Sequencing (MUTECT)
=================================
:doc:`mutect_modules`
	Modules included in the cancer (paired tumor/normal) whole genome sequencing pipeline. 
	
* FASTQC

* BWA-MEM

* MUTECT

.. rst-class:: html-toggle
ChIP-seq (MACS)
======================
:doc:`ChIPseq_modules_MACS`
	Modules included in the ChIP-seq MACS pipeline.

* FASTQC

* Homer ChIP Trim

* Bowtie

* MACS

.. rst-class:: html-toggle	
ChIP-seq (HOMER)
====================
:doc:`ChIPseq_modules_HOMER`
	Modules included in the ChIP-seq HOMER pipeline.

* FASTQC

* Homer ChIP Trim

* Bowtie

* Homer Read Density

* Homer Peaks

* Homer Peak Track

* Homer Annotate Peaks

* Homer Find Motifs

.. rst-class:: html-toggle
Breast Cancer Personalized Genomics Report- RNAseq
================================
:doc:`RNAseq_cancer_modules`
	Modules included in the RNAseq Cancer pipeline.
	
* FASTQC

* STAR

* RSEQC

* Fusion Catcher

* BWA/SNPiR

* Filter Variants

* HTseq

* Intogen

* `OncoRep Cancer Report`_

.. rst-class:: html-toggle
TCGA Reanalysis Pipeline - RNAseq
=========================
:doc:`RNAseq_cancer_modules_TCGA`
	Modules included in the RNAseq Cancer pipeline.

* TCGA Download (GeneTorrent)

* FASTQC

* STAR

* RSEQC

* Fusion Catcher

* BWA/SNPiR

* Filter Variants

* HTseq

* Intogen

* `OncoRep Cancer Report`_

.. rst-class:: html-toggle
TCGA Reanalysis Pipeline - RNAseq Counts
=========================
:doc:`RNAseq_counts_modules_TCGA`
	Modules included in the RNAseq counts pipeline for TCGA reanalysis.

* TCGA Download (GeneTorrent)

* FASTQC

* STAR

* HTSEQ

* Report


.. rst-class:: html-toggle
miRNAseq Counts (Anders 2013)
=========================
:doc:`miRNAseq_counts_modules`
	Modules included in the miRNAseq counts pipeline.

* Cutadapt

* FASTQ Length Filter

* FASTQC

* STAR

* HTSEQ

* Report

.. rst-class:: html-toggle
miRNAseq (Tuxedo)
=========================
:doc:`miRNAseq_Tuxedo_Modules`
	Modules included in the miRNAseq Tuxedo pipeline.

* Cutadapt

* FASTQ Length Filter

* TopHat

* Cufflinks

* Cuffmerge

* Cuffmergetocompare

* Cuffdiff

* R Summary Report

.. _OncoRep Cancer Report: https://bitbucket.org/sulab/oncorep

.. rst-class:: html-toggle
All Available Modules
=========================
:doc:`all_modules`