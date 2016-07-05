=======================
Third Party Software Dependencies
=======================
Omics Pipe is dependent upon several third-party software packages.  Before running Omics Pipe, please install all of the
required tools for the pipeline you will be running (see below) as `Modules`_ on your local cluster.  If you are running the AWS
distribution, all third party software is already installed.

.. _Modules: http://modules.sourceforge.net/

.. rst-class:: html-toggle	
R Packages Needed
===================

In R, you can cut and copy this to install all required packages::

	install.packages(c("bibtex",	"AnnotationDbi", "cluster", "cummeRbund", "data.table", "DBI", "DESeq2", "devtools",	"dplyr",	"gdata",
	"ggplot2", "graphite", "igraph", "KEGGREST","knitr", "knitrBootstrap",	"lattice", "locfit",	"pamr",	"pander",	"pathview",
	"plyr","RColorBrewer","Rcpp",	"RcppArmadillo", "RCurl", "ReactomePA", "RefManageR","RJSONIO","RSQLite",
	"stringr","survival", "XML", "xtable",	"yaml"))

.. rst-class:: html-toggle	
RNA-seq (Tuxedo)
===================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	TOPHAT:			`<http://tophat.cbcb.umd.edu/>`_
	CUFFLINKS:			`<http://cufflinks.cbcb.umd.edu/>`_
	====================== ===================================================
	
.. rst-class:: html-toggle
RNA-seq (Anders 2013)
=======================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	STAR:			`<http://code.google.com/p/rna-star/>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	HTSEQ:			`<http://www-huber.embl.de/users/anders/HTSeq/doc/index.html>`_
	====================== ===================================================
	
.. rst-class:: html-toggle
Whole Exome Sequencing (GATK)
============================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	BWA:			`<http://bio-bwa.sourceforge.net/>`_
	PICARD:			`<http://picard.sourceforge.net/>`_
	GATK:			`<https://www.broadinstitute.org/gatk/download>`_	
	====================== ===================================================
			

.. rst-class:: html-toggle
Whole Genome Sequencing (GATK)
=================================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	BWA:			`<http://bio-bwa.sourceforge.net/>`_
	PICARD:			`<http://picard.sourceforge.net/>`_
	GATK:			`<https://www.broadinstitute.org/gatk/download>`_	
	====================== ===================================================
	

.. rst-class:: html-toggle
Whole Genome Sequencing (MUTECT)
=================================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	BWA:			`<http://bio-bwa.sourceforge.net/>`_
	MUTECT:			http://www.broadinstitute.org/cancer/cga/mutect	
	====================== ===================================================
		

.. rst-class:: html-toggle
ChIP-seq (MACS)
======================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	HOMER:			`<http://homer.salk.edu/homer/ngs/index.html>`_
	BOWTIE:			`<http://bowtie-bio.sourceforge.net/index.shtml>`_
	MACS:			`<http://liulab.dfci.harvard.edu/MACS/>`_
	BEDTOOLS:			`<https://github.com/arq5x/bedtools2>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	====================== ===================================================
	

.. rst-class:: html-toggle	
ChIP-seq (HOMER)
====================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	HOMER:			`<http://homer.salk.edu/homer/ngs/index.html>`_
	BOWTIE:			`<http://bowtie-bio.sourceforge.net/index.shtml>`_
	BEDTOOLS:			`<https://github.com/arq5x/bedtools2>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	====================== ===================================================
	

.. rst-class:: html-toggle
Breast Cancer Personalized Genomics Report- RNAseq
================================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	STAR:			`<http://code.google.com/p/rna-star/>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	HTSEQ:			`<http://www-huber.embl.de/users/anders/HTSeq/doc/index.html>`_
	RSEQC:			`<http://rseqc.sourceforge.net/>`_
	PICARD:			`<http://picard.sourceforge.net/>`_
	GATK:			`<https://www.broadinstitute.org/gatk/download>`_	
	FusionCatcher:			`<https://code.google.com/p/fusioncatcher/>`_
	Oncofuse:			`<http://www.unav.es/genetica/oncofuse.html>`_
	BWA:			`<http://bio-bwa.sourceforge.net/>`_
	DNANEXUS SAMTOOLS:			`<https://github.com/dnanexus/samtools>`_
	BEDTOOLS:			`<https://github.com/arq5x/bedtools2>`_
	BLAT:			`<https://genome.ucsc.edu/FAQ/FAQblat.html#blat3>`_
	SNPiR:			`<http://lilab.stanford.edu/SNPiR/>`_
	SNPEFF:			`<http://snpeff.sourceforge.net/>`_
	SNPSIFT:			`<http://snpeff.sourceforge.net/SnpSift.html>`_
	VCFTOOLS:			`<http://vcftools.sourceforge.net/>`_	
	====================== ===================================================
	

.. rst-class:: html-toggle
TCGA Reanalysis Pipeline - RNAseq
=========================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	GeneTorrent:			`<https://cghub.ucsc.edu/docs/user/software.html>`_
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	STAR:			`<http://code.google.com/p/rna-star/>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	HTSEQ:			`<http://www-huber.embl.de/users/anders/HTSeq/doc/index.html>`_
	RSEQC:			`<http://rseqc.sourceforge.net/>`_
	PICARD:			`<http://picard.sourceforge.net/>`_
	GATK:			`<https://www.broadinstitute.org/gatk/download>`_	
	FusionCatcher:			`<https://code.google.com/p/fusioncatcher/>`_
	Oncofuse:			`<http://www.unav.es/genetica/oncofuse.html>`_
	BWA:			`<http://bio-bwa.sourceforge.net/>`_
	DNANEXUS SAMTOOLS:			`<https://github.com/dnanexus/samtools>`_
	BEDTOOLS:			`<https://github.com/arq5x/bedtools2>`_
	BLAT:			`<https://genome.ucsc.edu/FAQ/FAQblat.html#blat3>`_
	SNPiR:			`<http://lilab.stanford.edu/SNPiR/>`_
	SNPEFF:			`<http://snpeff.sourceforge.net/>`_
	SNPSIFT:			`<http://snpeff.sourceforge.net/SnpSift.html>`_
	VCFTOOLS:			`<http://vcftools.sourceforge.net/>`_	
	====================== ===================================================

		
	
.. rst-class:: html-toggle
TCGA Reanalysis Pipeline - RNAseq Counts
=========================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	GeneTorrent:			`<https://cghub.ucsc.edu/docs/user/software.html>`_
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	STAR:			`<http://code.google.com/p/rna-star/>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	HTSEQ:			`<http://www-huber.embl.de/users/anders/HTSeq/doc/index.html>`_
	====================== ===================================================
	


.. rst-class:: html-toggle
miRNAseq Counts (Anders 2013)
=========================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	CutAdapt:			`<http://code.google.com/p/cutadapt/>`_
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	STAR:			`<http://code.google.com/p/rna-star/>`_
	SAMTOOLS:			`<http://samtools.sourceforge.net/>`_
	HTSEQ:			`<http://www-huber.embl.de/users/anders/HTSeq/doc/index.html>`_
	====================== ===================================================
		

.. rst-class:: html-toggle
miRNAseq (Tuxedo)
=========================

	.. rst-class:: html-plain-table
	
	====================== ===================================================
	CutAdapt:			`<http://code.google.com/p/cutadapt/>`_
	FASTQC:			`<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
	TOPHAT:			`<http://tophat.cbcb.umd.edu/>`_
	CUFFLINKS:			`<http://cufflinks.cbcb.umd.edu/>`_
	====================== ===================================================
	
	