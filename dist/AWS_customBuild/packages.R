# Install R packages needed for OmicsPipe

source("~/.Rprofile")

#install.packages(c("bibtex",
	"cluster",
	"data.table",
	"DBI",
	"devtools",
	"dplyr",
	"gdata",
	"ggplot2",
	"igraph",
	"knitcitations",
	"knitr",
	"knitrBootstrap",
	"lattice",
	"locfit",
	"pamr",
	"pander",
	"plyr",
	"RColorBrewer",
	"Rcpp",
	"RcppArmadillo",
	"RCurl",
	"RefManageR",
	"RJSONIO",
	"RSQLite",
	"stringr",
	"survival",
	"XML",
	"xtable",
	"yaml"
	)
)

require(devtools)
install_github("rCharts", "ramnathv")

#source("http://bioconductor.org/biocLite.R")
#biocLite(ask=FALSE)
#biocLite(c("annotate",
	"AnnotationDbi",
	"Biobase",
	"BiocGenerics",
	"cummeRbund",
	"DESeq",
	"DESeq2",
	"edgeR",
	"GenomicRanges",
	"graph",
	"graphite",
	"IRanges",
	"KEGGgraph",
	"KEGGREST",
	"limma",
	"org.Hs.eg.db",
	"pathview",
	"ReactomePA",
	"SPIA",
	"XVector"
	),
	ask=FALSE
) 