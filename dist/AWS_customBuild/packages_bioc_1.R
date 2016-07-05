source("~/.Rprofile")

source("http://bioconductor.org/biocLite.R")
biocLite(ask=FALSE)
biocLite(c("annotate",
	"AnnotationDbi",
	"Biobase",
	"BiocGenerics",
	"cummeRbund",
	"DESeq",
	"DESeq2",
	"edgeR",
	"GenomicRanges",
	"graph",
	"graphite"
	),
	ask=FALSE
)