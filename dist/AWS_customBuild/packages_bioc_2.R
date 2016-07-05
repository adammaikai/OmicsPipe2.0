source("~/.Rprofile")

source("http://bioconductor.org/biocLite.R")
biocLite(ask=FALSE)
biocLite(c("IRanges",
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