library(plyr)
library(mygene)

args <- commandArgs(TRUE)
prefix = args[1] ###/data/storage/patients
sample = args[2]

.collapse <- function (...) {
     paste(unlist(list(...)), sep = ",", collapse = ",")
 }

act <- read.csv("/data/database/druggability/cancer_genes.txt", sep="\t", stringsAsFactors = F)
## annotate CNV output, intersect with Actionable Genes curated by Tobias
annotateCNV <- function(prefix, sample){
  cnv <- read.csv(paste0(prefix, "/cnv/cnvkit-medexome/", sample, "-T1-DNA.bed"), sep = "\t", stringsAsFactors = FALSE)
  cnv$query <- paste0("hg19.", cnv[[1]], ":", cnv[[2]], "-", cnv[[3]])
  cnv$symbol <- sapply(cnv$query, function(i) .collapse(subset(query(i, species="human")$hits, taxid=="9606")$symbol[[1]]))
  final <- arrange(subset(cnv, symbol %in% act$Symbol), -abs(seg_mean))
  write.table(final, file=paste0(prefix, "/cnv/cnvkit-medexome/", sample, "_cnv_annotations.txt"), sep="\t", row.names=FALSE, quote=FALSE)
}
annotateCNV(prefix, sample)
