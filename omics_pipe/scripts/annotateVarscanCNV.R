library(mygene)

sample = args[1]

act <- read.csv("/data/database/druggability/actionable_avera.tsv", sep="\t", stringsAsFactors = F)
## annotate Varscan CNV output, intersect with Actionable Genes curated by Tobias
annotateVarscanCNV <- function(sample){
  cnv <- read.csv(paste0("~/AWS/storage/patients/vcf/", sample, "/", sample, ".varscan_cnv.events.outfile"), sep = "\t", stringsAsFactors = FALSE)
  cnv$query <- paste0(cnv$chrom, ":", cnv$chr_start, "-", cnv$chr_stop)
  cnv$symbol <- sapply(cnv$query, function(i) query(i, species="human")$hits$symbol[[1]])
  final <- arrange(subset(cnv, symbol %in% act$Symbol), -abs(seg_mean))
  write.table(final, file=paste0("/data/storage/patients/cnv/", sample, "/", sample, "_varscan_cnv_annotations.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  }