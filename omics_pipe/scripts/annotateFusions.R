library(data.table)
## Fusion annotation
samples <- list.files("/data/storage/patients/mojo/")


## Dienstman database
DIENST <- read.csv("/data/database/starfusion/dienstmanV14_fusions.csv", sep=",", stringsAsFactors = F)
## known oncogenic fusions
TICdb <- read.csv("/data/database/starfusion/allseqs_TICdb.txt", header=F, sep="\t", stringsAsFactors=F)
TICdb$Gene <- paste0(TICdb$V1, "--", TICdb$V2)
setnames(TICdb, old=c("V3"), new=c("PMID"))
TICdb <- TICdb[c("Gene", "PMID")]

annotateFusions <- function(sample){
  
  star_fusion <- paste0("/data/storage/patients/alignments/", sample, "/star-fusion.fusion_candidates.final")
  SF <- data.frame(read.table(star_fusion, header=T,comment.char="", check.names=F), check.names=F)
  setnames(SF, old="#fusion_name", new="Gene")
  
  mojo <- paste0("/data/storage/patients/mojo/", sample, "/", sample, "/",  sample, ".fusions")
  MOJO <- read.csv(mojo, header=F, stringsAsFactors=F, sep="\t")
  MOJO <- MOJO[c(6:17)]
  MOJO$Gene <- paste(MOJO$V6, MOJO$V12, sep="--")
  
  fusions <- merge(MOJO, SF, by="Gene")
  fusions$Actionable <- ifelse(fusions$Gene %in% TICdb$Gene, yes = "known", no = "unknown")
#  filtered <- merge(fusions, TICdb, by="Gene")
  write.table(fusions, file=paste0("/data/storage/patients/mojo/", sample, "/", sample, "/",  sample, "_fusion_annotations.txt"),sep = "\t", quote=F, row.names=F)
  fusions
  }

lapply(samples, function(i) try(annotateFusions(i), silent=T))


