library(myvariant)
library(magrittr)
library(IRanges)
library(plyr)
library(httr)
library(jsonlite)
library(myvariant)
library(VariantAnnotation)
library(data.table)
source("~/.virtualenvs/pm/omics_pipe/omics_pipe/scripts/annotateIndels.R")

cancer_genes <- read.csv("/data/database/druggability/cancer_genes.txt", stringsAsFactors = FALSE, sep="\t")

.collapse <- function (...) {
  paste(unlist(list(...)), sep = ",", collapse = ",")
}

somaticAnnotation <- function(vcf.path, do_filter=TRUE){
  ## MyVariant.info annotations
  vcf <- readVcf(vcf.path, genome="hg19")
  snp <- vcf[isSNV(vcf)]
  hgvs <- formatHgvs(snp, "snp")
  annos <- getVariants(hgvs, fields=c("cadd.gene.prot.protpos", "cadd.oaa", "cadd.naa", "dbsnp.rsid", "cadd.consequence", 
                                      #"dbnsfp.aa.pos", "dbnsfp.aa.ref", "dbnsfp.aa.alt", 
                                      "snpeff.ann.hgvs_p",
                                      # "cadd.gene.prot.protpos", "cadd.oaa", "cadd.naa",
                                      "cadd.gene.genename",
                                      "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af",
                                      "dbnsfp.1000gp3.af", "cadd.phred",
                                      "dbnsfp.polyphen2.hdiv.rankscore", "dbnsfp.polyphen2.hdiv.pred", 
                                      "dbnsfp.mutationtaster.converted_rankscore", "dbnsfp.mutationtaster.pred"
  ))
  annos$Position <- paste0(seqnames(snp), ":", start(snp))
  # Coverage by Depth
  dp <- data.frame(geno(snp)$DP)
  row.names(dp) <- NULL
  # Allelic Depth
  ad <- data.frame(geno(snp)$AD)
  row.names(ad) <- NULL
  if (grepl("varscan", vcf.path)){
    setnames(dp, old=c(names(dp)[grepl("TUMOR", names(dp))], names(dp)[grepl("NORMAL", names(dp))]), 
                 new=c("TUMOR.DP", "NORMAL.DP"))
    setnames(ad, old=c(names(ad)[grepl("TUMOR", names(ad))], names(ad)[grepl("NORMAL", names(ad))]), 
                 new=c("TUMOR.AD", "NORMAL.AD"))
    coverage <- cbind(ad, dp)
    coverage$NORMAL.ALT.AF <- round(coverage$NORMAL.AD/coverage$NORMAL.DP, 2)
    coverage$TUMOR.ALT.AF <- round(coverage$TUMOR.AD/coverage$TUMOR.DP, 2)
    coverage$VariantCaller <- "Varscan2 Somatic"
   }
  if (grepl("mutect", vcf.path)){
    setnames(dp, old=c(names(dp)[grepl(".T", names(dp))], names(dp)[grepl(".B", names(dp))]), 
                 new=c("TUMOR.DP", "NORMAL.DP"))
    setnames(ad, old=c(names(ad)[grepl(".T", names(ad))], names(ad)[grepl(".B", names(ad))]), 
                 new=c("TUMOR.AD", "NORMAL.AD"))
    coverage <- cbind(ad, dp)
    naf <- lapply(coverage$NORMAL.AD, function(i) as.vector(i[[2]]/(i[[1]] + i[[2]])))
    names(naf) <- NULL
    taf <- lapply(coverage$TUMOR.AD, function(i) as.vector(i[[2]]/(i[[1]] + i[[2]])))
    names(taf) <- NULL
    coverage$NORMAL.ALT.AF <- lapply(naf, round, 2)
    coverage$TUMOR.ALT.AF <- lapply(taf, round, 2)
    coverage$VariantCaller <- "MuTect"
   }
  annos <- cbind(annos, coverage)
  ## filter consequence
  if(do_filter){
    annos <- subset(annos, cadd.consequence %in% c("STOP_GAINED","STOP_LOST", 
                                                   "NON_SYNONYMOUS", "SPLICE_SITE", 
                                                   "CANONICAL_SPLICE", "REGULATORY"))
    annos <- arrange(data.frame(subset(annos, is.na(exac.af) | exac.af < 0.05)), -dbnsfp.mutationtaster.converted_rankscore)
  }
  annos <- data.frame(annos)
  setnames(annos, 
           old = c("query", "cadd.naa", "cadd.oaa", "dbsnp.rsid", "cadd.consequence", "snpeff.ann",
                   "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af", "cadd.phred",
                   "dbnsfp.polyphen2.hdiv.rankscore", "dbnsfp.polyphen2.hdiv.pred",
                   "dbnsfp.mutationtaster.converted_rankscore", "dbnsfp.mutationtaster.pred"), 
           new = c("Variant", "Ref.AA", "Alt.AA", "dbSNP rsid", "Consequence", "Amino Acid",
                   "COSMIC ID", "COSMIC Tumor Site", "ExAC AF", "CADD Score",
                   "Polyphen-2 Score", "Polyphen-2 Prediction", "MutationTaster Score", "MutationTaster Prediction"))
  #names(annos)[names(annos) %in% c("cadd.gene.genename", "cadd.gene")] <- "Gene"
  annos <- DataFrame(annos) ##for some reason have to do this to eliminate the following columns
  annos[c("X_id", "notfound", "X_score", "cadd._license")] <- NULL
  annos$Gene <- sapply(annos$cadd.gene, function(i) i$genename)
  annos$`Amino Acid Position` <- lapply(annos$cadd.gene, function(i) .collapse(i$prot))
  annos <- lapply(annos, function(i) sapply(i, .collapse))
  ## merge annotations
  annos <- data.frame(annos)
  annotations <- merge(annos, cancer_genes, by.x="Gene", by.y="symbol")
  ## write file
  if(!nrow(annotations==0)){
  annotations <- data.frame(sapply(annotations, as.character), stringsAsFactors = FALSE)
  }
  annotations <- rbind.fill(annotations, annotateIndels(vcf.path))
  annotations[is.na(annotations)] <- ""
  annotations
}
args <- commandArgs(TRUE)
prefix <- args[1]
varscan <- args[2]
mutect <- args[3]
filter <- args[4]
som <- do.call(rbind.fill, lapply(c(varscan, mutect), function(i) somaticAnnotation(i, do_filter=filter)))
df <- som[order(som$Gene), ]

write.table(df, paste(prefix, "_somatic_annotations.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

