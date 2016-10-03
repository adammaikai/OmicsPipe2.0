library(myvariant)
library(magrittr)
library(IRanges)
library(plyr)
library(httr)
library(jsonlite)
library(myvariant)
library(VariantAnnotation)
library(data.table)
source("~/.virtualenvs/op2/OmicsPipe2.0/omics_pipe/scripts/annotateIndels.R")

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
                                      "snpeff.ann.hgvs_p",
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
    fs <- read.csv(gsub("vcf.gz", "snp.strand.dp.csv", vcf.path), stringsAsFactors = FALSE, header = FALSE)
    names(fs) <- c("ref_fwd", "ref_rev", "var_fwd", "var_rev")
    fs$p_value <- sapply(seq(1, nrow(fs)), function(i) fisher.test(matrix(c(fs$ref_fwd[i], fs$ref_rev[i], fs$var_fwd[i], fs$var_rev[i]), nrow = 2), workspace = 2000000000)[[1]])
    annos$FS_bias_p_adjust <- p.adjust(fs$p_value, method = "BH")
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
    annos <- arrange(data.frame(subset(annos, is.null(unlist(exac.af)) | as.numeric(unlist(exac.af)) < 0.05)), -dbnsfp.mutationtaster.converted_rankscore)
    #annos <- arrange(data.frame(subset(annos, is.na(exac.af) | exac.af < 0.05)), -dbnsfp.mutationtaster.converted_rankscore)
  }
  annos <- data.frame(annos)
  oldnames <- c("query", "cadd.naa", "cadd.oaa", "dbsnp.rsid", "cadd.consequence",
                     "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af", "cadd.phred",
                     "dbnsfp.polyphen2.hdiv.pred",
                     "dbnsfp.mutationtaster.pred")
  newnames <- c("Variant", "Ref.AA", "Alt.AA", "dbSNP rsid", "Consequence",
                     "COSMIC ID", "COSMIC Tumor Site", "ExAC AF", "CADD Score",
                     "Polyphen-2 Prediction", "MutationTaster Prediction")
  setnames(annos, old = intersect(oldnames, names(annos)), new = newnames[oldnames %in% names(annos)])
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
  #annotations[is.na(annotations)] <- ""
  annotations
}
args <- commandArgs(TRUE)
varscan <- args[1]
mutect <- args[2]
filter <- args[3]
som <- do.call(rbind.fill, lapply(c(varscan, mutect), function(i) somaticAnnotation(i, do_filter=filter)))
df <- som[order(som$Gene), ]

write.table(df, gsub("_varscan_somatic.vcf.gz", "_somatic_annotations.tsv", varscan), sep="\t", row.names=FALSE, quote=FALSE)
