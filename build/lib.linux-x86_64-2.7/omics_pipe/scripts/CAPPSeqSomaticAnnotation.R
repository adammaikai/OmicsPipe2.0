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
  ## filter consequence
  if(do_filter){
    annos <- subset(annos, cadd.consequence %in% c("STOP_GAINED","STOP_LOST", 
                                                   "NON_SYNONYMOUS", "SPLICE_SITE", 
                                                   "CANONICAL_SPLICE", "REGULATORY"))
    annos <- arrange(data.frame(subset(annos, is.na(exac.af) | exac.af < 0.05)), -dbnsfp.mutationtaster.converted_rankscore)
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

  annos <- DataFrame(annos) ##for some reason have to do this to eliminate the following columns  
  annos[c("X_id", "notfound", "X_score", "cadd._license")] <- NULL
  if("cadd.gene.genename" %in% names(annos)){Gene <- annos$cadd.gene.genename}else{Gene <- sapply(annos$cadd.gene, function(i) i$genename)}

  annos$Gene <- Gene
  if("cadd.gene.prot.protpos" %in% names(annos)){aaprot <- annos$cadd.gene.prot.protpos}else{aaprot <- lapply(annos$cadd.gene, function(i) .collapse(i$prot))}
  annos$`Amino Acid Position` <- aaprot
  annos <- DataFrame(annos) ##for some reason have to do this to eliminate the following columns
  annos[c("X_id", "notfound", "X_score", "cadd._license")] <- NULL
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
varscan <- args[1]
filter <- args[2]
som <- somaticAnnotation(varscan, do_filter=filter)
df <- som[order(som$Gene), ]

write.table(df, gsub("vcf.gz", "somatic_annotations.txt", varscan), sep="\t", row.names=FALSE, quote=FALSE)

