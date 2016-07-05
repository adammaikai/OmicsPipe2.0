library(myvariant)
library(magrittr)
library(IRanges)
library(plyr)
library(httr)
library(jsonlite)
library(myvariant)
library(VariantAnnotation)
#source("~/avera/repos/variant-analysis/annotateIndels.R")

args <- commandArgs(TRUE)
vcfPath <- args[1]

.collapse <- function (...) {
    paste(unlist(list(...)), sep = ",", collapse = ",")
}

annotateVcf <- function(vcf.path, do_filter=TRUE){
  ## intogen file
  intogen <- read.csv("/data/database/druggability/Mutational_drivers_per_tumor_type.tsv", sep="\t", comment.char="#")
  names(intogen) <- c("Gene", "tumor_type")
  intogen <- aggregate(tumor_type ~ Gene, intogen, .collapse)
  
  ## panCancer file
  panCancer <- read.csv("/data/database/druggability/NanostringPanCancerGeneList.csv", sep=",", comment.char="#")
  names(panCancer) <- c("Gene", "reported")
  panCancer <- subset(aggregate(reported ~ Gene, panCancer, .collapse), !grepl("reference", reported) & !grepl("control", reported))
  
  ## intogen druggability
  druggability <- read.csv("/data/database/druggability/Protein_Drug_Interactions.tsv", sep="\t", comment.char="#")
  druggability <- druggability[, c("geneHGNCsymbol", "Drug_name")]
  druggability <- aggregate(Drug_name ~ geneHGNCsymbol, druggability, .collapse)
  names(druggability) <- c("Gene", "drug_name")
  
  ## MyVariant.info annotations
  vcf <- readVcf(vcf.path, genome="hg19")
  vcf <- vcf[isSNV(vcf) | isIndel(vcf)]
  hgvs <- formatHgvs(vcf)
  annos <- getVariants(hgvs, fields=c("dbsnp.rsid", "cadd.consequence", "dbnsfp.aa.pos", "dbnsfp.aa.ref",
                                      "dbnsfp.aa.alt", "cadd.gene.genename",
                                      "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af",
                                      "dbnsfp.1000gp1.af", "cadd.phred", "dbnsfp.sift.converted_rankscore", "dbnsfp.sift.pred",
                                      "dbnsfp.polyphen2.hdiv.rankscore", "dbnsfp.polyphen2.hdiv.pred", "dbnsfp.mutationtaster.converted_rankscore", 
                                      "dbnsfp.mutationtaster.pred", "dbnsfp.mutationassessor.rankscore", "dbnsfp.mutationassessor.pred", 
                                      "dbnsfp.lrt.converted_rankscore", "dbnsfp.lrt.pred", "dbnsfp.metasvm.rankscore", "dbnsfp.metasvm.pred"))
  names(annos)[names(annos) %in% "query"] <- "Variant"
  names(annos)[names(annos) %in% c("cadd.gene.genename", "cadd.gene")] <- "Gene"
  annos[c("notfound", "X_id", "cadd._license")] <- NULL
  row.names(vcf) <- NULL
  annos <- DataFrame(lapply(annos, function(i) sapply(i, .collapse)))
  annotationsIntogen <- merge(annos, intogen, all.x=TRUE)
  annotationsPanCancer <- merge(annotationsIntogen, panCancer, all.x=TRUE)
  annotationsDruggability <- merge(annotationsPanCancer, druggability, all.x=TRUE, sort=TRUE)
  annotations <- data.frame(sapply(annotationsDruggability, as.character), stringsAsFactors = FALSE)
  
  newInfo <- DataFrame(Number=1, sapply(annotations, class), names(annotations))
  names(newInfo) <- c("Number", "Type", "Description")
  info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
  ## merge annotations
  info(vcf)$Variant <- hgvs
  info(vcf) <- merge(annos, info(vcf), all.x=TRUE)
  if(do_filter){
    ## filter consequence
    vcf <- vcf[info(vcf)$cadd.consequence %in% c("STOP_GAINED","STOP_LOST", 
                                                 "NON_SYNONYMOUS", "SPLICE_SITE", 
                                                 "CANONICAL_SPLICE", "REGULATORY")]
    vcf <- vcf[is.na(info(vcf)$exac.af) | info(vcf)$exac.af < 0.1]
    vcf <- vcf[!is.na(info(vcf)$reported)]
  }
  writeVcf(vcf, gsub(".vcf", ".annotated.vcf", vcf.path))
  vcf
}

annotateVcf(vcfPath)
