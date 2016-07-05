library(myvariant)
library(magrittr)
library(IRanges)
library(plyr)
library(httr)
library(jsonlite)
library(myvariant)
library(VariantAnnotation)
library(data.table)

args <- commandArgs(TRUE)
annotateIndelsScript <- args[1]
source(annotateIndelsScript)
vcfPath <- args[2]
filter <- args[3]


.collapse <- function (...) {
     paste(unlist(list(...)), sep = ",", collapse = ",")
 }

annotateVariantsFromVcf <- function(vcf.path, do_filter=FALSE){
    ## Avera in-house gene panel
    cancer_genes <- read.csv("/data/database/druggability/cancer_genes.txt", stringsAsFactors = FALSE, sep="\t")
    
    ## MyVariant.info annotations
    vcf <- readVcf(vcf.path, genome="hg19")
    snp <- vcf[isSNV(vcf)]
    hgvs <- formatHgvs(snp, "snp")
    annos <- getVariants(hgvs, fields=c("cadd.gene.prot.protpos", "cadd.oaa", "cadd.naa", "dbsnp.rsid", "cadd.consequence", 
                                        "cadd.gene.genename",
                                        "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af", "cadd.phred", "dbnsfp.polyphen2.hdiv.pred", 
                                        "dbnsfp.mutationtaster.pred"
    ))
    annos$Position <- paste0(seqnames(snp), ":", start(snp))
      dp <- geno(snp)$DP
      row.names(dp) <- NULL
      ad <- geno(snp)$AD
      row.names(ad) <- NULL
      coverage <- cbind(data.frame(ad), dp)
      names(coverage) <- c("AD", "DP")
    ## filter consequence
    if(do_filter){
      annos <- subset(annos, cadd.consequence %in% c("STOP_GAINED","STOP_LOST", 
                                                   "NON_SYNONYMOUS", "SPLICE_SITE", 
                                                   "CANONICAL_SPLICE", "REGULATORY"))
      annos <- data.frame(subset(annos, is.na(exac.af) | exac.af < 0.05))}
    annos <- data.frame(annos)
    setnames(annos, 
             old = c("query", "cadd.naa", "cadd.oaa", "dbsnp.rsid", "cadd.consequence",
                     "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af", "cadd.phred",
                     "dbnsfp.polyphen2.hdiv.pred",
                     "dbnsfp.mutationtaster.pred"), 
             new = c("Variant", "Ref.AA", "Alt.AA", "dbSNP rsid", "Consequence",
                     "COSMIC ID", "COSMIC Tumor Site", "ExAC AF", "CADD Score",
                     "Polyphen-2 Prediction", "MutationTaster Prediction"))
    #names(annos)[names(annos) %in% c("cadd.gene.genename", "cadd.gene")] <- "Gene"
    annos <- DataFrame(annos) ##for some reason have to do this to eliminate the following columns
    annos[c("X_id", "notfound", "X_score", "cadd._license")] <- NULL
    annos$Gene <- sapply(annos$cadd.gene, function(i) i$genename)
    annos$`Amino Acid Position` <- lapply(annos$cadd.gene, function(i) .collapse(i$prot))
    annos <- lapply(annos, function(i) sapply(i, .collapse))
    annos <- data.frame(annos)
    ## merge annotations
    annos$in_cancer <- ifelse(annos$Gene %in% cancer_genes$symbol, yes="yes", no="no")
    #annotations <- merge(annos, cancer_genes, by.x="Gene", by.y="symbol")
    annotations <- rbind.fill(annos, annotateIndels(vcf.path))
    #annotations[is.na(annotations)] <- ""
    write.table(annotations, gsub(".vcf.gz", ".annotated.txt", vcf.path), sep="\t", row.names=FALSE, quote=FALSE)
    #annotations
}

annotateVariantsFromVcf(vcfPath, do_filter=filter)
