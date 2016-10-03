library(myvariant)
library(mygene)
library(magrittr)
library(IRanges)
library(plyr)
library(httr)
library(jsonlite)
library(myvariant)
library(VariantAnnotation)
library(data.table)

args <- commandArgs(TRUE)
vcfPath <- args[1]
filter <- args[2]


.collapse <- function (...) {
     paste(unlist(list(...)), sep = ",", collapse = ",")
 }

annotateIndels <- function(vcf.path){

  cancer_genes <- read.csv("/data/database/druggability/cancer_genes.txt", stringsAsFactors = FALSE, sep="\t")

  vcf.object <- readVcf(vcf.path, genome="hg19")
  indel <- vcf.object[isIndel(vcf.object)]
  if (dim(indel)[1] == 0){
    return(data.frame())
  }
  loc <- paste("hg19.", seqnames(indel), ":", start(indel), "-", end(indel), sep="")
  Variant <- formatHgvs(indel, variant_type=c("insertion", "deletion"))
  hits <- lapply(loc, function(i) query(i, species="human")$hits$symbol)
  annotations <- DataFrame(Position=paste0(seqnames(indel), ":", start(indel)), Variant=Variant, Gene=sapply(hits, .collapse))
  annotations <- data.frame(merge(annotations, cancer_genes, by.x="Gene", by.y="symbol", all.x=TRUE))
  annotations
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
    ## filter consequence
    if(do_filter){
      annos <- subset(annos, cadd.consequence %in% c("STOP_GAINED","STOP_LOST", 
                                                   "NON_SYNONYMOUS", "SPLICE_SITE", 
                                                   "CANONICAL_SPLICE", "REGULATORY"))
      annos <- data.frame(subset(annos, is.na(exac.af) | exac.af < 0.05))}
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
    if("cadd.gene.genename" %in% names(annos)){Gene <- annos$cadd.gene.genename}else{Gene <- sapply(annos$cadd.gene, function(i) i$genename)}

    annos$Gene <- Gene
    if("cadd.gene.prot.protpos" %in% names(annos)){aaprot <- annos$cadd.gene.prot.protpos}else{aaprot <- lapply(annos$cadd.gene, function(i) .collapse(i$prot))}
    annos$`Amino Acid Position` <- aaprot
    annos <- lapply(annos, function(i) sapply(i, .collapse))
    annos <- data.frame(annos)
    ## merge annotations
    annos$in_cancer <- ifelse(annos$Gene %in% cancer_genes$symbol, yes="yes", no="no")
    annotations <- merge(annos, cancer_genes, by.x="Gene", by.y="symbol", all.x=TRUE)
    annotations <- rbind.fill(annos, annotateIndels(vcf.path))
    write.table(annotations, gsub(".vcf.gz", ".annotated.txt", vcf.path), sep="\t", row.names=FALSE, quote=FALSE)
}

annotateVariantsFromVcf(vcfPath, do_filter=filter)
