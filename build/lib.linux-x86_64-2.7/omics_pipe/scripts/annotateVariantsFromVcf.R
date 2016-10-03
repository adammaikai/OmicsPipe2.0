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
    
    ## Read in VCF
    vcf <- readVcf(vcf.path, genome="hg19")
    snp <- vcf[isSNV(vcf)]

    ## Filter out variants with strand bias
#    fs <- read.csv(gsub("vcf.gz", "snp.strand.dp.csv", vcf.path), stringsAsFactors = FALSE)
#    names(fs) <- c("ref_fwd", "ref_rev", "var_fwd", "var_rev")
#    fs$p_value <- sapply(seq(1, nrow(fs)), function(i) fisher.test(matrix(c(fs$ref_fwd[i], fs$ref_rev[i], fs$var_fwd[i], fs$var_rev[i]), nrow = 2), workspace = 2000000000)[[1]])
#    fs$p_adjust <- p.adjust(fs$p_value, method = "BH")

    ## MyVariant.info annotations
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
#   annos$FS_bias_adj_p_value <- fs$p_adjust
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
    #annotations[is.na(annotations)] <- ""
    write.table(annotations, gsub(".vcf.gz", ".annotated.tsv", vcf.path), sep="\t", row.names=FALSE, quote=FALSE)
    #annotations
}

annotateVariantsFromVcf(vcfPath, do_filter=filter)
