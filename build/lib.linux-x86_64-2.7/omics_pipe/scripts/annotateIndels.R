library(VariantAnnotation)
library(myvariant)
library(mygene)
library(data.table)

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
  
  dp <- data.frame(geno(indel)$DP)
  row.names(dp) <- NULL
  # Allelic Depth
  ad <- data.frame(geno(indel)$AD)
  row.names(ad) <- NULL
  p_adjust <- c()
  if (grepl("varscan", vcf.path)){
  ## Flag Strand Bias Variants
    fs <- read.csv(gsub("vcf.gz", "indel.strand.dp.csv", vcf.path), stringsAsFactors = FALSE, header = FALSE)
    names(fs) <- c("ref_fwd", "ref_rev", "var_fwd", "var_rev")
    fs$p_value <- sapply(seq(1, nrow(fs)), function(i) fisher.test(matrix(c(fs$ref_fwd[i], fs$ref_rev[i], fs$var_fwd[i], fs$var_rev[i]), nrow = 2), workspace = 2000000000)[[1]])
    p_adjust <- p.adjust(fs$p_value, method = "BH")
  ## Read Depth
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
  } else{
    coverage <- cbind(ad, dp)
    names(coverage) <- c("AD", "DP")
  }
  if (length(p_adjust) == 0){
  p_adjust <- rep(0, nrow(length(hits)))
  }
  annotations <- DataFrame(Position=paste0(seqnames(indel), ":", start(indel)), Variant=Variant, Gene=sapply(hits, .collapse), FS_bias_p_adjust=p_adjust)
  annotations <- cbind(annotations, coverage)
 # annotations$FS_bias_adj_p_value <- fs$p_adjust
  annotations <- data.frame(merge(annotations, cancer_genes, by.x="Gene", by.y="symbol", all.x=TRUE))
#  annotations <- sapply(annotations, as.character)
  annotations
}
