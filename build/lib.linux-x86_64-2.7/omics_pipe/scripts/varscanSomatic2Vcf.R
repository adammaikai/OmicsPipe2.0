library(VariantAnnotation)
library(plyr)

args <- commandArgs(TRUE)
varscanVcf <- args[1]
outVcf <- args[2]

varscan2Vcf <- function(varscan.vcf, out.vcf) {
  varscan<- read.csv(varscan.vcf, header = TRUE, sep="\t", 
            colClasses=c("character", "integer", "character", "character", "integer", "integer", "character", "character", "integer",
                          "integer", "character", "character", "character", "numeric", "numeric", "integer", "integer", 
                           "integer", "integer", "integer", "integer", "integer", "integer"))
  snps <- subset(varscan, !grepl("\\-|\\+", var))
  snps[c("REF", "ALT")] <- snps[c("ref", "var")]
  dels <- subset(varscan, grepl("\\-", var))
  dels$REF <- gsub("\\-", "", paste(dels$ref, dels$var, sep=""))
  dels$ALT <- dels$ref
  ins <- subset(varscan, grepl("\\+", var))
  ins$REF <- ins$ref
  ins$ALT <- gsub("\\+", "", paste(ins$ref, ins$var, sep=""))
  varscan <- do.call(rbind.fill, list(snps, dels, ins))
  vr <- VRanges(seqnames=Rle(varscan$chrom), ranges=IRanges(varscan$position, varscan$position), 
              ref=varscan$REF, alt=varscan$ALT, refDepth=varscan$tumor_reads1, altDepth=varscan$tumor_reads2, 
              totalDepth = (varscan$tumor_reads1 + varscan$tumor_reads2), AF=as.numeric(gsub("%", "", varscan$tumor_var_freq))/100,
              sampleNames=gsub(".vcf", "", (lapply(strsplit(varscan.vcf, "/"), tail, 1))))
  info <- DataFrame(varscan[c(5:23)])
  vcf <- as(vr, "VCF")
  suppressWarnings(info(vcf) <- info)
  row.names(vcf) <- NULL
  newInfo <- DataFrame(Number=1, sapply(info(vcf), class), names(info(vcf)))
  names(newInfo) <- c("Number", "Type", "Description")
  info(header(vcf)) <- newInfo
  writeVcf(vcf, out.vcf)
  vcf
}

varscan2Vcf(varscanVcf, outVcf)
