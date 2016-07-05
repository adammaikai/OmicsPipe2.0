library(VariantAnnotation)
library(plyr)

args <- commandArgs(TRUE)
varscanVcf <- args[1]
outVcf <- args[2]

varscan2Vcf <- function(varscan.vcf, out.vcf) {
  varscan<- read.csv(varscan.vcf, stringsAsFactors = FALSE, header = TRUE, sep="\t")
  varscan[c("Position", "Reads1", "Reads2")] <- as(c(varscan$Position, varscan$Reads1, varscan$Reads2), "integer")
  snps <- subset(varscan, !grepl("\\-|\\+", VarAllele))
  snps[c("ref", "alt")] <- snps[c("Ref", "VarAllele")]
  dels <- subset(varscan, grepl("\\-", VarAllele))
  dels$ref <- gsub("\\-", "", paste(dels$Ref, dels$VarAllele, sep=""))
  dels$alt <- dels$Ref
  ins <- subset(varscan, grepl("\\+", VarAllele))
  ins$ref <- ins$Ref
  ins$alt <- gsub("\\+", "", paste(ins$Ref, ins$VarAllele, sep=""))
  varscan <- do.call(rbind.fill, list(snps, dels, ins))
  vr <- VRanges(seqnames=Rle(varscan$Chrom), ranges=IRanges(varscan$Position, varscan$Position), 
              ref=varscan$ref, alt=varscan$alt, refDepth=varscan$Reads1, altDepth=varscan$Reads2, 
              totalDepth = (varscan$Reads1 + varscan$Reads2), AF=as.numeric(gsub("%", "", varscan$VarFreq))/100,
              sampleNames=gsub(".vcf", "", (lapply(strsplit(varscan.vcf, "/"), tail, 1))))
  info <- DataFrame(varscan[c("Cons", "VarFreq", "Strands1", "Strands2", "Qual1", "Qual2", "Pvalue", 
                                          "MapQual1", "MapQual2", "Reads1Plus", "Reads1Minus", "Reads2Plus", "Reads2Minus")])
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
