library(VariantAnnotation)

args <- commandArgs(TRUE)
vcf <- args[1]

filter_annovar <- function(annovarVcf){
  anno <- readVcf(annovarVcf, "hg19")
  nonsyn <- anno[!(unlist(info(anno)$ExonicFunc.knownGene) %in% c("synonymous_SNV", "."))]
  vcf <- nonsyn[unlist(info(nonsyn)$ExAC_ALL) < 0.05]
  writeVcf(vcf, filename = gsub(".hg19_multianno.vcf", "hg19_multianno.filtered.vcf", annovarVcf))
}

filter_annovar(vcf)