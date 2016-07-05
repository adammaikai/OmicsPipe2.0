args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
maf <- args[2]
file <- args[3]
skip <- args[4]

maf_file <- as.numeric(read.csv2(maf, sep='\t', header=TRUE, stringsAsFactors=FALSE)[,1])
variant_file <- read.csv2(file, sep='\t', skip=skip, header=FALSE)

zz <- which(maf_file >= 0.01)
v <- variant_file[zz,]

write.table(v, file=paste(dir, '/maf_filt.vcf', sep=''), sep='\t', quote=FALSE, col.names=F, row.names=F)