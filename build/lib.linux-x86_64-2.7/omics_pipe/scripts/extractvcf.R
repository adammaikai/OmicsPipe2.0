args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
skip <- args[2]

x <- read.csv2(paste(dir, '/raw_variants.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed', sep=''), sep='\t', header=FALSE)
y <- read.csv2(paste(dir, '/raw_variants.vcf', sep=''), sep='\t', skip=skip, header=FALSE)

colnames(x)[c(1,3)] <- c('chr', 'loc')
colnames(y)[c(1,2)] <- c('chr', 'loc')

extract <- y[match(x$V3, y$V2),]
extract <- merge(x,y,by=c('chr','loc'))
extract <- extract[,c(1,2,8:15)]
write.table(extract, file=paste(dir, '/red_variants.vcf', sep=''), sep='\t', quote=FALSE, col.names=F, row.names=F)