## ---- data ----
### this function scales the new sample based on stored loggeomean data 
### from the ref cohort
addNewSampleDESeq <- function(x, loggeomeansRef) {
  sf <- exp(median((log(x) - loggeomeansRef)[is.finite(loggeomeansRef)]))
  norm <- scale(x, center=FALSE, scale=sf) 
  return(norm)
}

### gene list (cancer type relevant genes)
geneList <- read.csv2(config$GENELIST, header=FALSE, stringsAsFactors=FALSE, sep='\t')$V2

### get the patients raw read counts 
patientCounts <- read.csv2(paste(config$HTSEQ_RESULTS, '/', patientID, '_counts.txt', sep=''), header=F, sep="\t", na.strings="", as.is=T, row.names=1)
ntemp <- rownames(patientCounts)
htseqStats <- patientCounts[(dim(patientCounts)[1]-4):(dim(patientCounts)[1]), ] # save htseq stats
names(htseqStats) <- ntemp[(dim(patientCounts)[1]-4):(dim(patientCounts)[1])]
patientCounts <- patientCounts[1:(dim(patientCounts)[1]-5), ] # remove the last 5 lines, they hold no genes
names(patientCounts) <- ntemp[1:length(patientCounts)]

### load reference cohort
# TCGA BRCA
if (config$TUMOR_TYPE=="BRCA") {
  # deseq normalized TCGA reference cohort  
  load(paste(config$R_SOURCE_PATH, '/data/deseq.tcga_brca.Rdata', sep=''))
  reference <- referenceBRCA$reference_norm # set the reference data to a more common variable that is used across all modules
  reference_count <- referenceBRCA$reference_raw_count
  reference_sf <- referenceBRCA$reference_sf
  group <- referenceBRCA$group # define groups within the reference cohort
  
  rownames(reference) <- unlist(lookUp(rownames(reference), 'org.Hs.eg', 'SYMBOL'))
  rownames(reference_count) <- unlist(lookUp(rownames(reference_count), 'org.Hs.eg', 'SYMBOL'))
  
  # log geometric mean value from TCGA cohort for add on preprocessing
  load(paste(config$R_SOURCE_PATH, '/data/loggeoameansBRCA.Rdata', sep='')) 
  
  # reduce patient genes to reference genes and vice versa
  # normalize new patient based on TCGA data
  # transform the data on log2 scale
  uniqueGenes <- intersect(rownames(reference), names(patientCounts))
  anno <- toTable(org.Hs.egSYMBOL)
  anno <- anno[match(uniqueGenes, anno$symbol),]

  patientCounts <- patientCounts[uniqueGenes]
  reference <- reference[uniqueGenes, ]
  reference_count <- reference_count[uniqueGenes, ]
  patientCountsNorm <- log2(addNewSampleDESeq(as.matrix(patientCounts), loggeomeansRef[anno$gene_id])+1)
}
# TCGA HNSC
if (config$TUMOR_TYPE=="HNSC") {
  # load ...
}
#TCGA BLCA
if (config$TUMOR_TYPE=="BLCA") {
  # load ...
}