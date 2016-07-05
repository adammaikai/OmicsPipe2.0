## ---- GeneTable ---
geneTable <- function(patientCountsNorm, reference, group, diff=NULL, glist=NULL, relExpress=NULL, anno) {
  Patient <- round(as.numeric(patientCountsNorm),1)
  if(any(group=='Healthy')) {
    Normal <- round(as.numeric(apply(reference[, which(group=="Healthy")], 1, median)),1)
    NormalSD <- round(as.numeric(apply(reference[, which(group=="Healthy")], 1, sd)),1)
  } else {
    if (any(group=='Normal')) {
      Normal <- round(as.numeric(apply(reference[, which(group=="Normal")], 1, median)),1)
      NormalSD <- round(as.numeric(apply(reference[, which(group=="Normal")], 1, sd)),1)
    } else {
      Normal <- NA
      NormalSD <- NA
    }
  }

  Ref <- round(as.numeric(apply(reference[, which(group==config$TUMOR_TYPE)], 1, median)),1)
  RefSD <- round(as.numeric(apply(reference[, which(group==config$TUMOR_TYPE)], 1, sd)),1)
  
  if(!is.null(relExpress)) {
    relExpressed <- rep('', length(rownames(patientCountsNorm)))
    relExpressed[which(relExpress$Pcalls=='P' | relExpress$Pcalls=='M')] <- 'expressed'
    #relExpressed[which(relExpress$Pcalls!='P')] <- 'not reliable expressed'
  } else {
    relExpressed <- rep('', length(rownames(patientCountsNorm)))
  }
  
  if(!is.null(diff)) {
    diffPatient <- rep('', length(rownames(patientCountsNorm)))
    diffPatient[which(rownames(patientCountsNorm) %in% names(diff$diffupAll))] <- 'upregulated'
    diffPatient[which(rownames(patientCountsNorm) %in% names(diff$diffdownAll))] <- 'downregulated'
    
    geneDiffRank <- numeric(length=length(rownames(patientCountsNorm)))
    geneDiffRank[which(rownames(patientCountsNorm) %in% names(diff$rankUp))] <- na.omit(as.numeric(diff$rankUp[match(rownames(patientCountsNorm), names(diff$rankUp))]))
    geneDiffRank[which(rownames(patientCountsNorm) %in% names(diff$rankDown))] <- na.omit(as.numeric(diff$rankDown[match(rownames(patientCountsNorm), names(diff$rankDown))]))
    geneDiffRank[which(geneDiffRank==0)] <- NA
    
    pvalPat <- NULL
    pvalPat[which(rownames(patientCountsNorm) %in% names(diff$FCdiff))] <- na.omit(diff$padj[match(rownames(patientCountsNorm), names(diff$FCdiff))])
    pvalPat[is.na(pvalPat)] <- ''
  } else {
    diffPatient <- rep('', length(rownames(patientCountsNorm)))
    geneDiffRank <- rep(NA, length=length(rownames(patientCountsNorm)))
    pvalPat <- rep(NA, length=length(rownames(patientCountsNorm)))
  }

  # gene list 
  if(!is.null(glist)) {
    # convert the entrez ids in the list to symbols
    # note: since tcga dataset holds less genes, reduce input 
    # list by the one that do not match up, this will be resolved in future versions...
    glist <- as.vector(na.omit(anno$symbol[match(glist, anno$gene_id)]))
    geneList <- rep(NA, length=length(rownames(patientCountsNorm)))
    geneList[which(rownames(patientCountsNorm) %in% glist)] <- 'in_gene_list'
  } else {
    geneList <- rep(NA, length=length(rownames(patientCountsNorm)))
  }

  # calc the percentile rank for patient gene value within reference data
  percRankPat <- rank(Patient)/length(Patient)
  percRankRef <- rank(Ref)/length(Ref)
  
  allGeneTable <- data.table(Gene = rownames(patientCountsNorm),
                             Diff_expressed = diffPatient,
                             Reliable_epxressed = relExpressed,
                             Rank = geneDiffRank,
                             perRank = percRankPat,
                             Patient = Patient,
                             Normal = Normal,
                             Normal_sd = NormalSD,
                             d_normal = Patient - Normal,
                             #TNBC = as.numeric(apply(tnbcGeneList[, match(meta$IDR[which(meta$Tissue.Type=='TNBC')], colnames(tnbcGeneList))], 1, median)),
                             Ref = Ref,
                             perRankRef = percRankRef,
                             Ref_sd = RefSD,
                             d_ref = Patient - Ref,
                             padj = pvalPat,
                             List_Gene = geneList
  )
  allGeneTable <- transform(allGeneTable, 
                            Gene = paste('<a href = ', shQuote(paste('http://www.ncbi.nlm.nih.gov/gene/?term=', rownames(patientCountsNorm), sep='')), 'target="_blank" >', allGeneTable$Gene, '</a>')
  )
  allGeneTable <- transform(allGeneTable, 
                            GeneWiki = paste('<a href = ', shQuote(paste('http://en.wikipedia.org/wiki/', rownames(patientCountsNorm), sep='')), 'target="_blank" >', 'Link', '</a>')
  )
  
  return(allGeneTable)
}

allGeneTable <- geneTable(patientCountsNorm, reference, group, diff=diffExpNBinom, glist=geneList, relExpress=patientREC, anno)
#kable(allGeneTable, format='html', table.attr = 'id=\"gep_table\"')

# export table to json (aaData format)
dir.create(paste(config$REPORT_RESULTS, '/', patientID, '/json', sep=''), showWarnings = FALSE) #supress warning if dir already exists

require(RJSONIO)
ll <- list(aaData = lapply(seq(nrow(allGeneTable)), function(x) as.character(allGeneTable[x,])))
cat(toJSON(ll, pretty=FALSE), file=paste(paste(config$REPORT_RESULTS, '/', patientID, '/json/gepall.txt',sep='')))
