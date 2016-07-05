## ---- diffExpressionBasic ----
# all 'diff' expressed genes between patient and normal 
# diff. meaning 3 sd above normal.. 
# this is a "basic" diff. expression calling
# if there is a healthy reference, use this, otherwise use matched tumor normals

# TO-DO: 
# make desiccion wether to include marginal expressed genes into diff. expression
# calling a configurable parameter in the config file

diffExpressionBasic <- function(patientCountNorm, reference, group, patientREC=NULL) {
  if(any(group=='Healthy')) {
    medianNormal <- apply(reference[, which(group=='Healthy')], 1, median)  
    sdNormal <- apply(reference[, which(group=='Healthy')], 1, sd)
  } else {
    medianNormal <- apply(reference[, which(group=='Normal')], 1, median)  
    sdNormal <- apply(reference[, which(group=='Normal')], 1, sd) 
  }
  
  # reduce differential expression calling to the list of reliable expressed genes if available..
  if(!is.null('patientREC')) {
    relExpressed <- names(patientREC$Pcalls)[which(patientREC$Pcalls=='P' | patientREC$Pcalls=='M')]
  } else {
    relExpressed <- rownames(patientCountsNorm)
  }
  
  diffupAll <- which(patientCountsNorm[match(relExpressed, rownames(patientCountsNorm)),1] > 
                       (medianNormal[match(relExpressed, names(medianNormal))] + 3 * sdNormal[match(relExpressed, names(sdNormal))]))
  diffdownAll <- which(patientCountsNorm < (medianNormal - 3 * sdNormal))
  names(diffdownAll) <- rownames(patientCountsNorm)[diffdownAll]
  
  FCdiffx <- patientCountsNorm - medianNormal
  FCdiff <- FCdiffx[match(c(names(diffupAll), names(diffdownAll)), rownames(FCdiffx))]
  names(FCdiff) <- c(names(diffupAll), names(diffdownAll))
  
  # rank the diff. expressed genes
  rankUp <- 1:(length(diffupAll))
  names(rankUp) <- names(rev(sort((FCdiff[match(names(diffupAll), names(FCdiff))]))))
  rankDown <- 1:(length(diffdownAll))
  names(rankDown) <- names(sort((FCdiff[match(names(diffdownAll), names(FCdiff))])))
  
  return(list(diffupAll=diffupAll,
              diffdownAll=diffdownAll,
              FCdiff=FCdiff,
              rankUp=rankUp,
              rankDown=rankDown,
              relExpressed=relExpressed))
}

diffExpBasic <- diffExpressionBasic(patientCountsNorm, reference, group, patientREC=patientREC)

## ---- diffExpressionNBinom ----
diffNBinom <- function(referece_count, patientCounts, filter=NULL) {
  library(DESeq2)
  
  if(any(group=='Healthy')) {
    dds <- DESeqDataSetFromMatrix(countData = cbind(reference_count[,group %in% 'Healthy'], patientCounts),
                                  colData = data.frame(condition=rep(c('Healthy', 'Patient'), c(length(which(group=='Healthy')), 1))),
                                  design = ~ condition)
    sizeFactors(dds) <- c(reference_sf[match(colnames(cbind(reference_count[,group %in% 'Healthy'], patientCounts))[1:length(which(group=='Healthy'))], names(reference_sf))],
                          exp(median((log(patientCounts) - loggeomeansRef)[is.finite(loggeomeansRef)]))) # use the size factors from the global dataset
    dds <- estimateDispersions(dds)
    if(!is.null(filter)) {
      dds <- dds[filter,]
    }
    dds <- nbinomWaldTest(dds)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = cbind(reference_count[,group %in% 'Normal'], patientCounts),
                                  colData = data.frame(condition=rep(c('Normal', 'Patient'), c(length(which(group=='Normal')), 1))),
                                  design = ~ condition)
    sizeFactors(dds) <- c(reference_sf[match(colnames(cbind(reference_count[,group %in% 'Normal'], patientCounts))[1:length(which(group=='Normal'))], names(reference_sf))],
                          exp(median((log(patientCounts) - loggeomeansRef)[is.finite(loggeomeansRef)]))) # use the size factors from the global dataset
    dds <- estimateDispersions(dds)
    if(!is.null(filter)) {
      dds <- dds[filter,]
    }
    dds <- nbinomWaldTest(dds)
  }

  res <- results(dds)
  res <- as.data.frame(res)
  res <- cbind(id=rownames(res), res)
  resSig <- arrange(res[which(res$padj < 0.05), ], padj)
  # resSig$log2FC <- log2(resSig$baseMeanA+1) - log2(resSig$baseMeanB+1)
  
  diffupAll <- which(resSig$log2FoldChange>0)
  names(diffupAll) <- resSig$id[which(resSig$log2FoldChange>0)]
  diffdownAll <- which(resSig$log2FoldChange<0)
  names(diffdownAll) <- resSig$id[which(resSig$log2FoldChange<0)]
  
  FCdiff <- resSig$log2FoldChange
  names(FCdiff) <- resSig$id
  
  # rank the diff. expressed genes
  rankUp <- 1:(length(diffupAll))
  names(rankUp) <- names(rev(sort((FCdiff[match(names(diffupAll), names(FCdiff))]))))
  rankDown <- 1:(length(diffdownAll))
  names(rankDown) <- names(sort((FCdiff[match(names(diffdownAll), names(FCdiff))])))
  
  return(list(diffupAll=diffupAll,
              diffdownAll=diffdownAll,
              FCdiff=FCdiff,
              rankUp=rankUp,
              rankDown=rankDown,
              relExpressed=filter,
              padj=resSig$padj))
}

diffExpNBinom <- diffNBinom(reference_count, patientCounts, filter=names(which(patientREC$Pcalls=='P' | patientREC$Pcalls=='M')))

## ---- diffExpPlot ----
diffExpPlot <- function(diff) {
  diffExpTable <- data.frame(status=c('upregulated',
                                      'downregulated',
                                      'unchanged'),
                             Freq=c(length(diff$diffupAll),
                                    length(diff$diffdownAll),
                                    length(diff$relExpressed)-length(diff$diffupAll)-length(diff$diffdownAll)
                             )
  )
  diffExpPlot <- rCharts:::Highcharts$new()
  diffExpPlot$chart(type = "column")
  diffExpPlot$title(text = "Differentially Expressed Genes")
  diffExpPlot$xAxis(categories = diffExpTable$status)
  diffExpPlot$yAxis(title = list(text = "number of genes"))
  diffExpPlot$data(diffExpTable)
  diffExpPlot$params$width = 800
  diffExpPlot$params$height = 400
  diffExpPlot
}
diffExpPlot(diffExpNBinom)



