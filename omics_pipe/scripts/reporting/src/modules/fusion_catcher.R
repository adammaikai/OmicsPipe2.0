## ---- fusionCatcher ----
fusionTable <- function() {
  if(!is.null(config$FUSIONCATCHER_RESULTS)) {
    fusions <- read.csv2(paste(config$FUSIONCATCHER_RESULTS, '/', patientID, '/final-list_candidate-fusion-genes.txt', sep=''), header=T, sep="\t", na.strings="", as.is=T, stringsAsFactor=FALSE)
    if(dim(fusions)[1]==0) {
      print('No fusions detected')
    } else {
      oncofuse <- read.csv2(paste(config$FUSIONCATCHER_RESULTS, '/', patientID, '/oncofuse_res.txt', sep=''), header=T, sep="\t", na.strings="", as.is=T, stringsAsFactor=FALSE) 
      fusionTable <- data.table(fusions[,-c(11,12)])
      colnames(fusionTable) <- c('fusion_Partner_5',
                                 'fusion_Partner_3',
                                 'fusion_description',
                                 '#_common_mapping_reads',
                                 'spanning_pairs',
                                 'spanning_uq_reads',
                                 'longest_anchor',
                                 'method',
                                 'fusion_point_5',
                                 'fusion_point_3',
                                 'exon_id_5',
                                 'exon_id_3',
                                 'sequence'
      ) 
      # add results from oncofuse prediction
      fusionID <- paste(fusions$Gene_1_symbol.5end_fusion_partner., fusions$Gene_2_symbol.3end_fusion_partner., sep='')
      oncofuseID <- paste(oncofuse$X5_FPG_GENE_NAME, oncofuse$X3_FPG_GENE_NAME, sep='')
      
      fusionTable <- cbind(fusionTable, driver_prob=NA, expression_gain=NA)
      fusionTable$driver_prob <- round(as.numeric(as.vector(oncofuse$DRIVER_PROB[match(fusionID, oncofuseID)])),2)
      fusionTable$expression_gain <- round(as.numeric(as.vector(oncofuse$EXPRESSION_GAIN[match(fusionID, oncofuseID)])),2)
      
      #export as csv
      write.table(fusionTable, paste(config$REPORT_RESULTS, '/', patientID, '/fusions.csv', sep=''), row.names=FALSE, sep='\t')
      
      kable(fusionTable, format='html', table.attr = 'id=\"fusion_table\"')
    }
  } else {
    print('No FusionCatcher results available')
  }
}
fusionTable()