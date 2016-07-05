## ---- DoG ------
pickDrug <- function(gene, dbFDA, regulation='up') {
  if(regulation=='up') {
    dic <- c('inhibitor', 'antagonist', 'antibody', 'suppressor')
  }
  if(regulation=='down') {
    dic <- c('agonist', 'potentiator')
  }
  dbFDAsub <- dbFDA[which(dbFDA$action %in% dic), ]
  res <- which(dbFDAsub$ENTREZ %in% gene)   
  if(length(res!=0)) {
    ids <- dbFDAsub$dbid[res]
    return(ids)
  } else {
    return(NA)
  }
}

DoG <- function(reference, diff) {
  library(AnnotationDbi) # for unlist2 function
  library(data.table)
  
  source(paste(config$R_SOURCE_PATH, '/src/modules/KEGGParse.R', sep=''))
  load(paste(config$R_SOURCE_PATH, '/data/DoG/dbFDA.Rdata', sep='')) # load dbFDA data frame
  dbFDA <- dbFDA3
  
  # head(dbFDA)
  #       action partner    dbid               names  FDA KEGGDrug  HGNC ENTREZ
  # 1  inhibitor      54 DB00001           Lepirudin TRUE   D06880  3535   2147
  # 2 antagonist     844 DB00002           Cetuximab TRUE   D03455  3236   1956
  # 3     binder     724 DB00004 Denileukin diftitox TRUE   D03682  6008   3559
  # 4    agonist     717 DB00004 Denileukin diftitox TRUE   D03682  6009   3560
  # 5   antibody     777 DB00005          Etanercept TRUE   D00742 11892   7124
  # 6  inhibitor      54 DB00006         Bivalirudin TRUE   D03136  3535   2147
  
  xx <- anno$gene_id[match(rownames(reference), anno$symbol)]
  
  allupdrugs <- unlist2(sapply(xx[match(names(diff$diffupAll), rownames(reference))], pickDrug, dbFDA=dbFDA))
  allupdrugs <- allupdrugs[-which(is.na(allupdrugs))]
  
  alldowndrugs <- unlist2(sapply(xx[match(names(diff$diffdownAll), rownames(reference))], pickDrug, dbFDA=dbFDA, regulation='down'))
  alldowndrugs <- alldowndrugs[-which(is.na(alldowndrugs))]
  
  eDrug <- c(allupdrugs, alldowndrugs)
  
  evaluationTable <- data.table(Drug=as.vector(eDrug),
                                KEGG=as.vector(dbFDA$KEGGDrug[match(eDrug, dbFDA$dbid)]),
                                Drug_Function=dbFDA$action[match(eDrug, dbFDA$dbid)],
                                Target=rownames(reference)[match(names(eDrug), xx)],
                                Target_FC=as.vector(round(diff$FCdiff[match(rownames(reference)[match(names(eDrug), xx)], names(diff$FCdiff))],2)),
                                Categorie=dbFDA$categorie[match(eDrug, dbFDA$dbid)]
  )
  
  # export as #csv
  evaluationTableExp <- evaluationTable
  evaluationTableExp$Name <- dbFDA$names[match(evaluationTable$Drug, dbFDA$dbid)]
  write.table(evaluationTableExp, paste(config$REPORT_RESULTS, '/', patientID, '/dog.csv', sep=''), row.names=FALSE, sep='\t')
  
  # create the html links
  evaluationTable <- transform(evaluationTable, 
                               Drug = paste('<a href = ', shQuote(paste('http://www.drugbank.ca/drugs/', evaluationTable$Drug, sep='')), 'target="_blank" >', dbFDA$names[match(evaluationTable$Drug, dbFDA$dbid)], '</a>'))
  evaluationTable <- transform(evaluationTable, 
                               KEGG = paste('<a href = ', shQuote(paste('http://www.genome.jp/dbget-bin/www_bget?dr:', evaluationTable$KEGG, sep='')), 'target="_blank" >', evaluationTable$KEGG, '</a>'))
  evaluationTable <- transform(evaluationTable, 
                               Target = paste('<a href = ', shQuote(paste('http://www.ncbi.nlm.nih.gov/gene/?term=', evaluationTable$Target, sep='')), 'target="_blank" >', evaluationTable$Target, '</a>'))
  kable(evaluationTable, format='html', table.attr = 'id=\"dog_table\"')
}
DoG(reference, diffExpNBinom)
