## ---- IntogenConsequences ----
intogenConsequences <- function() {
  if(!is.null(config$INTOGEN_RESULTS)) {
    intogenConsequences <- read.csv2(paste(config$INTOGEN_RESULTS, '/', patientID, '/consequences.tsv', sep=''), header=T, sep="\t", na.strings="", as.is=T)
    
    intogenConsequencesTable <- as.data.table(subset(intogenConsequences, IMPACT_CLASS=='medium' | IMPACT_CLASS=='high'))
    intogenConsequencesTable <- intogenConsequencesTable[,c(2:9,13,16,19,22:24),with=F]
    # dTable(intogenConsequencesTable, sPaginationType='full_numbers', iDisplayLength=25, sScrollX='100%', width='100%')  
    kable(intogenConsequencesTable, format='html', table.attr = 'id=\"intogenConsequences_table\"')
  } else {
    print("No intogen results available")
  }
}

intogenConsequences()


## ---- IntogenVarGenes ----
intogenVarGenes <- function() {
  if(!is.null(config$INTOGEN_RESULTS)) {
    intogenVarGenes <- read.csv2(paste(config$INTOGEN_RESULTS, '/', patientID, '/variant_genes.tsv', sep=''), header=T, sep="\t", na.strings="", as.is=T)
    
    intogenVarGenesTables <- as.data.table(subset(intogenVarGenes, VAR_IMPACT_DESC=='medium' | VAR_IMPACT_DESC=='high'))
    intogenVarGenesTables <- intogenVarGenesTables[,c(2:9,13:16),with=F]
    intogenVarGenesTables[,PROTEIN_CHANGES:=gsub(',', ', ', intogenVarGenesTables[,PROTEIN_CHANGES])]
    #dTable(intogenVarGenesTables, sPaginationType='full_numbers', iDisplayLength=25, sScrollX='100%', width='100%')
    kable(intogenVarGenesTables, format='html', table.attr = 'id=\"intogenVarGenes_table\"')
  } else {
    print("No intogen results available")
  }
}

intogenVarGenes()