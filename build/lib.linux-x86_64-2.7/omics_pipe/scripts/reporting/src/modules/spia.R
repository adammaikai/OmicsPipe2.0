## ---- spiaPW ----
require(SPIA)
require(graphite)

spiaPW <- function(diff, reference) {
  up <- names(diff$FCdiff)[which(diff$FCdiff > 1)]
  down <- names(diff$FCdiff)[which(diff$FCdiff < -1)]
  
  spiaIn <- diff$FCdiff[c(up,down)]
  # background: all diff exprssed genes + all expressed genes in the patient sample
  spiaBackground <- unique(c(names(spiaIn), rownames(reference)[match(names(which(patientCountsNorm[,1]!=0)), rownames(reference))]))
  
  # convert symbols to entrez id 
  names(spiaIn) <- anno$gene_id[match(names(spiaIn), anno$symbol)]
  spiaBackground <- anno$gene_id[match(spiaBackground, anno$symbol)]
  
  # move PW stuff to working dir
  system(paste('cp', paste(config$R_SOURCE_PATH, '/data/SPIA/*.RData', sep=''), getwd(), sep=' '))
  
  kegg <-runSPIA(de=spiaIn, all=spiaBackground, organism="hsa", nB=2000, plots=FALSE, verbose=FALSE, "keggEx")
  biocarta <-runSPIA(de=spiaIn, all=spiaBackground, organism="hsa", nB=2000, plots=FALSE, verbose=FALSE, "biocartaEx")
  nc <-runSPIA(de=spiaIn, all=spiaBackground, organism="hsa", nB=2000, plots=FALSE, verbose=FALSE, "nciEx")
  react <-runSPIA(de=spiaIn, all=spiaBackground, organism="hsa", nB=2000, plots=FALSE, verbose=FALSE, "reactomeEx")
  
  kegg.sig <- kegg[which(kegg$pGFdr <= 0.05), ]
  biocarta.sig <- biocarta[which(biocarta$pGFdr <= 0.05), ]
  nc.sig <- nc[which(nc$pGFdr <= 0.05), ]
  react.sig <- react[which(react$pGFdr <= 0.05), ]
  
  # remove PW stuff
  system(paste('rm ', getwd(), '/*.RData', sep=''))
  
  return(list(kegg=kegg.sig,
              biocarta=biocarta.sig,
              nci=nc.sig,
              reactome=react.sig
              )
         )
}

spiaRES <- spiaPW(diffExpNBinom, reference)

