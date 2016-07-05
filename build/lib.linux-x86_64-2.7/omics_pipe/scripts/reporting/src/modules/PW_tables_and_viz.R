## ---- PWFunctions ----
suppressPackageStartupMessages(suppressWarnings(library(tcltk))) # silend load tcltk to supress no Display warning on headless nodes

my.netplot <- function (g, vertex.label.font = 2, vertex.label.color = "#666666", 
          vertex.label.cex = 1.5, layout = layout.fruchterman.reingold, 
          foldChange = NULL, fixed = TRUE, col.bin = 10, legend.x = 1, 
          legend.y = 1) 
{
  if (fixed) {
    plot.igraph(g, vertex.label.font = vertex.label.font, 
                vertex.label.color = vertex.label.color, vertex.label.cex = vertex.label.cex, 
                vertex.frame.color = V(g)$color, layout = layout)
    if (!is.null(foldChange)) {
      fc <- foldChange
      lbs <- hist(fc, breaks = col.bin - 1, plot = FALSE)$breaks
      col.legend <- DOSE:::get.col.scale(lbs)
      x <- seq(from = legend.x, by = 0.03, length.out = col.bin)
      y <- rep(legend.y, col.bin)
      points(x, y, pch = 15, col = col.legend, cex = 2)
      idx <- c(1, seq(4, col.bin - 1, by = 3), col.bin)
      text(x = x[idx], y = rep(legend.y - 0.05, length(idx)), 
           label = lbs[idx], cex = 0.8)
      text(x = mean(x), y = legend.y + 0.05, labels = "Fold Change", 
           cex = 0.8, font = 2)
    }
  }
  else {
    plot.igraph(g, vertex.label.font = vertex.label.font, vertex.label.color = vertex.label.color, 
           vertex.label.cex = vertex.label.cex, vertex.frame.color = V(g)$color, 
           layout = layout)
  }
}

my.viewPathway <- function (pathName, organism = "human", readable = TRUE, foldChange = NULL, db = 'reactome', ...) 
{
  require(graphite)
  require(igraph)
  db <- eval(parse(text = db))
  p <- db[[pathName]]
  if (is.null(p)) {
    plot.new() # instead of breaking here, just create an empty white plot..
   return('pathway not found ...')
  }
  if (organism != "human") {
    stop("the specific organism is not supported yet...")
  }
  if (readable) {
    p <- convertIdentifiers(p, "symbol")
    if (!is.null(foldChange)) {
      gn <- DOSE:::EXTID2NAME(names(foldChange), organism)
      names(foldChange) <- gn
    }
  }
  else {
    if (!is.null(foldChange)) {
      p <- convertIdentifiers(p, "entrez")
    }
  }
  g <- pathwayGraph(p)
  gg <- igraph.from.graphNEL(g)
  gg <- DOSE:::setting.graph.attributes(gg)
  if (!is.null(foldChange)) {
    gg <- DOSE:::scaleNodeColor(gg, foldChange)
  }
  #DOSE:::netplot(gg, foldChange = foldChange, ...)
  my.netplot(gg, foldChange = foldChange, ...)
}

## ---- KEGGTable ----
keggViz <- function(pw, geneList, patientID, anno, dir) {
  require(pathview)
  require(KEGGREST)
  
  fc <- geneList
  names(fc) <- anno$gene_id[match(names(fc), anno$symbol)]
  
  xx <- keggList("pathway")
  if(!is.na(match(pw, xx))) {
    keggID <- gsub('path:map', '', names(xx[match(pw, xx)]))
    
    pathview(gene.data = fc, 
             pathway.id = keggID, 
             species = "hsa", 
             out.suffix = patientID, 
             kegg.native = T, 
             discrete = list(gene = T, cpd = F), 
             limit = list(gene = 15, cpd = 1), 
             na.col = "gray",
             kegg.dir=paste(dir, '/', patientID, '/figure', sep='')
    )
    
    # as the .png will be in the workdir, move it to figures..
    system(paste('mv', paste('hsa', keggID, '.png', sep=''), paste(getwd(), '/figure', sep='')))
    
    return(paste('hsa', keggID, '.png', sep=''))
  } else {
    return(NA)
  }
}

keggTable <- function(kegg.sig, geneList, patientID, dir) {
  if(dim(kegg.sig)[1]!=0) {
    keggTable <- as.data.table(kegg.sig)
    keggTable <- cbind(keggTable, Image=NA)
    
    # create the kegg pw images and get the filenames
    keggVizFiles <- sapply(spiaRES$kegg$Name, function(x) {
      keggViz(x, geneList, patientID, anno, config$REPORT_RESULTS)
    })
    
    keggTable <- transform(keggTable, 
                           Image = paste('<a href = ', 
                                         paste('figure/', as.vector(keggVizFiles), sep=''), 
                                         ' target="_blank"> <img height=50 width=50 src=', 
                                         paste('figure/', as.vector(keggVizFiles), sep=''), 
                                         '></a>', 
                                         sep=''
                                         )
    )
    kable(keggTable, format='html', table.attr = 'id=\"kegg_table\"')
  } else print('No sig. KEGG pathways')  
}
keggTable(spiaRES$kegg, diffExpNBinom$FCdiff, patientID, config$REPORT_RESULTS)

## ---- BiocartaPlots ----
biocartaViz <- function(pw, geneList, anno) {
  require(ReactomePA)
  fc <- geneList
  names(fc) <- anno$gene_id[match(names(fc), anno$symbol)]
  my.viewPathway(pw, readable=TRUE, foldChange=fc, db='biocarta')  
}

# create the biocartz pw images
if(length(spiaRES$biocarta$Name)>0) {
  sapply(spiaRES$biocarta$Name, function(x) {
    biocartaViz(x, diffExpNBinom$FCdiff, anno)
  })
}

## ---- BiocartaTable ----
biocartaTable <- function(biocarta.sig) {
  if(dim(biocarta.sig)[1]!=0) {
    biocartaTable <- as.data.table(biocarta.sig)
    biocartaTable <- cbind(biocartaTable, Image=NA)
    biocartaTable <- transform(biocartaTable, 
                          Image = paste('<a href = figure/BiocartaPlots', seq(1:length(spiaRES$biocarta$Name)), '.png target="_blank"', 
                                        '> <img height=50 width=50 src=figure/BiocartaPlots', seq(1:length(spiaRES$biocarta$Name)), 
                                        '.png>', '</a>', sep='')
    )
    kable(biocartaTable, format='html', table.attr = 'id=\"biocarta_table\"')
  } else print('No sig. Biocarta pathways')
}
biocartaTable(spiaRES$biocarta)

## ---- NCIPlots ----
nciViz <- function(pw, geneList, anno) {
  require(ReactomePA)
  fc <- geneList
  names(fc) <- anno$gene_id[match(names(fc), anno$symbol)]
  my.viewPathway(pw, readable=TRUE, foldChange=fc, db='nci')  
}

# create the nci pw images
if(length(spiaRES$nci$Name)>0) {
  sapply(spiaRES$nci$Name, function(x) {
    nciViz(x, diffExpNBinom$FCdiff, anno)
  })
}

## ---- NCITable ----
nciTable <- function(nci.sig) {
  if(dim(nci.sig)[1]!=0) {
    nciTable <- as.data.table(nci.sig)
    nciTable <- cbind(nciTable, Image=NA)
    nciTable <- transform(nciTable, 
                          Image = paste('<a href = figure/NCIPlots', seq(1:length(spiaRES$nci$Name)), '.png target="_blank"', '> <img height=50 width=50 src=figure/NCIPlots', seq(1:length(spiaRES$nci$Name)), '.png>', '</a>', sep='')
    )
    kable(nciTable, format='html', table.attr = 'id=\"nci_table\"')
  } else print('No sig. NCI pathways')  
}
nciTable(spiaRES$nci)

## ---- ReactomePlots ----
reactomeViz <- function(pw, geneList, anno) {
  require(ReactomePA)
  fc <- geneList
  names(fc) <- anno$gene_id[match(names(fc), anno$symbol)]
  my.viewPathway(pw, readable=TRUE, foldChange=fc, db='reactome')  
}

# create the reactome pw images
if(length(spiaRES$reactome$Name)>0) {
  sapply(spiaRES$reactome$Name, function(x) {
    reactomeViz(x, diffExpNBinom$FCdiff, anno)
  })
}

## ---- ReactomeTable ----
reactomeTable <- function(react.sig) {
  if(dim(react.sig)[1]!=0) {
    reactTable <- as.data.table(react.sig)
    reactTable <- cbind(reactTable, Image=NA)
    reactTable <- transform(reactTable, 
                            Image = paste('<a href = figure/ReactomePlots', seq(1:length(spiaRES$reactome$Name)), '.png target="_blank"', '> <img height=50 width=50 src=figure/ReactomePlots', seq(1:length(spiaRES$reactome$Name)), '.png>', '</a>', sep='')
    )
    kable(reactTable, format='html', table.attr = 'id=\"reactome_table\"')
  } else print('No sig. Reactome pathways')  
}
reactomeTable(spiaRES$reactome)

