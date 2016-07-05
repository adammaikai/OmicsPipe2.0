## Functions taken from the suppl. material from
## A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis
## http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.long

require(limma)
require(edgeR)
require(DESeq)
require(RColorBrewer)
require(cluster)

## -------------------------------------------------------
## getSizeFactor
## Calculate size factor for the following approaches :
## - Reads count normalization
## - Upper Quartile normalization
## - Median normalization
## - DESeq normalization
## - edgeR (TMM) normalization
##
## group : describes the condition and is only needed for DESeq package
##
## return : a list of size factor
## --------------------------------------------------------
getSizeFactor <- function(geneCount, group){
  
  ##Reads counts normalization
  N <- colSums(geneCount)
  TC <- N/mean(N)
  ##Upper Quartile normalization
  UQ <- apply(geneCount, 2,  quantile, 0.75)
  UQ <- UQ/mean(UQ)
  ##Median normalization
  Med <- apply(geneCount, 2, median)
  Med <- Med/mean(Med)
  ## DESeq normalization
  cds <- newCountDataSet(geneCount, group)
  cds <- estimateSizeFactors(cds)
  deseq <- sizeFactors(cds)
  ## edgeR normalization
  f <- calcNormFactors(as.matrix(geneCount),method="TMM")
  f.scale <- N*f / mean(N*f)
  
  return(list(TC=TC,UQ=UQ,Med=Med,DESeq=deseq,TMM=f.scale))
}


rpkm <- function(counts, trsLength){
  rpkmCounts <- apply(counts,2,function(s){
    lengthsInKb <- trsLength/1000
    millionMapped <- sum(s)/1e+06
    rpm <- s/millionMapped
    rpkm <- rpm/lengthsInKb
    rpkm
  })
  rpkmCounts
}


## -------------------------------------------------------
## getNormData
## Calculate normalized data from raw count data
##
## group : describes the condition
## info : describes the features
##
## return : a list of normalized factor
## --------------------------------------------------------
getNormData <- function(geneCount,group,trsLength=NA, addRaw=TRUE){
  sf <- getSizeFactor(geneCount, group)
  
  ##When size factor is available
  geneCountNorm <- lapply(sf,function(x){
    scale(geneCount, center=FALSE, scale=x)
  })
  names(geneCountNorm) <- names(sf)
  
  ##Other cases
  geneCountNorm[["FQ"]]<- limma::normalizeQuantiles(geneCount)
  
  ##rpkm
  if (!is.na(trsLength) && length(trsLength)>0){
    geneCountNorm[["RPKM"]]<- rpkm(geneCount, trsLength)
  }
  
  if (addRaw)
    geneCountNorm[["RawCount"]]<- geneCount
  
  return(geneCountNorm)
}

## -------------------------------------------------------
## diffAnDESseq
## Differential analysis using DESeq
##
## --------------------------------------------------------


diffAnDESseq <- function(geneCount, group, trsLength=NA, round=FALSE, condA=NA, condB=NA){
  sf <- getSizeFactor(geneCount, group)
  normData <- getNormData(geneCount, group, trsLength)
  
  if (is.na(condA)) condA <- levels(group)[1]
  if (is.na(condB)) condB <- levels(group)[2]
  
  if (!round){
    ## if size factor available
    res1 <- lapply(sf, function(sf){
      cds <- newCountDataSet(geneCount, group)
      pData(cds)$sizeFactor <- sf
      cds <- estimateDispersions( cds, method="pooled" )
      nbinomTest(cds,condA=condA,condB=condB)
    })
    
    ## DESeq differential expression when size factor are not available -> round function
    res2 <- lapply(normData[setdiff(names(normData), names(sf))], function(normCount){
      cds <- newCountDataSet(round(normCount), group)
      pData(cds)$sizeFactor <- rep(1,ncol(normCount))
      cds <- estimateDispersions( cds, method="pooled" )
      nbinomTest(cds,condA=condA, condB=condB)
    })
    
    res <- c(res1, res2)
  }else{
    res <- lapply(normData, function(normCount){ ##+rpkm
      cds <- newCountDataSet(round(normCount), group)
      pData(cds)$sizeFactor <- rep(1,ncol(normCount))
      cds <- estimateDispersions( cds, method="pooled" )  
      rep(1,ncol(normCount))
      nbinomTest(cds,condA=condA, condB=condB)
    })
  }
}

## -------------------------------------------------------
## normBoxplot
## --------------------------------------------------------
normBoxplot <- function(norm){
  norm.df <- as.data.frame(sapply(norm, rbind))
  sel.col <- brewer.pal(length(norm), "Set2")
  par(font.lab=2, mar=c(8,4,1,1))
  boxplot(log2(1+norm.df), 
          outline=FALSE, pch=16, las=2,frame=FALSE, 
          ylab="Log2 counts",
          col=as.vector(sapply(sel.col, rep.int, times=ncol(norm.df)/length(norm))))
}
## -------------------------------------------------------
## varIntra
## --------------------------------------------------------
varIntra <- function(norm, group){
  var.intra <- list()
  groups <- unique(group)
  N <- vector("numeric", length(groups))
  par(font.lab=2, font.axis=1, mar=c(4,4,1,1), mfrow=c(length(groups),1))
  for (g in 1:length(groups)){
    var.intra[[g]] <- list()
    for (j in 1:length(norm)){
      tab <- norm[[j]][,group == groups[g]]
      tab <- tab[rowSums(tab) > 0,]
      N[g] <- max(N[g], nrow(tab))
      var.intra[[g]][[j]] <- apply(tab, 1, sd)/rowMeans(tab)
    }
    names(var.intra[[g]]) <- names(norm)
    
    mat <- matrix(NA, nrow=N[g], ncol=length(norm))
    for (i in 1:length(var.intra[[g]])) mat[1:length(var.intra[[g]][[i]]),i] <- var.intra[[g]][[i]]
    boxplot(mat, names=names(var.intra[[g]]), main=paste("Intra-group Variance - ",groups[g], sep=""), frame=FALSE, cex.main=.7, cex.axis=.7, col="gray90", las=2)
  }
  names(var.intra) <- groups
  par(mfrow=c(1,1))
}
## -------------------------------------------------------
## hclustDE
## --------------------------------------------------------

hclustDE <- function(DEG){
  uDEG<-unique(unlist(DEG))
  DEGmat<-matrix(0, nrow=length(uDEG), ncol=length(DEG))
  rownames(DEGmat)<-uDEG
  colnames(DEGmat)<-names(DEG)
  for (i in 1:length(DEG)) DEGmat[DEG[[i]],i]<-1
  d <- dist(t(DEGmat),method="binary")
  par(font.lab=2, mar=c(4,4,1,1))
  plot(as.dendrogram(agnes(d,method="ward", diss=TRUE)),edgePar=list(lwd=2) ,xlab="",sub="",ylab="Ward-Jaccard Distances", cex.lab=.8)
}
## -------------------------------------------------------
## MAPlot
## --------------------------------------------------------

intersect2 <- function(...) {
  args <- list(...)
  emptyin <- unlist(lapply(args,is.null))
  if (any(emptyin))
    args <- args[which(!emptyin)]
  nargs <- length(args) 
  if(nargs <= 1) {
    if(nargs == 1 && is.list(args[[1]])) {
      do.call("intersect2", args[[1]])
    } else {
      stop("cannot evaluate intersection fewer than 2 arguments")
    }
  } else if(nargs == 2) {
    intersect(args[[1]], args[[2]])
  } else {
    intersect(args[[1]], intersect2(args[-1]))
  }
}

##NS - internal function
setdiff2 <- function(...) {
  args <- list(...)
  emptyin <- unlist(lapply(args,is.null))
  if (any(emptyin))
    args <- args[which(!emptyin)]
  nargs <- length(args)
  if(nargs <= 1) {
    if(nargs == 1 && is.list(args[[1]])) {
      do.call("setdiff2", args[[1]])
    } else {
      stop("cannot evaluate intersection fewer than 2 arguments")
    }
  } else if(nargs == 2) {
    setdiff(args[[1]], args[[2]])
  } else {
    setdiff(setdiff2(args[-length(args)]), args[[length(args)]])
  }
}


MAplotDE <- function(geneCount, group, condA, condB, DEG, showComm=FALSE, showSpe=TRUE){
  leg <- lfill <- character()
  
  dA=log2(1+apply(geneCount[,which(group==condA)],1,mean))
  dB=log2(1+apply(geneCount[,which(group==condB)],1,mean))
  M=dB-dA
  A=(dA+dB)/2
  
  par(font.lab=2, mar=c(5,4,1,1))
  plot(x=0, y=0, type="n", frame=FALSE, ylim=c(min(M),max(M)), xlim=c(0,max(A)), xlab="Mean Expression Level [A]", ylab="Expression Level Ratio [M]", cex.lab=.8, cex.axis=.7)
  if (showComm){
    comm <- intersect2(DEG)
    points(A[comm], M[comm], col="gray38", pch=20, cex=.7)
    leg <- "Common"
    lfill <- "gray38"
  }
  if (showSpe){
    sel.col <- brewer.pal(length(DEG), "Set2")
    for (i in names(DEG)){
      iind <- which(names(DEG)==i)
      spe <- setdiff2(DEG[c(iind,setdiff(1:length(DEG),iind))])
      points(A[spe], M[spe], col=sel.col[iind], pch=20, cex=.7)
    }
    leg <- c(leg, names(DEG))
    lfill <- c(lfill, sel.col)
  }
  
  if (showComm || showSpe) legend(x="topright",fill=lfill, legend=leg, cex=.7)
}