## ---- reliableGEP ----
# TO-DO:
# make tightCutoff & looseCutoff & background measure a configurable parameter in the yaml config 

reliableRNAExp <- function(x, reference, looseCutoff=0.02, tightCutoff=0.01, background='cero') {
  # Check cutoff ranges and relationship
  if ((tightCutoff > 1.0)|(tightCutoff < 0.0) | (looseCutoff > 1.0)|(looseCutoff < 0.0)) {
    stop("\nAborting: cutoffs must be between 0.0 and 1.0\n\n")
  }
  if (tightCutoff > looseCutoff){ 
    stop("\nAborting: tightCutoff must be lower than looseCutoff\n\n")
  }
  if (background!='median' & background!='cero') {
    stop('\nAborting: background measure must be either median or cero\n\n')
  }
  
  if(background=='median') {
    medianRef <- apply(reference, 1, median)
    background <- names(medianRef)[which(medianRef==0)]
  } 
  if(background=='cero') {
    ceroRef <- apply(reference, 1, sum)
    background <- names(ceroRef)[which(ceroRef==0)]
  }

  # Get all exprs, then get just the background subset
  allExprs <- x
  backExprs <- allExprs[background,]
  allExprs<- as.matrix(allExprs)
  backExprs<-as.matrix(backExprs)
  # Create cutoff function
  cutoff_fcn <- function(x) {
    if (x <= tightCutoff) return("P") 
    else if (x > looseCutoff) return("A")
    else return("M")
  }
  xBack <- backExprs[, 1]
  xBack <- sort(xBack, decreasing = TRUE) # sort the Background exprs, decreasing
  yBack <- seq(0,1,1/(length(xBack)-1)) # generate the y axis for the survivor fcn.
  interp <- approxfun(xBack, yBack, yleft=1, yright=0) # create interpolation fcn.
  Pvals <- sapply(allExprs[, 1], interp) # calculate column of P-values
  Pcalls <- sapply(Pvals, cutoff_fcn) # create P/A/M column
  # compute intensity values at the two cutoffs:
  myX <- backExprs
  myY <- seq(0,1,1/(length(myX[,1])-1))
  myX[,1] <- sort(myX[,1], decreasing = TRUE)
  revInterp <- approxfun(myY, myX[,1], yleft=1, yright=0)
  revTight <- revInterp(tightCutoff)
  revLoose <- revInterp(looseCutoff)
  # print out intensity values at the two cutoffs:  
  cat("\nIntensities at cutoff P-values of ", looseCutoff," and ", tightCutoff, ":\n")
  cat(colnames(allExprs), "\t", format.pval(revLoose,digits=3), "\t\t", format.pval(revTight,digits=3), "\n")
  cat("\n")
  return(list(Pcalls=Pcalls, Pvals=Pvals, revTight=revTight, revLoose=revLoose, background=background))  
}
patientREC <- reliableRNAExp(patientCountsNorm, reference)

## ---- reliableGEPDensityPlot ----
densPlot <- function(patientCountsNormEntrez, backgroundGenes, looseCutoff=0.02, tightCutoff=0.01) {
  # create x and y ranges for background intensities
  xBack <- sort(patientCountsNormEntrez[backgroundGenes,1], decreasing=TRUE)
  yBack <- seq(0,1,1/(length(xBack)-1))
  # Plot expression densities of all probesets, then just NSMPs
  plot(density(patientCountsNormEntrez[,1]),
       col="blue",
       xlim = c(1,18),
       ylim = c(0,0.25),
       main = "Expression density: background genes vs. all, and background survivor curve",
       xlab = "Log2(Intensity)",
       ylab = "Probability density")
  lines(density(xBack), col=6)
  # interpolate over the background exprs to draw survivor function
  interp <- approxfun(xBack, yBack, yleft=1, yright=0) 
  x = xBack
  curve(interp(x),add=TRUE, lwd=2)	 #add it to the plot
  # reverse interpolate to get intensity values at p-value cutoffs
  revInterp <- approxfun(yBack, xBack, yleft=1, yright=0)
  rev01=revInterp(tightCutoff)
  rev02=revInterp(looseCutoff)
  # Pinpoint the x-y locations
  points(rev01, tightCutoff, pch=21, cex=2, lwd=2,col=1)
  points(rev02, looseCutoff, pch=21, cex=2, lwd=2,col=1)
  # Draw horiz. lines & labels for both Pval cutoffs: tightCutoff, looseCutoff:
  abline(h = tightCutoff, col = 1, lty = 2)
  abline(h = looseCutoff, col = 1, lty = 2)
  text(2.4, tightCutoff, pos=3, offset=0.1, cex=.9, as.character(tightCutoff))
  text(1.7, looseCutoff, pos=3, offset=0.1, cex=.9, as.character(looseCutoff))
  # vertical lines & labels for interpolated intensities at cutoffs
  revTight=revInterp(tightCutoff)
  revLoose=revInterp(looseCutoff)
  abline(v = revTight, col = 1, lty = 2)
  abline(v = revLoose, col = 1, lty = 2)
  text(revLoose, .25, pos=2, offset=0.1, cex=.9, format.pval(revLoose,digits=3))
  text(revTight, .20, pos=2, offset=0.1, cex=.9, format.pval(revTight,digits=3))
  text(revLoose, .25, pos=4, cex=.8, "Log(intensity)")
  text(revTight, .20, pos=4, cex=.8, "Log(intensity)")
  lines(density(patientCountsNormEntrez[,1][patientREC$Pcalls=="P"], bw=.1, n=512),col=2, lty=2, lwd=1)
  lines(density(patientCountsNormEntrez[,1][patientREC$Pcalls=="A"], bw=.1, n=512),col=3, lty=2, lwd=1)
  legend(11,.25, c("background genes exprs, survivor fcn.","background genes exprs, density","all genes, density","'reliable expressed' genes, density","'not so reliable expressed' genes, density"),
         col = c(1,6,4,2,3), lty=c(1,1,1,2,2), lwd=c(2,1,1,1,1), cex=.75,
         text.col= "darkgreen",
         bg='gray90')
}
# plot calling of reliable expressed genes
densPlot(patientCountsNorm, patientREC$background)

## ---- reliableGEPTable ----
relTable <- data.frame(Expression=c('not reliable expressed', 'marginal expressed', 'reliable expressed'),
                       Freq=as.vector(table(patientREC$Pcalls)))
relTable <- relTable[c(3,2,1),]

## ---- reliableGEPPlot ----
relTablePlot <- rCharts:::Highcharts$new()
relTablePlot$chart(type = "column")
relTablePlot$title(text = "Reliable Expressed Genes")
relTablePlot$xAxis(categories = relTable$Expression)
relTablePlot$yAxis(title = list(text = "number of genes"))
relTablePlot$data(relTable)
relTablePlot$params$width = 800
relTablePlot$params$height = 400
relTablePlot

