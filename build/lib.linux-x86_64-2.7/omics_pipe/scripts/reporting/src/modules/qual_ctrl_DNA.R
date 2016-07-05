## ---- fastQC ----
#fastqcReportP1 <- paste(config$QC_PATH, '/', patientID, '_1/', patientID, '_1_fastqc/fastqc_data.txt',  sep="")
#fastqcReportP2 <- paste(config$QC_PATH, '/', patientID, '_2/', patientID, '_2_fastqc/fastqc_data.txt',  sep="")
fastqcReportP1 <- paste(config$QC_PATH, '/', patientID, '/fastqc/fastqc_data.txt',  sep="")
fastqcReportP2 <- paste(config$QC_PATH, '/', patientID, '/fastqc/fastqc_data.txt',  sep="")

files <- c(fastqcReportP1,fastqcReportP2)
fileName <- vector()
basicStats <- vector()
fileType <- vector()
encoding <- vector()
totalSeq <- numeric()
filteredSeq <- numeric()
seqLength <- numeric()
percentGC <- numeric()
url <- vector()
for(i in 1:length(files)) {
  xx <- read.delim(files[i], header=F, stringsAsFactors=FALSE)$V2[1:10]
  fileName[i] <- as.vector(xx[4])
  basicStats[i] <- as.vector(xx[2])
  fileType[i] <- as.vector(xx[5])
  encoding[i] <- as.vector(xx[6] )
  totalSeq[i] <- as.numeric(xx[7])
  filteredSeq[i] <- as.numeric(xx[8])
  seqLength[i] <- as.numeric(xx[9])
  percentGC[i] <- as.numeric(xx[10])
}
fastqcDf <- data.frame(file=fileName,
                       basic_statistics=basicStats,
                       file_type=fileType,
                       encoding=encoding,
                       total_seq=totalSeq,
                       filtered_seq=filteredSeq,
                       seq_length=seqLength,
                       percent_gc=percentGC,
#                       link_to_report=c(paste(config$QC_PATH, '/', patientID, "_1_fastqc/fastqc_report.html",  sep=""), 
#                                        paste(config$QC_PATH, '/', patientID, "_2_fastqc/fastqc_report.html",  sep="")
                       link_to_report=c(paste(config$QC_PATH, '/', patientID, "/fastqc/fastqc_report.html",  sep=""), 
                                        paste(config$QC_PATH, '/', patientID, "/fastqc/fastqc_report.html",  sep="")
 
                                        )
                      )
fastqcDf <- transform(fastqcDf, link_to_report = paste('<a href = ', shQuote(link_to_report), 'target="_blank" > FastQC report </a>'))
fastqcDf <- t(fastqcDf)
fastqcDf <- rbind(test=c('Pair_1', 'Pair_2'), fastqcDf)
rownames(fastqcDf)[1] <- ''
kable(fastqcDf, format='html', table.attr = 'id=\"fastqc_table\", class=\"stdTable\"')

## ---- starQCTable ----
starlogpath <- paste(config$STAR_RESULTS, '/', patientID, "/Log.final.out",  sep="")
starlog <- read.delim(starlogpath, header=F, stringsAsFactors=FALSE)
starlog$V1 <- gsub('^[ ]+', '', starlog$V1)
starlog$V1 <- gsub('\\|', '', starlog$V1)

starDf <- data.frame(Metric=starlog$V1,
                     Value=starlog$V2,
                     stringsAsFactors=FALSE
                     )
starDf <- rbind(c('Metric', 'Value'), starDf)
colnames(starDf) <- NULL # remove colnames as this works better with the current .css layout, thead is not defined in this yet

kable(starDf, format='html', table.attr = 'id=\"star_table\", class=\"stdTable\"')

## ---- starQCPie ----
starDfPie <- data.frame(stats=c('unique mapped reads', 
                                'multi mapped reads',
                                'too many loci',
                                'unmapped reads'),
                        Freq=c(as.numeric(gsub('%', '', starlog$V2[9])),
                               as.numeric(gsub('%', '', starlog$V2[24])),
                               as.numeric(gsub('%', '', starlog$V2[26])),
                               as.numeric(gsub('%', '', starlog$V2[29]))
                               )
                        )
starPie <- nPlot(x = "stats", y = "Freq", data = starDfPie, type = "pieChart")
starPie$params$width <- '600'
starPie$params$height <- '250'
starPie

## ---- htseqQCTable ----
htseqCountDf <- data.frame(Metric=names(htseqStats),
                           Value=as.numeric(htseqStats),
                           stringsAsFactors=FALSE
                           )
htseqCountDf <- rbind(c('Metric', 'Value'), htseqCountDf)
colnames(htseqCountDf) <- NULL
kable(htseqCountDf, format='html', table.attr = 'id=\"htseq_table\", class=\"stdTable\"')

## ---- htseqQCPie ----
htseqCountDfPie <- data.frame(stats=as.matrix(htseqCountDf)[c(2:6),1],
                              Freq=as.matrix(htseqCountDf)[c(2:6),2])
htseqPie <- nPlot(x = "stats", y = "Freq", data = htseqCountDfPie, type = "pieChart")
htseqPie$params$width <- '600'
htseqPie$params$height <- '250'
htseqPie

## ---- RNASeqMetricsTable ----
rnaSeqMetrics <- function(file) {
  dat <- read.csv2(file, sep='\t', skip=6, stringsAsFactor=FALSE)[1,]
  rnaSeqMetricsTable <- data.frame(Metric=names(dat),
                                   Value=as.numeric(dat),
                                   stringsAsFactors=FALSE
                                   )
  rnaSeqMetricsTable$Value <- format(rnaSeqMetricsTable$Value, scientific=FALSE, digits=2)
  rnaSeqMetricsTable$Value <- gsub('\\.00', '', rnaSeqMetricsTable$Value)
  return(rnaSeqMetricsTable)
}
rnaSeqMetricsTable <- rnaSeqMetrics(paste(config$QC_PATH, '/', patientID, '/rnaseqmetrics', sep=''))
rnaSeqMetricsTable <- rbind(c('Metric', 'Value'), rnaSeqMetricsTable)
colnames(rnaSeqMetricsTable) <- NULL
kable(rnaSeqMetricsTable, format='html', table.attr = 'id=\"rnaseqmetrics_table\", class=\"stdTable\"')

## ---- RNASeqMetricsPie ----
rnaSeqMetricsDfPie <- rnaSeqMetricsTable[3:7,]
colnames(rnaSeqMetricsDfPie) <- c('stats', 'Freq')
rnaSeqMetricsDfPie$Freq <- as.numeric(rnaSeqMetricsDfPie$Freq)
rnaSeqMetricsDfPie$stats <- as.vector(rnaSeqMetricsDfPie$stats)
rnaSeqMetricsPie <- nPlot(x = "stats", y = "Freq", data = rnaSeqMetricsDfPie, type = "pieChart")
rnaSeqMetricsPie$params$width <- '600'
rnaSeqMetricsPie$params$height <- '250'
rnaSeqMetricsPie

## ---- InsertSizeTable ----
insertData <- read.csv2(paste(config$QC_PATH, '/', patientID, '/insertSize.txt', sep=''), sep='\t', skip=6, stringsAsFactor=FALSE)
insertData <- insertData[1,1:18]
insertSizeDf <- data.frame(Metric=colnames(insertData),
                           Value=as.vector(as.matrix(insertData)),
                           stringsAsFactors=FALSE)
insertSizeDf <- rbind(c('Metric', 'Value'), insertSizeDf)
colnames(insertSizeDf) <- NULL
kable(insertSizeDf, format='html', table.attr = 'id=\"insertSize_table\", class=\"stdTable\"')

## ---- InsertSizeImage ----
# convert the pdf histogram to .png
system(paste('convert', 
             paste(config$QC_PATH, '/', patientID, '/insertSizeHist.pdf', sep=''), 
             paste(config$QC_PATH, '/', patientID, '/insertSizeHist.png', sep='')
             ),
       wait=TRUE
       )
cat(paste('<img src=', shQuote(paste(config$QC_PATH, '/', patientID, '/insertSizeHist.png', sep='')), '>', sep=''))


