## ---- setup ----
library(knitr)
library(knitrBootstrap) # https://github.com/jimhester/knitrBootstrap
library(rCharts) # https://github.com/ramnathv/rCharts/
library(DESeq)
library(org.Hs.eg.db)
library(data.table)
library(plyr)
library(dplyr)
library(gdata)
library(annotate)
library(XML)
library(xtable)
library(pander)
library(KEGGREST)
library(knitcitations)
library(reshape)# https://github.com/cboettig/knitcitations

source(paste(config$R_SOURCE_PATH, '/src/core/functions.R', sep=''))

#cite_options(tooltip=TRUE) # use bootstrap to create mouse over citations