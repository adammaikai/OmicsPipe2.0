## ---- DPS ----
dps <- function(diff) {
  # export genes with FC > 1 to txt
  up <- names(diff$FCdiff)[which(diff$FCdiff > 1)]
  down <- names(diff$FCdiff)[which(diff$FCdiff < -1)]
  write.table(up, paste(config$REPORT_RESULTS, '/', patientID, '/up.txt', sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(down, paste(config$REPORT_RESULTS, '/', patientID, '/down.txt', sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  if(exists('/opt/applications/dps/1.3.1111/DPS.jar')) {
    system('cp /opt/applications/dps/1.3.1111/DPS.jar .') # copy DPS into working dir (only seems to run within the dir that holds diff exprssed genes)
    system('java -jar DPS.jar down.txt up.txt drugs') # run DPS (-r option does not work, so swtich input files)
    Sys.sleep(60) # give DPS some time to create the output..
    system('rm DPS.jar') # delete DPS again
    
    drugTable <- data.table(read.csv2(paste(config$REPORT_RESULTS, '/', patientID, '/drugs.csv', sep=''), sep=',', header=TRUE))
    
    dTable(drugTable, sPaginationType='full_numbers', iDisplayLength=25, sScrollX='100%', sScrollXInner='100%', width='100%')
  } else {
    print('DPS could not be found on the system')
  }
}

dps(diffExpNBinom)