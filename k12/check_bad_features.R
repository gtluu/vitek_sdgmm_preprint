register(SnowParam(workers=(detectCores()-2), progressbar=TRUE, log=TRUE, stop.on.error=FALSE, logdir=getwd()))
set.seed(233)
sdgmm <- spatialDGMM(combined, r=1, k=2, init="kmeans", iter.max=100, tol=1e11, groups=run(combined))

check_bad_features <- function(dataset, r, k, compare) {
  bad_features_vector <- c()
  for (i in 1:(detectCores()-2)) {
    filename <- paste0(getwd(), "/BPJOB.task", as.character(i), ".log")
    logfile <- readChar(filename, file.info(filename)$size)
    logstring <- str_split(logfile,
                           '############### LOG OUTPUT ###############')[[1]]
    logstring <- logstring[length(logstring)]
    logstring <- str_split(logstring, 'feature = ')[[1]]
    logstring <- logstring[grepl('ERROR', logstring)]
    for (errorline in logstring) {
      bad_features_vector <- c(bad_features_vector,
                               as.numeric(str_split(errorline, ' ')[[1]][1]))
    }
  }
  if (length(bad_features_vector) != 0) {
    for (bad_feature in sort(bad_features_vector, decreasing=TRUE)) {
      dataset <- MSImagingExperiment(imageData=rbind(imageData(dataset)[1:(bad_feature-1),],
                                                     imageData(dataset)[(bad_feature+1):length(mz(dataset)),]),
                                     featureData=rbind(featureData(dataset)[1:(bad_feature-1),],
                                                       featureData(dataset)[(bad_feature+1):length(mz(dataset)),]),
                                     pixelData=pixelData(dataset),
                                     metadata=metadata(dataset),
                                     processing=processingData(dataset),
                                     centroided=TRUE)
    }
    register(SnowParam(workers=(detectCores()-2), progressbar=TRUE))
    if (compare == FALSE) {
      return(spatialDGMM(dataset, r=r, k=k, init="kmeans", iter.max=100, tol=1e11))
    } else {
      return(spatialDGMM(dataset, r=r, k=k, init="kmeans", iter.max=100, tol=1e11, groups=pixelData(dataset)$sample))
    }
  }
}

sdgmm <- check_bad_features(combined, r=1, k=2, compare=TRUE)