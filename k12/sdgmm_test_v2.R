setwd('D:/code/vitek_sdgmm_preprint/data')

library(Cardinal)
library(stringr)
library(fuzzyjoin)
library(dplyr)
library(ggplot2)

register(SnowParam(workers=(detectCores()-1), progressbar=TRUE))

edit_metadata <- function(dataset, spot_file, roi_list) {
  spots <- as.data.frame(read.table(spot_file, sep=" ", col.names=c("x", "y", "spot", "region")))
  spots$x <- as.integer(str_split(str_split(spots$spot, "X", simplify=TRUE)[,2], "Y", simplify=TRUE)[,1])
  spots$y <- as.integer(str_split(spots$spot, "Y", simplify=TRUE)[,2])
  spots$region <- as.character(spots$region)
  for (roi in roi_list) {
    spots$region[spots$region == roi[1]] <- roi[1]
    spots$region2[spots$region == roi[1]] <- roi[2]
  }
  pd_df <- as.data.frame(pixelData(dataset))
  pd_df <- difference_inner_join(pd_df, spots, by=c('x','y'), max_dist=0)
  pixelData(dataset)$region <- as.factor(pd_df$region)
  pixelData(dataset)$region2 <- as.factor(pd_df$region2)
  run(dataset) <- as.factor(paste(run(dataset), pixelData(dataset)$region, sep='_'))
  return(dataset)
}

preprocess_data <- function(dataset) {
  preprocessed <- normalize(dataset, method='tic')
  preprocessed <- smoothSignal(preprocessed, method='sgolay')
  preprocessed <- reduceBaseline(preprocessed, method='median')
  preprocessed <- process(preprocessed)
  peaks <- peakPick(dataset, method='simple', SNR=10)
  peaks <- peakAlign(peaks, tolerance=0.2, units='mz')
  peaks <- peakFilter(peaks, freq.min=0.5)
  peaks <- process(peaks)
  return(list(preprocessed, peaks, peaks))
}

generate_ref_peak_list <- function(list_of_datasets) {
  peaks <- list()
  for (i in 1:length(list_of_datasets)) {
    peaks[[i]] <- data.frame('peaks'=mz(list_of_datasets[[i]][['peaks']]))
  }
  df <- data.frame('peaks'=c())
  for (i in 1:length(peaks)) {
    df <- rbind(df, peaks[[i]])
  }
  df$peaks <- sort(df$peaks)
  df2 <- difference_inner_join(df, df, by='peaks', max_dist=0.2)
  colnames(df2) <- c('peaks', 'average')
  df3 <- aggregate(df2, by=list(df2$peaks), FUN=mean)
  return(sort(unique(df3$average)))
}

peak_process_reference <- function(list_of_datasets, ref_peaks) {
  for (i in 1:length(list_of_datasets)) {
    peaks <- peakPick(list_of_datasets[[i]][['preprocessed']], method='simple', SNR=6)
    peaks <- peakAlign(peaks, tolerance=0.2, units='mz', ref=ref_peaks)
    list_of_datasets[[i]][['processed_reference']] <- process(peaks)
  }
  return(list_of_datasets)
}

check_bad_features <- function(dataset, r=1, k=2) {
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
  return(bad_features_vector)
}

ims1 <- readMSIData('k12_lb.imzML')
ims2 <- readMSIData('mini_ims.imzML')

spots1 <- 'k12_lb.txt'
spots2 <- 'mini_ims.txt'

rois1 <- list(c('lb_1', 'lb'),
              c('lb_2', 'lb'),
              c('k12_1', 'k12'),
              c('k12_2', 'k12'))

rois2 <- list(c('lb_1', 'lb'),
              c('lb_2', 'lb'),
              c('k12_1', 'k12'),
              c('k12_2', 'k12'))

unproc1 <- edit_metadata(ims1, spots1, rois1)
unproc2 <- edit_metadata(ims2, spots2, rois2)

tmp1 <- preprocess_data(unproc1)
preproc1 <- tmp1[[1]]
peaks1 <- tmp1[[2]]
proc1 <- tmp1[[3]]
tmp2 <- preprocess_data(unproc2)
preproc2 <- tmp2[[1]]
peaks2 <- tmp2[[2]]
proc2 <- tmp2[[3]]

data1 <- list('unprocessed'=unproc1, 'preprocessed'=preproc1, 'peaks'=peaks1, 'processed'=proc1)
data2 <- list('unprocessed'=unproc2, 'preprocessed'=preproc2, 'peaks'=peaks2, 'processed'=proc2)

save.image("C:/Users/gordon/Downloads/sdgmm/sdgmm_test_v2.RData")

fuzzy_join <- difference_inner_join(data.frame('peaks'=mz(data1[['peaks']])), data.frame('peaks'=mz(data2[['peaks']])), by='peaks', max_dist=0.2)
fuzzy_join$mean <- rowMeans(fuzzy_join, na.rm=TRUE)
ref_peaks <- unique(fuzzy_join$mean)

ims_data_list <- peak_process_reference(list(data1, data2), ref_peaks)

save.image("C:/Users/gordon/Downloads/sdgmm/sdgmm_test_v2.RData")

combined <- Cardinal::combine(ims_data_list[[1]][['processed_reference']],
                              ims_data_list[[2]][['processed_reference']])
list_of_datasets <- list(ims_data_list[[1]][['processed_reference']], ims_data_list[[2]][['processed_reference']])

register(SnowParam(workers=(detectCores()-2), progressbar=TRUE, log=TRUE, stop.on.error=FALSE, logdir=getwd()))
set.seed(233)
sdgmm <- spatialDGMM(combined, r=1, k=2, groups=run(combined))
bad_features_vector <- check_bad_features(combined, r=1, k=2)
bad_features_vector <- c(11, 12)

if (is.na(bad_features_vector) != TRUE) {
  bad_mz_vector <- mz(combined)
  for (bad_feature in sort(bad_features_vector, decreasing=TRUE)) {
    combined <- combined[features(combined, mz < bad_mz_vector[bad_feature] - 0.0001 | mz > bad_mz_vector[bad_feature] + 0.0001),]
  }
}

register(SnowParam(workers=(detectCores()-2), progressbar=TRUE))
set.seed(233)
sdgmm <- spatialDGMM(combined, r=1, k=2, groups=run(combined))

segtest <- segmentationTest(sdgmm, ~ region2, clasControl='Ymax')
summary(segtest)
topFeatures(segtest, n=32)

save.image("C:/Users/gordon/Downloads/sdgmm/sdgmm_test_v2.RData")
