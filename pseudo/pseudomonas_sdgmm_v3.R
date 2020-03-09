setwd("G:/vitek_sdgmm_preprint")
library(Cardinal)
library(stringr)

################################################################################
# functions                                                                    #
################################################################################

# load in .imzML file and "Spot List" text file exported from Bruker flexImaging
load_datasets <- function(vector_of_datasets) {
  listo <- list()
  count <- 1
  for (i in vector_of_datasets) {
    imzml <- readMSIData(paste0(i, ".imzML"), mass.range=c(1,2000),
                         resolution=0.2, units="mz")
    spots <- as.data.frame(read.table(paste0(i, ".txt"), sep=" ",
                                      col.names=c("x", "y", "spot", "region")))
    listo[[count]] <- list(imzml, spots)
    count <- count + 1
  }
  return(listo)
}

# correct x and y coordinates based on "Spot List" text file
# update run names and add column for "condition" in object metadata
# level_list param == list of vectors containing region name and condition
process_metadata <- function(dataset, spot_table, level_list) {
  spot_table$x <- as.integer(str_split(str_split(spot_table$spot, "X",
                                                 simplify=TRUE)[,2], "Y",
                                       simplify=TRUE)[,1])
  spot_table$y <- as.integer(str_split(spot_table$spot, "Y", simplify=TRUE)[,2])
  for (i in level_list) {
    tmp_table <- spot_table[which(spot_table$region==i[1]),]
    tmp_table$condition <- factor(i[2])
    if (exists("new_spot_table")) {
      new_spot_table <- rbind(new_spot_table, tmp_table)
    } else {
      new_spot_table <- tmp_table
    }
  }
  new_spot_table$region <- as.factor(new_spot_table$region)
  new_spot_table$run <- paste(run(dataset)[1], new_spot_table$region, sep="_")
  pd_tmp <- merge(new_spot_table, coord(pixelData(dataset)), by=c("x", "y"))
  pixelData(dataset) <- PositionDataFrame(coord=pd_tmp[,1:2], run=pd_tmp$run,
                                          region=pd_tmp$region,
                                          condition=pd_tmp$condition)
  return(dataset)
}

# preprocess dataset using the following settings
# note: both "tic" and "rms" normalization were attempted
preprocess_datasets <- function(dataset) {
  preprocessed <- normalize(dataset, method="tic")
  preprocessed <- smoothSignal(preprocessed, method="sgolay")
  preprocessed <- reduceBaseline(preprocessed, method="median")
  preprocessed <- process(preprocessed)
  peaks <- peakPick(preprocessed, method="simple")
  peaks <- peakAlign(peaks, tolerance=0.2, units="mz")
  peaks <- peakFilter(peaks, freq.min=0.05)
  peaks <- process(peaks)
  bin <- peakBin(preprocessed, ref=mz(peaks), type="height")
  processed <- process(bin)
  return(processed)
}

################################################################################
# main()                                                                       #
################################################################################
# load dataset
dataset_files <- c("G:/vitek_sdgmm_preprint/PA14_TLCA_2 DAYS_2016-06-09")
datasets <- load_datasets(dataset_files)

datasets[[1]][[3]] <- list(c("PA14", "meoh"),
                           c("PA14_250uM_TLCA", "tlca"),
                           c("No_TLCA_Control", "meoh"),
						   c("250uM_TLCA_Control", "tlca"))

# edit metadata
pseudo_1 <- process_metadata(datasets[[1]][[1]], datasets[[1]][[2]],
                             datasets[[1]][[3]])

# preprocess data
pseudo_1_proc <- preprocess_datasets(pseudo_1)

# run sdgmm
# note: k=4 and k=5 attempts required significant amounts of time and RAM
pseudo_1_sdgmm <- spatialDGMM(pseudo_1_proc, r=1, k=3, method="adaptive",
                              init="kmeans", iter.max=1000, tol=1e-11)
summary(pseudo_1_sdgmm)

# run segmentation test
pseudo_1_segtest <- segmentationTest(pseudo_1_sdgmm, ~ condition,
                                     classControl="Ymax")
summary(pseudo_1_segtest)
topFeatures(pseudo_1_segtest, p.adjust="fdr", AdjP<0.05, n=1000)

# save segmentation test top features to dataframe and convert to vector
pseudo_1_segtest_features <- as.data.frame(topFeatures(pseudo_1_segtest,
                                           p.adjust="fdr", AdjP<0.05, n=1000))
pseudo_1_features_vec <- as.vector(pseudo_1_segtest_features$feature)
# view ion images for all top features IDed by segmentation test
image(pseudo_1_segtest, model=list(feature=pseudo_1_features_vec),
      values="mapping", layout=c(ceiling(length(pseudo_1_features_vec)/4), 4))

# view ion images from processed dataset for all top features IDed by
# segmentation test
pseudo_1_mz_vec <- as.vector(pseudo_1_segtest_features$mz)
image(pseudo_1_proc, mz=pseudo_1_mz_vec,
      layout=c(ceiling(length(pseudo_1_mz_vec)/4), 4),
	  contrast.enhance="histogram", smooth.image="gaussian")

################################################################################
# main2()                                                                      #
################################################################################
# subset for region: Pseudomonas aeruginosa w/ TLCA media supplementation
pseudo_1_roi <- pseudo_1_proc[,pixels(pseudo_1_proc,
                run=="PA14_TLCA_2 DAYS_2016-06-09_PA14_250uM_TLCA")]

# run sdgmm
pseudo_1_roi_sdgmm <- spatialDGMM(pseudo_1_roi, r=1, k=2, method="adaptive",
                                  init="kmeans", iter.max=1000,tol=1e-11)
summary(pseudo_1_roi_sdgmm)

# visualize various features
image(pseudo_1_roi_sdgmm, model=list(feature=c(302,303,304,320,321,322)),
      contrast.enhance="histogram", smooth.image="gaussian")
image(pseudo_1_roi_sdgmm, contrast.enhance="histogram",
      smooth.image="gaussian", layout=c(165,3))
