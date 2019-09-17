#' Align features between batches
#'
#' Aligns features with m/z or retention time shift between batches using alignBatches from batchCorr package
#' by Carl Brunius. See more details in the help file and the original paper.
#'
#' @param object_na a MetaboSet object with missing values as NA
#' @param object_fill a similar MetaboSet object with imputed values
#' (used to compute distances between features, can contain missing values as well)
#' @param batch character, column name of pData with batch labels
#' @param mz,rt column names of m/z and retention time columns in fData
#' @param mzdiff,rtdiff the windows for m/z and retention time for aligning features
#' @param plot_folder path to the location where the plots should be saved, if NULL, no plots are saved
#'
#' @return a MetaboSet object with the aligned features
align_batches_brunius <- function(object_na, object_fill, batch, mz, rt, mzdiff, rtdiff, plot_folder = NULL) {

  # Set working directory for plotting (the bathCorr functions saves plots in the current working directory...)
  if (!is.null(plot_folder)) {
    old_wd <- getwd()
    setwd(plot_folder)
    report <- TRUE
  } else {
    report <- FALSE
  }


  # Extract peak mz and rt information
  pInfo <- as.matrix(fData(object_na)[, c(mz, rt)])
  colnames(pInfo) <- c("mz", "rt")

  # Align batches based on the QCs
  aligned <- batchCorr::alignBatches(peakInfo = pInfo, PeakTabNoFill = t(exprs(object_na)), PeakTabFilled = t(exprs(object_fill)),
                          batches = pData(object_na)[, batch], sampleGroups = object_na$QC, selectGroup = "QC",
                          mzdiff = mzdiff, rtdiff = rtdiff, report = report)

  # Reset working directory
  if (!is.null(plot_folder)) {
    setwd(old_wd)
  }

  # Attach aligned features
  exprs(object_fill) <- t(aligned$PTalign)
  object_fill
}


#' Normalize batches
#'
#' Normalize bathces by either reference samples of population median.
#' Uses normalizeBatches function from batchCorr package by Carl Brunius.
#'
#' @param object a MetaboSet object
#' @param batch,group character, column names of pData with batch labels and group labels
#' @param ref_label the label of the reference group i.e. the group that is constant through batches
#'
#' @return list, the object with normalized features and information on which features were corrected by ref samples in each batch.
normalize_batches_brunius <- function(object, batch, group = group_col(object), ref_label, ...) {

  normData <- batchCorr::normalizeBatches(peakTable = t(exors(object)), batches = pData(object)[, batch],
                               sampleGroup = pData(object)[, group], refGroup = ref_label, ...)

  exprs(object) <- normData$peakTable
  return(list(object = object, ref_corrected = normData$refCorrected))
}







