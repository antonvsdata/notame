#' Batch correction
#'
#' "Basic" batch correction by median? from BatchCorrMetabolomics::doBC
#'
#' @param object a MetaboSet object
#' @param batch the column name for batch labels
#' @param ref the column name for reference sample labels
#' @param ref_label the label for reference samples
#' @param ... other parameters pased to doBC
#'
#' @return a MaetaboSet object with the corrected abundances
#'
#' @export
dobc <- function(object, batch, ref, ref_label, ...) {
  ref_idx <- which(pData(object)[, ref] == ref_label)
  seq_idx <- object$Injection_order
  batch_idx <- pData(object)[, batch]

  batch_corrected <- foreach::foreach(feature = featureNames(object), .combine = rbind) %dopar% {
    tmp <- BatchCorrMetabolomics::doBC(Xvec = exprs(object)[feature, ],
         ref.idx = ref_idx,
         batch.idx = batch_idx,
         seq.idx = seq_idx,
         minBsamp = 1,
         method = "lm")
    matrix(tmp, nrow = 1, dimnames = list(feature, names(tmp)))
  }

  exprs(object) <- batch_corrected

  object
}

#' Remove Unwanted Variation
#'
#' An interface for the RUVs method in RUVSeq package.
#'
#' @param object a MetaboSet object
#' @param batch the column name for batch labels
#' @param ref the column name for reference sample labels
#' @param ref_label the label for reference samples
#' @param k The number of factors of unwanted variation to be estimated from the data.
#' @param ... other parameters passed to RUVSeq::RUVs
#'
#' @return a MetaboSet object with the normalized data
#'
#' @export
ruvs_qc <- function(object, batch, ref, ref_label, k = 3, ...) {

  # Transform data to pseudo counts for RUVs
  exprs(object)[exprs(object) == 0] <- 1
  exprs(object) <- round(exprs(object))

  # Create list of the replicate samples
  # Cross-batch reference samples and QC samples of each batch
  replicates <- list(which(pData(object)[, ref] == ref_label))
  bathces <- pData(object)[, batch]
  for (b in unique(bathces)) {
    batch_qcs <- which(object[batches == b, ]$QC == "QC")
    replicates <- c(replicates, list(batch_qcs))
  }
  # Pad each replicate vector with -1 and transform to matrix
  max_len <- max(sapply(replicates, length))
  scIdx <- matrix(-1, nrow = length(replicates), ncol = max_len)
  for (i in seq_along(replicates)) {
    scIdx[i, seq_along(replicates[[i]])] <- replicates[[i]]
  }

  ruv_results <- RUVSeq::RUVs(x = exprs(object), cIdx = featureNames(object),
                              k = k, scIdx = scIdx, ...)

  exprs(object) <- ruv_results$normalizedCounts
  pData(object) <- cbind(pData(object), ruv_results$W)
  object
}
