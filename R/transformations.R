
#' Mark specified values as missing
#'
#' Replaces all values in the exprs part that equal the specified value with NA.
#' For example, vendor software might use 0 or 1 to signal a missing value,
#' which is not understood by R.
#'
#' @param object a MetaboSet object
#' @param value the value to be converted to NA
#'
#' @return MetaboSet object as the one supplied, with missing values correctly set to NA
#'
#' @export
mark_nas <- function(object, value) {
  ex <- exprs(object)
  ex[ex == value] <- NA
  exprs(object) <- ex
  object
}



#' Replace NA values in pheno data columns with "QC"
#'
#' Loops through specified columns in pheno data and replaces all NA values with "QC". Also transforms
#' the column to factor and sets "QC" as the first level.
#'
#' @param data a data frame such as pheno_data returned by \code{\link{read_from_excel}}
#' @param cols the columns to fix
#'
#' @return a data frame with the column modified
#'
#' @export
mark_qcs <- function(data, cols) {

  for (col in cols) {
    data[col] <- as.character(data[, col])
    data[is.na(data[, col]), col] <- "QC"
    data[col] <- factor(data[, col])
    data[col] <- relevel(data[, col], ref = "QC")
  }
  data
}


#' Drop QC samples
#'
#' @param object a MetaboSet object
#'
#' @return MetaboSet object as the one supplied, without QC samples
#'
#' @importFrom Biobase pData "pData<-"
#'
#' @export
drop_qcs <- function(object) {
  object <- object[, object$QC != "QC"]
  pData(object) <- droplevels(pData(object))
  object
}


#' Drop flagged features
#'
#' Removes all features that have been flagged by quality control functions.
#' Only features that do not have a flag (Flag == NA) are retained.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be retained? Mainly used by internal functions
#'
#' @export
drop_flagged <- function(object, all_features = FALSE) {
  if (!all_features) {
    object <- object[is.na(flag(object)), ]
  }
  object
}

#' Partially replace exprs field with new values
#'
#' Replaces a subset of data in exprs of an object by new values.
#' Used after an operation such as imputation is computed only for a subset of features or samples
#' such as only good-quality features that have not been flagged.
#' This function is mainly used internally, but can be useful.
#'
#' @param object a MetaboSet object
#' @param y matrix containing new values to be merged into exprs
#'
#' @export
merge_exprs <- function(object, y) {
  # Colnames and rownames should be found in the object
  if (!all(colnames(y) %in% colnames(exprs(object))) | is.null(colnames(y))) {
    stop("Column names of y do not match column names of exprs(object)")
  }
  if(!all(rownames(y) %in% rownames(exprs(object))) | is.null(rownames(y))) {
    stop("Row names of y do not match row names of exprs(object)")
  }

  exprs_tmp <- exprs(object)
  exprs_tmp[rownames(y), colnames(y)] <- y
  exprs(object) <- exprs_tmp
  object
}

#' Impute missing values using random forest
#'
#' Impute the missing values in the exprs part of the object using a
#' random forest. The estimated error in the imputation is logged.
#' It is recommended to set the seed number for reproducibility
#' (it is called random forest for a reason).
#' This a wrapper around \code{missForest::missForest}.
#' Use parallelize = "variables" to run in parallel for faster testing.
#' NOTE: running in parallel prevents user from setting a seed number.
#' \strong{CITATION:} When using this function, cite the \code{missForest} package
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before imputation.
#' @param ... passed to MissForest function
#'
#' @return MetaboSet object as the one supplied, with missing values imputed.
#'
#' @seealso \code{\link[missForest]{missForest}} for detail about the algorithm and the parameters
#'
#' @export
impute_rf <- function(object, all_features = FALSE, ...) {
  # Start log
  log_text(paste("\nStarting random forest imputation at", Sys.time()))
  # Drop flagged features
  dropped <- drop_flagged(object, all_features)

  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("missForest package not found")
  }

  # Impute missing values
  mf <- missForest::missForest(xmis = t(exprs(dropped)), ...)
  imputed <- t(mf$ximp)
  # Log imputation error
  log_text(paste0("Out-of-bag error in random forest imputation: ",
                 round(mf$OOBerror, digits = 3), "\n"))
  # Assign imputed data to the droppped
  rownames(imputed) <- rownames(exprs(dropped))
  colnames(imputed) <- colnames(exprs(dropped))
  # Attach imputed abundances to object
  object <- merge_exprs(object, imputed)
  log_text(paste("Random forest imputation finished at", Sys.time()))
  object
}

#' Inverse-rank normalization
#'
#' Applies inverse rank normalization to all features to approximate
#' a normal distribution.
#'
#' @param object a MetaboSet object
#'
#' @return MetaboSet object as the one supplied, with normalized features
#'
#' @export
inverse_normalize <- function(object) {

  exprs(object) <- exprs(object) %>%
    apply(1, function(x) {
      qnorm((rank(x, na.last="keep")-0.5) / sum(!is.na(x)))}) %>%
    t()
  object
}
