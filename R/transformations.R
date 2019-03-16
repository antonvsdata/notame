
#' Mark specified values as missing
#'
#' Replaces all values in the exprs part that equal the specified value with NA.
#' For example, vendor softwares might use 0 or 1 to signal a missing value,
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

#' Impute missing values using random forest
#'
#' Impute the missing values in the exprs part of the object using a
#' random forest. The estimated error in the imputation is logged.
#' It is recommended to ste the seed number for reproducibility
#' (it is called random forest for a reason).
#' This a wrapper around missForest::missForest.
#' Use parallelize = "variables" to run in parallel for faster testing.
#' NOTE: running in parallel prevents user from setting a seed number.
#'
#' @param object a MetaboSet object
#' @param ... passed to MissForest function
#'
#' @return MetaboSet object as the one supplied, with missing values imputed.
#'
#' @seealso \code{\link[missForest]{missForest}} for detail about the algorithm and the parameters
#'
#' @export
impute_rf <- function(object, ...) {

  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("missForest package not found")
  }

  # Impute missing values
  mf <- missForest::missForest(xmis = t(exprs(object)))
  imputed <- t(mf$ximp)
  # Log imputation error
  log_text(paste0("\nOut-of-bag error in random forest imputation: ",
                 round(mf$OOBerror, digits = 3), "\n"))
  # Assign imputed data to the object
  rownames(imputed) <- rownames(exprs(object))
  colnames(imputed) <- colnames(exprs(object))
  exprs(object) <- imputed
  object
}

#' Inverse-rank normalization
#'
#' Applies inverse rank normalization to all features to approxiamte
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
