
#' Summary statistics of finite elements
#'
#' These functions first remove non-finite and missing values, then compute the summary statistic in question.
#' They are helper functions used for computing quality measurements.
#'
#' @param x a numeric vector.
#' @name finite_helpers
NULL

#' @rdname finite_helpers
finite_sd <- function(x) {
  sd(x[is.finite(x)], na.rm = TRUE)
}

#' @rdname finite_helpers
finite_mean <- function(x) {
  mean(x[is.finite(x)], na.rm = TRUE)
}

#' @rdname finite_helpers
finite_median <- function(x) {
  median(x[is.finite(x)], na.rm = TRUE)
}

#' @rdname finite_helpers
finite_mad <- function(x) {
  mad(x[is.finite(x)], center = finite_median(x), na.rm = TRUE)
}


prop_na <- function(x) {
  sum(is.na(x)) / length(x)
}


prop_found <- function(x) {
  sum(!is.na(x)) / length(x)
}

# Join a dataframe to fData
join_fdata <- function(object, dframe) {
  fData(object) <- dplyr::left_join(fData(object),
                                    dframe,
                                    by = "Feature_ID")
  rownames(fData(object)) <- fData(object)$Feature_ID
  object
}


looks_numeric <- function(x) {
  suppressWarnings(all(!is.na(as.numeric(x))))
}

best_class <- function(x) {
  x <- as.character(x)
  if (looks_numeric(x)) {
    as.numeric(x)
  } else if (length(unique(x)) == length(x)) {
    as.character(x)
  } else {
    as.factor(x)
  }
}


best_classes <- function(x) {
  as.data.frame(lapply(x, best_class), stringsAsFactors = FALSE)
}
