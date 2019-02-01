
#' Summary statistics of finite elements
#'
#' These functions first remove non-finite and missing values, then compute the summary statistic in question.
#' They are helper functions used for computing quality measurements.
#'
#' @param x a numeric vector.
#' @name finite_helpers
NULL

#' @rdname finite_helpers
finite_sd <- function(x){
  sd(x[is.finite(x)], na.rm = TRUE)
}

#' @rdname finite_helpers
finite_mean <- function(x){
  mean(x[is.finite(x)], na.rm = TRUE)
}

#' @rdname finite_helpers
finite_median <- function(x){
  median(x[is.finite(x)], na.rm = TRUE)
}

#' @rdname finite_helpers
finite_mad <- function(x){
  mad(x[is.finite(x)], center = finite_median(x), na.rm = TRUE)
}
