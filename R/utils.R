

# Set default color scales on load
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.amp <- list(
    amp.color_scale_c = ggplot2::scale_color_viridis_c(),
    amp.color_scale_d = ggplot2::scale_color_brewer(palette = "Set1")
  )
  toset <- !(names(op.amp) %in% names(op))
  if(any(toset)) options(op.amp[toset])

  invisible()
}


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
  mad(x[is.finite(x)], center = median(x[is.finite(x)], na.rm = TRUE), na.rm = TRUE)
}


# Defaults for NULL values
`%||%` <- function(a, b) if (is.null(a) | is.na(a)) b else a

# Proportion of NA values in a vector
prop_na <- function(x) {
  sum(is.na(x)) / length(x)
}

# Proportion of non-missing values in a vector
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

# Replace all instances of value in exprs with NA
mark_nas <- function(object, value) {
  ex <- exprs(object)
  ex[ex == value] <- NA
  exprs(object) <- ex
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
