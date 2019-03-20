

# Set default color scales on load
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.amp <- list(
    amp.logging = FALSE,
    amp.log_file = NULL,
    amp.color_scale_con = ggplot2::scale_color_viridis_c(),
    amp.color_scale_dis = ggplot2::scale_color_brewer(palette = "Set1"),
    amp.fill_scale_con = ggplot2::scale_fill_viridis_c(),
    amp.fill_scale_dis = ggplot2::scale_fill_brewer(palette = "Set1"),
    amp.fill_scale_div = ggplot2::scale_fill_distiller(palette = "RdBu"),
    amp.shape_scale = ggplot2::scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 13))
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

#' @rdname finite_helpers
finite_quantile <- function(x, ...) {
  quantile(x[is.finite(x)], na.rm = TRUE, ...)
}



# Defaults for NULL values
`%||%` <- function(a, b) {
  suppressWarnings(if (is.null(a)){
    b
  } else if (is.na(a)){
    b
  } else {
    a
  })
}

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

best_class <- function(x) {
  x <- type.convert(as.character(x), as.is = TRUE)
  if (class(x) == "numeric") {
    x <- x
  } else if (length(unique(x)) < length(x)/4) {
    x <- as.factor(x)
  } else if (is.integer(x)) {
    x <- as.numeric(x)
  } else {
    x <- as.character(x)
  }
  x
}


best_classes <- function(x) {
  as.data.frame(lapply(x, best_class), stringsAsFactors = FALSE)
}
