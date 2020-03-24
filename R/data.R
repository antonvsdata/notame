#' Toy data set
#'
#' Contains imaginary data used in testing the package functions.
#' For more realistic data, see \code{\link{merged_sample}}
#' This dataset includes multiple observations from same subjects,
#' sampled at two timepoints and divided to two groups.
#' The dataset has 30 samples and 20 features.
"example_set"


#' Sample dataset
#'
#' Sample data from samples in two groups, at two timepoints.
#' Contains 83 features from four analytical modes. The analytical modes are
#' available as separate MetaboSet objects and also as a merged object.
#'
#' @format
#' \describe{
#'   \item{hilic_neg_sample}{features from HILIC column with negative ionization}
#'   \item{hilic_pos_sample}{features from HILIC column with positive ionization}
#'   \item{rp_neg_sample}{features from RP column with negative ionization}
#'   \item{rp_pos_sample}{features from RP column with positive ionization}
#'   \item{merged_sample}{all features merged together}
#' }
"merged_sample"

#' @rdname merged_sample
"hilic_neg_sample"

#' @rdname merged_sample
"hilic_pos_sample"

#' @rdname merged_sample
"rp_neg_sample"

#' @rdname merged_sample
"rp_pos_sample"
