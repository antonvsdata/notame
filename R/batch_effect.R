
#' Bhattacharyya distance between bathces in PCA space
#'
#' Computes Bhattacharyya distance between all pairs of batches after
#' projecting the samples into PCA space with pcaMethods::pca
#'
#' @param object a MetaboSet object
#' @param batch column name of pData givinh the batch labels
#'
#' @return matrix of Bhattacharyya distances between batches
pca_bhattacharyya_dist <- function(object, batch, ...) {
  # PCA to 2 dimenstions
  pca_res <- pcaMethods::pca(object, ...)
  pca_scores <- pcaMethods::scores(pca_res)

  # Split to batches
  batches <- list()
  for (b in unique(pData(object)[, batch])) {
    batches[[b]] <- pca_scores[pData(object)[, batch] == b, ]
  }

  # Compute means and covariance matrices for Bhattacharyya distance
  muarray <- sapply(batches, colMeans)
  sigmaarray <- array(sapply(batches, cov), dim = c(2, 2, 3))

  fpc::bhattacharyya.matrix(muarray,sigmaarray,ipairs="all", misclassification.bound = FALSE)
}


pooled_variance <- function(x, group) {
  # Remove missing values
  group <- group[!is.na(x)]
  x <- x[!is.na(x)]
  # Split to groups
  group_list <- split(x, group)
  n_1 <- sapply(group_list, length) - 1 # n - 1
  # Pooled variance
  sum(n_1 * sapply(group_list, var)) / sum(n_1)
}

between_variance <- function(x, group) {
  # Remove missing values
  group <- group[!is.na(x)]
  x <- x[!is.na(x)]
  # Split to groups
  group_list <- split(x, group)
  n <- sapply(group_list, length)
  means <- sapply(group_list, mean)
  k_1 <- length(unique(group)) - 1
  # Between group variance formula
  sum(n * (means - mean(x))^2) / k_1
}

repeatability <- function(x, group) {
  pv <- pooled_variance(x, group)
  bv <- between_variance(x, group)
  bv / (bv + pv)
}


#' Compute repeatability measures
#'
#' Computes repeatability for each feature with the following formula:
#' \deqn{\frac{\sigma^2_{between}}{\sigma^2_{between} + \sigma^2_{within}}}
#' The repeatability ranges from 0 to 1.
#'
#'
#' @param object a MetaboSet object
#' @param group column name of pData givinh the group labels
#'
#' @return data frame with one row per feature with the repeatability measure
perform_repeatability <- function(object, group) {

  group <- pData(object)[, group]
  repeatabilities <- foreach::foreach(feature = featureNames(object), .combine = rbind) %dopar% {
    result_row <- data.frame(Feature_ID = feature,
                             Repeatability = repeatability(exprs(object)[feature, ], group))
  }
  repeatabilities
}

