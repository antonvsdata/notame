
#' Fit Random Forest
#'
#' Fits a random forest, where given response column in pheno data is predicted using the features. Can be used
#' both for classification and regression. For more information, see the documentation of \code{randomForest::randomForest}.
#' After fitting the random forest, use rf_importance as a shortcut for getting the feature importance in random forest prediction.
#'
#' @param object a MetaboSet object
#' @param response character, column name of phenoData giving response
#' @param all_features should all features be included in the model? if FALSE, flagged features are left out
#' @param importance Should importance of features be assessed?
#' @param ... other parameters passed to \code{randomForest::randomForest}
#'
#' @return An object of class randomForest
#'
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{rf_importance}}
#'
#' @examples
#' rf <- fit_rf(example_set, response = "Group")
#' rf
#' rf_importance(rf)
#'
#' @export
fit_rf <- function(object, response, all_features = FALSE, importance = TRUE, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Package \"randomForest\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!all_features) {
    object <- drop_flagged(object)
  }

  rf <- randomForest::randomForest(x = t(exprs(object)), y = pData(object)[, response], importance = importance, ...)

  rf
}


#' Feature importance in random forest
#'
#' Extracts feature importance in random forest in a nice format.
#'
#' @param rf An object of class randomForest
#'
#' @return data frame of feature importance
#'
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{fit_rf}}
#'
#' @examples
#' rf <- fit_rf(example_set, response = "Group")
#' rf
#' rf_importance(rf)
#'
#' @export
importance_rf <- function(rf) {
  # Extract metrics and feature ID
  df <- data.frame(Feature_ID = rownames(rf$importance),
                   as.data.frame(rf$importance),
                   stringsAsFactors = FALSE, check.names = FALSE)
  df
}


#' PLS-DA
#'
#' A simple wrapper for fitting a PLS-DA model using plsda function of the mixOmics package.
#' The object can then be passed to many of the mixOmics functions for prediction,
#' performance evaluation etc. To run initial optimization for a PLS-DA model,
#' try \code{\link{mixomics_plsda_optimize}}.
#'
#' @param object a MetaboSet object
#' @param y character, column name of the grouping variable to predict
#' @param ... any parameters passed to \code{mixOmics::plsda}
#'
#' @return an object of class "plsdsa"
#'
#' @examples
#' plsda_res <- mixomics_plsda(merged_sample, y = "Group")
#'
#' @export
mixomics_plsda <- function(object, y, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # Extract X and Y matrices
  X <- t(exprs(object))
  Y <- pData(object)[, y]
  # Y needs to be a factor, this ensures the levels are right
  if (class(Y) != "factor") {
    Y <- as.factor(Y)
    warning(paste(y, "is not encoded as a factor, converted to factor with levels:", paste(levels(Y), collapse = ", ")))
  }
  log_text("Fitting PLS-DA")
  mixOmics::plsda(X, Y, ...)
}


#' PLS-DA
#'
#' A wrapper for fitting a PLS-DA model using plsda function of the mixOmics package.
#' Automatically evaluates performance with different number of components and
#' chooses the optimal number of components.
#'
#' @param object a MetaboSet object
#' @param y character, column name of the grouping variable to predict
#' @param ncomp_max numeric, the maximum number of components to try
#' @param ... any parameters passed to \code{mixOmics::plsda}
#'
#' @return an object of class "plsdsa"
#'
#' @examples
#' plsda_res <- mixomics_plsda_optimize(merged_sample, y = "Group", ncomp_max = 5)
#'
#' @export
mixomics_plsda_optimize <- function(object, y, ncomp_max, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  plsda_res <- mixomics_plsda(object = object, y = y, ncomp_max, ...)

  log_text("Evaluating PLS-DA performance")
  perf_plsda <- mixOmics::perf(plsda_res, validation = "Mfold", folds = 5, auc = TRUE, nrepeat = 5)

  plot(perf_plsda, col = mixOmics::color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

  # Find the distance metric with minimum BER
  ber <- perf_plsda$error.rate$BER
  inds <- which(ber == min(ber), arr.ind = TRUE)[1,]
  dist_met <- colnames(ber)[inds[2]]
  # Find the optimal number of components
  ncomp_opt <- perf_plsda$choice.ncomp["BER", dist_met]
  log_text(paste("Choosing a PLS-DA model with", ncomp_opt, "components using", dist_met))

  mixomics_plsda(object = object, y = y, ncomp = ncomp_opt, ...)
}


#' sPLS-DA
#'
#' A wrapper for fitting an sPLS-DA model using splsda function of the mixOmics package.
#' Automatically evaluates performance with different number of components and
#' different number of features per component, then chooses the optimal number of components and
#' optimal number of features for each component.
#'
#' @param object a MetaboSet object
#' @param y character, column name of the grouping variable to predict
#' @param ncomp_max numeric, the maximum number of components to try
#' @param dist the distance metric to use, one of "max.dist", "mahalanobis.dist", "centroids.dist"
#' @param n_features the number of features to try for each component
#' @param ... any parameters passed to \code{mixOmics::plsda}
#'
#' @return an object of class "splsdsa"
#'
#' @examples
#' plsda_res <- mixomics_splsda_optimize(merged_sample, dist = "max.dist", y = "Group", ncomp_max = 5)
#'
#' @export
mixomics_splsda_optimize <- function(object, y, ncomp_max, dist,
                                     n_features = c(1:10, seq(20, 300, 10)), ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("Package \"mixOmics\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  # Extract X and Y matrices
  X <- t(exprs(object))
  Y <- pData(object)[, y]
  # Y needs to be a factor, this ensures the levels are right
  if (class(Y) != "factor") {
    Y <- as.factor(Y)
    warning(paste(y, "is not encoded as a factor, converted to factor with levels:", paste(levels(Y), collapse = ", ")))
  }

  log_text("Tuning sPLS-DA")
  tuned_splsda <- mixOmics::tune.splsda(X, Y, ncomp = ncomp_max,
                                        validation = "Mfold", folds = 5,
                                        dist = dist,
                                        measure = "BER", nrepeat = 5,
                                        test.keepX = n_features)

  plot(tuned_splsda)

  ncomp_opt <- tuned_splsda$choice.ncomp$ncomp
  keep_x <- tuned_splsda$choice.keepX[1:ncomp_opt]
  log_text(paste("Final model has", ncomp_opt, "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))

  splsda_final <- mixOmics::splsda(X, Y, ncomp = ncomp_opt, keepX = keep_x)

  if (ncomp_opt > 1) {
    background <- mixOmics::background.predict(splsda_final, comp.predicted=2, dist = dist)
    mixOmics::plotIndiv(splsda_final, comp = 1:2, group = Y, ind.names = FALSE,
                        title = paste("Final sPLS-DA model,", dist, "prediction areas"), legend = TRUE, background = background)
    mixOmics::plotIndiv(splsda_final, comp = 1:2, group = Y, ind.names = FALSE,
                        title = paste("Final sPLS-DA model,", dist), legend = TRUE, ellipse = TRUE)
  }

  splsda_final
}
