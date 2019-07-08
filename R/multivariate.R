
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
fit_rf <- function(object, response, importance = TRUE, all_features = FALSE, ...) {

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
  # Choose importance metrics based on RF type
  if (rf$type == "classification") {
    cols <- c("MeanDecreaseAccuracy", "MeanDecreaseGini")
  } else {
    cols <- c("%IncMSE", "IncNodePurity")
  }
  # Extract metrics and feature ID
  df <- data.frame(Feature_ID = rownames(rf$importance),
                   as.data.frame(rf$importance)[cols],
                   stringsAsFactors = FALSE, check.names = FALSE)
  df
}
