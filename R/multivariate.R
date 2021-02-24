


# ----------- Random Forest ----------------

#' Fit Random Forest
#'
#' Fits a random forest, where given response column in pheno data is predicted using the features. Can be used
#' both for classification and regression. For more information, see the documentation of \code{randomForest::randomForest}.
#' After fitting the random forest, use rf_importance as a shortcut for getting the feature importance in random forest prediction.
#'
#' @param object a MetaboSet object
#' @param y character, column name of phenoData giving the dependent variable of the model
#' @param all_features logical, should all features be included in the model? if FALSE, flagged features are left out
#' @param covariates character, column names of pData to use as covariates in the model, in addition to
#' molecular features
#' @param importance Should importance of features be assessed?
#' @param ... other parameters passed to \code{randomForest::randomForest}
#'
#' @return An object of class randomForest
#'
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{importance_rf}}
#'
#' @examples
#' rf <- fit_rf(example_set, y = "Group")
#' rf
#' importance_rf(rf)
#'
#' @export
fit_rf <- function(object, y, all_features = FALSE, covariates = NULL, importance = TRUE, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Package \"randomForest\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  add_citation("randomFOrest package was used to fit random forest models:", citation("randomForest"))

  object <- drop_flagged(object, all_features = all_features)

  X <- combined_data(object)[, c(featureNames(object), covariates)]
  rf <- randomForest::randomForest(x = X, y = pData(object)[, y], importance = importance, ...)

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
#' rf <- fit_rf(example_set, y = "Group")
#' rf
#' importance_rf(rf)
#'
#' @export
importance_rf <- function(rf) {
  # Extract metrics and feature ID
  df <- data.frame(Feature_ID = rownames(rf$importance),
                   as.data.frame(rf$importance),
                   stringsAsFactors = FALSE, check.names = FALSE)
  df
}


# ------------------ mixOmics PLS ---------------------------

#' A helper function for extracting predictor matrix with covariates
#'
#' @param object a MetaboSet object
#' @param covariates character, column names of pData to use as covariates in the model, in addition to
#' molecular features
get_x <- function(object, covariates) {
  # Convert covariates to numeric
  if (any(!sapply(pData(object)[, covariates], looks_numeric))) {
    stop("All covariates should be convertable to numeric")
  }
  pData(object)[covariates] <- lapply(pData(object)[covariates], as.numeric)

  # Extract X
  X <- combined_data(object)[, c(featureNames(object), covariates)]
  X
}

#' Plot points in PLS space
#'
#' A helper function for \code{mixomics_pls} and \code{mixomics_spls}.
#'
#' @param model a PLS or sPLS model
#' @param Y the Y matrix
#' @param y the name of the y variable
#' @param title plot title
plot_pls <- function(model, Y, y, title) {
  # Extract scores and add y variable
  scores <- data.frame(model$variates$X[, 1:2])
  colnames(scores) <- c("X1", "X2")
  scores[, y[1]] <- Y[, 1]
  # Explained variance as percentage
  var_exp <- 100 * model$explained_variance$X[1:2] %>% round(digits = 3)
  p <- ggplot(scores, aes_string(x = "X1", y = "X2", color = y)) +
    geom_point() +
    getOption("notame.color_scale_con") +
    theme_minimal() +
    labs(x = paste("X1:", var_exp[1], "%"),
         y = paste("X2:", var_exp[2], "%"),
         title = title)
  print(p)
}


#' PLS
#'
#' Simple wrappers for fitting a PLS model using pls function of the mixOmics package.
#' The object can then be passed to many of the mixOmics functions for prediction,
#' performance evaluation etc. Also plot a scores plot of the first two components.
#' \itemize{
#' \item{\code{mixomics_pls} A simple PLS model with set number of components and all features}
#' \item{\code{mixomics_pls_optimize} Test different numbers of components, choose the one with minimal mean square error}
#' \item{\code{mixomics_spls_optimize} sPLS model: Test different numbers of components and features,
#' choose the one with minimal mean square error}
#' }
#'
#' @param object a MetaboSet object
#' @param y character vector, column names of the grouping variable to predict
#' @param ncomp number of X components
#' @param folds the number of folds to use in k-fold cross validation
#' @param nrepeat the number of times to repeat the cross validation. Lower this for faster testing.
#' @param plot_scores logical, if TRUE, a scatter plot with the first two PLS-components as x and y-axis will
#' be drawn, colored by the Y-variable. Only really makes sense if y is a single variable
#' @param all_features logical, should all features be included in the model? if FALSE, flagged features are left out
#' @param covariates character, column names of pData to use as covariates in the model, in addition to
#' molecular features
#' @param n_features the number of features to try for each component
#' @param ... any parameters passed to \code{mixOmics::pls} or \code{mixOmics::spls}
#'
#' @return an object of class "mixo_pls" or "mixo_spls"
#'
#' @examples
#' \dontrun{
#' pls_res <- mixomics_pls(merged_sample, y = "Injection_order", ncomp = 3)
#' pls_opt <- mixomics_pls_optimize(merged_sample, y = "Injection_order", ncomp = 3)
#' pls_res <- mixomics_spls_optimize(merged_sample, y = "Injection_order", ncomp = 3,
#'                                   n_features <- c(1:10, 12, 15, 20))
#' }
#' @name pls
#' @seealso \code{\link[mixOmics]{pls}}, \code{\link[mixOmics]{perf}},
#' \code{\link[mixOmics]{spls}}, \code{\link[mixOmics]{tune.spls}}
NULL

#' @rdname pls
#' @export
mixomics_pls <- function(object, y, ncomp, plot_scores = TRUE, all_features = FALSE,
                         covariates = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  add_citation("mixOmics package was used to fit PLS models:", citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)

  X <- get_x(object, covariates)
  Y <- pData(object)[y]

  log_text("Fitting PLS")
  pls_model <- mixOmics::pls(X, Y, ncomp = ncomp, ...)

  if (plot_scores & ncomp > 1) {
    plot_pls(pls_model, Y, y, title = "PLS: first 2 components and the Y variable")
  }
  pls_model
}

#' @rdname pls
#'
#' @export
mixomics_pls_optimize <- function(object, y, ncomp, folds = 5, nrepeat = 50, plot_scores = TRUE,
                                  all_features = FALSE, covariates = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  add_citation("mixOmics package was used to fit PLS models:", citation("mixOmics"))

  pls_res <- mixomics_pls(object = object, y = y, ncomp = ncomp, plot_scores = FALSE,
                          all_features = all_features, covariates = covariates, ...)

  log_text("Evaluating PLS performance")
  perf_pls <- mixOmics::perf(pls_res, validation = "Mfold", folds = folds, nrepeat = nrepeat)

  # Plot Mean Square Error
  p1 <- ggplot(data.frame(ncomp = seq_len(ncomp),
                          MSEP = as.vector(perf_pls$MSEP)), aes(x = ncomp, y = MSEP)) +
    geom_line() +
    labs(color = NULL, title = "Mean Square Error") +
    theme_bw() +
    scale_x_continuous(breaks = seq_len(ncomp)) +
    theme(panel.grid.minor.x = element_blank())

  # Plot R2 and Q2
  plot_data <- data.frame(R2 = as.vector(perf_pls$R2),
                          Q2 = as.vector(perf_pls$Q2),
                          ncomp = seq_len(ncomp)) %>%
    tidyr::gather(key = "key", value = "value", -ncomp)

  p2 <- ggplot(plot_data, aes(x = ncomp, y = value, color = key)) +
    geom_line() +
    labs(color = NULL, title = "R2 and Q2") +
    theme_bw() +
    getOption("notame.color_scale_dis") +
    scale_x_continuous(breaks = seq_len(ncomp)) +
    theme(panel.grid.minor.x = element_blank())

  if (requireNamespace("cowplot", quietly = TRUE)) {
      print(cowplot::plot_grid(p1, p2, nrow = 1))
  } else {
    print(p1)
    print(p2)
  }


  # Find the optimal number of components
  ncomp_opt <- which(perf_pls$MSEP == min(perf_pls$MSEP))[1]
  log_text(paste0("Choosing a PLS model with ", ncomp_opt, " component(s) based on the minimal MSE\n",
                 "Take a look at the plot and make sure this is the correct number of components"))

  mixomics_pls(object = object, y = y, ncomp = ncomp_opt, plot_scores = plot_scores,
               covariates = covariates, ...)
}

#' @rdname pls
#'
#' @export
mixomics_spls_optimize <- function(object, y, ncomp,
                                   n_features = c(1:10, seq(20, 300, 10)), folds = 5, nrepeat = 50,
                                   plot_scores = TRUE, all_features = FALSE, covariates = NULL, ...) {

  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  add_citation("mixOmics package was used to fit PLS models:", citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)

  X <- get_x(object, covariates)
  Y <- pData(object)[y]

  # Test different number of components and features with cross validation
  log_text("Tuning sPLS")
  tuned_spls <- mixOmics::tune.spls(X, Y, ncomp = ncomp,
                            test.keepX = n_features,
                            validation = "Mfold", folds = folds,
                            nrepeat =  nrepeat,
                            measure = 'MAE')
  # Plot error for each component with different number of features
  plot(tuned_spls)
  title("Performance of sPLS models")
  # Choose optimal numbers of components and features
  ncomp_opt <- tuned_spls$choice.ncomp$ncomp
  keep_x <- tuned_spls$choice.keepX[1:ncomp_opt]
  log_text(paste("Final model has", ncomp_opt, "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))
  # Fit the final model
  spls_final <- mixOmics::spls(X, Y, ncomp = ncomp_opt, keepX = keep_x, ...)
  # SCatter plot of points in PLS space
  if (ncomp_opt > 1) {
    plot_pls(spls_final, Y, y, title = "Final sPLS model: first 2 components and the Y variable")
  }

  spls_final
}

plot_plsda <- function(model, Y, title, dist = "max.dist") {
  background <- mixOmics::background.predict(model, comp.predicted = 2, dist = dist)
  mixOmics::plotIndiv(model, comp = 1:2, group = Y, ind.names = FALSE,
                      title = paste(title), legend = TRUE, ellipse = TRUE)
  mixOmics::plotIndiv(model, comp = 1:2, group = Y, ind.names = FALSE,
                      title = paste(title, "prediction areas"), legend = TRUE, background = background)
}


#' PLS-DA
#'
#' A simple wrapper for fitting a PLS-DA model using plsda function of the mixOmics package.
#' The object can then be passed to many of the mixOmics functions for prediction,
#' performance evaluation etc.
#' \itemize{
#' \item{\code{mixomics_plsda} A simple PLS-DA model with set number of components and all features}
#' \item{\code{mixomics_plsda_optimize} Test different numbers of components, choose the one with minimal balanced error rate}
#' \item{\code{mixomics_splsda_optimize} sPLS-DA model: Test different numbers of components and features,
#' choose the one with minimal balanced error rate}
#' }
#'
#' @param object a MetaboSet object
#' @param y character, column name of the grouping variable to predict
#' @param ncomp the number of X components
#' @param folds the number of folds to use in k-fold cross validation
#' @param nrepeat the number of times to repeat the cross validation. Lower this for faster testing.
#' @param n_features the number of features to try for each component
#' @param dist the distance metric to use, one of "max.dist", "mahalanobis.dist", "centroids.dist".
#' use \code{\link{mixomics_plsda_optimize}} to find the best distance metric
#' @param plot_scores logical, if TRUE, a scatter plot with the first two PLS-components as x and y-axis will
#' be drawn, with both prediction surface and ellipses
#' @param all_features logical, should all features be included in the model? if FALSE, flagged features are left out
#' @param covariates character, column names of pData to use as covariates in the model, in addition to
#' molecular features
#' @param ... any parameters passed to \code{mixOmics::plsda}
#'
#' @return an object of class "mixo_plsda"
#'
#' @examples
#' \dontrun{
#' noqc <- drop_qcs(merged_sample)
#' plsda_res <- mixomics_plsda(noqc, y = "Group", ncomp = 2)
#' set.seed(38)
#' plsda_opt <- mixomics_plsda_optimize(noqc, y = "Group", ncomp = 3)
#' set.seed(38)
#' splsda_opt <- mixomics_splsda_optimize(noqc, y = "Group", dist = "max.dist", ncomp = 2,
#'                                       n_features <- c(1:10, 12, 15, 20))
#' }
#' @name pls_da
#' @seealso \code{\link[mixOmics]{plsda}}, \code{\link[mixOmics]{perf}},
#' \code{\link[mixOmics]{splsda}}, \code{\link[mixOmics]{tune.splsda}}
NULL

#' @rdname pls_da
#' @export
mixomics_plsda <- function(object, y, ncomp, plot_scores = TRUE, all_features = FALSE,
                           covariates = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  add_citation("mixOmics package was used to fit PLS models:", citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)

  X <- get_x(object, covariates)
  Y <- pData(object)[, y]
  # Y needs to be a factor, this ensures the levels are right
  if (class(Y) != "factor") {
    Y <- as.factor(Y)
    warning(paste(y, "is not encoded as a factor, converted to factor with levels:", paste(levels(Y), collapse = ", ")))
  }
  log_text("Fitting PLS-DA")
  plsda_model <- mixOmics::plsda(X, Y, ncomp = ncomp,...)
  if (plot_scores & ncomp > 1) {
    plot_plsda(plsda_model, Y, title = "PLS-DA")
  }

  plsda_model
}


#' @rdname pls_da
#' @export
mixomics_plsda_optimize <- function(object, y, ncomp, folds = 5, nrepeat = 50, plot_scores = TRUE,
                                    all_features = FALSE, covariates = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  add_citation("mixOmics package was used to fit PLS models:", citation("mixOmics"))
  plsda_res <- mixomics_plsda(object = object, y = y, ncomp = ncomp, plot_scores = FALSE,
                              all_features = all_features, covariates = covariates, ...)

  log_text("Evaluating PLS-DA performance")
  perf_plsda <- mixOmics::perf(plsda_res, validation = "Mfold", folds = 5, auc = TRUE, nrepeat = 50)

  plot(perf_plsda, col = mixOmics::color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
  title("Performance of PLS-DA models")
  # Find the distance metric with minimum BER
  ber <- perf_plsda$error.rate$BER
  inds <- which(ber == min(ber), arr.ind = TRUE)[1,]
  dist_met <- colnames(ber)[inds[2]]
  # Find the optimal number of components
  ncomp_opt <- perf_plsda$choice.ncomp["BER", dist_met]
  log_text(paste("Choosing a PLS-DA model with", ncomp_opt, "components using", dist_met))

  mixomics_plsda(object = object, y = y, ncomp = ncomp_opt, plot_scores = plot_scores, ...)
}


#' @rdname pls_da
#' @export
mixomics_splsda_optimize <- function(object, y, ncomp, dist,
                                     n_features = c(1:10, seq(20, 300, 10)),
                                     folds = 5, nrepeat = 50,
                                     plot_scores = TRUE,
                                     all_features = FALSE,
                                     covariates = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("Package \"mixOmics\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  add_citation("mixOmics package was used to fit PLS models:", citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)

  X <- get_x(object, covariates)
  Y <- pData(object)[, y]
  # Y needs to be a factor, this ensures the levels are right
  if (class(Y) != "factor") {
    Y <- as.factor(Y)
    warning(paste(y, "is not encoded as a factor, converted to factor with levels:", paste(levels(Y), collapse = ", ")))
  }
  # Test different components and numbers of features with cross validation
  log_text("Tuning sPLS-DA")
  tuned_splsda <- mixOmics::tune.splsda(X, Y, ncomp = ncomp,
                                        validation = "Mfold", folds = folds,
                                        dist = dist,
                                        measure = "BER", nrepeat = nrepeat,
                                        test.keepX = n_features)
  # Plot error rate of different components as a function of number of features
  plot(tuned_splsda)
  title("Performance of sPLS-DA models")
  # Choose optimal numbers of components and features
  ncomp_opt <- tuned_splsda$choice.ncomp$ncomp
  keep_x <- tuned_splsda$choice.keepX[1:ncomp_opt]
  log_text(paste("Final model has", ncomp_opt, "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))
  # Fit the final model
  splsda_final <- mixOmics::splsda(X, Y, ncomp = ncomp_opt, keepX = keep_x)
  # Scatterplots with prediction surface and ellipses
  if (plot_scores & ncomp_opt > 1) {
    plot_plsda(splsda_final, Y, title = "Final sPLS-DA model", dist = dist)
  }

  splsda_final
}



# ------------------- MUVR --------------------------------

#' MUVR
#'
#' A wrapper around the MUVR algorithm from the MUVR package. For more information
#' about the algorithm, visit https://gitlab.com/CarlBrunius/MUVR.
#'
#' @param object a MetaboSet object
#' @param y character, column name in pData of the target variable to predict
#' @param id character, column name in pData of the subject ID variable in case of repeated measurements
#' @param multi_level logical, whether multi-level modeling should be applied, see Details
#' @param multi_level_var character, column name in pData of the variable for splitting the data in multi-level modeling
#' @param all_features logical, should all features be included in the model? if FALSE, flagged features are left out
#' @param covariates,static_covariates character, column names of pData to use as covariates in the model, in addition to
#' molecular features. For multi-level moddels, the change in \code{covariates} is computed, while
#' \code{static_covariates} are taken from the first time point. \code{static_covariates} are ignored for
#' non-multi-level models.
#' @param nRep Number of repetitions of double CV, parameter of MUVR
#' @param nOuter Number of outer CV loop segments, parameter of MUVR
#' @param nInner Number of inner CV loop segments, parameter of MUVR
#' @param varRatio Ratio of variables to include in subsequent inner loop iteration,
#'  parameter of MUVR
#' @param method Multivariate method. Supports 'PLS' and 'RF', parameter of MUVR
#' @param ... other parameters to \code{MUVR::MUVR}
#'
#' @details For example, sex should be entered as a static covariate, since the change in sex
#' is zero for all individuals, so computing the change and using that as a covariate does not make sense.
#'
#' @examples
#' \dontrun{
#' # Simple model, only 1 repetition for a quick example
#' rf_model <- muvr_analysis(drop_qcs(merged_sample), y = "Group", nRep = 1, method = "RF")
#'
#' # PLS on multilevel variable
#' pls_model <- muvr_analysis(drop_qcs(example_set), multi_level = TRUE,
#' id = "Subject_ID", multi_level_var = "Time")
#' }
#'
#' @seealso \code{\link[MUVR]{MUVR}}
#'
#' @export
muvr_analysis <- function(object, y = NULL, id = NULL, multi_level = FALSE, multi_level_var = NULL,
                          covariates = NULL, static_covariates = NULL, all_features = FALSE,
                          nRep = 5, nOuter = 6, nInner = nOuter - 1,
                          varRatio = 0.75, method = c("PLS", "RF"), ...) {

  if (!requireNamespace("MUVR", quietly = TRUE)) {
    stop("Package \"MUVR\" needed for this function to work. Please install it from
         https://gitlab.com/CarlBrunius/MUVR",
         call. = FALSE)
  }
  add_citation("MUVR package was used to fit multivariate models with variable selection:", citation("MUVR"))

  # MUVR can only use numeric input
  classes <- sapply(pData(object)[, c(covariates, static_covariates)], class)
  if (length(classes) && any(classes != "numeric")) {
    stop("MUVR can only deal with numeric inputs, please transform all covariates to numeric",
         call. = FALSE)
  }

  object <- drop_flagged(object, all_features = all_features)

  if (any(!sapply(pData(object)[, covariates], looks_numeric))) {
    stop("All covariates should be convertable to numeric")
  }
  pData(object)[covariates] <- lapply(pData(object)[covariates], as.numeric)

  # Classic MUVR
  if (!multi_level) {
    if (is.null(y)) {
      stop("y variable needs to be defined unless doing multi-level modeling")
    }
    X <- combined_data(object)[, c(featureNames(object), covariates)]
    Y <- pData(object)[, y]

    # Independent samples
    if (is.null(id)) {
      muvr_model <- MUVR::MUVR(X = X, Y = Y, nRep = nRep, nOuter = nOuter, nInner = nInner,
                               varRatio = varRatio, method = method, ...)
    } else {
      # Multiple measurements
      ID <- pData(object)[, id]
      muvr_model <- MUVR::MUVR(X = X, Y = Y, ID = ID, nRep = nRep, nOuter = nOuter, nInner = nInner,
                               varRatio = varRatio, method = method, ...)

    }
  } else { # Multi-level analysis
    if (is.null(id) | is.null(multi_level_var)) {
      stop("id and multi_level_var needed for multi-level modeling")
    }
    # Check that multi_level_var has only 2 unique values
    ml_var <- pData(object)[, multi_level_var] <- as.factor(pData(object)[, multi_level_var])
    if (length(levels(ml_var)) != 2) {
      stop("The multilevel variable should have exactly 2 unique values")
    } else {
      cat(paste("Computing effect matrix according to", multi_level_var, ":",
                levels(ml_var)[2], "-", levels(ml_var)[1]))
    }

    # Compute effect matrix with covariates
    cd <- combined_data(object)
    cd <- cd[order(cd[, id]), ]
    X1 <- cd[cd[, multi_level_var] == levels(ml_var)[1], c(featureNames(object), covariates)]
    X2 <- cd[cd[, multi_level_var] == levels(ml_var)[2], c(featureNames(object), covariates)]
    X <- X2 - X1
    # Add static covariates, where we don't want to compute change, such as sex
    X[, static_covariates] <- cd[cd[, multi_level_var] == levels(ml_var)[1], static_covariates]
    rownames(X) <- unique(cd[, id])

    # Modeling
    if (!is.null(y)) { # Compare change of multi_level_var between levels of y
      Y <- cd[cd[, multi_level_var] == levels(ml_var)[1], y]
      muvr_model <- MUVR::MUVR(X = X, Y = Y, nRep = nRep, nOuter = nOuter, nInner = nInner,
                               varRatio = varRatio, method = method, ...)
    } else { # Compare levels of multi_level_var
      muvr_model <- MUVR::MUVR(X = X, ML = TRUE, nRep = nRep, nOuter = nOuter, nInner = nInner,
                               varRatio = varRatio, method = method, ...)
    }

  }
  # Plot performance
  MUVR::plotVAL(muvr_model)

  muvr_model
}
