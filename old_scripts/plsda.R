devtools::load_all("~/amp")


mixomics_plsda <- function(object, y, ...) {
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


mixomics_plsda_optimize <- function(object, y, ncomp_max, ...) {
  
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

xd <- mixomics_splsda_optimize(object = merged_sample, y = "Group", ncomp_max = 3, dist = "mahalanobis.dist")

mixomics_splsda_optimize <- function(object, y, ncomp_max, dist,
                                     n_features = c(1:10, seq(20, 300, 10)), ...) {
  
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