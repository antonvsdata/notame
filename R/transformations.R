

impute_rf <- function(object, ...) {

  # Extract transposed metabolite data
  data <- t(exprs(object))

  # Impute missing values
  mf <- missForest(xmis = data.responses, parallelize = parallelize)
  imputed <- t(mf$ximp)
  # Log imputation error
  log_text(paste("\n Out-of-bag error in random forest imputation:",
                 round(mf$OOBerror, digits = 3), "\n"))
  # Assign imputed data to the object
  rownames(imputed) <- rownames(exprs(object))
  colnames(imputed) <- colnames(exprs(object))
  exprs(object) <- imputed
  object
}


inverse_normalize <- function(object) {

  exprs(object) <- exprs(object) %>%
    apply(1, function(x) {
      qnorm((rank(x, na.last="keep")-0.5) / sum(!is.na(x)))}) %>%
    t()
  object
}
