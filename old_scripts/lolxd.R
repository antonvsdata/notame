
devtools::document()

res <- read_from_excel(file = "~/amp/inst/extdata/split_data.xlsx",
                       sheet = 2,
                       corner_row = 4, corner_column = "G",
                       split_by = c("Column"))


lcms <- construct_MetaboSet(assay_data = res$assay_data,
                            pheno_data = res$pheno_data,
                            feature_data = res$feature_data)

lol <- lcms[[1]]

group_col(lol) <- "Group"

lol <- flag_detection(lol)

fData(lol)

group(lol)
group(lol) <- "LOL"
group(lol) <- "Group"

?Biobase::esApply

?spline


object <- lol


qc <- object[, object$QC == "QC"]
qc_order <- qc$Injection_order
qc_data <- exprs(qc)

full_order <- object$Injection_order
full_data <- exprs(object)

comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

dc_data <- foreach::foreach(i = seq_len(nrow(object)), .combine = comb) %dopar% {

  # Spline cannot be fitted if there are les than 4 QC values
  if (sum(!is.na(qc_data[i, ])) < 4) {
    return(list(dc_row = matrix(NA_real_, nrow = 1, ncol = ncol(full_data)),
                predicted = matrix(NA_real_, nrow = 1, ncol = ncol(full_data))))
  }

  fit <- smooth.spline(x = qc_order, y = qc_order)

}







splinem <- function(x) {
  if (sum(!is.na(x)) >= 4){
    smooth.spline(x = qc$Injection, y = x[!is.na(x)])
  } else {
    NA
  }
}



fit <- Biobase::esApply(qc, 1, splinem)


predictem <- function(object, x) {
  if (!is.na(object)) {
    predict(object, x = x)$y
  } else {
    NA
  }
}

lapply(fit, predictem, x = lol$Injection) %>% as.data.frame()
