
#' Mark specified values as missing
#'
#' Replaces all values in the exprs part that equal the specified value with NA.
#' For example, vendor software might use 0 or 1 to signal a missing value,
#' which is not understood by R.
#'
#' @param object a MetaboSet object
#' @param value the value to be converted to NA
#'
#' @return MetaboSet object as the one supplied, with missing values correctly set to NA
#'
#' @examples
#' nas_marked <- mark_nas(merged_sample, value = 0)
#'
#' @export
mark_nas <- function(object, value) {
  ex <- exprs(object)
  ex[ex == value] <- NA
  exprs(object) <- ex
  object
}

#' Transform the MS/MS output to publication ready
#'
#' Change the MS/MS output from MS-DIAL format to publication-ready format.
#' Original spectra is sorted according to abundance percentage and clarified. See the example below.
#'
#' Original MS/MS spectra from MS-DIAL:
#' m/z:Raw Abundance
#'
#' 23.193:254 26.13899:5 27.50986:25 55.01603:82 70.1914:16 73.03017:941 73.07685:13 73.13951:120
#'
#' Spectra after transformation:
#' m/z  (Abundance)
#'
#' 73.03 (100), 23.193 (27), 73.14 (12.8), 55.016 (8.7)
#'
#' @param object a MetaboSet object
#' @param peak_num maximum number of peak that is kept (Recommended: 4-10)
#' @param min_abund minimum relative abundance to be kept (Recommended: 1-5)
#' @param deci_num maximum number of decimals to m/z value (Recommended: >2)
#'
#' @return MetaboSet object as the one supplied, with publication-ready MS/MS peak information
#'
#' @examples
#' fixed_MSMS_peaks <- fix_MSMS(imputed)
#'
#' @export
fix_MSMS <- function(object, peak_num = 10,
                     min_abund = 5, deci_num = 3) {
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package \"stringr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  spec <- fData(object)$MS.MS.spectrum
  to_metab <- NULL

  for (i in 1:length(spec)) {
    # Check the feature spectra and skip if it doesn't exist
    if(is.na(spec[i])){
      to_metab[i] <- NA
      next()
    }
    # Transform format
    spectrum <- spec[i]
    spectrum2 <- t(stringr::str_split(spectrum, pattern = " ", simplify = T))
    spectrum3 <- as.data.frame(stringr::str_split(spectrum2, pattern = ":", simplify = T))
    spectrum3 <- as.data.frame(lapply(spectrum3,as.numeric))
    spectrum3 <- spectrum3[order(spectrum3$V2,decreasing = T),]

    # Leave n most intense fragment peaks or all peaks if number of peaks < n
    ifelse(nrow(spectrum3) > peak_num, num <- peak_num, num <- nrow(spectrum3))
    spectrum4 <- spectrum3[c(1:num),]

    # Round the m/z of fragments to n decimals and calculate the relative intensity (%)
    spectrum4$V1 <- round(spectrum4$V1,digits = deci_num)
    spectrum4$relative <- round(spectrum4$V2 / max(spectrum3$V2)*100, digits = 1)

    # Remove fragment peaks with relative intensity less than n% (recommended: 1-5)
    spectrum5 <- spectrum4[c(spectrum4$relative > min_abund),]

    # Finalize format and write results
    to_metab[i] <- paste(paste0(spectrum5$V1, " (", spectrum5$relative, ")"), collapse = ", ")
  }

  fData(object)$MS_MS_Spectrum_clean <- to_metab
  object
}

#' Replace NA values in pheno data columns with "QC"
#'
#' Loops through specified columns in pheno data and replaces all NA values with "QC". Also transforms
#' the column to factor and sets "QC" as the last level.
#'
#' @param data a data frame such as pheno_data returned by \code{\link{read_from_excel}}
#' @param cols the columns to fix
#'
#' @return a data frame with the column modified
#'
#' @examples
#' # Remove QC labels first
#' pheno_data <- pData(merged_sample)
#' pheno_data$Group[pheno_data$Group == "QC"] <- NA
#' head(pheno_data$Group, 20)
#' pheno_data <- mark_qcs(pheno_data, cols = "Group")
#' head(pheno_data$Group, 20)
#'
#' @export
mark_qcs <- function(data, cols) {

  for (col in cols) {
    data[col] <- as.character(data[, col])
    data[is.na(data[, col]), col] <- "QC"
    data[col] <- factor(data[, col])
    good_levels <- levels(data[, col])
    good_levels <- c(good_levels[good_levels != "QC"], "QC")
    data[col] <- factor(as.character(data[, col]), levels = good_levels)
  }
  data
}


#' Drop QC samples
#'
#' @param object a MetaboSet object
#'
#' @return MetaboSet object as the one supplied, without QC samples
#'
#' @examples
#' dim(merged_sample)
#' noqc <- drop_qcs(merged_sample)
#' dim(noqc)
#'
#' @importFrom Biobase pData "pData<-"
#'
#' @export
drop_qcs <- function(object) {
  object <- object[, object$QC != "QC"]
  pData(object) <- droplevels(pData(object))
  object
}


#' Drop flagged features
#'
#' Removes all features that have been flagged by quality control functions.
#' Only features that do not have a flag (Flag == NA) are retained.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be retained? Mainly used by internal functions
#'
#' @examples
#' dim(merged_sample)
#' flagged <- flag_quality(merged_sample)
#' noflags <- drop_flagged(flagged)
#' dim(noflags)
#'
#' @export
drop_flagged <- function(object, all_features = FALSE) {
  if (!all_features) {
    object <- object[is.na(flag(object)), ]
  }
  object
}

#' Partially replace exprs field with new values
#'
#' Replaces a subset of data in exprs of an object by new values.
#' Used after an operation such as imputation is computed only for a subset of features or samples
#' such as only good-quality features that have not been flagged.
#' This function is mainly used internally, but can be useful.
#'
#' @param object a MetaboSet object
#' @param y matrix containing new values to be merged into exprs
#'
#' @return a MetaboSet object with the new exprs values
#'
#' @examples
#' ex_set <- merged_sample[1:5, 1:5]
#' exprs(ex_set)
#' # Create a matrix of replacment values for rows 1, 3, 5 and columns 1, 3, 4
#' replacement <- matrix(1:9, ncol = 3,
#'                       dimnames = list(featureNames(ex_set)[c(1, 3, 5)],
#'                                       sampleNames(ex_set)[c(1, 3, 4)]))
#' replacement
#' merged <- merge_exprs(ex_set, replacement)
#' exprs(merged)
#'
#' @export
merge_exprs <- function(object, y) {
  # Colnames and rownames should be found in the object
  if (!all(colnames(y) %in% colnames(exprs(object))) | is.null(colnames(y))) {
    stop("Column names of y do not match column names of exprs(object)")
  }
  if(!all(rownames(y) %in% rownames(exprs(object))) | is.null(rownames(y))) {
    stop("Row names of y do not match row names of exprs(object)")
  }

  exprs_tmp <- exprs(object)
  exprs_tmp[rownames(y), colnames(y)] <- y
  exprs(object) <- exprs_tmp
  object
}

#' Impute missing values using random forest
#'
#' Impute the missing values in the exprs part of the object using a
#' random forest. The estimated error in the imputation is logged.
#' It is recommended to set the seed number for reproducibility
#' (it is called random forest for a reason).
#' This a wrapper around \code{missForest::missForest}.
#' Use parallelize = "variables" to run in parallel for faster testing.
#' NOTE: running in parallel prevents user from setting a seed number.
#' \strong{CITATION:} When using this function, cite the \code{missForest} package
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before imputation.
#' @param ... passed to MissForest function
#'
#' @return MetaboSet object as the one supplied, with missing values imputed.
#'
#' @examples
#' missing <- mark_nas(example_set, 0)
#' set.seed(38)
#' imputed <- impute_rf(missing)
#'
#' @seealso \code{\link[missForest]{missForest}} for detail about the algorithm and the parameters
#'
#' @export
impute_rf <- function(object, all_features = FALSE, ...) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
      stop("Package \"missForest\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  add_citation("missForest package was used for random forest imputation:", citation("missForest"))
  # Start log
  log_text(paste("\nStarting random forest imputation at", Sys.time()))
  # Drop flagged features
  dropped <- drop_flagged(object, all_features)

  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("missForest package not found")
  }

  # Impute missing values
  mf <- missForest::missForest(xmis = t(exprs(dropped)), ...)
  imputed <- t(mf$ximp)
  # Log imputation error
  log_text(paste0("Out-of-bag error in random forest imputation: ",
                 round(mf$OOBerror, digits = 3)))
  # Assign imputed data to the droppped
  rownames(imputed) <- rownames(exprs(dropped))
  colnames(imputed) <- colnames(exprs(dropped))
  # Attach imputed abundances to object
  object <- merge_exprs(object, imputed)
  log_text(paste("Random forest imputation finished at", Sys.time(), "\n"))
  object
}

#' Simple imputation
#'
#' Impute missing values using a simple imputation strategy. All missing values
#' of a feature are imputed with the same value. It is possible
#' to only impute features with a large number of missing values this way. This can be useful
#' for using this function before random forest imputation to speed things up.
#' The imputation strategies available are:
#' \itemize{
#' \item a numeric value: impute all missing values in all features with the same value, e.g. 1
#' \item "mean": impute missing values of a feature with the mean of observed values of that feature
#' \item "median": impute missing values of a feature with the median of observed values of that feature
#' \item "min": impute missing values of a feature with the minimum observed value of that feature
#' \item "half_min": impute missing values of a feature with half the minimum observed value of that feature
#' \item "small_random": impute missing values of a feature with random numbers between 0 and the
#' minimum of that feature (uniform distribution, remember to set the seed number!).
#' }
#'
#' @param object a MetaboSet object
#' @param value the value used for imputation, either a numeric or one of "min", "half_min", "small_random",
#' see above
#' @param na_limit only impute features with the proportion of NAs over this limit. For example, if
#' \code{na_limit = 0.5}, only features with at least half of the values missing are imputed.
#'
#' @examples
#' missing <- mark_nas(merged_sample, 0)
#' imputed <- impute_simple(missing, value = "min")
#'
#' @export
impute_simple <- function(object, value, na_limit = 0) {

  imp <- exprs(object)
  nas <- apply(imp, 1, prop_na)
  imp <- imp[nas > na_limit, ]
  if (nrow(imp) == 0) {
    warning("none of the features satisfy the NA limit, returning the original object")
    return(object)
  }

  # Replace all missing values with the given constant
  if (is.numeric(value)) {
    imp[is.na(imp)] <- value
  } else if (value == "mean") {
    imp <- t(apply(imp, 1, function(x){
      x[is.na(x)] <- finite_mean(x)
      x
    }))
  } else if (value == "median") {
    imp <- t(apply(imp, 1, function(x){
      x[is.na(x)] <- finite_median(x)
      x
    }))
  } else if (value == "min") {
    imp <- t(apply(imp, 1, function(x){
      x[is.na(x)] <- finite_min(x)
      x
    }))
  } else if (value == "half_min") {
    imp <- t(apply(imp, 1, function(x){
      x[is.na(x)] <- finite_min(x) / 2
      x
    }))
  } else if (value == "small_random") {
    imp <- t(apply(imp, 1, function(x){
      x[is.na(x)] <- runif(n = sum(is.na(x)), min = 0, max = finite_min(x))
      x
    }))
  } else {
    stop('value should be a numeric value or one of "min", "half_min", "small_random"')
  }

  obj <- merge_exprs(object, imp)
  obj
}

#' Inverse-rank normalization
#'
#' Applies inverse rank normalization to all features to approximate
#' a normal distribution.
#'
#' @param object a MetaboSet object
#'
#' @return MetaboSet object as the one supplied, with normalized features
#'
#' @examples
#' normalized <- inverse_normalize(merged_sample)
#'
#' @export
inverse_normalize <- function(object) {

  exprs(object) <- exprs(object) %>%
    apply(1, function(x) {
      qnorm((rank(x, na.last="keep")-0.5) / sum(!is.na(x)))}) %>%
    t()
  object
}


#' A report of flagged features
#'
#' Computes the number of features at each stage of flagging for each mode.
#'
#' @param object a MetaboSet object
#'
#' @return data frame with the number of features at each stage of flagging
#'
#' @examples
#' flagged <- merged_sample %>%
#'  mark_nas(0) %>%
#'  flag_detection() %>%
#'  flag_quality()
#' flag_report(flagged)
#'
#' @export
flag_report <- function(object) {

  splits <- sort(unique(fData(object)$Split))
  report <- data.frame()
  flag(object)[is.na(flag(object))] <- "Kept"
  for (split in splits) {
    tmp <- object[fData(object)$Split == split, ]
    report_row <- flag(tmp) %>% table %>% as.matrix() %>% t()
    report_row <- data.frame(Split = split, report_row)
    if (is.null(report_row$Kept)) {
      report_row$Kept <- 0
    }
    report_row$Total <- nrow(tmp)
    report_row$Flagged <- report_row$Total - report_row$Kept
    report <- suppressWarnings(dplyr::bind_rows(report, report_row))
  }
  report
}


# ---------- Logarithms ----------

#' Logarithm
#'
#' Log-transforms the exprs part of a MetaboSet object. Shortcust log2 and log10 also implemented.
#' For more information, see \code{\link{log}}
#'
#' @param x a MetaboSet object
#' @param base the base of the logarithm
#'
#' @export
setMethod("log", "MetaboSet", function(x, base = exp(1)) {
  exprs(x) <- log(exprs(x), base = base)
  x
})

#' @rdname log-MetaboSet-method
#' @export
setMethod("log2", "MetaboSet", function(x) {
  exprs(x) <- log2(exprs(x))
  x
})

#' @rdname log-MetaboSet-method
#' @export
setMethod("log10", "MetaboSet", function(x) {
  exprs(x) <- log10(exprs(x))
  x
})

# scale
setGeneric("scale")

#' Scale exprs data
#'
#' Applies the base R function scale to transposed exprs matrix. See ?scale for details
#'
#' @param x a MetaboSet object
#' @param center,scale as in base scale function
#'
#' @export
setMethod("scale", "MetaboSet", function(x, center = TRUE, scale = TRUE) {
  exprs(x) <- t(scale(t(exprs(x)), center = center, scale = scale))
  x
})

#' Exponential function
#'
#' Apply the exponential function to feature abundances (exprs)
#'
#' @param object a MetaboSet object
#' @param base base of the exponential
#'
#' @return a MetaboSet object with altered feature abundances
#'
#' @export
setGeneric("exponential", signature = "object",
           function(object, base = exp(1)) standardGeneric("exponential"))


#' @export
setMethod("exponential", c(object = "MetaboSet"),
          function(object, base = exp(1)) {
            exprs(object) <- base^exprs(object)
            object
          })
