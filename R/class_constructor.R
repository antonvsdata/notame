
check_pheno_data <- function(x, id_prefix) {

  # Check that Injection order is included
  if (!"Injection_order" %in% colnames(x)) {
    stop('"Injection_order" not found for the samples')
  }
  # Injection order should be unique
  if (length(unique(x$Injection_order)) != nrow(x)) {
    stop("Injection_order is not unique")
  }
  # If Sample_ID is not provided explicitly, it will be created
  if (!"Sample_ID" %in% colnames(x)) {
    x$Sample_ID <- paste0(id_prefix, x$Injection_order)
  } else {
    # Add a running index to all "QC" identifiers
    x$Sample_ID <- as.character(x$Sample_ID)
    x$Sample_ID[x$Sample_ID == "QC"] <- paste0("QC_", seq_len(sum(x$Sample_ID == "QC")))
    # After this, the Sample IDs should be unique
    if (length(unique(x$Sample_ID)) != nrow(x)) {
      stop("Sample_ID is not unique")
    }
  }

  x <- best_classes(x)
  rownames(x) <- x$Sample_ID
  x <- as.data.frame(dplyr::select(x, Sample_ID, dplyr::everything()))
  x
}


check_position <- function(x, cc, cr) {
  condition <- (is.na(x[cr - 1, cc - 1])) &
    (is.numeric(type.convert(x[cr + 1, cc + 1]))) &
    (!is.na(x[cr, cc]))
  if (!condition) {
    stop("The corner row and column coordinates seem to be incorrect!")
  }

}

#' @importFrom magrittr "%>%"
read_from_excel <- function(file, sheet, corner_row, corner_column, id_prefix = "ID_", split_by = NULL, name = NULL) {

  if (is.null(split_by) & is.null(name)) {
    stop("Etiher namr or split_by needs to be defined, see documentation")
  } else if ((!is.null(split_by)) & (!is.null(name))) {
    stop("Only define split_by OR name, see documentation")
  }

  dada <- openxlsx::read.xlsx(file, sheet, colNames = FALSE)

  # Define excel column order A-Z, AA - ZZ
  combinations <- expand.grid(LETTERS, LETTERS)
  excel_columns <- c(LETTERS, paste0(combinations$Var2, combinations$Var1))

  # Column can be given as a character
  cc <- ifelse(is.character(corner_column),
               which(excel_columns == corner_column),
               corner_column)
  cr <- corner_row

  check_position(dada, cc, cr)

  # Extract sample information
  pheno_data <- as.data.frame(t(dada[1:cr, (cc+1):ncol(dada)]), stringsAsFactors = FALSE)
  colnames(pheno_data) <- gsub(" ", "_", c(dada[1:(cr-1), cc], "Datafile"))

  pheno_data <- check_pheno_data(x = pheno_data, id_prefix = id_prefix)

  # Exctract feature information
  feature_data <- dada[(cr+1):nrow(dada), 1:cc]
  colnames(feature_data) <- dada[cr, 1:cc]

  # If the file only contains one mode, add the mode name as Split column
  if (!is.null(name)) {
    feature_data$Split <- name
    split_by <- "Split"
  } else { # Multiple modes in the file, create Split column to separate modes
    feature_data <- feature_data %>%
      tidyr::unite("Split", split_by, remove = FALSE)
  }

  # Create feature ID if necessary
  if (!"Feature_ID" %in% colnames(feature_data)){
    feature_data <- name_features(feature_data = feature_data)
  }
  # Reorganise columns and add Flag column
  feature_data <- feature_data %>%
    dplyr::select(Feature_ID, Split, dplyr::everything()) %>%
    dplyr:: mutate(Flag = NA_character_)
  rownames(feature_data) <- feature_data$Feature_ID

  # Extract LC-MS measurements as matrix
  assay_data <- dada[(cr+1):nrow(dada), (cc+1):ncol(dada)] %>%
    apply(2, as.numeric)
  rownames(assay_data) <- rownames(feature_data)
  colnames(assay_data) <- rownames(pheno_data)

  return(list(assay_data = assay_data, pheno_data = pheno_data, feature_data = feature_data))
}
# Combines mode name, mass and retention time to create a Feature ID
#' @importFrom magrittr "%>%"
name_features <- function(feature_data) {

  # Find mass and retention time columns
  mz_tags <- c("mass", "average mz", "average.mz")
  rt_tags <-  c("retention time", "retentiontime", "average rt[(]min[)]",
                "^rt$")

  mz_col <- NULL
  for (tag in mz_tags) {
    hits <- grepl(tag, tolower(colnames(feature_data)))
    if (any(hits)) {
      mz_col <- colnames(feature_data)[which(hits)[1]]
      break
    }
  }
  rt_col <- NULL
  for (tag in rt_tags) {
    hits <- grepl(tag, tolower(colnames(feature_data)))
    if (any(hits)) {
      rt_col <- colnames(feature_data)[which(hits)[1]]
      break
    }
  }

  if (is.null(mz_col)){
    stop(paste0("No mass to charge ratio column found - should match one of:\n",
                paste(mz_tags, collapse = ", "), " (not case-sensitive)"))
  }
  if (is.null(rt_col)){
    stop(paste0("No retention time column found - should match one of:\n",
                paste(rt_tags, collapse = ", "), " (not case-sensitive)"))
  }

  # Concatenate rounded mass and retention time
  round_mz <- as.numeric(feature_data[, mz_col]) %>% round(digits = 4) %>%
    as.character() %>% gsub("[.]", "_", .)
  round_rt <- as.numeric(feature_data[, rt_col]) %>% round(digits = 4) %>%
    as.character() %>% gsub("[.]", "_", .)
  feature_data$Feature_ID <- paste0(round_mz, "a", round_rt)

  # Add the split columns (usually column, mode and possibly tissue)

  feature_data <- feature_data %>%
    tidyr::unite("Feature_ID", c("Split", "Feature_ID"), remove = FALSE)

  feature_data
}

#' @import methods
#' @importClassesFrom Biobase ExpressionSet
MetaboSet <- setClass("MetaboSet",
                      slots = c(group_col = "character",
                                time_col = "character",
                                subject_col = "character",
                                predicted = "matrix"),
                      contains = "ExpressionSet")

setValidity("MetaboSet",
            function(object) {
              if (!is.na(object@group_col) & !object@group_col %in% colnames(object@phenoData@data)) {
                paste0("Column '", object@group_col, "' not found in pheno data")
              } else if (!is.na(object@time_col) & !object@time_col %in% colnames(object@phenoData@data)) {
                paste("Column", object@time_col, "not found in pheno data")
              } else if (!is.na(object@subject_col) & !object@subject_col %in% colnames(object@phenoData@data)) {
                paste("Column", object@subject_col, "not found in pheno data")
              } else {
                TRUE
              }
            })

construct_MetaboSet <- function(assay_data, pheno_data, feature_data,
                                group_col = NA_character_, time_col = NA_character_,
                                subject_col = NA_character_) {

  pheno_data <- Biobase::AnnotatedDataFrame(data=pheno_data)

  # Split the data by the Split column of feature data
  parts <- unique(feature_data$Split)
  obj_list <- list()
  for (part in parts) {
    fd_tmp <- Biobase::AnnotatedDataFrame(data= feature_data[feature_data$Split == part, ])
    ad_tmp <- assay_data[fd_tmp$Feature_ID,]
    obj_list[[part]] <- MetaboSet(exprs = ad_tmp,
                        phenoData = pheno_data,
                        featureData = fd_tmp,
                        group_col = group_col,
                        time_col = time_col,
                        subject_col = subject_col,
                        predicted = matrix(NA_real_, nrow = nrow(ad_tmp),
                                           ncol = ncol(ad_tmp),
                                           dimnames = dimnames(ad_tmp)))
  }

  obj_list
}


# ------------ Accessors and Replacers -----------------

setGeneric("combined_data", signature = "object",
           function(object) standardGeneric("combined_data"))

#' @importFrom Biobase exprs pData
setMethod("combined_data", c(object = "MetaboSet"),
          function(object) {
            cbind(pData(object), t(exprs(object)))
          })


# stage
setGeneric("stage", signature = "object",
           function(object) standardGeneric("stage"))

setMethod("stage", "MetaboSet",
          function(object) object@stage)

setGeneric("stage<-", signature = "object",
           function(object, value) standardGeneric("stage<-"))

setMethod("stage<-", "MetaboSet",
          function(object, value) {
            object@stage <- value
            if (validObject(object)) {
              return(object)
            }
          })

# group
setGeneric("group_col", signature = "object",
           function(object) standardGeneric("group_col"))

setMethod("group_col", "MetaboSet",
          function(object) object@group_col)

setGeneric("group_col<-", signature = "object",
           function(object, value) standardGeneric("group_col<-"))

setMethod("group_col<-", "MetaboSet",
          function(object, value) {
            object@group_col <- value
            if (validObject(object)) {
              return(object)
            }
          })

# time
setGeneric("time_col", signature = "object",
           function(object) standardGeneric("time_col"))

setMethod("time_col", "MetaboSet",
          function(object) object@time_col)

setGeneric("time_col<-", signature = "object",
           function(object, value) standardGeneric("time_col<-"))

setMethod("time_col<-", "MetaboSet",
          function(object, value) {
            object@time_col <- value
            if (validObject(object)) {
              return(object)
            }
          })

# subject ID
setGeneric("subject_col", signature = "object",
           function(object) standardGeneric("subject_col"))

setMethod("subject_col", "MetaboSet",
          function(object) object@subject_col)

setGeneric("subject_col<-", signature = "object",
           function(object, value) standardGeneric("subject_col<-"))

setMethod("subject_col<-", "MetaboSet",
          function(object, value) {
            object@subject_col <- value
            if (validObject(object)) {
              return(object)
            }
          })

# predicted values from spline regression
setGeneric("predicted", signature = "object",
           function(object) standardGeneric("predicted"))

setMethod("predicted", "MetaboSet",
          function(object) object@predicted)

setGeneric("predicted<-", signature = "object",
           function(object, value) standardGeneric("predicted<-"))

setMethod("predicted<-", "MetaboSet",
          function(object, value) {
            object@predicted <- value
            if (validObject(object)) {
              return(object)
            }
          })
















