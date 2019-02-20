
#' @importFrom magrittr "%>%"
read_from_excel <- function(file, corner_row, corner_column, split_by) {
  dada <- openxlsx::read.xlsx(file, colNames = FALSE)

  cc <- ifelse(is.character(corner_column),
               which(LETTERS == corner_column),
               corner_column)
  cr <- corner_row

  # Extract sample information
  pheno_data <- as.data.frame(t(dada[1:cr, (cc+1):ncol(dada)]), stringsAsFactors = FALSE)
  colnames(pheno_data) <- gsub(" ", "_", c(dada[1:(cr-1), cc], "Datafile"))
  if ("Sample_ID" %in% colnames(pheno_data)) {
    rownames(pheno_data <- pheno_data$Sample_ID)
  } else {
    rownames(pheno_data) <- paste0("ID_", 1:nrow(pheno_data))
  }


  # Exctract feature information
  feature_data <- dada[(cr+1):nrow(dada), 1:cc]
  colnames(feature_data) <- dada[cr, 1:cc]

  # Create feature ID if necessary
  if (!"Feature_ID" %in% colnames(feature_data)){
    feature_data <- name_features(feature_data = feature_data,
                                  split_by = split_by)
  }
  feature_data <- feature_data %>%
    tidyr::unite("Split", split_by, remove = FALSE) %>%
    dplyr:: mutate(Flag = NA_character_) %>%
    dplyr::select(Feature_ID, Split, dplyr::everything())

  feature_data <- name_features(feature_data, split_by)
  rownames(feature_data) <- feature_data$Feature_ID

  # Extract LC-MS measurements as matrix
  assay_data <- dada[(cr+1):nrow(dada), (cc+1):ncol(dada)] %>%
    apply(2, as.numeric)
  rownames(assay_data) <- rownames(feature_data)
  colnames(assay_data) <- rownames(pheno_data)

  return(list(assay_data = assay_data, pheno_data = pheno_data, feature_data = feature_data))
}

#' @importFrom magrittr "%>%"
name_features <- function(feature_data, split_by) {

  # Find mass and retention time columns
  mz_tags <- c("mass", "average mz", "average.mz")
  rt_tags <-  c("retention time", "retentiontime", "average rt(min)",
                "^rt$")

  mzCol <- NULL
  for (tag in mz_tags) {
    hits <- grepl(tag, tolower(colnames(feature_data)))
    if (any(hits)) {
      mz_col <- colnames(feature_data)[which(hits)[1]]
      break
    }
  }
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
    tidyr::unite("Feature_ID", c(split_by, "Feature_ID"), remove = FALSE)

  feature_data
}

#' @import methods
#' @importClassesFrom Biobase ExpressionSet
MetaboSet <- setClass("MetaboSet",
                      slots = c(stage = "character",
                                group = "character",
                                time = "character",
                                subject_id = "character"),
                      contains = "ExpressionSet")

setValidity("MetaboSet",
            function(object) {
              if (!is.na(object@group) & !object@group %in% colnames(object@phenoData@data)) {
                paste0("Column '", object@group, "' not found in pheno data")
              } else if (!is.na(object@time) & !object@time %in% colnames(object@phenoData@data)) {
                paste("Column", object@time, "not found in pheno data")
              } else if (!is.na(object@subject_id) & !object@subject_id %in% colnames(object@phenoData@data)) {
                paste("Column", object@subject_id, "not found in pheno data")
              } else {
                TRUE
              }
            })


construct_MetaboSet <- function(assay_data, pheno_data, feature_data,
                                group_col = NA_character_, time_col = NA_character_,
                                subject_id_col = NA_character_) {

  pheno_data <- new("AnnotatedDataFrame",
                    data=pheno_data)

  # Split the data by the Split column of feature data
  parts <- unique(feature_data$Split)
  obj_list <- list()
  for (part in parts) {
    fd_tmp <- new("AnnotatedDataFrame",
                  data= feature_data[feature_data$Split == part, ])
    ad_tmp <- assay_data[fd_tmp$Feature_ID,]
    obj_list[[part]] <- MetaboSet(exprs = ad_tmp,
                                  phenoData = pheno_data,
                                  featureData = fd_tmp,
                                  stage = "Original",
                                  group = group_col,
                                  time = time_col,
                                  subject_id = subject_id_col)
  }

  obj_list
}



setGeneric("lcms_data", signature = "object",
           function(object) standardGeneric("lcms_data"))

#' @importFrom Biobase exprs pData
setMethod("lcms_data", c(object = "MetaboSet"),
          function(object) {
            cbind(pData(object), t(exprs(object)))
          })

# ------------ Accessors and Replacers -----------------

# group
setGeneric("group", signature = "object",
           function(object) standardGeneric("group"))

setMethod("group", "MetaboSet",
          function(object) object@group)

setGeneric("group<-", signature = "object",
           function(object, value) standardGeneric("group<-"))

setMethod("group<-", "MetaboSet",
          function(object, value) {
            object@group <- value
            if (validObject(object)) {
              return(object)
            }
          })

# time
setGeneric("time", signature = "object",
           function(object) standardGeneric("time"))

setMethod("time", "MetaboSet",
          function(object) object@time)

setGeneric("time<-", signature = "object",
           function(object, value) standardGeneric("time<-"))

setMethod("time<-", "MetaboSet",
          function(object, value) {
            object@time <- value
            if (validObject(object)) {
              return(object)
            }
          })

# subject_id
setGeneric("subject_id", signature = "object",
           function(object) standardGeneric("subject_id"))

setMethod("subject_id", "MetaboSet",
          function(object) object@subject_id)

setGeneric("subject_id<-", signature = "object",
           function(object, value) standardGeneric("subject_id<-"))

setMethod("subject_id<-", "MetaboSet",
          function(object, value) {
            object@subject_id <- value
            if (validObject(object)) {
              return(object)
            }
          })



















