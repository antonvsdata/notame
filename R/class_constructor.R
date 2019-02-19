
#' @importFrom magrittr "%>%"
read_from_excel <- function(file, corner_row, corner_column, split_by){
  dada <- openxlsx::read.xlsx(file, colNames = FALSE)

  cc <- corner_column
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
name_features <- function(feature_data, split_by){

  # Find mass and retention time columns
  mz_tags <- c("mass", "average mz", "average.mz")
  rt_tags <-  c("retention time", "retentiontime", "average rt(min)",
                "^rt$")

  mz_col <- NULL
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


MetaboSet <- setClass("LCMSSet",
                      slots = c(stage = "character",
                                group = "character",
                                time = "character",
                                subject_id = "character"),
                      contains = "ExpressionSet")


construct_MetaboSet <- function(assay_data, pheno_data, feature_data,
                                group_col = NULL, time_col = NULL, subject_id_col = NULL) {

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


































