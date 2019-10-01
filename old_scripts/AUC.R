#' Area under curve
#'
#' Compute area under curve for each subject and feature.
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param ci_level the confidence level used in constructing the confidence intervals
#' for regression coefficients
#' @param ... additional parameters passed to lm
#'
#' @return a data frame with one row per feature, with all the
#' relevant statistics of the linear model as columns
#'
#'
#' @examples
#' 
#'
#' @seealso \code{\link[DescTools]{AUC}}
#'
#' @export
perform_auc <- function(object, time = time_col(object), subject = subject_col(object),
                        group = group_col(object)) {
  
  # Start log
  log_text(paste("\nStarting AUC computation at", Sys.time()))
  
  data <- combined_data(object)
  
  # Create new pheno data, only one row per subject and group
  pheno_data <- data[, c(subject, group)] %>%
    dplyr::distinct() %>%
    tidyr::unite("Sample_ID", subject, group, remove = FALSE)
  rownames(pheno_data) <- pheno_data$Sample_ID
  
  # AUCs
  features <- featureNames(object)
  aucs <- foreach::foreach(i = seq_along(features), .combine = rbind) %dopar% {
    feature <- features[i]
    result_row <- rep(NA_real_, nrow(pheno_data))
    # Compute AUC for each subject in each group
    tryCatch({
      for(j in seq_len(nrow(pheno_data))){
        subset_idx <- data[, subject] == pheno_data[j, subject] & data[, group] == pheno_data[j, group]
        result_row[j] <- PK::auc(time = as.numeric(data[subset_idx, time]),
                                  conc = data[subset_idx, feature], design = "complete")$est[1]
      }
    })
    
    matrix(result_row, nrow = 1, dimnames = list(feature, pheno_data$Sample_ID))
  }
  
  # Construct new MetaboSet object (with all modes together)
  new_object <- construct_MetaboSet(exprs = aucs, feature_data = fData(object),
                                    pheno_data = pheno_data, group_col = group,
                                    subject_col = subject) %>%
    merge_metabosets()
  
  new_object
  
}

load("~/ville_gradu/merged.RData")

auced <- compute_auc(imputed)
subject_col(imputed)
View(pData(imputed))

pheno_data <- pData(imputed)

breads <- unique(pheno_data$Bread)
n_time <- length(levels(pheno_data$Time))

subject_id <- c()
for (bread in breads) {
  subject_id <- c(subject_id, rep(paste(bread, 1:6), n_time, sep = "_"))
}
pheno_data$Subject_ID <- subject_id
pData(imputed) <- pheno_data
subject_col(imputed) <- "Subject_ID"


auced <- perform_auc(imputed)




