## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(notame)
library(doParallel)

## ------------------------------------------------------------------------
ppath <- "~/test_project/"

## ------------------------------------------------------------------------
init_log(log_file = paste0(ppath, "log.txt"))
# Check logging state
log_state()

## ------------------------------------------------------------------------
data <- read_from_excel(file = system.file("extdata", "sample_data_whole.xlsx", package = "notame"), sheet = 1,
                        corner_row = 4, corner_column = "X",
                        split_by = c("Column", "Ion mode"))

## ---- out.width = "600px", echo=FALSE------------------------------------
knitr::include_graphics("Data_input.png")

## ------------------------------------------------------------------------
names(data)
sapply(data, class)
sapply(data, dim)

## ------------------------------------------------------------------------
modes <- construct_metabosets(exprs = data$exprs, pheno_data = data$pheno_data,
                             feature_data = data$feature_data,
                             group_col = "Group")

## ------------------------------------------------------------------------
names(modes)
sapply(modes, class)

## ------------------------------------------------------------------------
# Initialize empty list for processed objects
processed <- list()
for (i in seq_along(modes)) {
  # PREPROCESSING STEPS
}

## ------------------------------------------------------------------------
i <- 1
name <- names(modes)[i]
mode <- modes[[i]]

## ------------------------------------------------------------------------
# Set all zero abundances to NA
mode <- mark_nas(mode, value = 0)

## ------------------------------------------------------------------------
mode <- flag_detection(mode, qc_limit = 0.7, group_limit = 0.8)

## ---- eval=FALSE---------------------------------------------------------
#  visualizations(mode, prefix = paste0(ppath, "figures/", name, "_ORIG"))

## ---- include = FALSE----------------------------------------------------
corrected <- correct_drift(mode)

## ----eval = FALSE--------------------------------------------------------
#  corrected <- correct_drift(mode)
#  visualizations(corrected, prefix = paste0(ppath, "figures/", name, "_DRIFT"))

## ------------------------------------------------------------------------
fData(corrected)$DC_note

## ---- include = FALSE----------------------------------------------------
corrected <- flag_quality(corrected)
processed[[i]] <- corrected

## ---- eval = FALSE-------------------------------------------------------
#  corrected <- flag_quality(corrected)
#  processed[[i]] <- corrected
#  
#  visualizations(corrected, prefix = paste0(ppath, "figures/", name, "_CLEANED"))

## ------------------------------------------------------------------------
# Initialize empty list for processed objects
processed <- list()
for (i in seq_along(modes)) {
  name <- names(modes)[i]
  mode <- modes[[i]]
  # Set all zero abundances to NA
  mode <- mark_nas(mode, value = 0)
  mode <- flag_detection(mode, qc_limit = 0.7, group_limit = 0.8)
  # visualizations(mode, prefix = paste0(ppath, "figures/", name, "_ORIG"))
  corrected <- correct_drift(mode)
  # visualizations(corrected, prefix = paste0(ppath, "figures/", name, "_DRIFT"))
  
  corrected <- corrected %>% assess_quality() %>% flag_quality()
  processed[[i]] <- corrected

  # visualizations(corrected, prefix = paste0(ppath, "figures/", name, "_CLEANED"))
}

## ---- include = FALSE----------------------------------------------------
merged <- merge_metabosets(processed)

## ---- eval = FALSE-------------------------------------------------------
#  merged <- merge_metabosets(processed)
#  
#  visualizations(merged, prefix = paste0(ppath, "figures/_FULL"))

## ---- include = FALSE----------------------------------------------------
merged_no_qc <- drop_qcs(merged)


## ---- eval = FALSE-------------------------------------------------------
#  merged_no_qc <- drop_qcs(merged)
#  
#  visualizations(merged_no_qc, prefix = paste0(ppath, "figures/FULL_NO_QC"))

## ------------------------------------------------------------------------
#Set seed number for reproducibility
set.seed(38)
imputed <- impute_rf(merged_no_qc)

## ------------------------------------------------------------------------
imputed <- impute_rf(imputed, all_features = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  visualizations(imputed, prefix = paste0(ppath, "figures/FULL_IMPUTED"))

## ---- eval = FALSE-------------------------------------------------------
#  save(imputed, file = paste0(ppath, "full_data.RData"))

## ------------------------------------------------------------------------
anova_results <- perform_oneway_anova(imputed, formula_char = "Feature ~ Group")

## ------------------------------------------------------------------------
top_features <- anova_results$Feature_ID[which(anova_results$ANOVA_P_FDR < 0.2)]

top_index <- fData(imputed)$Feature_ID %in% top_features

pairwise_results <- perform_pairwise_t_test(imputed[top_index, ], group = "Group")

## ---- eval = FALSE-------------------------------------------------------
#  combined_results <- dplyr::left_join(anova_results, pairwise_results)
#  imputed <- join_fData(imputed, combined_results)
#  
#  write_to_excel(imputed, file = paste0(ppath, "results.xlsx"))

## ------------------------------------------------------------------------
finish_log()

