## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(amp)
library(doParallel)

## ------------------------------------------------------------------------
path <- "~/test_project/"

## ------------------------------------------------------------------------
init_log(log_file = paste0(path, "log.txt"))
# Check logging state
log_state()

## ------------------------------------------------------------------------
data <- read_from_excel(file = system.file("extdata", "sample_data_bread.xlsx", package = "amp"), sheet = 1,
                        corner_row = 5, corner_column = "X",
                        split_by = c("Column", "Ion mode"))

## ---- out.width = "600px"------------------------------------------------
knitr::include_graphics("Data_input.png")

## ------------------------------------------------------------------------
names(data)
sapply(data, class)
sapply(data, dim)

## ------------------------------------------------------------------------
modes <- construct_MetaboSet(exprs = data$exprs, pheno_data = data$pheno_data,
                             feature_data = data$feature_data,
                             group_col = "Bread")

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

## ------------------------------------------------------------------------
visualizations(mode, prefix = paste0(path, "figures/", name, "_ORIG"))

## ------------------------------------------------------------------------
corrected <- correct_drift(mode)
visualizations(corrected, prefix = paste0(path, "figures/", name, "_DRIFT"))

## ------------------------------------------------------------------------
fData(corrected)$DC_note

## ------------------------------------------------------------------------
corrected <- corrected %>% assess_quality() %>% flag_quality()
processed[[i]] <- corrected

visualizations(corrected, prefix = paste0(path, "figures/", name, "_CLEANED"))

## ------------------------------------------------------------------------
# Initialize empty list for processed objects
processed <- list()
for (i in seq_along(modes)) {
  name <- names(modes)[i]
  mode <- modes[[i]]
  # Set all zero abundances to NA
  mode <- mark_nas(mode, value = 0)
  mode <- flag_detection(mode, qc_limit = 0.7, group_limit = 0.8)
  # visualizations(mode, prefix = paste0(path, "figures/", name, "_ORIG"))
  corrected <- correct_drift(mode)
  # visualizations(corrected, prefix = paste0(path, "figures/", name, "_DRIFT"))
  
  corrected <- corrected %>% assess_quality() %>% flag_quality()
  processed[[i]] <- corrected

  # visualizations(corrected, prefix = paste0(path, "figures/", name, "_CLEANED"))
}

## ------------------------------------------------------------------------
merged <- merge_metabosets(processed)

visualizations(merged, prefix = paste0(path, "figures/", name, "_FULL"))

## ------------------------------------------------------------------------
merged <- drop_qcs(merged)

visualizations(mode, prefix = paste0(path, "figures/", name, "_FULL_NO_QC"))

## ------------------------------------------------------------------------
#Set seed number for reproducibility
set.seed(38)
imputed <- impute_rf(merged)

## ------------------------------------------------------------------------
imputed <- impute_rf(imputed, all_features = TRUE)

## ------------------------------------------------------------------------
visualizations(imputed, prefix = paste0(path, "figures/", name, "_FULL_IMPUTED"))

## ------------------------------------------------------------------------
save(imputed, file = paste0(path, "full_data.RData"))

## ------------------------------------------------------------------------
anova_results <- perform_oneway_anova(imputed, formula_char = "Feature ~ Bread")

## ------------------------------------------------------------------------
top_features <- anova_results$Feature_ID[which(anova_results$ANOVA_P_FDR < 0.2)]

top_index <- fData(imputed)$Feature_ID %in% top_features
# Setting the group column is redundant, but made for clarity
pairwise_results <- perform_pairwise_t_test(imputed[top_index, ], group = "Bread")

## ------------------------------------------------------------------------
results(imputed) <- dplyr::left_join(anova_results, pairwise_results)

write_to_excel(imputed, file = paste0(path, "results.xlsx"))

## ------------------------------------------------------------------------
finish_log()

