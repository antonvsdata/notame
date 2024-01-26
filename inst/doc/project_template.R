## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  library(notame)
#
#  # Define the project path
#  ppath <- ""
#
#  # Start logging
#  init_log(log_file = paste0(ppath, "results/log_preprocessing.txt"))
#
#  # Read data from Excel file
#  data <- read_from_excel(file = paste0(ppath, "data/"),
#                          split_by = c("RP_HILIC", "POS_NEG"))
#
#  # Modify pheno_data if necessary (set factor levels etc.)
#
#  # Create MetaboSet objects for each analytical mode
#  objects <- construct_metabosets(exprs = data$exprs,
#                                  pheno_data = data$pheno_data,
#                                  feature_data = data$feature_data,
#                                  group_col = "status")
#
#
#  # Preprocessing by mode
#  processed <- list()
#  for (name in names(objects)) {
#    log_text(paste("Processing", name))
#    object <- objects[[name]]
#    log_text(paste(nrow(object), "features,", ncol(object), "samples\n"))
#
#    # Mark NAs correctly and flag features with low detection
#    object <- mark_nas(object, 0)
#    detected <- flag_detection(object, qc_limit = 0.7, group_limit = 0.5)
#
#    # Visualizations with original data
#    visualizations(detected, prefix = paste0(ppath, "results/figures/", name, "_ORIG"))
#
#    # Drift correction with cubic spline in log space
#    corrected <- correct_drift(detected)
#
#    # Flag low-quality features
#    flagged <- flag_quality(corrected)
#
#    # Visualizations after quality control
#    visualizations(flagged, prefix = paste0(ppath, "results/figures/", name, "_CLEANED"))
#
#    # Remove QC samples before imputation
#    noqc <- drop_qcs(flagged)
#
#    # Random forest imputation
#    # First impute qood quality features
#    set.seed(38)
#    imputed <- impute_rf(noqc)
#    # Then impute the rest of the features
#    imputed <- impute_rf(imputed, all_features = TRUE)
#
#    # Save the processed object to list
#    processed[[name]] <- imputed
#  }
#
#  # Merge analytical modes
#  merged <- merge_metabosets(processed)
#  log_text(paste("Merged analytical modes together, the merged object has", nrow(merged), "features and", ncol(merged), "samples."))
#
#  # Visualize complete dataset
#  visualizations(merged, prefix = paste0(ppath, "results/figures/FULL"))
#
#  # Save processed objects, useful if you want to run statistics later
#  save(pre_drift_correction, processed, merged, file = paste0(ppath, "results/preprocessed.RData"))
#  log_text("Saved preprocessed objects to results subfolder")
#
#  # Statistics, naturally depends a lot on the project, but this is a general workflow
#  lm_results <- perform_lm(merged, formula_char = "Feature ~ Group")
#
#  # Add results to object
#  # Remove intercept columns
#  lm_results <- dplyr::select(lm_results, -contains("Intercept"))
#  with_results <- join_fData(merged, lm_results)
#
#  # Write results to Excel
#  write_to_excel(with_results, file = paste0(ppath, "results/results.xlsx"))
#
#  # Draw plots of the most interesting features
#  # Select interesting features based on statistical relevance
#  interesting <- lm_results$Feature_ID[lm_results$GroupB_P_FDR < 0.05] # NOTE: just an example
#
#  # Save relevant visualizations of all the interesting features
#
#  # Group box plots as an example
#  save_group_boxplots(merged[interesting, ], color = "Group",
#                      file = paste0(ppath, "results/figures/group_boxplots.pdf"))
#
#  # Finish logging and save session information
#  finish_log()
