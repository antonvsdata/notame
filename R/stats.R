
#' Summary statistics
#'
#' Computes summary statistics for each feature, possibly grouped by a factor.
#' The statistics include mean, standard deviation (sd), median,
#' median absolute deviation (mad), minimum (min), maximum (max)
#' as well as 25% and 75% quantiles (Q25 & Q75).
#'
#' @param object a MetaboSet object
#' @param grouping_cols character vector, the columns by which grouping should be done. Use \code{NA}
#' to compute statistics without grouping.
#'
#' @examples
#' # Group by "Group"
#' sum_stats <- summary_statistics(example_set)
#' # Group by Group and Time
#'  sum_stats <- summary_statistics(example_set, grouping_cols = c("Group", "Time"))
#' # No Grouping
#' sum_stats <- summary_statistics(example_set, grouping_cols = NA)
#'
#' @return a data frame with the summary statistics
#'
#' @export
summary_statistics <- function(object, grouping_cols = NA) {

  data <- combined_data(object)
  features <- Biobase::featureNames(object)
  # Possible grouping

  statistics <- foreach::foreach(i = seq_along(features), .combine = rbind,
                                 .export = c("finite_sd", "finite_mad", "finite_mean",
                                             "finite_max", "finite_min",
                                             "finite_median", "finite_quantile")) %dopar% {
    feature <- features[i]
    f_levels <- data[, feature]
    if (is.na(grouping_cols)[1]) {
      groups <- rep(1, nrow(data))
      group_names <- ""
    } else {
      # Single grouping column
      if (length(grouping_cols) == 1) {
        groups <- data[, grouping_cols]
        if (class(groups) == "factor") {
          group_names <- levels(groups)
        } else {
          group_names <- unique(groups)
        }
      } else {
        groups <- rep("", nrow(data))
        for (grouping_col in grouping_cols) {
          tmp_group <- paste(grouping_col, data[, grouping_col], sep = "_")
          groups <- paste(groups, tmp_group, sep = "_")
        }
        groups <- as.factor(gsub("^_", "", groups))
        group_names <- levels(groups)
      }

    }

    # Define functions to use
    funs <- list(mean = finite_mean,
                 sd = finite_sd,
                 median = finite_median,
                 mad = finite_mad,
                 min = finite_min,
                 Q25 = function(x){finite_quantile(x, probs = 0.25)},
                 Q75 = function(x){finite_quantile(x, probs = 0.75)},
                 max = finite_max)
    # Initialize named vector for the results
    result_row <- rep(0, times = length(group_names) * length(funs))
    if (!is.na(grouping_cols[1])) {
      var_names <- expand.grid(names(funs), group_names)
      names(result_row) <- paste(var_names$Var2, var_names$Var1, sep = "_")
    } else {
      names(result_row) <- names(funs)
    }
    # Compute statistics
    for (fname in names(funs)) {
      tmp <- tapply(f_levels, groups, funs[[fname]])
      if (is.na(grouping_cols[1])) {
        result_row[fname] <- tmp[1]
      } else {
        result_row[paste(group_names, fname, sep = "_")] <- tmp
      }
    }
    # Combine as data frame
    result_row <- data.frame(Feature_ID = feature, as.list(result_row),
                             stringsAsFactors = FALSE)

    result_row
  }

  statistics
}

#' Statistics cleaning
#'
#' Uses regexp to remove unnecessary columns from statistics results data frame.
#' Can also rename columns effectively.
#'
#' @param df data frame, statistics results
#' @param remove list, should contain strings that are matching to unwanted columns
#' @param rename named list, names should contain matches that are replaced with values
#'
#' @examples
#' # Simple manipulation to linear model results
#' lm_results <- perform_lm(drop_qcs(example_set), formula_char = "Feature ~ Group + Time")
#' lm_results <- clean_stats_results(lm_results,
#' rename = c("GroupB" = "GroupB_vs_A", "Time2" = "Time2_vs_1"))
#'
#' @export
clean_stats_results <- function(
    df,
    remove = c("Intercept", "CI95", "Std_error", "t_value", "z_value", "R2"),
    rename = NULL) {
  df <- df[, !grepl(paste(remove, collapse = "|"), colnames(df))]
  if (!is.null(rename)) {
    for (name in names(rename)) {
      colnames(df) <- gsub(name, rename[name], colnames(df))
    }
  }

  df
}

cohens_d_fun <- function(object, group, id, time) {
  data <- combined_data(object)
  features <- Biobase::featureNames(object)
  group_levels <- levels(data[, group])
  time_levels <- NULL

  if (is.null(time)) {
    group1 <- data[which(data[, group] == group_levels[1]), ]
    group2 <- data[which(data[, group] == group_levels[2]), ]
    log_text(paste("Starting to compute Cohen's D between groups",
                   paste(rev(group_levels), collapse = " & ")
    ))
  } else {
    time_levels <- levels(data[, time])
    # Split to time points
    time1 <- data[which(data[, time] == time_levels[1]), ]
    time2 <- data[which(data[, time] == time_levels[2]), ]
    common_ids <- intersect(time1[, id], time2[, id])
    rownames(time1) <- time1[, id]
    rownames(time2) <- time2[, id]
    time1 <- time1[common_ids, ]
    time2 <- time2[common_ids, ]
    if (!identical(time1[, group], time2[, group])) {
      stop("Groups of subjects do not match between time points",
           call. = FALSE)
    }
    # Change between time points
    new_data <- time2[, features] - time1[, features]
    # Split to groups
    group1 <- new_data[which(time1[, group] == levels(time1[,group])[1]), ]
    group2 <- new_data[which(time1[, group] == levels(time1[,group])[2]), ]

    log_text(paste("Starting to compute Cohen's D between groups",
                   paste(rev(group_levels), collapse = " & "),
                   "from time change",
                   paste(rev(time_levels), collapse = " - ")
    ))
  }
  ds <-  foreach::foreach(i = seq_along(features), .combine = rbind) %dopar% {
    feature <- features[i]
    f1 <- group1[, feature]
    f2 <- group2[, feature]
    d <- data.frame(Feature_ID = feature,
                    Cohen_d = (finite_mean(f2) - finite_mean(f1)) /
                      sqrt((finite_sd(f1)^2 + finite_sd(f2)^2) / 2),
                    stringsAsFactors = FALSE)
  }

  rownames(ds) <- ds$Feature_ID

  if (is.null(time_levels)) {
    colnames(ds)[2] <- paste0(group_levels[2], "_vs_", group_levels[1], "_Cohen_d")
  } else {
    colnames(ds)[2] <- paste0( group_levels[2], "_vs_", group_levels[1],
                              "_", time_levels[2], "_minus_", time_levels[1],
                              "_Cohen_d"
    )
  }

  log_text("Cohen's D computed.")
  ds
}


#' Cohen's D
#'
#' Computes Cohen's D for each feature. If time and ID are supplied,
#' change between two time points is computed for each subject,
#' and Cohen's d is computed from the changes
#'
#' @param object a MetaboSet object
#' @param id character, name of the subject ID column
#' @param group character, name of the group column
#' @param time character, name of the time column
#'
#' @return data frame with Cohen's d for each feature
#'
#' @examples
#' d_results <- cohens_d(drop_qcs(example_set))
#' d_results_time <- cohens_d(drop_qcs(example_set),
#'                            time = "Time", id = "Subject_ID")
#'
#' @export
cohens_d <- function(object, group = group_col(object),
                     id = NULL, time = NULL) {
  res <- NULL
  # Check that both group and time are factors and have at least two levels
  for (column in c(group, time)) {
    if (is.null(column)) {
      next
    }
    if (!is.factor(pData(object)[, column])) {
      stop(paste0("Column ", column, " should be a factor!"))
    }
    if (length(levels(pData(object)[, column])) < 2) {
      stop(paste("Column", column, "should have at least two levels!"))
    }
  }
  group_combos <- combn(levels(pData(object)[, group]), 2)

  count_obs_geq_than <- function(x, n) {
    sum(x >= n)
  }

  if(is.null(time)) {
    for (i in seq_len(ncol(group_combos))) {
      object_split <- object[, which(
        pData(object)[, group] %in% c(group_combos[1, i], group_combos[2, i])
      )]
      pData(object_split) <- droplevels(pData(object_split))

      if (is.null(res)) {
        res <- cohens_d_fun(object_split, group, id, time)
        } else {
        res <- dplyr::full_join(res,
                                cohens_d_fun(object_split, group, id, time),
                                by = "Feature_ID"
        )
      }
    }
  } else {
    if (is.null(id)) {
      stop("Please specify id column.", call. = FALSE)
    }
    time_combos <- combn(levels(pData(object)[, time]), 2)
    for (i in seq_len(ncol(group_combos))) {
      for (j in seq_len(ncol(time_combos))) {
        object_split <- object[, which(
          pData(object)[, group] %in% group_combos[, i] &
            pData(object)[, time] %in% time_combos[, j]
        )]
        pData(object_split) <- droplevels(pData(object_split))
        # Check data is valid for Cohen's D
        group_table <- table(pData(object_split)[, c(id, group)])
        time_table <- table(pData(object_split)[, c(id, time)])
        column <- paste0("Cohen_d_", group_combos[1, i], "_", group_combos[2, i],
                         "_", time_combos[2, j], "_minus_", time_combos[1, j]
        )
        if (any(apply(group_table, 2, count_obs_geq_than, 2) < 2)) {
          warning(paste0("In ", column,
                         ": Groups don't have two observations of at least two subjects, skipping!"
          ))
          next
        }
        if (any(apply(time_table, 1, count_obs_geq_than, 2) != 0)) {
          warning(paste0("In ", column,
                         ": Same subject recorded more than once at same time, skipping!"
          ))
          next
        }
        if (any(apply(group_table, 1, count_obs_geq_than, 1) != 1)) {
          warning(paste0("In ", column,
                         ": Same subject recorded in two groups, skipping!"
          ))
          next
        }
        if (!all(apply(time_table, 1, count_obs_geq_than, 1) != 1)) {
          warning(paste("One or more subject(s) missing time points,",
                        column, "will be counted using common subjects in time points!"
          ))
        }

        if (is.null(res)) {
          res <- cohens_d_fun(object_split, group, id, time)
        } else {
          res <- dplyr::full_join(res,
                                  cohens_d_fun(object_split, group, id, time),
                                  by = "Feature_ID"
          )
        }
      }
    }
  }
  rownames(res) <- res$Feature_ID
  res
}

#' Fold change
#'
#' Computes fold change between each group for each feature.
#'
#' @param object a MetaboSet object
#' @param group character, name of the group column
#'
#' @return data frame with fold changes for each feature
#'
#' @examples
#' # Between groups
#' fc <- fold_change(example_set)
#' # Between time points
#' fc <- fold_change(example_set, group = "Time")
#'
#' @export
fold_change <- function(object, group = group_col(object)) {

  log_text("Starting to compute fold changes.")

  data <- combined_data(object)
  groups <- combn(levels(data[, group]), 2)

  features <- Biobase::featureNames(object)

  results_df <- foreach::foreach(i = seq_along(features), .combine = rbind,
                                 .export = "finite_mean") %dopar% {
    feature <- features[i]
    result_row <- rep(NA_real_, ncol(groups))
    # Calculate fold changes
    tryCatch({
      for(i in 1:ncol(groups)){
        group1 <- data[data[, group] == groups[1,i], feature]
        group2 <- data[data[, group] == groups[2,i], feature]
        result_row[i] <- finite_mean(group2)/finite_mean(group1)
      }
    })

    result_row
  }

  # Create comparison labels for result column names
  comp_labels <- groups %>%
    t() %>%
    as.data.frame() %>%
    tidyr::unite("Comparison", V2, V1, sep = "_vs_")
  comp_labels <- paste0(comp_labels[,1], "_FC")
  results_df <- data.frame(features, results_df, stringsAsFactors = FALSE)
  colnames(results_df) <- c("Feature_ID", comp_labels)
  rownames(results_df) <- results_df$Feature_ID

  log_text("Fold changes computed.")

  # Order the columns accordingly
  results_df[c("Feature_ID", comp_labels[order(comp_labels)])]
}



#' Perform correlation tests
#'
#' Performs a correlation test between two sets of variables. All the variables must be either
#' feature names or column names of pheno data (sample information).
#' There are two ways to use this function:
#' either provide a set of variables as \code{x}, and all correlations between
#' those variables are computed. Or
#' provide two distinct sets of variables \code{x, y} and correlations between each x variable
#' and each y variable are computed.
#'
#' @param object a MetaboSet object
#' @param x character vector, names of variables to be correlated
#' @param y character vector, either identical to x (the default) or a distinct set of variables
#' to be correlated agains x
#' @param id character, column name for subject IDs. If provided, the correlation will be computed
#' using the rmcorr package
#' @param object2 optional second MeatboSet object. If provided, x variables will be taken from object and
#' y variables will be taken from object2. Both objects should have the same number of samples.
#' @param fdr logical, whether p-values from the correlation test should be adjusted with FDR correction
#' @param all_pairs logical, whether all pairs between x and y should be tested.
#' If FALSE, x and y give the exact pairs of variables to test, and should have the same length.
#' @param duplicates logical, whether correlations should be dublicated. If \code{TRUE}, each correlation
#' will be included in the results twice, where the order of the variables (which is x and which is y)
#' is changed. Can be useful for e.g. plotting a heatmap of the results, see examples of
#' \code{\link{plot_effect_heatmap}}
#' @param ... other parameters passed to \code{\link{cor.test}}, such as method
#'
#' @return a data frame with the results of correlation tests: the pair of variables,
#' correlation coefficient and p-value
#'
#' @examples
#' # Correlations between all features
#' correlations <- perform_correlation_tests(example_set, x = featureNames(example_set))
#'
#' # Spearman Correlations between features and sample information variables
#' # Drop QCs and convert time to numeric
#' no_qc <- drop_qcs(example_set)
#' no_qc$Time <- as.numeric(no_qc$Time)
#' correlations <- perform_correlation_tests(no_qc, x = featureNames(example_set),
#'                                          y = c("Time", "Injection_order"), method = "spearman")
#'
#' # Correlations between variables from two distinct MetaboSets
#' cross_object_cor <- perform_correlation_tests(hilic_neg_sample,
#'                                               x = featureNames(hilic_neg_sample),
#'                                               object2 = hilic_pos_sample,
#'                                               y = featureNames(hilic_pos_sample),
#'                                               all_pairs = FALSE)
#' @seealso \code{\link{cor.test}}, \code{\link[rmcorr]{rmcorr}}
#'
#' @importFrom foreach %do%
#' @export
perform_correlation_tests <- function(object, x, y = x, id = NULL, object2 = NULL, fdr = TRUE,
                                      all_pairs = TRUE, duplicates = FALSE, ...) {

  log_text("Starting correlation tests.")

  data1 <- combined_data(object)

  if (!is.null(object2)) {
    if (ncol(object) != ncol(object2)) {
      stop("The objects have different numbers of samples")
    }
    data2 <- combined_data(object2)
  } else {
    data2 <- data1
  }

  # Checks for repeated measures correlation
  if (!is.null(id)) {
    if (!requireNamespace("rmcorr", quietly = TRUE)) {
      stop("Package \"rmcorr\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (!id %in% colnames(data1) || !id %in% colnames(data2)) {
      stop("id column not found", call. = FALSE)
    }
    if (!identical(data1[, id], data2[, id])) {
      stop("ids do not match between the two objects: make sure the subjects are in the same order!",
           call. = FALSE)
    }
    add_citation("rmcorr package was used to compute correlations with repeated measuremtns:", citation("rmcorr"))
  }

  # All x and y should be columns names of combined data
  not_found <- setdiff(x, colnames(data1))
  not_found <- c(not_found, setdiff(y, colnames(data2)))
  if (length(not_found)) {
    stop(paste("Following variables do not match to know variables in the object(s):",
               paste(not_found, collapse = ", ")))
  }

  if (all_pairs){
    # If the same variable is present in x and y, the correlation would be computed
    # twice. This makes sure only unique combinations of variables are treated.
    if (identical(x,y)) {
      var_pairs <- combn(x, 2) %>% t() %>% data.frame(stringsAsFactors = FALSE)
      colnames(var_pairs) <- c("x", "y")
      # Add correlations of all variables with themselves (useful for plotting)
      var_pairs <- rbind(var_pairs, data.frame(x = x, y = x, stringsAsFactors = FALSE))
    } else if (is.null(object2) & length(intersect(x, y))) {
      stop("Currently only identical x & y or completely separate x & y are supported for one object")
    } else {
      var_pairs <- expand.grid(x, y, stringsAsFactors = FALSE)
      colnames(var_pairs) <- c("x", "y")
    }
  } else {
    if (length(x) != length(y)) {
      stop("If all_pairs = FALSE, x and y should have the same length")
    }
    var_pairs <- data.frame(x = x, y = y, stringsAsFactors = FALSE)
  }


  # Compute correlations
  cor_results <- foreach::foreach(i = seq_len(nrow(var_pairs)), .combine = rbind) %dopar% {
    x_tmp = var_pairs$x[i]
    y_tmp = var_pairs$y[i]
    cor_tmp <- NULL
    tryCatch({
      if (is.null(id)) {
        cor_tmp <- cor.test(data1[, x_tmp], data2[, y_tmp], ...)
      } else {
        id_tmp <- data1[, id]
        df_tmp <- data.frame(id_var = id_tmp, x_var = data1[, x_tmp], y_var = data2[, y_tmp])
        cor_tmp <- rmcorr::rmcorr(participant = id_var,
                                  measure1 = x_var,
                                  measure2 = y_var,
                                  dataset = df_tmp)
        cor_tmp <- list(estimate = cor_tmp$r, p.value = cor_tmp$p)
      }
    }, error = function(e) cat(paste0(x_tmp, " vs ", y_tmp, ": ", e$message, "\n")))
    if (is.null(cor_tmp)) {
      cor_tmp <- list(estimate = NA,
                      p.value = NA)
    }
    data.frame(X = x_tmp, Y = y_tmp,
               Correlation_coefficient = cor_tmp$estimate,
               Correlation_P = cor_tmp$p.value,
               stringsAsFactors = FALSE)

  }

  if (duplicates) {
    cor_results_dup <- cor_results
    cor_results_dup$X <- cor_results$Y
    cor_results_dup$Y <- cor_results$X
    # Remove possible duplicated correlations of a variable with itself
    cor_results_dup <- dplyr::filter(cor_results_dup, X != Y)
    cor_results <- rbind(cor_results, cor_results_dup)
  }

  # FDR correction
  if (fdr) {
    flags <- rep(NA_character_, nrow(cor_results))
    cor_results <- adjust_p_values(cor_results, flags)
  }

  rownames(cor_results) <- seq_len(nrow(cor_results))

  log_text("Correlation tests performed.")

  cor_results
}

#' Area under curve
#'
#' Compute area under curve (AUC) for each subject and feature.
#' Creates a pseudo MetaboSet object, where the "samples" are subjects
#' (or subject/group combinations in case the same subjects are submitted to different treatments)
#' and the "abundances" are AUCs. This object can then be used to compute results of e.g. t-tests of
#' AUCs between groups.
#'
#' @param object a MetaboSet object
#' @param time,subject,group column names of pData(object), holding time, subject and group labels
#'
#' @return a pseudo MetaboSet object with the AUCs
#'
#' @examples
#' # Drop QC samples before computing AUCs
#' aucs <- perform_auc(drop_qcs(example_set))
#' # t-test with the AUCs
#' t_test_results <- perform_t_test(aucs, formula_char =  "Feature ~ Group")
#'
#' @seealso \code{\link[PK]{auc}}
#'
#' @export
perform_auc <- function(object, time = time_col(object), subject = subject_col(object),
                        group = group_col(object)) {
  if (!requireNamespace("PK", quietly = TRUE)) {
      stop("Package \"PK\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  add_citation("PK package was used to compute AUC:", citation("PK"))

  log_text("Starting AUC computation.")

  data <- combined_data(object)

  # Create new pheno data, only one row per subject and group
  pheno_data <- data[, c(subject, group)] %>%
    dplyr::distinct() %>%
    tidyr::unite("Sample_ID", subject, group, remove = FALSE)
  # QC and Injection_order columns to pass validObject check
  pheno_data$QC <- "Sample"
  pheno_data$Injection_order <- seq_len(nrow(pheno_data))
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
  new_object <- construct_metabosets(exprs = aucs, feature_data = fData(object),
                                    pheno_data = pheno_data, group_col = group,
                                    subject_col = subject) %>%
    merge_metabosets()

  log_text("AUC computation finished.")
  new_object
}

# Helper function for FDR correction
adjust_p_values <- function(x, flags) {
  p_cols <- colnames(x)[grepl("_P$", colnames(x))]
  for (p_col in p_cols) {
    p_values <- x[, p_col, drop = TRUE]
    p_values[!is.na(flags)] <- NA
    x <- tibble::add_column(.data = x,
                            FDR = p.adjust(p_values, method = "BH"),
                            .after = p_col)
    p_idx <- which(colnames(x) == p_col)
    colnames(x)[p_idx + 1] <- paste0(p_col, "_FDR")
  }
  x
}

# Helper function for filling missing rows in results files with NAs
# Some statistical tests may fail for some features, due to e.g. missing values.
fill_results <- function(results_df, features) {
  # Add NA rows for features where the test failed
  results_df <- results_df %>% dplyr::select(Feature_ID, dplyr::everything())
  missing_features <- setdiff(features, results_df$Feature_ID)
  fill_nas <- matrix(NA, nrow = length(missing_features), ncol = ncol(results_df) - 1) %>%
    as.data.frame()
  results_fill <- data.frame(Feature_ID = missing_features, fill_nas)
  rownames(results_fill) <- missing_features
  colnames(results_fill) <- colnames(results_df)
  results_df <- rbind(results_df, results_fill) %>% as.data.frame()
  rownames(results_df) <- results_df$Feature_ID
  # Set Feature ID to the original order
  results_df <- results_df[features, ]
  results_df
}

# Helper function for running a variety of simple statistical tests
perform_test <- function(object, formula_char, result_fun, all_features, fdr = TRUE, packages = NULL) {

  data <- combined_data(object)
  features <- Biobase::featureNames(object)

  results_df <- foreach::foreach(i = seq_along(features), .combine = dplyr::bind_rows,
                                 .packages = packages) %dopar% {
    feature <- features[i]
    # Replace "Feature" with the current feature name
    tmp_formula <- gsub("Feature", feature, formula_char)
    # Run test
    result_row <- result_fun(feature = feature, formula = as.formula(tmp_formula), data = data)
    # In case Feature is used as predictor, make the column names match
    if (!is.null(result_row)){
      colnames(result_row) <- gsub(feature, "Feature", colnames(result_row))
    }
    result_row
  }
  # Check that results actually contain something
  # If the tests are run on parallel, the error messages from failing tests are not visible
  if (nrow(results_df) == 0) {
    stop("All the test failed, to see the individual error messages run the tests withot parallelization.",
         call. = FALSE)
  }
  # Rows full of NA for features where the test failed
  results_df <- fill_results(results_df, features)

  # FDR correction
  if (fdr) {
    if (all_features) {
      flags <- rep(NA_character_, nrow(results_df))
    } else {
      flags <- flag(object)
    }
    results_df <- adjust_p_values(results_df, flags)
  }

  results_df
}

#' Linear models
#'
#' Fits a linear model separately for each feature. Returns all relevant
#' statistics.
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
#' @details The linear model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting, see the example
#'
#' @examples
#' # A simple example without QC samples
#' # Features predicted by Group and Time
#' lm_results <- perform_lm(drop_qcs(example_set), formula_char = "Feature ~ Group + Time")
#'
#' @seealso \code{\link[stats]{lm}}
#'
#' @export
perform_lm <- function(object, formula_char, all_features = FALSE, ci_level = 0.95, ...) {

  log_text("Starting linear regression.")

  lm_fun <- function(feature, formula, data) {
    # Try to fit the linear model
    fit <- NULL
    tryCatch({
      fit <- lm(formula, data = data, ...)
    }, error = function(e) cat(paste0(feature, ": ", e$message, "\n")))
    if(is.null(fit) | sum(!is.na(data[, feature])) < 2){
      result_row <- NULL
    } else {
      # Gather coefficients and CIs to one data frame row
      coefs <- summary(fit)$coefficients
      confints <- confint(fit, level = ci_level)
      coefs <- data.frame(Variable = rownames(coefs), coefs, stringsAsFactors = FALSE)
      confints <- data.frame(Variable = rownames(confints), confints, stringsAsFactors = FALSE)

      result_row <- dplyr::left_join(coefs,confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "t_value" ="t.value",
                      "P" = "Pr...t..", "LCI95" = "X2.5..", "UCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -Variable) %>%
        tidyr::unite("Column", Variable, Metric, sep="_") %>%
        tidyr::spread(Column, Value)
      # Add R2 statistics and feature ID
      result_row$R2 <- summary(fit)$r.squared
      result_row$Adj_R2 <- summary(fit)$adj.r.squared
      result_row$Feature_ID <- feature
    }
    result_row
  }

  results_df <- perform_test(object, formula_char, lm_fun, all_features)

  # Set a good column order
  variables <- gsub("_P$", "", colnames(results_df)[grep("P$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error", "t_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, variables, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", Var2, Var1)
  col_order <- c("Feature_ID", col_order$Column, c("R2", "Adj_R2"))


  log_text("Linear regression performed.")

  results_df[col_order]
}


#' Logistic regression
#'
#' Fits a logistic regression model separately for each feature. Returns all relevant
#' statistics.
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param ci_level the confidence level used in constructing the confidence intervals
#' for regression coefficients
#' @param ... additional parameters passed to glm
#'
#' @return a data frame with one row per feature, with all the
#' relevant statistics of the linear model as columns
#'
#' @details The logistic regression model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting, see the example
#'
#' @examples
#' # A simple example without QC samples
#' # Time predicted by features
#' logistic_results <- perform_logistic(drop_qcs(example_set),
#'                                      formula_char = "Time ~ Feature + Group ")
#'
#' @seealso \code{\link[stats]{glm}}
#'
#' @export
perform_logistic <- function(object, formula_char, all_features = FALSE, ci_level = 0.95, ...) {

  log_text("Starting logistic regression.")

  logistic_fun <- function(feature, formula, data) {
    # Try to fit the linear model
    fit <- NULL
    tryCatch({
      fit <- glm(formula, data = data, family = binomial(), ...)
    }, error = function(e) cat(paste0(feature, ": ", e$message, "\n")))
    if(is.null(fit) | sum(!is.na(data[, feature])) < 2){
      result_row <- NULL
    } else {
      # Gather coefficients and CIs to one data frame row
      coefs <- summary(fit)$coefficients
      suppressMessages(confints <- confint(fit, level = ci_level))
      coefs <- data.frame(Variable = rownames(coefs), coefs, stringsAsFactors = FALSE)
      confints <- data.frame(Variable = rownames(confints), confints, stringsAsFactors = FALSE)

      result_row <- dplyr::left_join(coefs,confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "z_value" ="z.value",
                      "P" = "Pr...z..", "LCI95" = "X2.5..", "UCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -Variable) %>%
        tidyr::unite("Column", Variable, Metric, sep="_") %>%
        tidyr::spread(Column, Value)
      result_row$Feature_ID <- feature
    }
    result_row
  }

  results_df <- perform_test(object, formula_char, logistic_fun, all_features)

  # Set a good column order
  variables <- gsub("_P$", "", colnames(results_df)[grep("P$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error", "z_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, variables, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", Var2, Var1)
  col_order <- c("Feature_ID", col_order$Column)
  results_df <- results_df[col_order]
  # Add odds ratios
  estimate_cols <- colnames(results_df)[grepl("_Estimate$", colnames(results_df))]
  for (estimate_col in estimate_cols) {
    estimate_values <- results_df[, estimate_col]
    results_df <- tibble::add_column(.data = results_df,
                            OR = exp(estimate_values),
                            .after = estimate_col)
    estimate_idx <- which(colnames(results_df) == estimate_col)
    colnames(results_df)[estimate_idx + 1] <- gsub("Estimate", "OR", estimate_col)
  }


  log_text("Logistic regression performed.")

  results_df
}


#' Linear mixed models
#'
#' Fits a linear mixed model separately for each feature. Returns all relevant
#' statistics.
#' \strong{CITATION:} When using this function, cite \code{lme4} and \code{lmerTest} packages
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' @param all_features should all features be included in FDR correction?
#' @param ci_level the confidence level used in constructing the confidence intervals
#' for regression coefficients
#' @param ci_method The method for calculating the confidence intervals, see documentation
#' of confint below
#' @param test_random logical, whether tests for the significance of the random effects
#' should be performed
#' @param ... additional parameters passed to lmer
#'
#' @return a data frame with one row per feature, with all the
#' relevant statistics of the linear mixed model as columns
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting, see the example
#'
#'
#' @examples
#' # A simple example without QC samples
#' # Features predicted by Group and Time as fixed effects with Subject ID as a random effect
#' \dontrun{
#' lmer_results <- perform_lmer(drop_qcs(example_set),
#'                              formula_char = "Feature ~ Group + Time + (1 | Subject_ID)",
#'                         ci_method = "Wald")
#' }
#' @seealso \code{\link[lmerTest]{lmer}} for model specification and
#' \code{\link[lme4]{confint.merMod}} for the computation of confidence intervals
#'
#' @export
perform_lmer <- function(object, formula_char, all_features = FALSE,  ci_level = 0.95,
                         ci_method = c("Wald", "profile", "boot"),
                         test_random = FALSE, ...) {

  log_text("Starting fitting linear mixed models.")

  if (!requireNamespace("lmerTest", quietly = TRUE)) {
      stop("Package \"lmerTest\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("MuMIn", quietly = TRUE)) {
      stop("Package \"MuMIn\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  add_citation("lme4 package was used to fit linear mixed models:", citation("lme4"))
  add_citation("lmerTest package was used for statistical tests in linear mixed models:", citation("lmerTest"))
  add_citation("MuMIn package was used to assess R2 values in linear mixed models:", citation("MuMIn"))

  # Check that ci_method is one of the accepted choices
  ci_method <- match.arg(ci_method)

  lmer_fun <- function(feature, formula, data) {
    # Set seed, needed for some of the CI methods
    set.seed(38)

    # Try to fit the linear model
    fit <- NULL
    # If fitting causes an error, a NULL row is returned
    result_row <- NULL
    tryCatch({
      fit <- lmerTest::lmer(formula, data = data, ...)
    },error = function(e) cat(paste0(feature, ": ", e$message, "\n")))
    if (!is.null(fit)) {
      # Extract model coefficients
      coefs <- summary(fit)$coefficients
      coefs <- data.frame(Variable = rownames(coefs), coefs, stringsAsFactors = FALSE)
      # Try to compute confidence intervals
      # If the computation fails, all CIs are NA
      confints <- data.frame(Variable = rownames(coefs), "X2.5.." = NA, "X97.5.." = NA)
      tryCatch({
        confints <- confint(fit, nsim = 1000, method = ci_method, oldNames = FALSE)
        confints <- data.frame(Variable = rownames(confints), confints, stringsAsFactors = FALSE)
      },error = function(e) cat(paste0(feature, ": ", e$message, "\n")))

      # Gather coefficients and CIs to one data frame row
      result_row <- dplyr::left_join(coefs,confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "t_value" ="t.value",
                      "P" = "Pr...t..", "LCI95" = "X2.5..", "LCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -Variable) %>%
        tidyr::unite("Column", Variable, Metric, sep="_") %>%
        tidyr::spread(Column, Value)
      # Add R2 statistics
      result_row$Marginal_R2 <- NA
      result_row$Conditional_R2 <- NA
      tryCatch({
        R2s <- suppressWarnings(MuMIn::r.squaredGLMM(fit))
        result_row$Marginal_R2 <- R2s[1]
        result_row$Conditional_R2 <- R2s[2]
      },error = function(e) cat(paste0(feature, ": ", e$message, "\n")))
      # Add Feature ID
      result_row$Feature_ID <- feature
    }

    # Add optional test results for the random effects
    if(test_random) {
      tryCatch({
        r_tests <- as.data.frame(ranova(fit))[-1,c(4,6)]
        r_tests$Variable <- rownames(r_tests) %>%
          gsub("[(]1 [|] ", "", .) %>% gsub("[)]", "", .)
        # Get confidence intervals for the standard deviations of the random effects
        confints$Variable <- confints$Variable  %>%
          gsub("sd_[(]Intercept[)][|]", "", .)
        # Get standard deviations of the random effects
        r_variances <- as.data.frame(summary(fit)$varcor)[c("grp", "sdcor")]
        # Join all the information together
        r_result_row <- dplyr::inner_join(r_variances, confints, by = c("grp" = "Variable")) %>%
          dplyr::left_join(r_tests,  by = c("grp" = "Variable")) %>%
          dplyr::rename(SD = sdcor, "LCI95" = "X2.5..", "UCI95" = "X97.5..", "P" = "Pr(>Chisq)") %>%
          tidyr::gather("Metric", "Value", -grp) %>%
          tidyr::unite("Column", grp, Metric, sep="_") %>%
          tidyr::spread(Column, Value)
        result_row <- cbind(result_row, r_result_row)
      },error = function(e) cat(paste0(feature, ": ", e$message, "\n")))
    }

    result_row
  }

  results_df <- perform_test(object, formula_char, lmer_fun, all_features, packages = "lmerTest")

  # Set a good column order
  fixed_effects <- gsub("_Estimate$", "", colnames(results_df)[grep("Estimate$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error", "t_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, fixed_effects, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", Var2, Var1)
  col_order <- c("Feature_ID", col_order$Column, c("Marginal_R2", "Conditional_R2"))

  if (test_random) {
    random_effects <- gsub("_SD$", "", colnames(results_df)[grep("SD$", colnames(results_df))])
    statistics <- c("SD", "LCI95", "UCI95", "LRT", "P", "P_FDR")
    random_effect_order <- expand.grid(statistics, random_effects, stringsAsFactors = FALSE) %>%
      tidyr::unite("Column", Var2, Var1)
    col_order <- c(col_order, random_effect_order$Column)
  }

  log_text("Linear mixed models fit.")

  results_df[col_order]
}

#' Test homoscedasticity
#'
#' Performs Bartlett's, Levene's and Fligner-Killeen tests for equality of variances
#' \strong{CITATION:} When using this function, cite the \code{car} package, which
#' provides the function for Levene's test. (The other tests are included in base R).
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' Defaults to "Feature ~ group_col(object)
#' @param all_features should all features be included in FDR correction?
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting. For example, if testing for equality of
#' variances in study groups, use "Feature ~ Group".
#'
#' @return data frame with the results
#'
#' @examples
#' perform_homoscedasticity_tests(example_set, formula_char = "Feature ~ Group")
#'
#' @seealso \code{\link{bartlett.test}}, \code{\link[car]{leveneTest}}, \code{\link{fligner.test}}
#'
#' @export
perform_homoscedasticity_tests <- function(object, formula_char, all_features = FALSE) {
  if (!requireNamespace("car", quietly = TRUE)) {
      stop("Package \"car\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  add_citation("car package was used for Levene's test of homoscedasticity:", citation("car"))

  log_text("Starting homoscedasticity tests.")

  homosced_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch({
      bartlett <- bartlett.test(formula = formula, data = data)
      levene <- car::leveneTest(y = formula, data = data)
      fligner <- fligner.test(formula = formula, data = data)

      result_row <- data.frame(Feature_ID = feature,
                               Bartlett_P = bartlett$p.value,
                               Levene_P = levene$`Pr(>F)`[1],
                               Fligner_P = fligner$p.value,
                               stringsAsFactors = FALSE)
    }, error = function(e) {cat(paste0(feature, ": ", e$message, "\n"))})

    result_row

  }

  results_df <- perform_test(object, formula_char, homosced_fun, all_features)

  log_text("Homoscedasticity tests performed.")

  results_df
}

#' Perform Kruskal-Wallis Rank Sum Tests
#'
#' Performs Kruskal-Wallis Rank Sum Test for equality
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' Defaults to "Feature ~ group_col(object)
#' @param all_features should all features be included in FDR correction?
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting. For example, if testing for equality of
#' means in study groups, use "Feature ~ Group".
#'
#' @return data frame with the results
#'
#' @seealso \code{\link{kruskal.test}}
#'
#' @examples
#' perform_kruskal_wallis(example_set, formula_char = "Feature ~ Group")
#'
#' @export
perform_kruskal_wallis <- function(object, formula_char, all_features = FALSE) {

  log_text("Starting Kruskal_wallis tests.")

  kruskal_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch({
      kruskal <- kruskal.test(formula = formula, data = data)

      result_row <- data.frame(Feature_ID = feature,
                               Kruskal_P = kruskal$p.value,
                               stringsAsFactors = FALSE)
    }, error = function(e) {cat(paste0(feature, ": ", e$message, "\n"))})

    result_row
  }

  results_df <- perform_test(object, formula_char, kruskal_fun, all_features)

  log_text("Kruskal_wallis tests performed.")

  results_df
}


#' Perform ANOVA
#'
#' Performs ANOVA with Welch's correction as default, to deal with heterogeneity of variances.
#' Can also perform classic ANOVA with assumption of equal variances.
#' Uses base R function \code{oneway.test}.
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' Defaults to "Feature ~ group_col(object)
#' @param all_features should all features be included in FDR correction?
#' @param ... other parameters to \code{\link{oneway.test}}
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting. For example, if testing for equality of
#' means in study groups, use "Feature ~ Group".
#'
#' @return data frame with the results
#'
#' @seealso \code{\link{oneway.test}}
#'
#' @examples
#' perform_oneway_anova(example_set, formula_char = "Feature ~ Group")
#'
#' @export
perform_oneway_anova <- function(object, formula_char, all_features = FALSE, ...) {

  log_text("Starting ANOVA tests.")

  anova_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch({
      anova_res <- oneway.test(formula = formula, data = data, ...)

      result_row <- data.frame(Feature_ID = feature,
                               ANOVA_P = anova_res$p.value,
                               stringsAsFactors = FALSE)
    }, error = function(e) {cat(paste0(feature, ": ", e$message, "\n"))})

    result_row

  }

  results_df <- perform_test(object, formula_char, anova_fun, all_features)

  log_text("ANOVA performed.")

  results_df
}

#' Perform t-tests
#'
#' Performs t-tests, the R default is Welch's t-test (unequal variances), use var.equal = TRUE
#' for Student's t-test
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' Defaults to "Feature ~ group_col(object)
#' @param all_features should all features be included in FDR correction?
#' @param ... additional parameters to t.test
#'
#' @details The model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting. For example, if testing for equality of
#' means in study groups, use "Feature ~ Group".
#'
#' @return data frame with the results
#'
#' @examples
#' t_test_results <- perform_t_test(drop_qcs(merged_sample), formula_char = "Feature ~ Group")
#'
#' @seealso \code{\link{t.test}}
#'
#' @export
perform_t_test <- function(object, formula_char, all_features = FALSE, ...) {
  message(paste0("Remember that t.test returns difference between group means",
                 "in different order than lm.\n",
                 "This function mimics this behavior, so the effect size is",
                 " mean of reference level minus mean of second level."))


  log_text("Starting t-tests.")

  t_fun <- function(feature, formula, data) {
    result_row <- NULL
    tryCatch({
      t_res <- t.test(formula = formula, data = data, ...)

      conf_level <- attr(t_res$conf.int, "conf.level") * 100

      result_row <- data.frame(Feature_ID = feature,
                               Mean1 = t_res$estimate[1],
                               Mean2 = t_res$estimate[2],
                               Mean_1_minus_2 = t_res$estimate[1] - t_res$estimate[2],
                               "Lower_CI_" = t_res$conf.int[1],
                               "Upper_CI_" = t_res$conf.int[2],
                               t_test_P = t_res$p.value,
                               stringsAsFactors = FALSE)
      colnames(result_row)[5:6] <- paste0(colnames(result_row)[5:6], conf_level)
    }, error = function(e) {cat(paste0(feature, ": ", e$message, "\n"))})

    result_row
  }

  results_df <- perform_test(object, formula_char, t_fun, all_features)

  log_text("t-tests performed.")

  results_df
}

#' Perform paired t-tests
#'
#' Performs paired t-tests between two groups.
#'
#' @param object a MetaboSet object
#' @param group character, column name of pData with the group information
#' @param id character, column name of pData with the identifiers for the pairs
#' @param all_features should all features be included in FDR correction?
#' @param ... additional parameters to t.test
#'
#' @return data frame with the results
#'
#' @examples
#' paired_t_results <- perform_paired_t_test(drop_qcs(example_set), group = "Time", id = "Subject_ID")
#'
#' @seealso \code{\link{t.test}}
#'
#' @export
perform_paired_t_test <- function(object, group, id, all_features = FALSE, ...) {

  log_text("Starting paired t-tests.")

  data <- combined_data(object)
  groups <- data[, group]
  if (class(groups) != "factor") {
    groups <- as.factor(groups)
  }
  if (length(levels(groups)) > 2) {
    warning("More than two groups detected, only using the first two", call. = FALSE)
  }

  # Split to groups
  group1 <- data[which(groups == levels(groups)[1]), ]
  group2 <- data[which(groups == levels(groups)[2]), ]
  # Keep only complete pairs, order by id
  common_ids <- intersect(group1[, id], group2[, id])
  group1 <- group1[group1[, id] %in% common_ids, ][order(common_ids), ]
  group2 <- group2[group2[, id] %in% common_ids, ][order(common_ids), ]
  log_text(paste("Found", length(common_ids), "complete pairs"))

  features <- featureNames(object)
  results_df <- foreach::foreach(i = seq_along(features), .combine = dplyr::bind_rows) %dopar% {
    feature <- features[i]
    result_row <- NULL
    tryCatch({
      t_res <- t.test(group1[, feature], group2[, feature], paired = TRUE, ...)

      conf_level <- attr(t_res$conf.int, "conf.level") * 100

      result_row <- data.frame(Feature_ID = feature,
                               Mean_diff = t_res$estimate[1],
                               "Lower_CI_" = t_res$conf.int[1],
                               "Upper_CI_" = t_res$conf.int[2],
                               t_test_P = t_res$p.value,
                               stringsAsFactors = FALSE)
      colnames(result_row)[3:4] <- paste0(colnames(result_row)[3:4], conf_level)
    }, error = function(e) {cat(paste0(feature, ": ", e$message, "\n"))})

    result_row
  }

  # Check that results actually contain something
  # If the tests are run on parallel, the error messages from failing tests are not visible
  if (nrow(results_df) == 0) {
    stop("All the test failed, to see the individual error messages run the tests withot parallelization.",
         call. = FALSE)
  }
  rownames(results_df) <- results_df$Feature_ID
  colnames(results_df)[2] <- paste0("Mean_diff_", levels(groups)[1],
                                    "_minus_", levels(groups)[2])
  # Rows full of NA for features where the test failed
  results_df <- fill_results(results_df, features)

  # FDR correction
  if (all_features) {
    flags <- rep(NA_character_, nrow(results_df))
  } else {
    flags <- flag(object)
  }
  results_df <- adjust_p_values(results_df, flags)

  log_text("Paired t-tests performed.")

  results_df
}

#' Perform pairwise t-tests
#'
#' Performs pairwise t-tests between all study groups. NOTE! Does not use formula interface
#'
#' @param object a MetaboSet object
#' @param group character, column name of phenoData giving the groups
#' @param all_features should all features be included in FDR correction?
#' @param ... other parameters passed to perform_t_test, and eventually to base R t.test
#'
#' @details P-values of each comparison are corrected separately from each other.
#'
#' @return data frame with the results
#'
#' @examples
#' #Including QCs as a study group for example
#' t_test_results <- perform_pairwise_t_test(merged_sample, group = "Group")
#'
#' @seealso \code{\link{perform_t_test}}, \code{\link{t.test}}
#'
#' @export
perform_pairwise_t_test <- function(object, group = group_col(object), all_features = FALSE, ...) {

  if (!is.factor(pData(object)[, group])) {
    stop("Group column should be a factor")
  }

  log_text("Starting pairwise t-tests.")

  if (class(pData(object)[, group]) == "factor") {
    groups <- levels(pData(object)[, group])
  } else {
    groups <- unique(pData(object)[, group])
  }
  combinations <- combn(groups, 2)

  for (i in seq_len(ncol(combinations))) {
    group1 <- as.character(combinations[1, i])
    group2 <- as.character(combinations[2, i])
    # Subset the pair of groups
    object_tmp <- object[, pData(object)[, group] %in% c(group1, group2)]
    pData(object_tmp) <- droplevels(pData(object_tmp))

    t_results <- perform_t_test(object_tmp, formula_char = paste("Feature ~", group), all_features)
    colnames(t_results) <- c("Feature_ID", paste0("Mean_", c(group1, group2)),
                          paste0("Mean_", group1, "_minus_", group2),
                          paste0(paste0(group1, "_", group2, "_"), colnames(t_results)[5:8]))

    if (i == 1) {
      results_df <- t_results
    } else {
      results_df <- dplyr::left_join(results_df, t_results,
                                     by = intersect(colnames(results_df), colnames(t_results)))
    }

  }

  log_text("Pairwise t-tests performed.")

  results_df
}
