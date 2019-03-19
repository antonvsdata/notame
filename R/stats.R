#' Cohen's D
#'
#' Computes Cohen's D for each feature
#'
#' @param object a MetaboSet object
#' @param id character, name of the subject ID column
#' @param group character, name of the group column
#' @param time character, name of the time column
#'
#' @return data frame with Cohen's d for each feature
#'
#' @export
cohens_d <- function(object, id = subject_col(object), group = group_col(object),
                     time = time_col(object)) {

  data = combined_data(object)
  data[time] <- ifelse(data[, time] == levels(data[, time])[1], "time1", "time2")
  data[group] <- ifelse(data[, group] == levels(data[, group])[1], "group1", "group2")

  features <- Biobase::featureNames(object)
  ds <- foreach::foreach(i = seq_along(features), .combine = rbind, .packages = c("dplyr", "tidyr")) %dopar% {
    feature <- features[i]
    tmp <- data[c(id, group, time, feature)]
    colnames(tmp) <- c("ID", "group", "time", "feature")
    tmp <- tmp %>%
      spread(time, feature) %>%
      mutate(diff = time2 - time1) %>%
      group_by(group) %>%
      dplyr::summarise(mean_diff = mean(diff, na.rm = TRUE), sd_diff = sd(diff, na.rm = TRUE))

    d <- data.frame(Feature_ID = feature, Cohen_d = (tmp$mean_diff[tmp$group == "group2"] - tmp$mean_diff[tmp$group == "group1"]) / mean(tmp$sd_diff))
    d
  }
  rownames(ds) <- ds$Feature_ID
  ds
}

#' Fold change
#'
#' Computes fold change between eeach group for each feature.
#'
#' @param object a MetaboSet object
#' @param group character, name of the group column
#'
#' @return data frame with fold changes for each feature
#'
#' @export
fold_change <- function(object, group = group_col(object)) {

  data <- combined_data(object)
  groups <- combn(levels(data[, group]), 2)

  features <- Biobase::featureNames(object)

  results <- foreach::foreach(i = seq_along(features), .combine = rbind) %dopar% {
    feature <- features[i]
    result_row <- rep(0, ncol(groups))
    # Calculate fold changes
    for(i in 1:ncol(groups)){
      group1 <- data[data[, group] == groups[1,i], feature]
      group2 <- data[data[, group] == groups[2,i], feature]
      result_row[i] <- mean(group2)/mean(group1)
    }
    result_row
  }

  # Create comparison labels for result column names
  comp_labels <- groups %>% t() %>% as.data.frame() %>% unite("Comparison", V2, V1, sep = "_vs_")
  comp_labels <- comp_labels[,1]
  results_df <- data.frame(features, results, stringsAsFactors = FALSE)
  colnames(results_df) <- c("Feature_ID", comp_labels)
  rownames(results_df) <- results_df$Feature_ID
  # Order the columns accordingly
  results_df[c("Feature_ID", comp_labels[order(comp_labels)])]
}


#' Linear models
#'
#' Fits a linear model separately for each feature. Returns all relevant
#' statistics.
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' @param ci_level the confidence level used in constructing the confidence intervals
#' for regression coefficients
#' @param ... additional parameters passed to lm function
#'
#' @return a data frame with one row per feature, with all the
#' relevant statistics of the linear model as columns
#'
#' @details The linear model is fit on combined_data(object). Thus, column names
#' in pData(object) can be specified. To make the formulas flexible, the word "Feature"
#' must be used to signal the role of the features in the formula. "Feature" will be replaced
#' by the actual Feature IDs during model fitting, see the example
#'
#' @example
#' # A simple example without QC samples
#' # Features predicted by Group and Time
#' perform_lm(example_set[, example_set$QC != "QC"], formula_char = "Feature ~ Group + Time")
#'
#' @seealso \code{\link[stats]{lm}}
perform_lm <- function(object, formula_char,  ci_level = 0.95, ...) {

  data <- combined_data(object)
  features <- Biobase::featureNames(object)

  results <- foreach::foreach(i = seq_along(features), .combine = rbind,
                              .packages = c("dplyr", "tidyr")) %dopar% {
    feature <- features[i]
    # Replace "Feature" with the current feature name
    tmp_formula <- gsub("Feature", feature, formula_char)

    # Try to fit the linear model
    fit <- NULL
    tryCatch({
      fit <- lm(as.formula(tmp_formula), data = data, ...)
    }, error = function(e) print(e$message))
    if(is.null(fit) | sum(!is.na(data[, feature])) < 2){
      result_row <- NULL
    } else {
      # Gather coefficients and CIs to one data frame row
      coefs <- summary(fit)$coefficients
      confints <- confint(fit, level = ci_level)
      coefs <- data.frame(Variable = rownames(coefs), coefs, stringsAsFactors = FALSE)
      confints <- data.frame(Variable = rownames(confints), confints, stringsAsFactors = FALSE)

      result_row <- dplyr::left_join(coefs,confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "t_value" ="t.value", "P" = "Pr...t..", "LCI95" = "X2.5..", "UCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -Variable) %>%
        tidyr::unite("Column", Variable, Metric, sep="_") %>%
        tidyr::spread(Column, Value)
      # Add R2 statistics
      result_row$R2 <- summary(fit)$r.squared
      result_row$Adj_R2 <- summary(fit)$adj.r.squared

    }
    result_row

  }

  # FDR correction per column
  p_cols <- colnames(results)[grep("P$", colnames(results))]
  for (p_col in p_cols) {
    results[paste0(p_col, "_FDR")] <- p.adjust(results[, p_col], method = "BH")
  }

  # Set a good column order
  variables <- gsub("_P$", "", p_cols)
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error", "t_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, variables, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", Var2, Var1)
  col_order <- c(col_order$Column, c("R2", "Adj_R2"))

  results[col_order]
}



