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
  results <- data.frame()
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
