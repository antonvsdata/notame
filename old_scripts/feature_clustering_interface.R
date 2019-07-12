


#' Cluster correltated features
#'
#'
#' Clusters features potentially originating from the same compound. Features with high
#' Pearson correlation coefficient and small retention time difference are linked together.
#' Then clusters are formed by setting a threshold for the relative degree that each node in a clusters needs
#' to fulfill. This is a wrapper around numerous functions that are based on the MATLAB code by David Broadhurst.
#'
#' @param object a MetaboSet object
#' @param rt_window the retention time window for potential links NOTE: use the samu unit as the retention time
#' @param corr_thresh the correlation threshold required for potential links between features
#' @param d_thresh the threshold for the relative degree required by each node
#' @param mz_col the column name in fData(object) that holds mass-to-charge ratios
#' @param rt_col the column name in fData(object) that holds retention times
#' @param plotting should plots be drawn for each cluster?
#' @param min_size_plotting the minimum number of features a cluster needs to have to be plotted
#' @param prefix the prefix to the files to be plotted
#'
#' @return a MetaboSet object, with the cluster ID added to results
#'
#' @examples
#' # The parameters are really weird because example data is imaginary
#' clustered <- cluster_features(example_set, rt_window = 1,)
cluster_features <- function(object, rt_window = 1/60, corr_thresh = 0.9, d_thresh = 0.8, mz_col= NULL, rt_col = NULL,
                             plotting = FALSE, min_size_plotting = 3, prefix = NULL) {

  # Start log
  log_text(paste("\nStarting feature clustering at", Sys.time()))

  # Find connections between features
  conn <- find_connections(data = as.data.frame(t(exprs(object))),
                           features = fData(object),
                           corr_thresh = corr_thresh,
                           rt_window = rt_window,
                           name_col = "Feature_ID",
                           mz_col = mz_col,
                           rt_col = rt_col)
  log_text(paste("Found", nrow(conn), "connections"))

  # Form clusters
  clusters <- find_clusters(conn, d_thresh)
  lens <- sapply(clusters, function(x){length(x$features)})
  log_text(paste("Found", sum(lens > 1), "clusters of 2 or more features"))

  features <- assign_cluster_id(clusters, fData(object), "Feature_ID")

  if (plotting) {
    features$MPA <- apply(exprs(object), 1, finite_median)
    visualize_clusters(data = as.data.frame(t(exprs(object))),
                       features = features, clusters = clusters,
                       min_size = min_size_plotting,
                       rt_window = rt_window, name_col = "Feature_ID",
                       mz_col = mz_col, rt_col = rt_col, file_path = prefix)
  }

  object <- join_results(object, features[c("Feature_ID", "Cluster_ID")])

  object
}

#' Assign Cluster ID to features
#'
#'
assign_cluster_id <- function(clusters, features, name_col) {

  features$Cluster_ID <- 0
  n_clust <- length(clusters)
  for (i in seq_len(n_clust)) {
    tmp_features <- clusters[[i]]$features
    print(tmp_features)
    features$Cluster_ID[features[, name_col] %in% tmp_features] <- i
  }
  # Assign cluster ID for leftalone features
  miss_idx <- features$Cluster_ID == 0
  features$Cluster_ID[miss_idx] <- seq(n_clust + 1, n_clust + sum(miss_idx))

  features
}
