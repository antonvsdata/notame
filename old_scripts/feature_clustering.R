
#' Find out which features are correlated within a specified retention time window
#' 
#' A part of the peak clustering algorithm. Iterates over all possible pairs of features
#' and records a connection between them if a) they have a Pearson correlation coefficient
#' higher than \code{corr_thresh} and b) their retention time difference is less than 
#' \code{rt_window}
#' 
#' @param data data frame with the abundaces (exprs(object))
#' @param feautres data frame with feature information
#' @param corr_thresh numeric, the threshhold of correlation to use in linking features
#' @param t_window numeric, the retention time window to use in linking features. NOTE you
#' need to use the same unit as in the retention time column
#' @param name_col character, name of the column in features that contains feature names
#' @param mz_col character, name of the column in features that contains mass-to-charge ratios
#' @param rt_col character, name of the column in features that contains retention times
#' 
#' @return a data frame of pairs of signals that are linked together
#' \begin{itemize}
#' \item x & y: indexes and names of the signals
#' \item cor: correlation coefficient
#' \item mz_diff & rt_diff: mass and retention time difference
#' \end{itemize}
#' 
#' @export 
find_connections <- function(data, features, corr_thresh = 0.9, rt_window = 1/60,
                             name_col, mz_col, rt_col) {
  
  D <- data[features[,name_col]]
  if (ncol(D) < 2) {
    stop("Need at least 2 signals to do any clustering!")
  }
  n <- nrow(features)
  connections <- foreach::foreach(i = seq_len(n-1), .combine = rbind) %dopar% {
    if (i %% 100 == 0){
      print(i)
    }
    connections_tmp <- data.frame()
    for (j in (i+1):n){
      rt_diff <- features[j, rt_col] - features[i, rt_col]
      cor_coef <- cor(D[, i], D[, j], use = "na.or.complete")
      if (!is.na(cor_coef)) {
        if (abs(rt_diff) < rt_window & cor_coef > corr_thresh){
          mz_diff <- features[j, mz_col] - features[i, mz_col]
          connections_tmp <- rbind(connections_tmp, data.frame(x = features[i, name_col], y = features[j, name_col],
                                                               cor = cor_coef, rt_diff = rt_diff, mz_diff = mz_diff))
        }
      }
    }
    connections_tmp
  }
  connections
}


#' Extract the densely connected clusters
#' 
#' First forms clusters of compounds that are linked together. THen the clusters are pruned
#' so that in the final clusters, each feature is linked to at least a set percentage
#' of the other features in th e cluster.
#'
#' @param connections data frame of pairs of signals that are linked together,
#' output of find_connections
#' @param d_thresh numeric, the minimum degree required for each signal in a cluster
#' expressed as a percentage of the maximum degree in the cluster
#' 
#' @return a lis of clusters, each a list of:
#' \begin{itemize}
#' \item features: character vector of the names of the features included in the cluster
#' \item graph: an igraph object of the cluster
#' \end{itemize}
#' 
#' @export
find_clusters <- function(connections, d_thresh = 0.8){
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("The igraph package is required for this function")
  }
  
  # Construct graph from the given edges
  g <- igraph::graph_from_edgelist(as.matrix(connections[1:2]), directed = FALSE)
  g <- igraph::set.edge.attribute(graph = g, name = "weight", value = connections$cor)
  
  # Initialize list of clusters
  clusters <- list()
  k <- 1
  
  # Repeatedly extract densely connected clusters from the graph
  while(length(igraph::V(g))){
    # Connected components of the remaining graph
    comp <- igraph::decompose(g)
    n_comp <- length(comp)
    cat(paste(n_comp, "components found\n\n"))
    
    # Only keep the densely connected part of each component (subgraph)
    clusters_tmp <- foreach(i = seq_len(n_comp), .combine = c) %dopar% {
      
      if (i %% 100 == 0) {
        cat(paste("Component", i, "/", n_comp,"\n"))
      }
      subg <- comp[[i]]
      
      n_nodes <- length(igraph::V(subg))
      d <- igraph::degree(subg)
      # The limit of the degree a node needs to be kept
      d_lim <- round(d_thresh * (n_nodes-1))
      
      # 
      if(n_nodes >= 3) {
        # Remove the node with the smallest degree until all nodes in the cluster have
        # a degree above the limit
        while(any(d < d_lim)){
          idx <- which(d == min(d))
          if (length(idx) > 1) {
            edgesums <- igraph::strength(subg, vids = igraph::V(subg)$name[idx])
            idx <- idx[which(edgesums == min(edgesums))[1]]
          }
          subg <- igraph::delete.vertices(subg, v = igraph::V(subg)[idx])
          d <- igraph::degree(subg)
          n_nodes <- n_nodes - 1
          d_lim <- round(d_thresh * (n_nodes-1))
        }
      }
      
      # Record the final cluster and remove the nodes from the main graph
      list(list(features = names(igraph::V(subg)),
                graph = subg))
      
    }
    
    for (j in seq_along(clusters_tmp)) {
      clusters[[k]] <- clusters_tmp[[j]]
      k <- k + 1
      subg <- clusters_tmp[[j]]$graph
      g <- igraph::delete.vertices(g, v = names(igraph::V(subg)))
    }
    
  }
  
  clusters
}

#' Extract information of the features in clusters
#'
#' The LC-MS data of the feature with largest median peak area is retained,
#' all the features in every cluster are recorded
#'
#' @param clusters list of cluster information, as returned by find_clusters
#' @param data data frame of the original LC-MS data
#' @param features data frame holding the feature information
#' @param name_col name_col character, name of the column in features that contains feature names
#' 
#' @return a list of two items:
#' \begin{itemize}
#' \item cdata: a new data frame with the combined LC-MS data
#' \item cfeatures: data frame, feature information per cluster
#' \end{itemize}
#' 
#' @export
pull_features <- function(clusters, data, features,
                          name_col){
  
  # Median peak area
  features$MPA <- sapply(data[features[, name_col]], median, na.rm = TRUE)
  
  cfeatures <- data.frame()
  sample_cols <- setdiff(colnames(data), features[, name_col])
  cdata <- data[sample_cols]
  handled_features <- c()
  
  n_clusters <- length(clusters)
  # Retain the strongest signal (MPA) from each cluster 
  for (i in seq_along(clusters)) {
    if (i %% 100 == 0) {
      cat(paste("Cluster", i, "/", n_clusters, "\n"))
    }
    
    cluster <- clusters[[i]]
    if (length(cluster$features) > 1) {
      features_tmp <- features[features[, name_col] %in% cluster$features, ]
      
      # Find the feature with maximal MPA
      max_mpa_idx <- which(features_tmp$MPA == max(features_tmp$MPA, na.rm = TRUE))[1]
      cluster_row <- features_tmp[max_mpa_idx, ]
      # Record all the features in the cluster
      cluster_row$Features <- paste(sort(cluster$features), collapse = ";")
      cluster_row$n_features <- length(cluster$features)
      # Create cluster ID
      cluster_row$Cluster_ID <- paste0("Cluster_", cluster_row[, name_col])
      cfeatures <- rbind(cfeatures, cluster_row)
      
      # Take the LC-MS data of the largest feature
      cdata_col <- data[features_tmp[max_mpa_idx, name_col]]
      colnames(cdata_col) <- cluster_row$Cluster_ID
      cdata <- cbind(cdata, cdata_col)
      
      handled_features <- c(handled_features, cluster$features)
    }
  }
  
  # Reorganise
  cfeatures <- dplyr::arrange(cfeatures, Cluster_ID)
  cdata <- cdata[c(sample_cols, cfeatures$Cluster_ID)]
  
  # All the features that were not in the clusters are retained unchanged
  missed_features <- features[!features[, name_col] %in% handled_features, ]
  missed_features$Features <- missed_features[, name_col]
  missed_features$n_features <- 1
  missed_features$Cluster_ID <- missed_features[, name_col]
  
  
  cfeatures <- rbind(cfeatures, missed_features)
  cdata <- cbind(cdata, data[missed_features[, name_col]])
  
  # Reorder columns
  cfeatures <- dplyr::select(cfeatures, "Cluster_ID", "n_features", "Features", name_col, dplyr::everything())
  rownames(cfeatures) <- 1:nrow(cfeatures)
  
  list(cdata = cdata, cfeatures = cfeatures)
}

# A helper function for plotting, scales the values in X
# between new min and max
rescale <- function(x, new_min, new_max) {
  (new_max - new_min) * (x - min(x)) / (max(x) - min(x)) + new_min
}


# UNFINISHED!!
plot_graph <- function(features, cluster, name_col, mz_col, rt_col) {
  
  # Ensure a correct order of the rows
  g <- cluster$graph
  vertices <- data.frame(Name = igraph::V(g)$name, stringsAsFactors = FALSE)
  colnames(vertices) <- name_col
  features_tmp <- dplyr::left_join(vertices, features, by = name_col)
  
  # Scaling of MPA to correct size
  # Square root to scale area, not radius
  size <- sqrt(rescale(features_tmp$MPA, new_min = 15^2, new_max = 40^2))
  
  if (length(igraph::E(g)) <= 20) {
    edge_labels <- as.character(round(igraph::E(g)$weight, digits = 2))
  } else {
    edge_labels <- NULL
  }
  
  
  g$palette <- RColorBrewer::brewer.pal(n = max(3, length(unique(igraph::degree(g)))), name = "Blues")
  plot(g, vertex.label = as.character(features_tmp[, mz_col]), vertex.size = size, vertex.label.dist = 0.1*size, vertex.label.degree = -pi/2,
       vertex.color = igraph::degree(g), edge.label = edge_labels)
}

plot_features <- function(features, cluster, name_col, mz_col, rt_col, rt_window) {
  
  features_tmp <- features[features[, name_col] %in% cluster$features, ]
  
  p1 <- ggplot(features_tmp, aes_string(mz_col, "MPA")) +
    geom_point(size = 3, color = "steelblue4") +
    geom_segment(aes_string(x = mz_col, yend = "MPA", xend = mz_col), y = 0, color = "steelblue4") +
    geom_label_repel(aes_string(label = mz_col), color = "steelblue4") +
    theme_minimal() +
    xlim(0.9*min(features_tmp[, mz_col], na.rm = TRUE),
         1.15*max(features_tmp[, mz_col], na.rm = FALSE)) +
    expand_limits(y=0) +
    labs(x = "Mass-to-charge ratio", y = "Median Peak Area")
  
  features_tmp$rtmin <- features_tmp[, rt_col] - rt_window
  features_tmp$rtmax <- features_tmp[, rt_col] + rt_window
  
  p2 <- ggplot(features_tmp, aes_string(rt_col, mz_col)) +
    geom_point(size = 3, color = "steelblue4") +
    geom_errorbarh(aes(xmin = rtmin, xmax = rtmax), color = "steelblue4") +
    theme_minimal() +
    labs(x = "Retention time", y = "Mass-to-charge ratio", title = "Retention time & tolerance")
  
  plot(p1)
  plot(p2)
}

plot_heatmaps <- function(data, features, cluster, name_col, mz_col, rt_col) {
  
  D <- data[, cluster$features]
  
  features_tmp <- features[features[, name_col] %in% cluster$features, ]
  
  n <- length(cluster$features)
  mz_rt <- data.frame()
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mz_rt <- rbind(mz_rt, data.frame(x = features_tmp[i, name_col],
                                       y = features_tmp[j, name_col],
                                       mz_diff = features_tmp[i, mz_col] - features_tmp[j, mz_col],
                                       rt_diff = features_tmp[i, rt_col] - features_tmp[j, rt_col],
                                       stringsAsFactors = FALSE))
    }
  }
  
  mz_ord <- features_tmp[, name_col][order(features_tmp[, mz_col])]
  mz_rt$x <- factor(mz_rt$x, levels = mz_ord)
  mz_rt$y <- factor(mz_rt$y, levels = rev(mz_ord))
  
  p1 <- ggplot(mz_rt, aes(x = x, y = y, fill = mz_diff)) +
    geom_tile(color = "grey80") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    scale_fill_gradient2()
  
  if (nrow(mz_rt) <= 10) {
    p <- p + geom_text(aes(label = round(mz_diff, digits = 2)))
  }
  
  # rt_ord <- features_tmp[, name_col][order(features_tmp[, rt_col])]
  # mz_rt$x <- factor(mz_rt$x, levels = rt_ord)
  # mz_rt$y <- factor(mz_rt$y, levels = rev(rt_ord))
  # 
  # ggplot(mz_rt, aes(x = x, y = y, fill = rt_diff)) +
  #   geom_tile(color = "grey80") +
  #   geom_text(aes(label = round(rt_diff, digits = 5))) +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #   scale_fill_gradient2()
  
  plot(p1)
}


#' Visualize clusters of fetaures
#' 
#' Draws multiple visualizations of each cluster, creating a separate file for each cluster
#' 
#' 
visualize_clusters <- function(data, features, clusters, min_size, rt_window, name_col, mz_col, rt_col, file_path) {
  
  for (i in seq_along(clusters)) {
    if (i %% 100 == 0) {
      print(paste(i, "/", length(clusters)))
    }
    cluster <- clusters[[i]]
    if (length(cluster$features) >= min_size) {
      features_tmp <- features[features[, name_col] %in% cluster$features, ]
      cluster_id <- features_tmp[, name_col][which(features_tmp$MPA == max(features_tmp$MPA, na.rm = TRUE))[1]]
      
      pdf(paste0(file_path, "Cluster_", i, "_", cluster_id, ".pdf"), width = 10, height = 10)
      plot_heatmaps(data, features, cluster, name_col, mz_col, rt_col)
      plot_features(features, cluster, name_col, mz_col, rt_col, rt_window)
      plot_graph(features, cluster, name_col, mz_col, rt_col)
      dev.off()
    }
  }
  
  
}