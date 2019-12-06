
#' Draw dendrograms
#'
#' Draws a dendrogram of a hierarchical clustering applied to the samples of an experiment
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param color character, name of the column used for coloring the sample labels
#' @param dist_method distance method used in clustering, see ?dist
#' @param clust_method method used in clustering, see ?hclust
#' @param center logical, should the data be centered?
#' @param scale scaling used, as in \code{pcaMethods::prep}. Default is "uv" for unit variance
#' @param title The plot title
#' @param subtitle The plot subtitle
#' @param color_scale the color scale as returned by a ggplot function
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_dendrogram(merged_sample)
#'
#' @seealso \code{\link{dist}} \code{\link{hclust}}
#'
#' @export
plot_dendrogram <- function(object, all_features = FALSE, color = group_col(object), dist_method = "euclidean", clust_method = "ward.D2",
                     center = TRUE, scale = "uv", title = "Dendrogram of hierarchical clustering",
                     subtitle = NULL, color_scale = NULL) {

  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
      stop("Package \"ggdendro\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  subtitle <- subtitle %||% paste("Distance method:", dist_method, "Clustering method:", clust_method)

  if (is.null(color_scale)) {
    if (class(pData(object)[, color]) == "numeric") {
      color_scale <- getOption("amp.color_scale_con")
    } else {
      color_scale <- getOption("amp.color_scale_dis")
    }
  }

  object <- pcaMethods::prep(object, center = center, scale = scale)

  d_data <- dist(t(exprs(object)), method = dist_method) %>%
    hclust(method = clust_method) %>%
    as.dendrogram() %>%
    ggdendro::dendro_data()

  labels <- ggdendro::label(d_data) %>%
    dplyr::mutate(label = as.character(label)) %>%
    dplyr::left_join(pData(object)[c("Sample_ID", color)], by = c("label" = "Sample_ID"))
  labels[, color] <- as.factor(labels[, color])

  p <- ggplot(ggdendro::segment(d_data)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data = labels, aes_string(x="x", y="y", label = "label", color = color),
              angle = 90, hjust = 1) +
    ggdendro::theme_dendro() +
    color_scale +
    labs(title = title, subtitle = subtitle)

  p

}
#' Draw heatmaps
#'
#' Draws a heatmap of the distances between the samples of an experiment,
#' the samples are ordered by hierarchical clustering.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param dist_method distance method used in clustering, see \code{\link{dist}}
#' @param clust_method clustering method used in clustering, see \code{\link{hclust}}
#' @param center logical, should the data  be centered?
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param group_bar logical, should a bar showing the groups be drawn under the heat map?
#' @param group character, name of the column used for coloring the group bar
#' @param title The plot title
#' @param subtitle The plot subtitle
#' @param fill_scale_con Continuous fill scale for the heatmap as returned by a ggplot function
#' @param fill_scale_dis Discrete fill scale for the group bar as returned by a ggplot function
#'
#' @return a ggplot object. If \code{group_bar} is \code{TRUE}, the plot will consist of multiple
#' parts and is harder to modify.
#'
#' @examples
#' plot_sample_heatmap(merged_sample)
#'
#' @seealso \code{\link{dist}} \code{\link{hclust}}
#'
#' @export
plot_sample_heatmap <- function(object, all_features = FALSE, dist_method = "euclidean", clust_method = "ward.D2",
                         center = TRUE, scale = "uv",
                         group_bar = TRUE, group = group_col(object),
                         title = "Heatmap of distances between samples",
                         subtitle = NULL, fill_scale_con = NULL, fill_scale_dis = NULL) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  # Default settings
  subtitle <- subtitle %||% paste("Distance method:", dist_method, "Clustering method:", clust_method)
  fill_scale_con <- fill_scale_con %||% getOption("amp.fill_scale_con")
  fill_scale_dis <- fill_scale_dis %||% getOption("amp.fill_scale_dis")

  object <- pcaMethods::prep(object, center = center, scale = scale)

  # Distances
  distances <- dist(t(exprs(object)), method = dist_method)
  # Hierarchical clustering for ordering
  hc <- hclust(distances, method = clust_method)
  hc_order <- hc$labels[hc$order]

  # From wide to long format for ggplot
  distances_df <- as.matrix(distances) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("X") %>%
    tidyr::gather("Y", "Distance", -X)
  # Set the correct order given by hclust
  distances_df$X <- factor(distances_df$X, levels = hc_order, ordered = TRUE)
  distances_df$Y <- factor(distances_df$Y, levels = rev(hc_order), ordered = TRUE)

  # Heatmap
  p <- ggplot(distances_df, aes(X, Y, fill = Distance)) +
    geom_tile(color = NA) +
    fill_scale_con +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.05)) +
    coord_fixed()
  # Group bar
  if (group_bar) {
    pheno_data <- pData(object)
    pheno_data$Sample_ID <- factor(pheno_data$Sample_ID, levels = hc_order, ordered = TRUE)

    gb <- ggplot(pheno_data, aes_string(x = "Sample_ID", y = 1, fill = group)) +
      geom_tile(color = "white") +
      theme_void() +
      fill_scale_dis

    p <- cowplot::plot_grid(p, gb, ncol = 1, align = "v", rel_heights = c(10/11,1/11))
  }

  p
}
