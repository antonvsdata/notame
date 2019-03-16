

plot_dendrogram <- function(object, color = group_col(object), dist_method = "euclidean", clust_method = "ward.D2",
                     center = TRUE, scale = "uv", title = "Dendrogram of hierarchical clustering",
                     subtitle = NULL, color_scale = NULL) {

  subtitle <- subtitle %||% paste("Distance method:", dist_method, "Clustering method:", clust_method)

  if (is.null(color_scale)) {
    if (class(pData(object)[, color]) == "numeric") {
      color_scale <- getOption("amp.color_scale_con")
    } else {
      color_scale <- getOption("amp.color_scale_dis")
    }
  }

  prepd <- pcaMethods::prep(object, center = center, scale = scale)

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
    geom_text(data = labels, aes_string(x="x", y="y", label = "label", color = color), angle = 90, hjust = 1) +
    ggdendro::theme_dendro() +
    color_scale +
    labs(title = title, subtitle = subtitle)

  p

}


plot_heatmap <- function(object, text_color = group_col(object), dist_method = "euclidean", clust_method = "ward.D2",
                         center = TRUE, scale = "uv", title = "Heatmap of distances between samples",
                         subtitle = NULL, fill_scale = NULL, color_scale = NULL) {

  subtitle <- subtitle %||% paste("Distance method:", dist_method, "Clustering method:", clust_method)
  fill_scale <- fill_scale %||% getOption("amp.fill_scale_con")
  color_scale <- color_scale %||% getOption("amp.color_scale_dis")

  distances <- dist(t(exprs(object)), method = dist_method)

  hc <- hclust(distances, method = clust_method)
  hc_order <- hc$labels[hc$order]


}
