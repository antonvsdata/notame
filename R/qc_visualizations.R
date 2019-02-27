

density_plot <- function(data, x, fill, fill_scale = NULL, title = NULL, subtitle = NULL, xlab = x, fill_lab = fill) {

  fill_scale <- fill_scale %||% getOption("amp.fill_scale_dis")

  # ggpubr::ggdensity(data, x, fill = fill, palette = c("#00AFBB", "#E7B800")) +
  #   #fill_scale +
  #   labs(title = title, subtitle = subtitle, x = xlab, fill = fill_lab)
  p <- ggplot(data, aes_string(x, fill = fill, color = NULL)) +
    geom_density(alpha = 0.2) +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, fill = fill_lab, color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())

  p
}

plot_dist_density <- function(object, dist_method = "euclidean", fill_scale = NULL,
                              title = NULL,
                              subtitle = NULL, xlab = x, fill_lab = fill) {

  title <- title %||% paste("Density plot of", dist_method, "distances between samples")

  qc_data <- t(exprs(object)[, object$QC == "QC"])
  sample_data <- t(exprs(object)[, object$QC != "QC"])

  qc_dist <- dist(qc_data, method = dist_method) %>% as.numeric()
  sample_dist <- dist(sample_data, method = dist_method) %>% as.numeric()
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  distances <- data.frame(dist = c(qc_dist, sample_dist), qc = qc)

  density_plot(distances, x = "dist", fill = "qc", fill_scale = fill_scale, xlab = "Distance", fill_lab = NULL,
               title = title, subtitle = subtitle)
}
