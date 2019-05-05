

density_plot <- function(data, x, fill, fill_scale = NULL, title = NULL, subtitle = NULL, xlab = x, fill_lab = fill) {

  fill_scale <- fill_scale %||% getOption("amp.fill_scale_dis")

  # ggpubr::ggdensity(data, x, fill = fill, palette = c("#00AFBB", "#E7B800")) +
  #   #fill_scale +
  #   labs(title = title, subtitle = subtitle, x = xlab, fill = fill_lab)
  p <- ggplot(data, aes_string(x, fill = fill)) +
    geom_density(alpha = 0.2, color = NA) +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, fill = fill_lab, color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())

  p
}

#' Distance density plot
#'
#' Plot density of distances between samples in QC samples and actual samples
#'
#' @param object a MetaboSet object
#' @param dist_method method for calculating the distances, passed to dist
#' @param fill_scale a scale for the fill of the density curves, as returned by a ggplot function
#' @param xlab the x
#'
#'
#' @seealso \code{\link[stats]{dist}}
#'
#' @export
plot_dist_density <- function(object, dist_method = "euclidean", fill_scale = NULL,
                              title = NULL, subtitle = NULL) {

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

#' Estimate the magnitude of drift
#'
#' Plots histograms of p-values from linear regression model, where each feature is predicted
#' by injection order alone. The expected uniform distribution is represented by a dashed red line.
#'
#' @param object A MetaboSet object
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{plot_p_histogram}}
#'
#' @export
plot_injection_lm <- function(object) {

  # Apply linear model to QC samples and biological samples separately
  lm_all <- perform_lm(object, "Feature ~ Injection_order")
  lm_sample <- perform_lm(object[, object$QC != "QC"], "Feature ~ Injection_order")
  lm_qc <- perform_lm(object[, object$QC == "QC"], "Feature ~ Injection_order")

  # Only interested in the p_values
  p_values <- list("All samples" = lm_all$Injection_order_P,
                   "Biological samples" = lm_sample$Injection_order_P,
                   "QC samples" = lm_qc$Injection_order_P)
  # Plotting
  plot_p_histogram(p_values)
}

#' Histogram of p-values
#'
#' Draws histograms of p-values with expected uniform distribution represented by a dashed red line
#'
#' @param p_values list, each element is a vector of p-values. The list names are used as plot titles
#'
#' @return A ggplot object.
#'
#' @export
plot_p_histogram <- function(p_values) {
  # Custom breaks for the x-axis
  breaks <- seq(0, 1, by = 0.05)

  # THree separate histograms
  plots <- list()
  for (i in seq_along(p_values)) {
    # Compute the position of the expected line
    finite_count <- sum(is.finite(p_values[[i]]))
    h_line <- finite_count/(length(breaks)-1)

    p <- ggplot(data.frame(P = p_values[[i]]), aes(P)) +
      geom_histogram(breaks = breaks, col = "grey50", fill = "grey80", size = 1) +
      geom_hline(yintercept = h_line, color="red", linetype = "dashed", size = 1) +
      labs(x="p-value", y="Frequency") +
      ggtitle(names(p_values)[i]) +
      theme_minimal() +
      theme(plot.title = element_text(face="bold", hjust=0.5))

    plots <- c(plots, list(p))
  }
  # Combine plots
  p <- cowplot::plot_grid(plotlist = plots, ncol = 1)

  p
}
