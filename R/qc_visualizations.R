
density_plot <- function(data, x, fill, fill_scale = NULL, color_scale = NULL,
                         title = NULL, subtitle = NULL,
                         xlab = x, fill_lab = fill) {

  fill_scale <- fill_scale %||% getOption("amp.fill_scale_dis")
  color_scale <- color_scale %||% getOption("amp.color_scale_dis")

  p <- ggplot(data, aes_string(x, fill = fill, color = fill)) +
    geom_density(alpha = 0.2) +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, fill = fill_lab, color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    color_scale

  p
}

#' Distance density plot
#'
#' Plot density of distances between samples in QC samples and actual samples
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param dist_method method for calculating the distances, passed to dist
#' @param center logical, should the data be centered?
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param color_scale a scale for the color of the edge of density curves, as returned by a ggplot function
#' @param fill_scale a scale for the fill of the density curves, as returned by a ggplot function
#' @param title the plot title
#' @param subtitle the plot subtitle
#'
#' @return a ggplot object
#'
#' @examples
#' plot_dist_density(merged_sample)
#' # Drift correction tightens QCs together
#' plot_dist_density(correct_drift(merged_sample))
#'
#' @seealso \code{\link[stats]{dist}}
#'
#' @export
plot_dist_density <- function(object, all_features = FALSE, dist_method = "euclidean",
                              center = TRUE, scale = "uv",
                              color_scale = NULL, fill_scale = NULL,
                              title = NULL, subtitle = NULL) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  title <- title %||% paste("Density plot of", dist_method, "distances between samples")

  object <- pcaMethods::prep(object, center = center, scale = scale)

  qc_data <- t(exprs(object)[, object$QC == "QC"])
  sample_data <- t(exprs(object)[, object$QC != "QC"])

  qc_dist <- dist(qc_data, method = dist_method) %>% as.numeric()
  sample_dist <- dist(sample_data, method = dist_method) %>% as.numeric()
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  distances <- data.frame(dist = c(qc_dist, sample_dist), qc = qc)

  density_plot(distances, x = "dist", fill = "qc", fill_scale = fill_scale,
               color_scale = color_scale,
               xlab = "Distance", fill_lab = NULL,
               title = title, subtitle = subtitle)
}

#' Estimate the magnitude of drift
#'
#' Plots histograms of p-values from linear regression model, where each feature is predicted
#' by injection order alone. The expected uniform distribution is represented by a dashed red line.
#'
#' @param object A MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{plot_p_histogram}}
#'
#' @examples
#' plot_injection_lm(merged_sample)
#'
#' @export
plot_injection_lm <- function(object, all_features = FALSE) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

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
#' @param p_values list or data frame, each element/column is a vector of p-values. The list names are used as plot titles
#' @param hline logical, whether a horizontal line representing uniform distribution should be plotted
#' @param combine logical, whether plots of individual p-value vectors should be combined into a single object. Set to FALSE if
#' you want to add other plots to the list before plotting
#'
#' @return if combine = TRUE, a ggplot object. Otherwise a list of ggplot objects
#'
#' @export
plot_p_histogram <- function(p_values, hline = TRUE, combine = TRUE, x_label = "p-value") {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  # Custom breaks for the x-axis
  breaks <- seq(0, 1, by = 0.05)

  # THree separate histograms
  plots <- list()
  for (i in seq_along(p_values)) {


    p <- ggplot(data.frame(P = p_values[[i]]), aes(P)) +
      geom_histogram(breaks = breaks, col = "grey50", fill = "grey80", size = 1) +
      labs(x = x_label, y = "Frequency") +
      ggtitle(names(p_values)[i]) +
      theme_minimal() +
      theme(plot.title = element_text(face="bold", hjust=0.5))

    if (hline) {
      # Compute the position of the expected line
      finite_count <- sum(is.finite(p_values[[i]]))
      h_line <- finite_count/(length(breaks)-1)
      p <- p +
        geom_hline(yintercept = h_line, color="red", linetype = "dashed", size = 1)
    }

    plots <- c(plots, list(p))
  }

  if (combine) {
    return(cowplot::plot_grid(plotlist = plots, ncol = 1))
  } else {
    return(plots)
  }
}


#' Plot quality metrics
#'
#' Plots distribution of each quality metric, and a distribution of the flags.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#'
#' @return a ggplot object
#'
#' @examples
#' plot_quality(example_set)
#'
#' @export
plot_quality <- function(object, all_features = FALSE) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
  }

  # Plot bar plot of flags
  flags <- flag(object)
  flags[is.na(flags)] <- "Good"
  flags <- factor(flags) %>% relevel(ref = "Good")

  fp <- ggplot(data.frame(flags), aes(x = flags)) +
    geom_bar(col = "grey50", fill = "grey80", size = 1) +
    scale_y_continuous(sec.axis = sec_axis(~.*100/length(flags), name = "Percentage")) +
    theme_minimal() +
    labs(x = "Flag")

  # Drop flagged features
  object <- drop_flagged(object, all_features = all_features)

  if (is.null(quality(object))) {
    cat("\nQuality metrics not found, computing them now\n")
    object <- assess_quality(object)
  }

  # Distribution of quality metrics
  qps <- plot_p_histogram(quality(object)[, -1], hline = FALSE, combine = FALSE, x_label = "")

  p <- cowplot::plot_grid(plotlist = c(qps, list(fp)), ncol = 1)
  p
}


#' Plot a boxplot for each sample
#'
#' Plots a boxplot of the distribution of the metabolite values for each sample. The boxplots can
#' be ordered and filled by any combination of columns in the pheno data. By default, order and
#' fill are both determined by the combination of group and time columns.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param order_by character vector, names of columns used to order the samples
#' @param fill_by character vector, names of columns used to fill the boxplots
#' @param title,subtitle character, title and subtitle of the plot
#' @param fill_scale a scale for the fill of the boxplots, as returned by a ggplot function
#' @param zoom_boxplot logical, whether outliers should be left outside the plot and only
#' the boxplots shown. Defaults to TRUE.
#'
#' @return a ggplot object
#'
#' @examples
#' plot_sample_boxplots(merged_sample, order_by = "Group", fill_by = "Group")
#'
#' @export
plot_sample_boxplots <- function(object, all_features = FALSE, order_by = NULL, fill_by = NULL,
                             title = "Boxplot of samples", subtitle = NULL,
                             fill_scale = NULL, zoom_boxplot = TRUE) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  order_by <- order_by %||% as.character(na.omit(c(group_col(object), time_col(object))))
  fill_by <- fill_by %||% as.character(na.omit(c(group_col(object), time_col(object))))
  fill_scale <- fill_scale %||% getOption("amp.fill_scale_dis")

  data <- combined_data(object)

  if (length(order_by) == 1) {
    data$order_by <- data[, order_by]
  } else {
    data <- tidyr::unite(data, "order_by", order_by, remove = FALSE)
  }

  if (length(fill_by) == 1) {
    data$fill_by <- data[, fill_by]
  } else {
    data <- tidyr::unite(data, "fill_by", fill_by, remove = FALSE)
  }

  data <- data %>%
    dplyr::arrange(order_by)

  data$Sample_ID <- factor(data$Sample_ID, levels = data$Sample_ID)

  data <- tidyr::gather(data, "Variable", "Value", rownames(exprs(object)))

  p <- ggplot(data, aes(x = Sample_ID, y = Value, fill = fill_by))

  ## Zooming outliers out of view
  if(zoom_boxplot){
    # compute lower and upper whiskers
    ylimits <- data %>%
      dplyr::group_by(Sample_ID) %>%
      dplyr::summarise(low = boxplot.stats(Value)$stats[1],
                       high = boxplot.stats(Value)$stats[5])


    ylimits <- c(0, max(ylimits$high))
    # scale y limits based on ylim1
    p <- p +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = ylimits)
    # add text to main title
    subtitle <- paste(subtitle, "(zoomed in boxplot: outliers out of view)", sep=" ")
  } else {
    p <- p +
      geom_boxplot()
  }

  p <- p +
    labs(x = paste(order_by, collapse = "_"),
         y = "Abundance of metabolites",
         fill = paste(fill_by, collapse = "_"),
         title = title, subtitle = subtitle) +
    theme_bw() +
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, vjust=0.3)) +
    fill_scale

  p
}

