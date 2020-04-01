

# ------- HELPER FUNCTIONS ----------------

# Helper function for computing PCA
pca_helper <- function(object, center, scale, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  res_pca <- pcaMethods::pca(object, scale = scale, center = center, ...)
  pca_scores <- pcaMethods::scores(res_pca) %>% as.data.frame()
  R2 <- res_pca@R2
  labels <- paste0(c("PC1", "PC2"), " (", scales::percent(R2), ")")

  return(list(pca_scores = pca_scores, labels = labels))
}

# Helper function for computing t-SNE
t_sne_helper <- function(object, center, scale, perplexity, pca_method, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
      stop("Package \"Rtsne\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  prepd <- pcaMethods::prep(object, center = center, scale = scale)

  if (sum(is.na(exprs(prepd))) > 0) {
    res_pca <- pcaMethods::pca(object, method = pca_method, nPcs = min(nrow(object), ncol(object), 50),
                               scale = "none", center = FALSE)
    pca_scores <- pcaMethods::scores(res_pca)
    res_tsne <- Rtsne::Rtsne(pca_scores, perplexity = perplexity, pca = FALSE, ...)
  } else {
    res_tsne <- Rtsne::Rtsne(t(exprs(prepd)), perplexity = perplexity, ...)
  }
  data.frame(res_tsne$Y)
}

# -------------- SCATTER PLOTS ---------------

#' PCA scatter plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} package
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param label character, name of the column used for point labels
#' @param density logical, whether to include density plots to both axes. The density curves will be split and colored by the 'color' variable.
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function. Set to NA to choose the appropriate scale based on the class of the coloring variable.
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param fill_scale the fill scale used for density curves. If a continuous variable is used as color, density curve will be colorless.
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return a ggplot object. If \code{density} is \code{TRUE}, the plot will consist of multiple
#' parts and is harder to modify.
#'
#' @examples
#' plot_pca(merged_sample, color = "Injection_order", shape = "Group")
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                     color = group_col(object), shape = NULL, label = NULL, density = FALSE,  title = "PCA",
                     subtitle = NULL, color_scale = NA,
                     shape_scale = getOption("notame.shape_scale"), fill_scale = getOption("notame.fill_scale_dis"), ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  shape <- shape %||% color

  pca_results <- pca_helper(object, center, scale, ...)
  pca_scores <- pca_results$pca_scores
  pca_scores[color] <- pData(object)[, color]
  pca_scores[shape] <- pData(object)[, shape]
  pca_scores[label] <- pData(object)[, label]

  scatter_plot(pca_scores, x = "PC1", y = "PC2", xlab = pca_results$labels[1], ylab = pca_results$labels[2],
               color = color, shape = shape, label = label, density = density, title = title,
               subtitle = subtitle, color_scale = color_scale, shape_scale = shape_scale,
               fill_scale = fill_scale)
}

#' t-SNE scatter plot
#'
#' Computes t-SNE into two dimensions and plots the map points.
#' In case there are missing values, PCA is performed using the nipals method of \code{pcaMethods::pca},
#' the  method can be changed to "ppca" if nipals fails.
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} and \code{Rtsne} packages
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param perplexity the perplexity used in t-SNE
#' @param pca_method the method used in PCA if there are missing values
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param label character, name of the column used for point labels
#' @param density logical, whether to include density plots to both axes. The density curves will be split and colored by the 'color' variable.
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function. Set to NA to choose the appropriate scale based on the class of the coloring variable.
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param fill_scale the fill scale used for density curves. If a continuous variable is used as color, density curve will be colorless.
#' @param ... additional arguments passed to \code{Rtsne::Rtsne}
#'
#' @return a ggplot object. If \code{density} is \code{TRUE}, the plot will consist of multiple
#' parts and is harder to modify.
#'
#' @examples
#' plot_tsne(merged_sample, color = "Time", shape = "Group")
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne <- function(object, all_features = FALSE, center = TRUE, scale = "uv", perplexity = 30,
                      pca_method = "nipals",
                      color = group_col(object), shape = NULL, label = NULL, density = FALSE, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), color_scale = NA,
                      shape_scale = getOption("notame.shape_scale"), fill_scale = getOption("notame.fill_scale_dis"), ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  shape <- shape %||% color

  # t-SNE
  tsne_scores <- t_sne_helper(object, center, scale, perplexity, pca_method, ...)
  # Add columns for plotting
  tsne_scores[color] <- pData(object)[, color]
  tsne_scores[shape] <- pData(object)[, shape]
  tsne_scores[label] <- pData(object)[, label]

  scatter_plot(tsne_scores, x = "X1", y = "X2", color = color, shape = shape, label = label,
               density = density, title = title, subtitle = subtitle,
               color_scale = color_scale, shape_scale = shape_scale, fill_scale = fill_scale)

}

scatter_plot <- function(data, x, y, color, shape, label = NULL, density = FALSE, fixed = TRUE, color_scale = NA,
                         shape_scale = NULL, fill_scale = NA, title = NULL, subtitle = NULL, xlab = x, ylab = y,
                         color_lab = color, shape_lab = shape) {

  if (!is.null(color_scale)) {
    if (is.na(color_scale)) {
      if (class(data[, color]) %in% c("numeric", "integer")) {
        color_scale <- getOption("notame.color_scale_con")
      } else {
        color_scale <- getOption("notame.color_scale_dis")
      }
    }
  }


  p <- ggplot(data, aes_string(x = x, y = y, color = color)) +
    theme_bw() +
    color_scale +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
         color = color_lab)
  if (fixed) {
    p <- p + coord_fixed() + theme(aspect.ratio=1)
  }

  if (class(data[, shape]) == "character") {
    data[shape] <- as.factor(data[, shape])
    warning(paste("Shape variable not given as a factor, converted to factor with levels",
                  paste(levels(data[, shape]), collapse = ", ")),
            call. = FALSE)
  }

  if (class(data[, shape]) == "factor") {
    if (length(levels(data[, shape])) <= 8){
      p <- p +
        geom_point(aes_string(shape = shape)) +
        shape_scale +
        labs(shape = shape_lab)

    } else if (is.null(shape_scale)) {
      cat("Only 8 distinct shapes currently available!")
      p <- p +
        geom_point()
    }
  } else {
    p <- p +
      geom_point()
  }

  # Add point labels
  if (!is.null(label)) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package \"ggrepel\" needed for this function to label the points. Please install it.",
             call. = FALSE)
    }
    p <- p +
      ggrepel::geom_text_repel(aes_string(label = label))
  }

  # Add density plots to top and right
  if (density) {
    if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to add density curves. Please install it.",
           call. = FALSE)
    }
    xdens <- cowplot::axis_canvas(p, axis = "x")+
      geom_density(data = data, aes_string(x = x, fill = color),
                   alpha = 0.7, size = 0.2) +
      fill_scale

    ydens <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE)+
      geom_density(data = data, aes_string(x = y, fill = color),
                   alpha = 0.7, size = 0.2)+
      coord_flip() +
      fill_scale

    p <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p <- cowplot::insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
    p <- cowplot::ggdraw(p)
  }

  p
}

#' PCA loadings plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the loadings of first principal components
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} package
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param npc1,npc2 number of top feature to plot for first and second principal component
#' @param title,subtitle the titles of the plot
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return a ggplot object.
#'
#' @examples
#' plot_pca_loadings(merged_sample, npc1 = 5, npc2 = 5)
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_loadings <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                              npc1 = 10, npc2 = 10,
                              title = "PCA loadings", subtitle = NULL, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("Package \"ggrepel\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  pca_res <- pcaMethods::pca(object, center = TRUE, scale = "uv", ...)

  loads <- pca_res@loadings %>%
    as.data.frame()
  loads$Feature_ID <- rownames(loads)

  features_pc1 <- loads$Feature_ID[order(abs(loads$PC1), decreasing = TRUE)][seq_len(npc1)]
  features_pc2 <- loads$Feature_ID[order(abs(loads$PC2), decreasing = TRUE)][seq_len(npc2)]

  loads <- loads[union(features_pc1, features_pc2),]

  ggplot(loads, aes(x = PC1, y = PC2, label = Feature_ID)) +
    geom_point() +
    ggrepel::geom_text_repel() +
    theme_bw() +
    labs(title = title, subtitle = subtitle)

}


# --------------- HEXBIN PLOTS --------------------


#' PCA hexbin plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components as hexagonal bins, where the value of the coloring
#' variable is summarised for each bin, by default as the mean of the values inside the bin.
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} package
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param fill character, name of the column used for coloring the hexagons
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param title,subtitle the titles of the plot
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_pca_hexbin(merged_sample)
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_hexbin <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                     fill = "Injection_order", summary_fun = "mean", bins = 10, title = "PCA",
                     subtitle = NULL, fill_scale = getOption("notame.fill_scale_con"), ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  pca_results <- pca_helper(object, center, scale, ...)
  pca_scores <- pca_results$pca_scores
  pca_scores[fill] <- pData(object)[, fill]

  hexbin_plot(data = pca_scores, x = "PC1", y = "PC2", xlab = pca_results$labels[1], ylab = pca_results$labels[2],
              fill = fill, summary_fun = summary_fun, bins = bins, fill_scale = fill_scale,
                          title = title, subtitle = subtitle)
}

#' t-SNE hexbin plot
#'
#' Computes t-SNE into two dimensions and plots the map as hexagonal bins, where the value of the coloring
#' variable is summarised for each bin, by default as the mean of the values inside the bin.
#' In case there are missing values, PCA is performed using the nipals method of \code{pcaMethods::pca},
#' the  method can be changed to "ppca" if niipals fails.
#'
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} and \code{Rtsne} packages
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param pca_method the method used in PCA if there are missing values
#' @param perplexity the perplexity used in t-SNE
#' @param fill character, name of the column used for coloring the hexagons
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param title,subtitle the titles of the plot
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param ... additional arguments passed to Rtsne::Rtsne
#'
#' @return
#' A ggplot object.
#'
#' @examples
#' plot_tsne_hexbin(merged_sample)
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne_hexbin <- function(object, all_features = FALSE, center = TRUE, scale = "uv", pca_method = "nipals", perplexity = 30,
                      fill = "Injection_order", summary_fun = "mean", bins = 10, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), fill_scale = getOption("notame.fill_scale_con"), ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  # t-SNE
  tsne_scores <- t_sne_helper(object, center, scale, perplexity, pca_method, ...)
  # Add columns for plotting
  tsne_scores[fill] <- pData(object)[, fill]

  hexbin_plot(tsne_scores, x = "X1", y = "X2",
               fill = fill, summary_fun = summary_fun, bins = bins,
               fill_scale = fill_scale,
               title = title, subtitle = subtitle)

}


hexbin_plot <- function(data, x, y, fill, summary_fun = "mean", bins = 10, fill_scale = NULL,
                         title = NULL, subtitle = NULL, xlab = x, ylab = y,
                         fill_lab = fill) {
  if (!requireNamespace("hexbin", quietly = TRUE)) {
      stop("Package \"hexbin\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop("Package \"Hmisc\" needed for this function to work. Please install it.",
           call. = FALSE)
  }

  p <- ggplot(data, aes_string(x = x, y = y, z = fill)) +
    stat_summary_hex(bins = bins, fun = summary_fun) +
    theme_bw() +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
         fill = fill_lab) +
    coord_fixed() +
    theme(aspect.ratio=1)
  p
}


# --------- ARROW PLOTS -----------

arrow_plot <- function(data, x, y, color, time, subject, alpha, arrow_style,
                       color_scale, title, subtitle, xlab, ylab) {

  data <- data[order(data[, subject], data[, time]), ]

  p <- ggplot(data, aes_string(x = x, y = y, color = color, group = subject)) +
    geom_path(arrow = arrow_style, alpha = alpha) +
    theme_bw() +
    color_scale +
    coord_fixed() +
    theme(aspect.ratio=1) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle)
  p
}

#' PCA plot with arrows
#'
#' Plots changes in PCA space according to time. All the observations of a single subject are connected
#' by an arrow ending at the last observation.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param color character, name of the column used for coloring the arrows
#' @param time character, name of the column containing timepoints
#' @param subject character, name of the column containing subject identifiers
#' @param alpha numeric, value for the alpha parameter of the arrows (transparency)
#' @param arrow_style a description of arrow heads, the size and angle can be modified, see \code{?arrow}
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return a ggplot object.
#'
#' @examples
#' plot_pca_arrows(drop_qcs(example_set))
#' # If the sample size is large, plot groups separately
#' plot_pca_arrows(drop_qcs(example_set)) +
#' facet_wrap(~ Group)
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_arrows <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                            color = group_col(object), time = time_col(object), subject = subject_col(object),
                            alpha = 0.6, arrow_style = arrow(), title = "PCA changes",
                            subtitle = NULL, color_scale = getOption("notame.color_scale_dis"), ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  pca_results <- pca_helper(object, center, scale, ...)
  pca_scores <- pca_results$pca_scores
  pca_scores[color] <- pData(object)[, color]
  pca_scores[time] <- pData(object)[, time]
  pca_scores[subject] <- pData(object)[, subject]

  arrow_plot(data = pca_scores, x = "PC1", y = "PC2", color = color, time = time, subject = subject,
             alpha = alpha, arrow_style = arrow_style, color_scale = color_scale,
             title = title, subtitle = subtitle,
             xlab = pca_results$labels[1], ylab = pca_results$labels[2])

}



#' t-SNE plot with arrows
#'
#' Computes t-SNE into two dimensions and plots changes according to time.
#' All the observations of a single subject are connected by an arrow ending at the last observation.
#' In case there are missing values, PCA is performed using the nipals method of \code{pcaMethods::pca},
#' the  method can be changed to "ppca" if nipals fails.
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} and \code{Rtsne} packages
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param perplexity the perplexity used in t-SNE
#' @param pca_method the method used in PCA if there are missing values
#' @param color character, name of the column used for coloring the points
#' @param time character, name of the column containing timepoints
#' @param subject character, name of the column containing subject identifiers
#' @param alpha numeric, value for the alpha parameter of the arrows (transparency)
#' @param arrow_style a description of arrow heads, the size and angle can be modified, see \code{?arrow}
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function
#' @param ... additional arguments passed to \code{Rtsne::Rtsne}
#'
#' @return a ggplot object. If \code{density} is \code{TRUE}, the plot will consist of multiple
#' parts and is harder to modify.
#'
#' @examples
#' plot_tsne_arrows(drop_qcs(example_set), perplexity = 5)
#' # If the sample size is large, plot groups separately
#' plot_tsne_arrows(drop_qcs(example_set), perplexity = 5) +
#' facet_wrap(~ Group)
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne_arrows <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                             perplexity = 30, pca_method = "nipals",
                             color = group_col(object), time = time_col(object), subject = subject_col(object),
                             alpha = 0.6, arrow_style = arrow(), title = "t-SNE changes",
                             subtitle = paste("Perplexity:", perplexity), color_scale = getOption("notame.color_scale_dis"), ...) {

  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  tsne_scores <- t_sne_helper(object, center, scale, perplexity, pca_method, ...)
  tsne_scores[color] <- pData(object)[, color]
  tsne_scores[time] <- pData(object)[, time]
  tsne_scores[subject] <- pData(object)[, subject]

  arrow_plot(data = tsne_scores, x = "X1", y = "X2", color = color, time = time, subject = subject,
             alpha = alpha, arrow_style = arrow_style, color_scale = color_scale,
             title = title, subtitle = subtitle,
             xlab = "X1", ylab = "X2")
}


# -------- VOLCANO PLOT -----------

minus_log10 <- scales::trans_new("minus_log19",
                                 transform = function(x) {-log10(x)},
                                 inverse = function(x) {10^{-x}})

#' Volcano plot
#'
#' Draws a volcano plot of effect size and p-values.
#'
#' @param data a data frame with the effect size and p-values
#' @param x,p the column names of effect size (x-axis) and p-values
#' @param p_fdr column name of FDR corrected p-values, used to draw a line showing the fdr-corrected significance level
#' @param color column name used to color the plots
#' @param p_breaks a numerical vector of the p_values to show on the y-axis
#' @param fdr_limit the significance level used in the experiment
#' @param log2_x logical, whether effect size should be plotted on a log2 axis
#' @param center_x_axis logical, whether x-axis should be centered. If \code{TRUE}, the "zero-effect" will
#' be on the middle of the plot. The "zero effect" is 0 if \code{log2_x = FALSE} and 1 if  \code{log2_x = TRUE}
#' @param x_lim numerical vector of length 2 for manually setting the x-axis limits
#' @param color_scale the color scale as returned by a ggplot function
#' @param title,subtitle the title and subtitle of the plot
#' @param ...  parameters passed to \code{\link[ggplot2]{geom_point}}, such as shape and alpha values. New aesthetics can
#' also be passed using \code{mapping = aes(...)}.
#'
#' @return a ggplot object
#'
#' @examples
#' # naturally, this looks messy as there are not enough p-values
#' lm_results <- perform_lm(drop_qcs(merged_sample), formula_char = "Feature ~ Group")
#' volcano_plot(data = lm_results, x = "GroupB_Estimate",
#'              p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'              fdr_limit = 0.1)
#'
#' @export
volcano_plot <- function(data, x, p, p_fdr = NULL, color = NULL,
                         p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
                         log2_x = FALSE, center_x_axis = TRUE, x_lim = NULL,
                         color_scale = getOption("notame.color_scale_con"),
                         title = "Volcano plot", subtitle = NULL, ...) {

  if (center_x_axis & !is.null(x_lim)) {
    warning("Manually setting x-axis limits overrides x-axis centering")
    center_x_axis <- FALSE
  }
  if (min(data[, p]) > max(p_breaks)) {
    warning("All the p-values are larger than the p-value breaks supplied. Consider using larger p_breaks for plotting")
  }

  pl <- ggplot(data, aes_string(x = x, y = p, color = color)) +
    geom_point(...) +
    color_scale +
    theme_bw() +
    labs(title = title, subtitle = subtitle) +
    theme(panel.grid.minor.y = element_blank(), #only show the specified p-values, which might be unevenly spaced
          axis.ticks.y = element_blank())

  if(!is.null(p_fdr)) {

    if (any(data[, p_fdr] < fdr_limit)) {
      # Add horizontal line with the FDR < 0.05 limit
      q_limit <- max(data[data[, p_fdr] < fdr_limit, p], na.rm = TRUE)
      # sec_axis writes e.g. "q < 0.05" on the right sifde of the plot
      pl <- pl +
        geom_hline(yintercept = q_limit, linetype = "dashed") +
        scale_y_continuous(trans = minus_log10, breaks = p_breaks, labels = as.character(p_breaks),
                           sec.axis = sec_axis(~., breaks = q_limit, labels = paste("q =", fdr_limit)))
    } else {
      warning("None of the FDR-adjusted p-values are below the significance level, not plotting the horizontal line",
              call. = FALSE)
      pl <- pl +
        scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                           labels = as.character(p_breaks))
    }

  } else {
    pl <- pl +
      scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                         labels = as.character(p_breaks))
  }

  if (log2_x) {
    if (center_x_axis) {
      x_lim <- max(abs(log2(data[, x])))
      x_lim <- c(2^(-x_lim), 2^x_lim)
    }
    pl <- pl +
      scale_x_continuous(trans = "log2", limits = x_lim)
  } else {
    if (center_x_axis) {
      x_lim <- max(abs(data[, x]))
      x_lim <- c(-x_lim, x_lim)
    }
    pl <- pl +
      scale_x_continuous(limits = x_lim)
  }

  pl
}


# ---------- MANHATTAN PLOT -------


#' Manhattan plot
#'
#' Draws a (directed) Manhattan plot of p-values and versus e.g. retention time or mass-to-charge ratio.
#' If effect size and direction is supplied, the -log10(p-value) on the y-axis will be multiplied
#' by the direction (sign) of the effect, so part of the points will "drop" from the p = 1 (-log10(p) = 0) line.
#' This results in a so-called directed Manhattan plot.
#'
#' @param data a data frame with the effect size and p-values
#' @param x,p the column names of effect size (x-axis) and p-values
#' @param effect column name of effect size (should have negative and positive values).
#' @param p_fdr column name of FDR corrected p-values, used to draw a line showing the fdr-corrected significance level
#' @param color column name used to color the plots
#' @param p_breaks a numerical vector of the p_values to show on the y-axis
#' @param fdr_limit the significance level used in the experiment
#' @param x_lim,y_lim numerical vectors of length 2 for manually setting the axis limits
#' @param color_scale the color scale as returned by a ggplot function
#' @param title,subtitle the title and subtitle of the plot
#' @param ...  parameters passed to \code{\link[ggplot2]{geom_point}}, such as shape and alpha values. New aesthetics can
#' also be passed using \code{mapping = aes(...)}.
#'
#' @return a ggplot object
#'
#' @examples
#' # naturally, this looks messy as there are not enough p-values
#' lm_results <- perform_lm(drop_qcs(merged_sample), formula_char = "Feature ~ Group")
#' lm_data <- dplyr::left_join(fData(merged_sample), lm_results)
#' # Traditional Manhattan plot
#' manhattan_plot(data = lm_data, x = "Average.Mz",
#'              p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'              fdr_limit = 0.1)
#' # Directed Manhattan plot
#' manhattan_plot(data = lm_data, x = "Average.Mz", effect = "GroupB_Estimate",
#'              p = "GroupB_P", p_fdr = "GroupB_P_FDR",
#'              fdr_limit = 0.1)
#'
#' @export
manhattan_plot <- function(data, x, p, effect = NULL, p_fdr = NULL, color = NULL,
                           p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
                           x_lim = NULL, y_lim = NULL,
                           color_scale = getOption("notame.color_scale_con"),
                           title = "Manhattan plot", subtitle = NULL, ...) {

  if (min(data[, p]) > max(p_breaks)) {
    warning("All the p-values are larger than the p-value breaks supplied. Consider using larger p_breaks for plotting")
  }

  if (!is.null(effect)) {
    data$y <- -log10(data[, p]) * sign(data[, effect])
    p_labels <- outer(c(-1,1), p_breaks) %>% as.vector() %>% as.character()
    p_breaks <- outer(c(-1,1), -log10(p_breaks)) %>% as.vector()
    p_labels <- p_labels[order(p_breaks)]
    p_breaks <- sort(p_breaks)
  } else {
    data$y <- -log10(data[, p])
    p_labels <- as.character(p_breaks)
    p_breaks <- -log10(p_breaks)
    p_labels <- p_labels[order(p_breaks)]
    p_breaks <- sort(p_breaks)
  }

  pl <- ggplot(data, aes_string(x = x, y = "y", color = color)) +
    geom_point(...) +
    color_scale +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_hline(yintercept = 0, color = "grey") +
    labs(title = title, subtitle = subtitle, y = "p-value")


  if(!is.null(p_fdr)) {

    if (any(data[, p_fdr] < fdr_limit)) {
      # Add horizontal line with the FDR < 0.05 limit
      q_limit <- max(data[data[, p_fdr] < fdr_limit, p], na.rm = TRUE)
      # sec_axis writes e.g. "q < 0.05" on the right sifde of the plot
      if (!is.null(effect)) {
        pl <- pl +
          geom_hline(yintercept = log10(q_limit), linetype = "dashed") +
          geom_hline(yintercept = -log10(q_limit), linetype = "dashed") +
          scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim,
                             sec.axis = sec_axis(~., breaks = c(log10(q_limit), -log10(q_limit)), labels = rep(paste("q =", fdr_limit), 2)))
      } else {
        pl <- pl +
          geom_hline(yintercept = -log10(q_limit), linetype = "dashed") +
          scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim,
                             sec.axis = sec_axis(~., breaks = -log10(q_limit), labels = paste("q =", fdr_limit)))
      }
    } else {
      warning("None of the FDR-adjusted p-values are below the significance level, not plotting the horizontal line",
              call. = FALSE)
      pl <- pl +
        scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim)
    }
  } else {
    pl <- pl +
      scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim)
  }

  pl

}
