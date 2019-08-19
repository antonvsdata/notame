
#' PCA scatter plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components
#' \strong{CITATION:} When using this function, cite the \code{pcaMethods} package
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param center logical, should the data be centered prior to PCA? (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param density logical, whether to include density plots to both axes
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param fill_scale the fill scale used for density curves
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return a ggplot object. If \code{density} is \code{TRUE}, the plot will consist of multiple
#' parts and is harder to modify.
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                     color = group_col(object), shape = NULL, density = FALSE,  title = "PCA",
                     subtitle = NULL, color_scale = NULL, shape_scale = NULL, fill_scale = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  shape <- shape %||% color

  res_pca <- pcaMethods::pca(object, scale = scale, center = center, ...)
  pca_scores <- pcaMethods::scores(res_pca) %>% as.data.frame()
  pca_scores[color] <- pData(object)[, color]
  pca_scores[shape] <- pData(object)[, shape]
  R2 <- res_pca@R2
  labels <- paste0(c("PC1", "PC2"), " (", scales::percent(R2), ")")


  scatter_plot(pca_scores, x = "PC1", y = "PC2", xlab = labels[1], ylab = labels[2],
               color = color, shape = shape, density = density, title = title,
               subtitle = subtitle, color_scale = color_scale, shape_scale = shape_scale,
               fill_scale = fill_scale)
}

#' t-SNE scatter plot
#'
#' Computes t-SNE into two dimensions and plots the map points.
#' In case there are missing values, PCA is performed using the nipals method of \code{pcaMethods::pca},
#' the  method can be changed to "ppca" if niipals fails.
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
#' @param density logical, whether to include density plots to both axes
#' @param title,subtitle the titles of the plot
#' @param color_scale the color scale as returned by a ggplot function
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param fill_scale the fill scale used for density curves
#' @param ... additional arguments passed to \code{Rtsne::Rtsne}
#'
#' @return a ggplot object. If \code{density} is \code{TRUE}, the plot will consist of multiple
#' parts and is harder to modify.
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne <- function(object, all_features = FALSE, center = TRUE, scale = "uv", perplexity = 30,
                      pca_method = "nipals",
                      color = group_col(object), shape = NULL, density = FALSE, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), color_scale = NULL,
                      shape_scale = NULL, fill_scale = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  shape <- shape %||% color
  prepd <- pcaMethods::prep(object, center = center, scale = scale)

  # If there are missing values, use ppca method from pcaMethods instead of usual PCA
  if (sum(is.na(exprs(prepd))) > 0) {
    res_pca <- pcaMethods::pca(object, method = pca_method, nPcs = min(nrow(object), ncol(object), 50), scale = "none", center = FALSE)
    pca_scores <- pcaMethods::scores(res_pca)
    res_tsne <- Rtsne::Rtsne(pca_scores, perplexity = perplexity, pca = FALSE, ...)
  } else {
    res_tsne <- Rtsne::Rtsne(t(exprs(prepd)), perplexity = perplexity, ...)
  }

  tsne_scores <- data.frame(res_tsne$Y)
  tsne_scores[color] <- pData(object)[, color]
  tsne_scores[shape] <- pData(object)[, shape]

  scatter_plot(tsne_scores, x = "X1", y = "X2", color = color, shape = shape,
               density = density, title = title, subtitle = subtitle,
               color_scale = color_scale, shape_scale = shape_scale, fill_scale = fill_scale)

}

scatter_plot <- function(data, x, y, color, shape, density = FALSE, fixed = TRUE, color_scale = NULL,
                         shape_scale = NULL, fill_scale = NULL, title = NULL, subtitle = NULL, xlab = x, ylab = y,
                         color_lab = color, shape_lab = shape) {

  if (is.null(color_scale)) {
    if (class(data[, color]) %in% c("numeric", "integer")) {
      color_scale <- getOption("amp.color_scale_con")
    } else {
      color_scale <- getOption("amp.color_scale_dis")
    }
  }

  if (is.null(fill_scale)) {
    if (class(data[, color]) %in% c("numeric", "integer")) {
      fill_scale <- getOption("amp.fill_scale_con")
    } else {
      fill_scale <- getOption("amp.fill_scale_dis")
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

  if (class(data[, shape]) == "factor") {
    if (length(levels(data[, shape])) <= 8){
      shape_scale <- shape_scale %||% getOption("amp.shape_scale")
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
  # Add density plots to top and right
  if (density) {
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
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_hexbin <- function(object, all_features = FALSE, center = TRUE, scale = "uv",
                     fill = "Injection_order", summary_fun = "mean", bins = 10, title = "PCA",
                     subtitle = NULL, fill_scale = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  res_pca <- pcaMethods::pca(object, scale = scale, center = center, ...)
  pca_scores <- pcaMethods::scores(res_pca) %>% as.data.frame()
  pca_scores[fill] <- pData(object)[, fill]
  R2 <- res_pca@R2
  labels <- paste0(c("PC1", "PC2"), " (", scales::percent(R2), ")")

  hexbin_plot(data = pca_scores, x = "PC1", y = "PC2", xlab = labels[1], ylab = labels[2],
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
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne_hexbin <- function(object, all_features = FALSE, center = TRUE, scale = "uv", pca_method = "nipals", perplexity = 30,
                      fill = "Injection_order", summary_fun = "mean", bins = 10, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), fill_scale = NULL, ...) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  prepd <- pcaMethods::prep(object, center = center, scale = scale)

  if (sum(is.na(exprs(prepd))) > 0) {
    res_pca <- pcaMethods::pca(object, method = pca_method, nPcs = min(nrow(object), ncol(object), 50), scale = "none", center = FALSE)
    pca_scores <- pcaMethods::scores(res_pca)
    res_tsne <- Rtsne::Rtsne(pca_scores, perplexity = perplexity, pca = FALSE, ...)
  } else {
    res_tsne <- Rtsne::Rtsne(t(exprs(prepd)), perplexity = perplexity, ...)
  }
  tsne_scores <- data.frame(res_tsne$Y)
  tsne_scores[fill] <- pData(object)[, fill]


  hexbin_plot(tsne_scores, x = "X1", y = "X2",
               fill = fill, summary_fun = summary_fun, bins = bins,
               fill_scale = fill_scale,
               title = title, subtitle = subtitle)

}


hexbin_plot <- function(data, x, y, fill, summary_fun = "mean", bins = 10, fill_scale = NULL,
                         title = NULL, subtitle = NULL, xlab = x, ylab = y,
                         fill_lab = fill) {

  fill_scale <- fill_scale %||% getOption("amp.fill_scale_con")

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


minus_log10 <- scales::trans_new("minus_log19",
                                 transform = function(x) {-log10(x)},
                                 inverse = function(x) {10^{-x}})

#' Volcano plot
#'
#' Draws a volcano plot of effect size and p-values.
#'
#' @param data a data frame with the effect size and p-values
#' @param effect,p the column names of effect size and p-values
#' @param log2_effect logical, whether effect size should be plotted on a log2 axis
#' @param center_x_axis logical, whether x-axis should be centered. If \code{TRUE}, the "zero-effect" will
#' be on the middle of the plot. The "zero effect" is 0 if \code{log2_effect = FALSE} and 1 if  \code{log2_effect = TRUE}
#' @param title,subtitle the title and subtitle of the plot
#' @param ...  parameters passed to \code{\link[ggplot2]{geom_point}}, such as shape and alpha values.
#'
#' @return a ggplot object
#'
#' @export
volcano_plot <- function(data, effect, p, log2_effect = FALSE, center_x_axis = TRUE,
                         title = "Volcano plot", subtitle = NA, ...) {

  p <- ggplot(data, aes_string(x = effect, y = p)) +
    geom_point(...) +
    theme_bw() +
    scale_y_continuous(trans = minus_log10) +
    labs(title = title, subtitle = subtitle)

  if (log2_effect) {
    if (center_x_axis) {
      x_lim <- max(abs(log2(data[, effect])))
      x_lim <- c(2^(-x_lim), 2^x_lim)
    } else {
      x_lim = NULL
    }
    p <- p +
      scale_x_continuous(trans = "log2", limits = x_lim)
  } else {
    if (center_x_axis) {
      x_lim <- max(abs(data[, effect]))
      x_lim <- c(-x_lim, x_lim)
    } else {
      x_lim = NULL
    }
    p <- p +
      scale_x_continuous(limits = x_lim)
  }

  p
}
