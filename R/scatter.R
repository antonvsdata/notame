
#' PCA scatter plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components
#'
#' @param object a MetaboSet object
#' @param center logical, should the data  be centered prior to PCA (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param method the method to use, see documentation in pcaMethods
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param density logical, whether to include density plots to both axes
#' @param color_scale the color scale as returned by a ggplot function
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return if \code{density} is \code{TRUE}, an object of class gTable (from the cowplot package).
#' NOTE this object must be drawn with ggdraw() instead of plot(). Otherwise, a ggplot object.
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca <- function(object, center = TRUE, scale = "uv", method = "ppca",
                     color = group_col(object), shape = NULL, density = FALSE,  title = "PCA",
                     subtitle = NULL, color_scale = NULL, shape_scale = NULL, ...) {

  res_pca <- pcaMethods::pca(object, method = method, scale = scale, center = center, ...)
  pca_scores <- pcaMethods::scores(res_pca) %>% as.data.frame()
  pca_scores[color] <- pData(object)[, color]
  R2 <- res_pca@R2
  labels <- paste0(c("PC1", "PC2"), " (", scales::percent(R2), ")")

  scatter_plot(pca_scores, x = "PC1", y = "PC2", xlab = labels[1], ylab = labels[2],
               color = color, shape = color, density = density, title = title,
               subtitle = subtitle, color_scale = color_scale, shape_scale = shape_scale)
}

#' t-SNE scatter plot
#'
#' Computes t-SNE into two dimensions and plots the map points.
#'
#' @param object a MetaboSet object
#' @param center logical, should the data  be centered prior to PCA (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param perplexity the perplexity used in t-SNE
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param density logical, whether to include density plots to both axes
#' @param color_scale the color scale as returned by a ggplot function
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param ... additional arguments passed to Rtsne::Rtsne
#'
#' @return if \code{density} is \code{TRUE}, an object of class gTable (from the cowplot package).
#' NOTE this object must be drawn with ggdraw() instead of plot(). Otherwise, a ggplot object.
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne <- function(object, center = TRUE, scale = "uv", perplexity = 30,
                      color = group_col(object), shape = NULL, density = FALSE, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), color_scale = NULL,
                      shape_scale = NULL, ...) {

  prepd <- pcaMethods::prep(object, center = center, scale = scale)
  res_tsne <- Rtsne::Rtsne(t(exprs(prepd)), perplexity = perplexity, ...)
  tsne_scores <- data.frame(res_tsne$Y)
  tsne_scores[color] <- pData(object)[, color]


  scatter_plot(tsne_scores, x = "X1", y = "X2", color = color, shape = color,
               density = density, title = title, subtitle = subtitle,
               color_scale = color_scale, shape_scale = shape_scale)

}

scatter_plot <- function(data, x, y, color, shape, density = FALSE, fixed = TRUE, color_scale = NULL,
                         shape_scale = NULL, title = NULL, subtitle = NULL, xlab = x, ylab = y,
                         color_lab = color, shape_lab = shape) {

  if (is.null(color_scale)) {
    if (class(data[, color]) %in% c("numeric", "integer")) {
      color_scale <- getOption("amp.color_scale_con")
    } else {
      color_scale <- getOption("amp.color_scale_dis")
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
                   alpha = 0.7, size = 0.2)

    ydens <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE)+
      geom_density(data = data, aes_string(x = y, fill = color),
                   alpha = 0.7, size = 0.2)+
      coord_flip()

    p <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p <- cowplot::insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
  }

  p
}

# --------------- HEXBIN PLOTS --------------------


#' PCA hexbin plot
#'
#' Computes PCA using one of the methods provided in the Bioconductor package
#' pcaMethods and plots the two first principal components
#'
#' @param object a MetaboSet object
#' @param center logical, should the data  be centered prior to PCA (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param method the method to use, see documentation in pcaMethods
#' @param fill character, name of the column used for coloring the hexagons
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param ... additional arguments passed to pcaMethods::pca
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link[pcaMethods]{pca}}
#'
#' @export
plot_pca_hexbin <- function(object, center = TRUE, scale = "uv", method = "ppca",
                     fill = "Injection_order", summary_fun = "mean", bins = 10, title = "PCA",
                     subtitle = NULL, fill_scale = NULL, ...) {

  res_pca <- pcaMethods::pca(object, method = method, scale = scale, center = center, ...)
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
#' Computes t-SNE into two dimensions and plots the map points.
#'
#' @param object a MetaboSet object
#' @param center logical, should the data  be centered prior to PCA (usually yes)
#' @param scale scaling used, as in pcaMethods::prep. Default is "uv" for unit variance
#' @param perplexity the perplexity used in t-SNE
#' @param color character, name of the column used for coloring the points
#' @param fill character, name of the column used for coloring the hexagons
#' @param fill_scale the fill scale as returned by a ggplot function
#' @param summary_fun the function used to compute the value for each hexagon
#' @param bins the number of bins in x and y axes
#' @param ... additional arguments passed to Rtsne::Rtsne
#'
#' A ggplot object.
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#'
#' @export
plot_tsne_hexbin <- function(object, center = TRUE, scale = "uv", perplexity = 30,
                      fill = "Injection_order", summary_fun = "mean", bins = 10, title = "t-SNE",
                      subtitle = paste("Perplexity:", perplexity), fill_scale = NULL, ...) {

  prepd <- pcaMethods::prep(object, center = center, scale = scale)
  res_tsne <- Rtsne::Rtsne(t(exprs(prepd)), perplexity = perplexity, ...)
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

