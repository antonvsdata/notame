

plot_pca <- function(object, color = group_col(object), method = "ppca", center = TRUE, scale = "uv",
                     title = "PCA", subtitle = NULL, color_scale = NULL, shape_scale = NULL, ...) {

  res_pca <- pcaMethods::pca(object, method = method, scale = scale, center = center, ...)
  pca_scores <- pcaMethods::scores(res_pca) %>% as.data.frame()
  pca_scores[color] <- pData(object)[, color]
  R2 <- res_pca@R2
  labels <- paste0(c("PC1", "PC2"), " (", scales::percent(R2), ")")

  scatter_plot(pca_scores, x = "PC1", y = "PC2", xlab = labels[1], ylab = labels[2], color = color, shape = color, title = title,
               subtitle = subtitle, color_scale = color_scale, shape_scale = shape_scale)
}


plot_tsne <- function(object, color = group_col(object), perplexity = 30, center = TRUE, scale = "uv",
                      title = "t-SNE", subtitle = NULL, color_scale = NULL, shape_scale = NULL, ...) {

  prepd <- pcaMethods::prep(object, center = center, scale = scale)
  res_tsne <- Rtsne::Rtsne(t(exprs(prepd)), perplexity = perplexity, ...)
  tsne_scores <- data.frame(res_tsne$Y)
  tsne_scores[color] <- pData(object)[, color]


  scatter_plot(tsne_scores, x = "X1", y = "X2", color = color, shape = color, title = title,
               subtitle = subtitle, color_scale = color_scale, shape_scale = shape_scale)

}

scatter_plot <- function(data, x, y, color, shape, fixed = TRUE, color_scale = NULL, shape_scale = NULL,
                         title = NULL, subtitle = NULL, xlab = x, ylab = y,
                         color_lab = color, shape_lab = shape) {

  if (is.null(color_scale)) {
    if (class(data[, color]) == "numeric") {
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

    } else {
      cat("Only 8 distinct shapes currently available!")
      p <- p +
        geom_point()
    }
  } else {
    p <- p +
      geom_point()
  }

  p
}
