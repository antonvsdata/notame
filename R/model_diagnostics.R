
#' Plot model diagnostics
#'
#' Plot model diagnostics for lm or lmer models fit separately for each feature.
#' One page of plots is saved per feature. Uses autoplot function from ggfortify
#'
#' @param object a MetaboSet object
#' @param formula_char character, the formula to be used in the linear model (see Details)
#' @param model_type the type of model, either "lm" or "lmer"
#' @param file the file where the plots will be saved
#' @param all_features should all features be included?
#' @param width,height the dimension of the plot
#' @param ... other parameters passed to lm/lmer and/or autoplot function from the ggfortify package
#'
#' @examples
#' save_lm_diagnostic_plots(example_set, formula_char = "Feature ~ Group", model_type = "lm", file = "test.pdf")
#'
#' # Add colouring by time point (don't mind the warnings)
#' save_lm_diagnostic_plots(example_set, formula_char = "Feature ~ Group",
#'                          model_type = "lm", file = "test.pdf", colour = "Group")
#'
#' # Linear mixed model
#' save_lm_diagnostic_plots(example_set, formula_char = "Feature ~ Group + (1|Subject_ID)",
#'                          model_type = "lmer", file = "test.pdf")
#'
#'
#' @export
save_lm_diagnostic_plots <- function(object, formula_char, model_type = c("lm", "lmer"), file, all_features = FALSE,
                                     width = 12, height = 18, ...) {

  require(ggfortify)
  if (model_type == "lmer") {
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package \"lme4\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
  }

  if (missing(model_type)) model_type <- "lm"
  model_type <- match.arg(model_type)
  # Start log
  log_text(paste("\nStarting linear model diagnostics at", Sys.time()))

  object <- drop_flagged(object, all_features = all_features)
  data <- combined_data(object)
  features <- Biobase::featureNames(object)

  pdf(file = file, width = width, height = height)

  for (i in seq_along(features)) {
    if (i %% 100 == 0){
      print(paste(i, "/", length(features)))
    }
    feature <- features[i]
    # Replace "Feature" with the current feature name
    tmp_formula <- gsub("Feature", feature, formula_char)

    fit <- NULL
    if (model_type == "lm") {
      tryCatch({
        fit <- lm(tmp_formula, data = data, ...)
      }, error = function(e) print(paste0(feature, ": ", e$message)))
    } else if (model_type == "lmer") {
      tryCatch({
        fit <- lme4::lmer(tmp_formula, data = data, ...)
      }, error = function(e) print(paste0(feature, ": ", e$message)))
    }

    if (!is.null(fit)) {
      tryCatch({
        p <- autoplot(object = fit, which = 1:6, data = data) +
          theme_bw()
        print(p)
      })
    }
  }
  dev.off()
  log_text(paste("Saved linear model diagnostic plots to:", file))

}


#' Flatten dataframe contains matrix
#'
#' tains list or matrix as column
#'
#' @param df \code{data.frame} to be flatten
flatten <- function(df) {
  ismatrix <- vapply(df, is.matrix, logical(1))
  if (any(ismatrix)) {
    return(data.frame(c(df[!ismatrix], do.call(data.frame, df[ismatrix])),
                      stringsAsFactors = FALSE))
  } else {
    return(df)
  }
}


#' Calcurate layout matrix for \code{ggmultiplot}
#'
#' @param nplots Number of plots
#' @param ncol Number of grid columns
#' @param nrow Number of grid rows
#' @return matrix
#' @examples
#' ggfortify:::get.layout(3, 2, 2)
get.layout <- function(nplots, ncol, nrow) {
  if (ncol == 0 && nrow == 0) {
    ncol <- 2
  } else if (ncol == 0 && nrow != 0) {
    ncol <- ceiling(nplots / nrow)
  }

  if (nrow == 0) {
    nrow <- ceiling(nplots / ncol)
  } else {
    nrow <- nrow
  }

  if (nrow * ncol < nplots) {
    message <- paste('nrow * ncol (', nrow, ' * ', ncol,
                     ') must be larger than number of plots', nplots)
    stop(message)
  }

  t(matrix(1:(ncol * nrow), ncol = nrow, nrow = ncol))
}


#' Attach label to \code{ggplot2::ggplot}
#'
#' @param p \code{ggplot2::ggplot} instance
#' @param data Data contains text label
#' @param x x coordinates for label
#' @param y y coordinates for label
#' @param label Logical value whether to display labels
#' @param label.label Column name used for label text
#' @param label.colour Colour for text labels
#' @param label.alpha Alpha for text labels
#' @param label.size Size for text labels
#' @param label.angle Angle for text labels
#' @param label.family Font family for text labels
#' @param label.fontface Fontface for text labels
#' @param label.lineheight Lineheight for text labels
#' @param label.hjust Horizontal adjustment for text labels
#' @param label.vjust Vertical adjustment for text labels
#' @param label.repel Logical flag indicating whether to use \code{ggrepel}, enabling this may take some time for plotting
#' @param label.show.legend Logical value indicating whether to show the legend of the text labels
#' @return ggplot
plot_label <- function(p, data, x = NULL, y = NULL, label = TRUE, label.label = 'rownames',
                       label.colour = NULL, label.alpha = NULL,
                       label.size = NULL, label.angle = NULL,
                       label.family = NULL, label.fontface = NULL,
                       label.lineheight = NULL,
                       label.hjust = NULL, label.vjust = NULL,
                       label.repel = FALSE, label.show.legend = NA) {

  if (!is.data.frame(data)) {
    stop(paste0('Unsupported class: ', class(data)))
  }

  if (!missing(label.colour) && !is.null(label.colour) && missing(label)) {
    # if flag is missing but colour is specified, turn flag to TRUE
    label <- TRUE
  }

  if (label || label.repel) {
    # user wants label if they enables repel
    if (is.null(label.colour)) {
      # NULL may be explicitly passed from parent functions
      label.colour <- '#000000'
    }
    if (label.repel && 'ggrepel' %in% rownames(installed.packages())) {
      textfunc <- ggrepel::geom_text_repel
    } else {
      textfunc <- ggplot2::geom_text
    }
    p <- p + geom_factory(textfunc, data, x = x, y = y,
                          label = label.label,
                          colour = label.colour, alpha = label.alpha,
                          size = label.size, angle = label.angle,
                          family = label.family, fontface = label.fontface,
                          lineheight = label.lineheight,
                          hjust = label.hjust, vjust = label.vjust,
                          show.legend = label.show.legend)
  }
  p
}


#' Factory function to control \code{ggplot2::geom_xxx} functions
#'
#' @param geomfunc \code{ggplot2::geom_xxx} function
#' @param data plotting data
#' @param ... other arguments passed to methods
#' @return proto
geom_factory <- function(geomfunc, data = NULL, ...) {
  mapping <- list()
  option <- list()

  columns <- colnames(data)
  for (key in names(list(...))) {
    value <- list(...)[[key]]
    if (is.null(value)) {
      # pass
    } else if (value %in% columns) {
      mapping[[key]] <- value
    } else {
      option[[key]] <- value
    }
  }
  if (!is.null(data)) {
    option[['data']] <- data
  }
  option[['mapping']] <- do.call(ggplot2::aes_string, mapping)
  return(do.call(geomfunc, option))
}





#' Autoplot \code{lme4::lmerMod}
#'
#' @param object \code{lme4::lmerMod} instance
#' @param which If a subset of the plots is required, specify a subset of the numbers 1:6.
#' @param data original dataset, if needed
#' @param colour line colour
#' @param size point size
#' @param linetype line type
#' @param alpha alpha
#' @param fill fill colour
#' @param shape point shape
#' @param label.n Number of points to be laeled in each plot, starting with the most extreme
#' @inheritParams plot_label
#' @param smooth.colour Line colour for smoother lines
#' @param smooth.linetype Line type for smoother lines
#' @param ad.colour Line colour for additional lines
#' @param ad.linetype Line type for additional lines
#' @param ad.size Fill colour for additional lines
#' @param nrow Number of facet/subplot rows
#' @param ncol Number of facet/subplot columns
#' @param ... other arguments passed to methods
#' @return ggplot
#' @examples
#' \dontrun{
#' autoplot(lm(Petal.Width ~ Petal.Length, data = iris))
#' autoplot(glm(Petal.Width ~ Petal.Length, data = iris), which = 1:6)
#' autoplot(lm(Petal.Width~Petal.Length, data = iris), data = iris, colour = 'Species')
#' }
#' @export
autoplot.lmerMod <- function(object, which = c(1:3, 5), data = NULL,
                        colour = '#444444', size = NULL, linetype = NULL,
                        alpha = NULL, fill = NULL, shape = NULL,
                        label = TRUE, label.label = '.label',
                        label.colour = '#000000', label.alpha = NULL,
                        label.size = NULL, label.angle = NULL,
                        label.family = NULL, label.fontface = NULL,
                        label.lineheight = NULL,
                        label.hjust = NULL, label.vjust = NULL, label.repel = FALSE,
                        label.n = 3,
                        smooth.colour = '#0000FF', smooth.linetype = 'solid',
                        ad.colour = '#888888', ad.linetype = 'dashed', ad.size = .2,
                        nrow = NULL, ncol = NULL, ...) {
  # initialization
  p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- NULL

  dropInf <- function(x, h) {
    if (any(isInf <- h >= 1)) {
      warning(gettextf("not plotting observations with leverage one:\n  %s",
                       paste(which(isInf), collapse = ", ")), call. = FALSE,
              domain = NA)
      x[isInf] <- NaN
    }
    x
  }

  show <- rep(FALSE, 6)
  show[which] <- TRUE

  if (is.null(data)) {
    # ggplot2::fortify can't handle NULL properly
    plot.data <- ggplot2::fortify(object)
  } else {
    plot.data <- ggplot2::fortify(object, data = data)
  }
  # Filling in missing metrics
  plot.data$.cooksd <- cooks.distance(object)
  plot.data$.stdresid <- plot.data$.resid / sd(plot.data$.resid)
  plot.data$.hat <- hatvalues(object)
  n <- nrow(plot.data)

  plot.data$.index <- 1:n
  plot.data$.label <- rownames(plot.data)

  is_glm <- inherits(object, "glm")
  r <- residuals(object)
  w <- weights(object)
  if (any(show[2L:6L])) {
    s <- if (inherits(object, "rlm")) {
      object$s
    } else if (is_glm) {
      sqrt(summary(object)$dispersion)
    } else {
      sqrt(stats::deviance(object, REML = FALSE) / stats::df.residual(object))
    }
    hii <- hatvalues(object)

    is_const_lev <- all(hii == 0) ||
      diff(hii) < 1e-10 * mean(hii, na.rm = TRUE)

    fs <- dplyr::select_if(plot.data,
                           function(x) is.character(x) | is.factor(x))
    fs[[".label"]] <- NULL
    if (is_const_lev & ncol(fs) > 0){
      plot.data$.nf <- stringr::str_wrap(interaction(fs, sep = ":"), width = 10)
    }

    if (any(show[2L:3L])) {
      plot.data$.wresid <- if (is.null(w)) {
        r
      } else {
        sqrt(w) * r
      }
      plot.data$.wstdresid <- plot.data$.wresid / (s * sqrt(1 - hii))
    }
    if (show[2L]) {
      ylim <- range(plot.data$.wstdresid, na.rm = TRUE)
      ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      qn <- stats::qqnorm(plot.data$.wstdresid, ylim = ylim, plot.it = FALSE)
      plot.data$.qqx <- qn$x
      plot.data$.qqy <- qn$y
    }
  }

  label.fitted <- ifelse(is_glm, 'Predicted values', 'Fitted values')
  label.y23 <- ifelse(is_glm, 'Std. deviance resid.', 'Standardized residuals')

  if (is.logical(shape) && !shape) {
    if (missing(label)) {
      # if label is missing and shape=FALSE, turn label to TRUE
      label <- TRUE
    }
    if (missing(label.n)) {
      label.n <- nrow(plot.data)
    }
  }

  # target or explanatory variables may be nested, e.g. binomial
  plot.data <- flatten(plot.data)

  if (label.n > 0L) {
    if (show[1L]) {
      r.data <- dplyr::arrange_(plot.data, 'dplyr::desc(abs(.resid))')
      r.data <- utils::head(r.data, label.n)
    }
    if (".wresid" %in% colnames(plot.data)) {
      wr.data <- dplyr::arrange_(plot.data, 'dplyr::desc(abs(.wresid))')
      wr.data <- utils::head(wr.data, label.n)
    }
    if (any(show[4L:6L])) {
      cd.data <- dplyr::arrange_(plot.data, 'dplyr::desc(abs(.cooksd))')
      cd.data <- utils::head(cd.data, label.n)
    }
  }

  .smooth <- function(x, y) {
    stats::lowess(x, y, f = 2 / 3, iter = 3)
  }

  .decorate.label <- function(p, data) {
    if (label & label.n > 0) {
      p <- plot_label(p = p, data = data,
                      label = label,
                      label.label = label.label,
                      label.colour = label.colour,
                      label.alpha = label.alpha,
                      label.size = label.size,
                      label.angle = label.angle,
                      label.family = label.family,
                      label.fontface = label.fontface,
                      label.lineheight = label.lineheight,
                      label.hjust = label.hjust,
                      label.vjust = label.vjust,
                      label.repel = label.repel)
    }
    p
  }

  .decorate.plot <- function(p, xlab = NULL, ylab = NULL, title = NULL) {
    p +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::ggtitle(title)
  }

  smoother_m <- ggplot2::aes_string(x = 'x', y = 'y')

  if (show[1L]) {
    t1 <- 'Residuals vs Fitted'
    mapping <- ggplot2::aes_string(x = '.fitted', y = '.resid')
    smoother <- .smooth(plot.data$.fitted, plot.data$.resid)
    smoother <- as.data.frame(smoother)
    p1 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p1 <- p1 + geom_factory(geom_point, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p1 <- p1 +
      ggplot2::geom_line(data = smoother, mapping = smoother_m,
                         colour = smooth.colour, linetype = smooth.linetype) +
      ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype,
                          size = ad.size, colour = ad.colour)
    p1 <- .decorate.label(p1, r.data)
    p1 <- .decorate.plot(p1, xlab = label.fitted, ylab = 'Residuals', title = t1)
  }

  if (show[2L]) {
    t2 <- 'Normal Q-Q'
    qprobs <- c(0.25, 0.75)

    qy <- stats::quantile(plot.data$.wstdresid, probs = qprobs, names = FALSE,
                          type = 7, na.rm = TRUE)
    qx <- stats::qnorm(qprobs)
    slope <- diff(qy) / diff(qx)
    int <- qy[1L] - slope * qx[1L]

    mapping <- ggplot2::aes_string(x = '.qqx', y = '.qqy')
    p2 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      # Do not use stat_qq here for labeling
      p2 <- p2 + geom_factory(geom_point, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p2 <- p2 + ggplot2::geom_abline(intercept=int, slope=slope,
                                    linetype = ad.linetype, size = ad.size,
                                    colour = ad.colour)
    p2 <- .decorate.label(p2, wr.data)
    p2 <- .decorate.plot(p2, xlab = 'Theoretical Quantiles',
                         ylab = label.y23, title = t2)
  }

  if (show[3L]) {
    t3 <- 'Scale-Location'
    mapping <- ggplot2::aes_string(x = '.fitted', y = 'sqrt(abs(.wstdresid))')
    smoother <- .smooth(plot.data$.fitted, sqrt(abs(plot.data$.wstdresid)))
    smoother <- as.data.frame(smoother)
    p3 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p3 <- p3 + geom_factory(geom_point, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p3 <- p3 + ggplot2::geom_line(data = smoother, mapping = smoother_m,
                                  colour = smooth.colour, linetype = smooth.linetype)
    p3 <- .decorate.label(p3, wr.data)
    label.y3 <- ifelse(is_glm, expression(sqrt(abs(`Std. deviance resid.`))),
                       expression(sqrt(abs(`Standardized residuals`))))
    p3 <- .decorate.plot(p3, xlab = label.fitted, ylab = label.y3,
                         title = t3)
  }

  if (show[4L]) {
    t4 <- "Cook's distance"
    mapping <- ggplot2::aes_string(x = '.index', y = '.cooksd',
                                   ymin = 0, ymax = '.cooksd')
    p4 <-  ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p4 <- p4 + geom_factory(geom_linerange, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p4 <- .decorate.label(p4, cd.data)
    p4 <- .decorate.plot(p4, xlab = 'Obs. Number',
                         ylab = "Cook's distance", title = t4)
  }

  if (show[5L]) {
    if (is_const_lev & ncol(fs) > 0){
      t5 <- 'Constant Leverage:\nResiduals vs Factor Levels'
      mapping <- ggplot2::aes_string(x = '.nf', y = '.stdresid')

      p5 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
      if (!is.logical(shape) || shape) {
        p5 <- p5 + geom_factory(geom_point, plot.data,
                                colour = colour, size = size, linetype = linetype,
                                alpha = alpha, fill = fill, shape = shape)
      }
      p5 <- p5 +
        ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype,
                            size = ad.size, colour = ad.colour) +
        ggplot2::expand_limits(x = 0)
      p5 <- .decorate.label(p5, cd.data)
      label.y5 <- ifelse(is_glm, 'Std. Pearson resid.', 'Standardized Residuals')
      p5 <- .decorate.plot(p5, xlab = 'Factor Level Combination',
                           ylab = label.y5, title = t5)
    } else {
      t5 <- 'Residuals vs Leverage'
      mapping <- ggplot2::aes_string(x = '.hat', y = '.stdresid')
      smoother <- .smooth(plot.data$.hat, plot.data$.stdresid)
      smoother <- as.data.frame(smoother)
      p5 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
      if (!is.logical(shape) || shape) {
        p5 <- p5 + geom_factory(geom_point, plot.data,
                                colour = colour, size = size, linetype = linetype,
                                alpha = alpha, fill = fill, shape = shape)
      }
      p5 <- p5 + ggplot2::geom_line(data = smoother, mapping = smoother_m,
                                    colour = smooth.colour, linetype = smooth.linetype) +
        ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype,
                            size = ad.size, colour = ad.colour) +
        ggplot2::expand_limits(x = 0)
      p5 <- .decorate.label(p5, cd.data)
      label.y5 <- ifelse(is_glm, 'Std. Pearson resid.', 'Standardized Residuals')
      p5 <- .decorate.plot(p5, xlab = 'Leverage', ylab = label.y5, title = t5)
    }
  }

  if (show[6L]) {
    t6 <- "Cook's dist vs Leverage"
    mapping <- ggplot2::aes_string(x = '.hat', y = '.cooksd')
    smoother <- .smooth(plot.data$.hat, plot.data$.cooksd)
    smoother <- as.data.frame(smoother)
    p6 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p6 <- p6 + geom_factory(geom_point, plot.data,
                              colour = colour, size = size, linetype = linetype,
                              alpha = alpha, fill = fill, shape = shape)
    }
    p6 <- p6 + ggplot2::geom_line(data = smoother, mapping = smoother_m,
                                  colour = smooth.colour, linetype = smooth.linetype) +
      ggplot2::expand_limits(x = 0, y = 0)
    p6 <- .decorate.label(p6, cd.data)
    p6 <- .decorate.plot(p6, xlab = 'Leverage', ylab = "Cook's distance",
                         title = t6)

    g <- dropInf(hii / (1 - hii), hii)
    p <- length(stats::coef(object))
    bval <- pretty(sqrt(p * plot.data$.cooksd / g), 5)
    for (i in seq_along(bval)) {
      bi2 <- bval[i] ^ 2
      p6 <- p6 + ggplot2::geom_abline(intercept=0, slope=bi2,
                                      linetype = ad.linetype, size = ad.size,
                                      colour = ad.colour)
    }
  }

  if (is.null(ncol)) { ncol <- 0 }
  if (is.null(nrow)) { nrow <- 0 }
  plot.list <- list(p1, p2, p3, p4, p5, p6)[which]
  new('ggmultiplot', plots = plot.list, nrow = nrow, ncol = ncol)
}


