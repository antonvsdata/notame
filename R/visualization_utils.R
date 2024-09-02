#' Save plot to chosen format
#'
#' DEPRECATED: Please use \code{\link[ggsave]{ggsave}} instead.
#'
#' @export
save_plot <- function(p, file, ...) {
  stop("This function is deprecated, please use ggsave instead.")
}

#' Write all relevant visualizations to pdf
#'
#' A wrapper around all the major visualization functions, used for visualizing data between
#' major steps of data preprocessing. Saves all visualizations as PDFs with a set prefix on filenames.
#'
#' @param object A MetaboSet object
#' @param prefix character, a file path prefix added to the file paths
#' @param format character, format in which the plots should be saved, DOES NOT support raster formats
#' @param perplexity perplexity for t-SNE plots
#' @param merge logical, whether the files should be merged to a single PDF, see Details
#' @param remove_singles logical, whether to remove single plot files after merging.
#' Only used if \code{merge = TRUE}
#'
#' @details If \code{merge} is \code{TRUE} and \code{format} id \code{pdf},
#' then a file containing all the visualizations named \code{prefix.pdf} will be created.
#' NOTE: on Windows this requires installation of pdftk (\url{https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/})
#' and on Linux you need to have pdfunite installed.
#' On MacOS, no external software is needed. Note that at least on Windows, prefix should be a path from the root,
#' so that the underlying system command will find the files.
#' The type of visualizations to be saved depends on the type of object.
#' Here is a comprehensive list of the visualizations:
#' \itemize{
#' \item Distribution of quality metrics and flags \code{\link{plot_quality}}
#' \item Boxplots of each sample in injection order \code{\link{plot_sample_boxplots}}
#' \item PCA scores plot of samples colored by injection order \code{\link{plot_pca}}
#' \item t-SNE plot of samples colored by injection order \code{\link{plot_tsne}}
#' \item If the object has over 60 samples, hexbin versions of the PCA and t-SNE plots above
#' \code{\link{plot_pca_hexbin}}, \code{\link{plot_tsne_hexbin}}
#' \item Dendrogram of samples ordered by hierarchical clustering, sample labels colored by group if present
#' \code{\link{plot_dendrogram}}
#' \item heat map of intersample distances, ordered by hierarchical clustering \code{\link{plot_sample_heatmap}}
#' \item If the object has QC samples: \itemize{
#' \item Density function of the intersample distances in both QCs and biological samples
#' \code{\link{plot_dist_density}}
#' \item Histograms of p-values from linear regression of features against injection order
#' in both QCs and biological samples \code{\link{plot_p_histogram}}}
#' \item If the object has a group column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by group \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' }
#' \item If the object has a time column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by time \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' \item Dendrogram of samples ordered by hierarchical clustering, sample labels colored by time point
#' \code{\link{plot_dendrogram}}
#' }
#' \item If the object has a group column OR a time column: \itemize{
#' \item Boxplots of samples ordered and colored by group and/or time \code{\link{plot_sample_boxplots}}
#' }
#' \item If the object has a group column AND a time column: \itemize{
#' \item PCA and tSNE plots with points shaped by group and colored by time
#' \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' }
#' \item If the object has a time column AND a subject column: \itemize{
#' \item PCA and tSNE plots with arrows connecting the samples of each subject in time point order
#' \code{\link{plot_pca_arrows}}, \code{\link{plot_tsne_arrows}}
#' }
#' }
#'
#' @export
visualizations <- function(object,
                           prefix,
                           format = "pdf",
                           perplexity = 30,
                           merge = FALSE,
                           remove_singles = FALSE) {
  # Keep all plots in a list
  plots <- list()
  # Helper function for handling errors and keeping track of file names
  save_name <- function(fun, name, ...) {
    p <- NULL
    tryCatch(
      {
        p <- fun(object, ...)
      },
      error = function(e) {
        futile.logger::flog.warn(
          paste0(
            "Error on plot ", name, ": ", geterrmessage()
          ),
          name = "notame"
        )
      }
    )

    if (!is.null(p)) {
      plots[[name]] <<- p
    }
  }

  if (sum(object$QC == "QC")) {
    save_name(
      fun = plot_dist_density,
      name = "density_plot"
    )
    save_name(plot_injection_lm, "lm_p_histograms")
  }
  # Quality metrics
  save_name(plot_quality, "quality_metrics")

  # Plots with injection order
  save_name(plot_sample_boxplots, "boxplots_injection",
    order_by = "Injection_order", fill_by = "QC"
  )
  set.seed(38)
  save_name(plot_pca, "PCA_injection", color = "Injection_order")
  set.seed(38)
  save_name(plot_tsne, "tSNE_injection",
    perplexity = perplexity,
    color = "Injection_order"
  )
  # Clustering
  save_name(plot_dendrogram, "dendrogram")
  save_name(plot_sample_heatmap, "heatmap_samples")

  # For large sets, plot hexbin plots
  if (ncol(object) > 60) {
    set.seed(38)
    save_name(plot_pca_hexbin, "PCA_hexbin")
    set.seed(38)
    save_name(plot_tsne_hexbin, "tSNE_hexbin", perplexity = perplexity)
  }

  # If not grouped, plot PCA and t-SNE on QC information
  if (is.na(group_col(object))) {
    group_col(object) <- "QC"
  }
  set.seed(38)
  save_name(plot_pca, "PCA_group")

  set.seed(38)
  save_name(plot_tsne, "tSNE_group", perplexity = perplexity)
  # Time point
  if (!is.na(time_col(object))) {
    set.seed(38)
    save_name(plot_pca, "PCA_time", color = time_col(object))

    set.seed(38)
    save_name(plot_tsne, "tSNE_time",
      color = time_col(object),
      perplexity = perplexity
    )
    save_name(plot_dendrogram, "dendrogram_time",
      color = time_col(object)
    )
  }
  # Time point OR group
  if (!is.na(group_col(object)) || !is.na(time_col(object))) {
    save_name(plot_sample_boxplots, "boxplots_group")
  }
  # Time point AND group
  if (!is.na(group_col(object)) && !is.na(time_col(object))) {
    set.seed(38)
    save_name(plot_pca, "PCA_group_time",
      color = time_col(object),
      shape = group_col(object)
    )

    set.seed(38)
    save_name(plot_tsne, "tSNE_group_time",
      color = time_col(object),
      shape = group_col(object),
      perplexity = perplexity
    )
  }
  # Multiple time points per subject
  if (!is.na(time_col(object)) &&
    !is.na(subject_col(object)) &&
    sum(object$QC == "QC") == 0) {
    set.seed(38)
    save_name(plot_pca_arrows, "PCA_arrows")

    set.seed(38)
    save_name(plot_tsne_arrows, "tSNE_arrows", perplexity = perplexity)
  }

  # Save plots
  if (merge) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("gridExtra is required for merging plots.")
    }
    log_text("Merging plots to a single PDF file")
    ggsave(paste0(prefix, ".pdf"),
      plot = gridExtra::marrangeGrob(grobs = plots, nrow = 1, ncol = 1, top = NULL),
      width = 12, height = 12
    )
  }
  if (!remove_singles) {
    log_text("Saving plots to separate files")
    for (name in names(plots)) {
      width <- 7
      height <- 7
      if (name %in% c("boxplots_injection", "dendrogram", "heatmap_samples", "boxplots_group", "dendrogram_time")) {
        width <- 15
      }
      if (name == "density_plot") height <- 6
      if (name == "heatmap_samples") height <- 16
      ggsave(paste0(prefix, "_", name, ".", format),
        plot = plots[[name]],
        width = width, height = height
      )
    }
  }
}
