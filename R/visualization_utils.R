
#' Save plot to PDF
#'
#' Saves the given plot to a PDF file
#'
#' @param p a ggplot object
#' @param file the file path
#' @param ... other arguments to pdf, like width and height
#'
#' @seealso \code{\link[grDevices]{pdf}}
#'
#' @export
save_plot <- function(p, file, ...) {

  pdf(file, ...)
  plot(p)
  dev.off()
  log_text(paste("Saved", file))
}

#' Write all relevant visualizations to pdf
#'
#' A wrapper around all the major visualization functions, used for visualizing data between
#' major steps of data preprocessing. Saves all visualizations as PDFs with a set prefix on filenames.
#'
#' @param object A MetaboSet object
#' @param prefix character, a file path prefix added to the file paths
#' @param perplexity perplexity for t-SNE plots
#' @param merge logical, whether the files should be merged to a single PDF, see Details
#'
#' @details If \code{merge} is \code{TRUE}, then a file containing all the visualizations
#' named \code{prefix.pdf} will be created. NOTE: on Windows this requires installation of pdftk
#' (\url{https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/}) and on Linux you need to have pdfunite installed.
#' Mac is not yet supported.
#'
#' @export
visualizations <- function(object, prefix, perplexity = 30, merge = FALSE) {

  # Helper function
  file_names <- ""
  save_name <- function(p, name, ...) {
    file_name <- paste0(prefix, "_", name, ".pdf")
    save_plot(p, file = file_name, ...)
    file_names <<- paste(file_names, file_name)
  }

  if (sum(object$QC == "QC")) {
    save_name(plot_dist_density(object), "density_plot", width = 8, height = 6)
    save_name(plot_injection_lm(object), "lm_p_histograms")
  }
  # Quality metrics
  save_name(plot_quality(object), "quality_metrics")

  # Boxplots
  save_name(plot_sample_boxplots(object), "boxplots_group", width = 15)
  save_name(plot_sample_boxplots(object, order_by = "Injection_order",
                                 fill_by = "QC"), "boxplots_injection", width = 15)
  # PCA
  set.seed(38)
  save_name(plot_pca(object), "PCA_group")
  set.seed(38)
  save_name(plot_pca(object, color = "Injection_order"), "PCA_injection")
  set.seed(38)
  save_name(plot_pca_hexbin(object), "PCA_hexbin")
  # t-SNE
  set.seed(38)
  save_name(plot_tsne(object, perplexity = perplexity), "tSNE_group")
  set.seed(38)
  save_name(plot_tsne(object, perplexity = perplexity, color = "Injection_order"), "tSNE_injection")
  set.seed(38)
  save_name(plot_tsne_hexbin(object, perplexity = perplexity), "tSNE_hexbin")
  # Clustering
  save_name(plot_dendrogram(object), "dendrogram", width = 15)
  save_name(plot_sample_heatmap(object), "heatmap_samples", width = 15, height = 16)


  if (merge) {
    prefix <- gsub("_$", "", prefix)
    merged_file <- paste0(prefix, ".pdf")
    os <- Sys.info()[["sysname"]]
    if (os == "Windows") {
      # Merge files
      system(paste("pdftk", file_names, "cat output", merged_file))
      log_text(paste("Merged plots to", merged_file))
    } else if (os == "Linux"){
      system(paste("pdfunite", file_names, merged_file))
      log_text(paste("Merged plots to", merged_file))
    } else {
      log_text("Unfortunately your operating system is not yet supported by the merging")
    }
  }
}
