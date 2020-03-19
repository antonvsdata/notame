
#' Save plot to PDF
#'
#' Saves the given plot to a PDF file. If an error occurs with the plot, an empty file is created.
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
  tryCatch({
    plot(p)
  }, error = function(e) {
    dev.off()
    stop(e$message, call. = FALSE)})
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
#' On MacOS, no external software is needed. Note that at least on Windows, prefix should be a path from the root,
#' so that the underlying system command will find the files.
#' The type of visualizations to be saved depends on the type of object.
#' Here is a comprehensive list of the visualizations:
#' \itemize{
#' \item Distribution of quality metrics and flags
#' \item Boxplots of each sample in injection order
#' \item PCA scores plot of samples colored by injection order
#' \item t-SNE plot of samples colored by injection order
#' \item If the object has over 60 samples, hexbin versions of the PCA and t-SNE plots above
#' \item Dendrogram of samples ordered by hierarchical clustering, sample labels colored by group if present
#' \item heat map of intersample distances, ordered by hierarchical clustering
#' \item If the object has QC samples: \itemize{
#' \item Density function of the intersample distances in both QCs and biological samples
#' \item Histograms of p-values from linear regression of features against injection order
#' in both QCs and biological samples}
#' \item If the object has a group column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by group
#' }
#' \item If the object has a time column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by time
#' \item Dendrogram of samples ordered by hierarchical clustering, sample labels colored by time point
#' }
#' \item If the object has a group column OR a time column: \itemize{
#' \item Boxplots of samples ordered and colored by group and/or time
#' }
#' \item If the object has a group column AND a time column: \itemize{
#' \item PCA and tSNE plots with points shaped by group and colored by time
#' }
#' \item If the object has a time column AND a subject column: \itemize{
#' \item PCA and tSNE plots with arrows connecting the samples of each subject in time point order
#' }
#' }
#'
#' @export
visualizations <- function(object, prefix, perplexity = 30, merge = FALSE) {

  # Record file names for merging
  file_names <- ""
  # Helper function for handling errors and keeping track of file names
  save_name <- function(fun, name, width = 7, height = 7, ...) {
    p <- NULL
    tryCatch({
      p <- fun(object, ...)},
      error = function(e) {
        cat(paste0("Error with plot named ", name, ":\n", e$message, "\n"))}
      )

    if (!is.null(p)) {
      file_name <- paste0(prefix, "_", name, ".pdf")
      save_plot(p, file = file_name, width = width, height = height)
      file_names <<- paste(file_names, file_name)
    }
  }

  if (sum(object$QC == "QC")) {
    save_name(fun = plot_dist_density, name = "density_plot", width = 8, height = 6)
    save_name(plot_injection_lm, "lm_p_histograms")
  }
  # Quality metrics
  save_name(plot_quality, "quality_metrics")

  # Plots with injection order
  save_name(plot_sample_boxplots, "boxplots_injection",
            order_by = "Injection_order", fill_by = "QC", width = 15)
  set.seed(38)
  save_name(plot_pca, "PCA_injection", color = "Injection_order")
  set.seed(38)
  save_name(plot_tsne, "tSNE_injection", perplexity = perplexity, color = "Injection_order")

  # Clustering
  save_name(plot_dendrogram, "dendrogram", width = 15)
  save_name(plot_sample_heatmap, "heatmap_samples", width = 15, height = 16)

  # For large sets, plot hexbin plots
  if (ncol(object) > 60) {
    set.seed(38)
    save_name(plot_pca_hexbin, "PCA_hexbin")
    set.seed(38)
    save_name(plot_tsne_hexbin, "tSNE_hexbin", perplexity = perplexity)
  }

  # If grouped
  if (!is.na(group_col(object))) {
    set.seed(38)
    save_name(plot_pca, "PCA_group")

    set.seed(38)
    save_name(plot_tsne, "tSNE_group", perplexity = perplexity)
  }
  # Time point
  if (!is.na(time_col(object))) {
    set.seed(38)
    save_name(plot_pca, "PCA_time", color = time_col(object))

    set.seed(38)
    save_name(plot_tsne, "tSNE_time", color = time_col(object), perplexity = perplexity)

    save_name(plot_dendrogram, "dendrogram_time", color = time_col(object),  width = 15)
  }
  # Time point OR group
  if (!is.na(group_col(object)) | !is.na(time_col(object))) {
    save_name(plot_sample_boxplots, "boxplots_group", width = 15)
  }
  # Time point AND group
  if (!is.na(group_col(object)) & !is.na(time_col(object))) {
    set.seed(38)
    save_name(plot_pca, "PCA_group_time", color = time_col(object), shape = group_col(object))

    set.seed(38)
    save_name(plot_tsne, "tSNE_group_time", color = time_col(object), shape = group_col(object), perplexity = perplexity)
  }
  # Multiple time points per subject
  if (!is.na(time_col(object)) & !is.na(subject_col(object)) & sum(object$QC == "QC") == 0) {
    set.seed(38)
    save_name(plot_pca_arrows, "PCA_arrows")

    set.seed(38)
    save_name(plot_tsne_arrows, "tSNE_arrows", perplexity = perplexity)
  }

  if (merge) {
    prefix <- gsub("_$", "", prefix)
    merged_file <- paste0(prefix, ".pdf")
    os <- Sys.info()[["sysname"]]
    output <- NULL
    if (os == "Windows") {
      # Merge files
      output <- system(paste("pdftk", file_names, "cat output", merged_file))
    } else if (os == "Linux"){
      output <- system(paste("pdfunite", file_names, merged_file))
    } else if (os == "Darwin") {
      output <- system(paste('"/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o',
                             merged_file, file_names))
    } else {
      log_text("Unfortunately your operating system is not yet supported by the merging")
      return()
    }
    if (length(output) & output != "0") {
      log_text(paste("Merging plots resulted in the following message:", paste0(output, collapse = " ")))
    } else {
      log_text(paste("Attempted merging plots to", merged_file))
    }
  }
}
