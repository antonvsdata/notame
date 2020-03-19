
#' Suprahex plots
#'
#' Plots supraHex plots of each sample using functions from the supraHex package.
#' See the supraHex paper and package vignette for more information.
#'
#' @param object a MetaboSet object
#' @param file filename (pdf)
#' @param width,height dimensions of the plot
#' @param all_features if FALSE, flagged features are droppped
#' @param sample_labels the column for labels of samples in the plot
#' @param grid_xdim,grid_ydim dimensions of the grid for the samples
#' @param title.xy position of the sample label relative to the supraHex
#' @param title.rotate rotation of the sample label in degrees
#' @param fontsize the fontsize for sample labels
#' @param colormap colormap for the hexagons
#' @param ... other parameters for supraHex::sPipeline
#'
#' @examples
#' \dontrun{
#' plot_sample_suprahex(merged_sample[, 1:20], file = "supra.pdf", xdim = 5, title.xy = c(0.35, 1),
#'                      width = 10, height = 10, grid_xdim = 7, grid_ydim = 7, sample_labels = "Group")
#' }
#' @seealso \code{\link[supraHex]{sPipeline}}, \code{\link[supraHex]{sCompReorder}}, \code{\link[supraHex]{visCompReorder}}
#'
#' @export
plot_sample_suprahex <- function(object, file, width = 16, height = 16, all_features = FALSE,
                                 sample_labels = "Sample_ID",
                                 grid_xdim = NULL, grid_ydim = NULL,
                                 title.xy = c(0.35,1), title.rotate = 0, fontsize = 10,
                                 colormap = "jet", ...) {
  if (!requireNamespace("supraHex", quietly = TRUE)) {
    stop("Package \"supraHex\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  object <- drop_flagged(object, all_features = all_features)
  data <- scale(exprs(object))
  colnames(data) <- pData(object)[,sample_labels]

  sMap <- supraHex::sPipeline(data = data, ...)
  sReorder <- supraHex::sCompReorder(sMap=sMap, xdim = grid_xdim, ydim = grid_ydim)

  pdf(file = file, width = width, height = height)
  supraHex::visCompReorder(sMap=sMap, sReorder=sReorder, newpage = FALSE, colormap = colormap,
                           title.xy = title.xy, title.rotate = title.rotate, height = height,
                           gp = grid::gpar(fontsize = fontsize))
  dev.off()
}
