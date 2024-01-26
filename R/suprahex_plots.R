#' Suprahex plots
#'
#' Plots supraHex plots of each sample using functions from the supraHex package.
#' See the supraHex paper and package vignette for more information.
#'
#' @param object a MetaboSet object
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
#' plot_sample_suprahex(merged_sample[, 1:20],
#'   xdim = 5, title.xy = c(0.35, 1),
#'   grid_xdim = 7, grid_ydim = 7, sample_labels = "Group"
#' )
#' }
#' @seealso \code{\link[supraHex]{sPipeline}},
#' \code{\link[supraHex]{sCompReorder}},
#'  \code{\link[supraHex]{visCompReorder}}
#'
#' @export
plot_sample_suprahex <- function(object, all_features = FALSE,
                                 sample_labels = "Sample_ID",
                                 grid_xdim = NULL, grid_ydim = NULL,
                                 title.xy = c(0.35, 1), title.rotate = 0, fontsize = 10, # nolint: object_name_linter.
                                 colormap = "jet", ...) {
  if (!requireNamespace("supraHex", quietly = TRUE)) {
    stop("Package \"supraHex\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  add_citation("supraHex package was used for suprahexagonal maps:", citation("supraHex"))

  object <- drop_flagged(object, all_features = all_features)
  data <- scale(exprs(object))
  colnames(data) <- pData(object)[, sample_labels]

  s_map <- supraHex::sPipeline(data = data, ...)
  s_reorder <- supraHex::sCompReorder(sMap = s_map, xdim = grid_xdim, ydim = grid_ydim)

  supraHex::visCompReorder(
    sMap = s_map, sReorder = s_reorder, newpage = FALSE, colormap = colormap,
    title.xy = title.xy, title.rotate = title.rotate, height = height,
    gp = grid::gpar(fontsize = fontsize)
  )
}
