

#' Save line plots with mean
#'
#' Plots the change in the feature abundaces as a function of e.g. time.
#' A line is drawn for each subject and a mean line is added.
#' One plot is drawn for each feature
#'
#' @param object a MetaboSet object
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param x character, name of the column to be used as x-axis
#' @param id character, name of the column containing subject IDs
#'
#' @export
save_subject_line_plots <- function(object, file, width = 8, height = 6,
                                    x = time_col(object), id = subject_col(object)) {

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (fname in Biobase::featureNames(object)) {

    p <- ggplot(data, aes_string(x = x, y = fname)) +
      geom_line(aes_string(group = id), color = "grey20", alpha = 0.35, size = 0.3) +
      stat_summary(aes(group = 1), fun.data = "mean_se",
                   geom = "line", size = 1.2, color = "red") +
      theme_bw()

    if (class(data[, x]) == "factor") {
      p <- p +
        scale_x_discrete(expand = c(0.05,0.05))
    }
    plot(p)

  }
  dev.off()

  log_text(paste("Saved line plots with mean line to:", file))
}
