

#' Save line plots with mean
#'
#' Plots the change in the feature abundaces as a function of e.g. time.
#' A line is drawn for each subject and a mean line is added.
#' A separate plot is drawn for each feature.
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

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- ggplot(data, aes_string(x = x, y = fname)) +
      geom_line(aes_string(group = id), color = "grey20", alpha = 0.35, size = 0.3) +
      stat_summary(aes(group = 1), fun.data = "mean_se",
                   geom = "line", size = 1.2, color = "red") +
      labs(title = fname, y = "Abundance") +
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

#' Save box plots of each feature by group
#'
#' Draws a boxplot of feature abundances in each group.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param group character, name of the column to be used as x-axis and color
#'
#' @export
save_group_boxplots <- function(object, file, width = 8, height = 6, group = group_col(object),
                                color_scale =  NULL) {

  color_scale <- color_scale %||% getOption("amp.color_scale_dis")

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- ggplot(data, aes_string(x = group, y = fname, color = group)) +
      geom_boxplot() +
      color_scale +
      labs(title = fname, y = "Abundance") +
      theme_bw()

    plot(p)

  }
  dev.off()

  log_text(paste("Saved group boxplots to:", file))

}

#' Save line plots with errorbars by group
#'
#' Plots the change in the feature abundaces as a function of e.g. time.
#' A line is drawn for each group and errorbars are added.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param x character, name of the column to be used as x-axis
#' @param group character, name of the column containing group information, used for coloring
#' @param fun.data passed to ggplot2::stat_summary and used for errorbars,
#' "A function that is given the complete data and should return a data frame with variables ymin, y, and ymax."
#' @param fun.ymin,fun.y,fun.ymax Alternative to fun.data, passed to ggplot2::stat_summary,
#' "supply three individual functions that are each passed a vector of x's and should return a single number"
#' @param position_dodge_amout numeric: how much the group mean points should dodge away from each other
#' @param color_scale the color scale as returned by a ggplot function
#'
#' @seealso \code{\link[ggplot2]{stat_summary}}
#'
#' @export
save_group_lineplots <- function(object, file, width = 8, height = 6,
                                 x = time_col(object), group = group_col(object),
                                 fun.data = "mean_cl_boot", fun.y = NULL,
                                 fun.ymin = NULL, fun.ymax = NULL, position_dodge_amount = 0.2,
                                 color_scale =  NULL) {

  color_scale <- color_scale %||% getOption("amp.color_scale_dis")

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- ggplot(data, aes_string(x = x, y = fname, group = group, color = group)) +
      # Errorbars with solid lines
      stat_summary(fun.data = fun.data,
                   geom = "errorbar", width = 0.5,
                   fun.y = fun.y, fun.ymin = fun.ymin, fun.ymax = fun.ymax,
                   position = position_dodge(position_dodge_amount)) +
      # Plot point to mean
      stat_summary(fun.data = fun.data,
                   geom = "point",
                   fun.y = fun.y, fun.ymin = fun.ymin, fun.ymax = fun.ymax,
                   position = position_dodge(position_dodge_amount), size = 4) +
      # Line from mean to mean between for example timepoints
      stat_summary(fun.data =fun.data,
                   geom = "line",
                   position = position_dodge(position_dodge_amount), size = 0.5,
                   fun.y = fun.y, fun.ymin = fun.ymin, fun.ymax = fun.ymax) +
      labs(title = fname, y = "Abundance") +
      color_scale +
      theme_bw()

    plot(p)

  }
  dev.off()

  log_text(paste("Saved line plots with mean line to:", file))
}


