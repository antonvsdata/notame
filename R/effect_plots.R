

#' Save line plots with mean
#'
#' Plots the change in the feature abundances as a function of e.g. time.
#' A line is drawn for each subject and a mean line is added.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param x character, name of the column to be used as x-axis
#' @param id character, name of the column containing subject IDs
#' @param color character, the column name to color the lines by (optional)
#' @param color_scale the color scale as returned by a ggplot function
#' @param facet character, the column name to facet by (optional, usually same as color)
#'
#' @examples
#' \dontrun{save_subject_line_plots(drop_qcs(example_set), file = "subject_lines.pdf")}
#'
#' @export
save_subject_line_plots <- function(object, all_features = FALSE, file, width = 8, height = 6,
                                    x = time_col(object), id = subject_col(object),
                                    color = NA, color_scale = getOption("notame.color_scale_dis"), facet = NULL) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  if (is.na(x)) {
    stop("The time column is missing")
  }
  if (is.na(id)) {
    stop("The subject column is missing")
  }

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- ggplot(data, aes_string(x = x, y = fname)) +
      labs(title = fname, y = "Abundance") +
      theme_bw()

    if (is.na(color)) {
      p <- p +
        geom_line(aes_string(group = id), color = "grey20", alpha = 0.35, size = 0.3) +
        stat_summary(aes(group = 1), fun.data = "mean_se",
                     geom = "line", size = 1.2, color = "red")
    } else {
      p <- p +
        geom_line(aes_string(group = id, color = color), alpha = 0.35, size = 0.3) +
        stat_summary(aes_string(group = color, color = color), fun.data = "mean_se",
                     geom = "line", size = 1.2) +
        color_scale
    }

    if (!is.null(facet)) {
      p <- p + facet_wrap(facets = facet)
    }

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
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param x character, name of the column to be used as x-axis
#' @param color character, name of the column to be used for coloring
#' @param color_scale the color scale as returned by a ggplot function
#'
#' @examples
#' \dontrun{
#' # Default boxplots by group
#' save_group_boxplots(drop_qcs(merged_sample), file = "boxplots.pdf")
#' # x and color can be a different variable
#' save_group_boxplots(drop_qcs(merged_sample), file = "boxplots_time.pdf",
#'                     x = "Time", color = "Group")
#' }
#'
#' @export
save_group_boxplots <- function(object, all_features = FALSE, file, width = 8, height = 6,
                                x = group_col(object), color = group_col(object),
                                color_scale =  getOption("notame.color_scale_dis")) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- ggplot(data, aes_string(x = x, y = fname, color = color)) +
      geom_boxplot(position = position_dodge(0.6), width = 0.5) +
      stat_summary(fun.data = mean_se,
                   geom = "point", shape = 18, size = 3,
                   position = position_dodge(0.6)) +
      color_scale +
      labs(title = fname, y = "Abundance") +
      theme_bw()

    plot(p)

  }
  dev.off()

  log_text(paste("Saved group boxplots to:", file))

}

#' Save beeswarm plots of each feature by group
#'
#' Draws a beeswarm plot of feature abundances in each group.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param x character, name of the column to be used as x-axis
#' @param add_boxplots logical, should boxplots be added to the figure?
#' @param color character, name of the column to be used for coloring
#' @param color_scale the color scale as returned by a ggplot function
#'
#' @examples
#' \dontrun{
#' # Default beeswarms by group
#' save_beeswarm_plots(drop_qcs(merged_sample), file = "beeswarms.pdf")
#' # x and color can be a different variable
#' save_beeswarm_plots(drop_qcs(merged_sample), file = "beeswarms_time.pdf",
#'                     x = "Time", color = "Group")
#' }
#'
#' @export
save_beeswarm_plots <- function(object, all_features = FALSE, file, width = 8, height = 6,
                                x = group_col(object),add_boxplots = FALSE,
                                color = group_col(object), color_scale =  NULL){
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  color_scale <- color_scale %||% getOption("notame.color_scale_dis")

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- ggplot(data, aes_string(x = x, y = fname, color = color))

    if (add_boxplots) {
      p <- p +
        geom_boxplot(position = position_dodge(0.6), width = 0.5, lwd = .3) +
        stat_boxplot(geom ='errorbar', width = 0.5, lwd = .3)
    }
    p <- p +
      ggbeeswarm::geom_beeswarm() +
      color_scale +
      labs(title = fname, y = "Abundance") +
      theme_bw()


    plot(p)

  }
  dev.off()

  log_text(paste("Saved beeswarm plots to:", file))
}

#' Save line plots with errorbars by group
#'
#' Plots the change in the feature abundances as a function of e.g. time.
#' A line is drawn for each group and error bars are added.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default), flagged features are removed before visualization.
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
#' @examples
#' \dontrun{save_group_lineplots(drop_qcs(merged_sample), file = "group_lineplots.pdf")}
#'
#' @export
save_group_lineplots <- function(object, all_features = FALSE, file, width = 8, height = 6,
                                 x = time_col(object), group = group_col(object),
                                 fun.data = "mean_cl_boot", fun.y = NULL,
                                 fun.ymin = NULL, fun.ymax = NULL, position_dodge_amount = 0.2,
                                 color_scale =  getOption("notame.color_scale_dis")) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  if (is.na(group)) {
    stop("The group column is missing")
  }
  if (is.na(x)) {
    stop("The time column is missing")
  }
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


