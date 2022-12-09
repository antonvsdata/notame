save_plot_helper <- function(plot, file, width, height, res = 300) {
  # don't ask how this works but it can handle multiple file extensions and files with no extension
  file_format = gsub("^(.*\\.|[^.]+)(?=[^.]*)", "", file, perl = TRUE)
  switch (file_format,
          "pdf" = {
            pdf(file = file,
                width = width,
                height = height)
          },
          "emf" = {
            if (!requireNamespace("devEMF", quietly = TRUE)) {
              stop("Package devEMF needed for this function to work. Please install it.",
                   call. = FALSE)
            }
            devEMF::emf(file = file,
                        width = width,
                        height = height)
          },
          "svg" = {
            svg(file = file,
                width = width,
                height = height)
          },
          "png" = {
            png(filename = file,
                width = width,
                height = height,
                units = "in",
                res = res)
          },
          "tiff" = {
            tiff(filename = file,
                 width = width,
                 height = height,
                 units = "in",
                 res = res,
                 compression = "lzw")
          },
          warning(paste0("File format '", file_format,
                         "' is not valid, saving failed"
          ))
  )
  print(plot)
  dev.off()
}



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
#' @param title,subtitle column names from fData to use as plot title and subtitle. Set to NULL for no title/subtitle
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
                                    title = "Feature_ID", subtitle = NULL,
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
      labs(title = fData(object)[i, title], subtitle = fData(object)[i, subtitle], y = "Abundance") +
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
#' @param title,subtitle column names from fData to use as plot title and subtitle. Set to NULL for no title/subtitle
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
                                title = "Feature_ID", subtitle = NULL,
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
      labs(title = fData(object)[i, title], subtitle = fData(object)[i, subtitle], y = "Abundance") +
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
#' @param title,subtitle column names from fData to use as plot title and subtitle. Set to NULL for no title/subtitle
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
                                title = "Feature_ID", subtitle = NULL,
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
      labs(title = fData(object)[i, title], subtitle = fData(object)[i, subtitle], y = "Abundance") +
      theme_bw()


    plot(p)

  }
  dev.off()

  log_text(paste("Saved beeswarm plots to:", file))
}

#' Save scatter plots of each feature against a set variable
#'
#' Draws a scatterplots with a feature on y-axis and another variable on x-axis.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param x character, name of the column to be used as x-axis
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param all_features logical, should all features be used? If FALSE
#' (the default), flagged features are removed before visualization.
#' @param color character, name of the column to be used for coloring
#' @param color_scale the color scale as returned by a ggplot function.
#' Set to NA to choose the appropriate scale based on the class of the coloring variable.
#' @param shape character, name of the column used for shape
#' @param title,subtitle column names from fData to use as plot title and subtitle. Set to NULL for no title/subtitle
#' @param shape_scale the shape scale as returned by a ggplot function
#'
#' @examples
#' \dontrun{
#' # Against injection order, colored by group
#' save_scatter_plots(object = merged_sample,
#'                            x = "Injection_order",
#'                            color = "Group",
#'                            file = "injection_scatter.pdf")
#' }
#'
#' @export
save_scatter_plots <- function(object, x = group_col(object), file,
                               width = 8, height = 6,
                               all_features = FALSE,
                               color = NULL, color_scale =  NA,
                               shape = NULL,
                               title = "Feature_ID", subtitle = NULL,
                               shape_scale = getOption("notame.shape_scale")) {

  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)

  pdf(file, width = width, height = height)

  data <- combined_data(object)

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- Biobase::featureNames(object)[i]

    p <- scatter_plot(data = data, x = x, y = fname,
                      color = color, color_scale = color_scale,
                      shape = shape, shape_scale = shape_scale,
                      title = fData(object)[i, title], subtitle = fData(object)[i, subtitle], ylab = "Abundance")

    plot(p)

  }
  dev.off()

  log_text(paste("Saved scatter plots to:", file))

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
#' @param title,subtitle column names from fData to use as plot title and subtitle. Set to NULL for no title/subtitle
#' @param fun.data passed to ggplot2::stat_summary and used for errorbars,
#' "A function that is given the complete data and should return a data frame with variables ymin, y, and ymax."
#' @param fun.ymin,fun.y,fun.ymax Alternative to fun.data, passed to ggplot2::stat_summary,
#' "supply three individual functions that are each passed a vector of x's and should return a single number"
#' @param position_dodge_amount numeric: how much the group mean points should dodge away from each other
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
                                 title = "Feature_ID", subtitle = NULL,
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
      labs(title = fData(object)[i, title], subtitle = fData(object)[i, subtitle], y = "Abundance") +
      color_scale +
      theme_bw()

    plot(p)

  }
  dev.off()

  log_text(paste("Saved line plots with mean line to:", file))
}


