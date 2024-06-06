#' Save plots of individual features
#'
#' Helper function for saving plots of individual features
#' to either one multi-page PDF or separate EMF figures
#'
#' @param plot_fun a function with arguments:
#' data frame from combined_data(object)
#' feature id
#' Should return a ggplot object for plotting
#' @param ... arguments to code{\link[ggplot2]{ggsave}}
save_feature_plots <- function(object, file_path, format,
                               title, subtitle, text_base_size,
                               plot_fun, ...) {
  if (is.null(file_path)) file_path <- getwd()
  if (endsWith(file_path, ".pdf") && format != "pdf") {
    message("Switching to PDF format based on file path")
    format <- "pdf"
  } else if (!endsWith(file_path, "/") && format != "pdf") {
    message("Adding an additional slash to file path to allow proper folder structure")
    file_path <- paste0(file_path, "/")
  }

  folder <- dirname(file_path)
  if (!file.exists(folder)) {
    message("Creating folder ", folder)
    dir.create(folder, recursive = TRUE)
  }
  if (format == "pdf") {
    pdf(file_path, ...)
  }

  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- featureNames(object)[i]
    name <- fData(object)[i, title]

    p <- plot_fun(object, fname)

    if (format != "pdf") {
      if (is.null(title)) {
        file <- paste0(file_path, fname, ".", format)
      } else {
        file <- paste0(file_path, gsub("[:;/]", "_", name), ".", format)
      }

      ggsave(file, plot = p, device = format, ...)
    } else {
      print(p)
    }
  }

  if (format == "pdf") {
    dev.off()
  }
}

#' Generate a list of plots
#'
#' Helper function for generating a list of feature-wise plots given a plot function
#'
#' @param object a MetaboSet object, should contain only features to be plotted
#' @param plot_fun function, a notame plot function
#' @return a list of ggplot objects
create_feature_plot_list <- function(object, plot_fun) {
  message("Just a remainder, creating a long list of plots takes a lot of memory!")
  plot_list <- vector("list", nrow(object))
  for (i in seq_len(nrow(object))) {
    if (i %% 500 == 0) {
      cat(paste0("Iteration ", i, "/", nrow(object), "\n"))
    }
    fname <- featureNames(object)[i]
    p <- plot_fun(object, fname)
    plot_list[[i]] <- p
  }

  plot_list
}

#' Save line plots with mean
#'
#' Plots the change in the feature abundances as a function of e.g. time.
#' A line is drawn for each subject and a mean line is added.
#' A separate plot is drawn and saved for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used?
#' If FALSE (the default), flagged features are removed before visualization.
#' @param save logical, if false, the plots are not saved but returned as a list
#' @param file_path character, a file path for PDF or prefix added to the file paths for other formats
#' @param format character, format in which the plots should be saved
#' @param x character, name of the column to be used as x-axis
#' @param id character, name of the column containing subject IDs
#' @param title,subtitle column names from fData to use as plot title/filename and subtitle.
#' Set to NULL for no title/subtitle, this creates running numbered filenames
#' @param color character, the column name to color the lines by (optional)
#' @param color_scale the color scale as returned by a ggplot function
#' @param facet character, the column name to facet by (optional, usually same as color)
#' @param text_base_size integer, base size for text in figures
#' @param line_width numeric, width of the lines
#' @param mean_line_width numeric, width of the mean line
#' @param title_line_length integer, maximum length of the title line in characters, passed to stringr::str_wrap
#' @param theme a ggplot theme to be added to the plot
#' @param ... arguments to code{\link[ggplot2]{ggsave}}
#'
#'
#' @examples
#' \dontrun{
#' save_subject_line_plots(drop_qcs(example_set),
#'   file_path = "./subject_line_plots.pdf",
#'   format = "pdf"
#' )
#' }
#' # Plot one feature
#' save_subject_line_plots(drop_qcs(example_set[1, ]), save = FALSE)
#' @export
save_subject_line_plots <- function(object,
                                    all_features = FALSE,
                                    save = TRUE,
                                    file_path = NULL,
                                    format = "emf",
                                    x = time_col(object),
                                    id = subject_col(object),
                                    title = "Feature_ID",
                                    subtitle = NULL,
                                    color = NA,
                                    color_scale = getOption("notame.color_scale_dis"),
                                    facet = NULL,
                                    text_base_size = 14,
                                    line_width = 0.3,
                                    mean_line_width = 1.2,
                                    title_line_length = 40,
                                    theme = theme_bw(base_size = text_base_size),
                                    ...) {
  if (is.na(x)) {
    stop("The time column is missing")
  }
  if (is.na(id)) {
    stop("The subject column is missing")
  }

  subject_line_fun <- function(object, fname) {
    data <- combined_data(object)

    p <- ggplot(data, aes(x = .data[[x]], y = .data[[fname]]))

    if (is.na(color)) {
      p <- p +
        geom_line(aes(group = .data[[id]]),
          color = "grey20",
          alpha = 0.35,
          size = line_width
        ) +
        stat_summary(aes(group = 1),
          fun.data = "mean_se",
          geom = "line",
          size = mean_line_width,
          color = color_scale$palette(1)[1]
        )
    } else {
      p <- p +
        geom_line(aes(group = .data[[id]], color = .data[[color]]),
          alpha = 0.35,
          size = line_width
        ) +
        stat_summary(aes(group = .data[[color]], color = .data[[color]]),
          fun.data = "mean_se",
          geom = "line",
          size = mean_line_width
        ) +
        color_scale
    }
    if (!is.null(facet)) {
      p <- p + facet_wrap(facets = facet)
    }
    if (class(data[, x]) == "factor") {
      p <- p +
        scale_x_discrete(expand = c(0.05, 0.05))
    }
    splitted_title <-
      p <- p +
      theme +
      labs(
        title = stringr::str_wrap(fData(object)[fname, title], title_line_length),
        subtitle = fData(object)[fname, subtitle],
        y = "Abundance"
      )
    p
  }

  object <- drop_flagged(object, all_features)
  if (save) {
    save_feature_plots(
      object, file_path, format,
      title, subtitle, text_base_size, subject_line_fun, ...
    )
    log_text(paste("Saved line plots with mean line to:", file_path))
  } else {
    return(create_feature_plot_list(object, subject_line_fun))
    log_text("Created a list of line plots with mean line")
  }
}

#' Save box plots of each feature by group
#'
#' Draws a boxplot of feature abundances in each group.
#' A separate plot is drawn and saved for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#' @param save logical, if false, the plots are not saved but returned as a list
#' @param file_path character, a file path for PDF or prefix added to the file paths for other formats
#' @param format character, format in which the plots should be saved
#' @param x character, name of the column to be used as x-axis
#' @param color character, name of the column to be used for coloring
#' @param title,subtitle column names from fData to use as plot title/filename and subtitle.
#' Set to NULL for no title/subtitle, this creates running numbered filenames
#' @param color_scale the color scale as returned by a ggplot function
#' @param text_base_size integer, base size for text in figures
#' @param box_width numeric, width of the boxes
#' @param line_width numeric, width of the lines
#' @param point_size numeric, size of the mean points
#' @param title_line_length integer, maximum length of the title line in characters, passed to stringr::str_wrap
#' @param theme a ggplot theme to be added to the plot
#' @param ... arguments to code{\link[ggplot2]{ggsave}}
#'
#' @examples
#' \dontrun{
#' # Default boxplots by group
#' save_group_boxplots(drop_qcs(merged_sample),
#'   file_path = "./group_boxplots.pdf",
#'   format = "pdf", title = NULL
#' )
#' # x and color can be a different variable
#' save_group_boxplots(drop_qcs(merged_sample)[1:10],
#'   file_path = "./time_boxplots/",
#'   format = "emf",
#'   x = "Time",
#'   color = "Group", title = NULL
#' )
#' }
#' # Plot one feature
#' save_group_boxplots(drop_qcs(merged_sample)[5, ], save = FALSE)
#' @export
save_group_boxplots <- function(object,
                                all_features = FALSE,
                                save = TRUE,
                                file_path = NULL,
                                format = "emf",
                                x = group_col(object),
                                color = group_col(object),
                                title = "Feature_ID",
                                subtitle = NULL,
                                color_scale = getOption("notame.color_scale_dis"),
                                text_base_size = 14,
                                box_width = 0.8,
                                line_width = 0.5,
                                point_size = 3,
                                title_line_length = 40,
                                theme = theme_bw(base_size = text_base_size),
                                ...) {
  boxplot_fun <- function(object, fname) {
    data <- combined_data(object)
    dodge_amount <- box_width + 0.05
    p <- ggplot(data, aes(x = .data[[x]], y = .data[[fname]], color = .data[[color]])) +
      geom_boxplot(position = position_dodge(dodge_amount), width = box_width, size = line_width) +
      stat_summary(
        fun.data = mean_se,
        geom = "point",
        shape = 18,
        size = point_size,
        position = position_dodge(dodge_amount)
      ) +
      color_scale +
      theme +
      labs(
        title = stringr::str_wrap(fData(object)[fname, title], title_line_length),
        subtitle = fData(object)[fname, subtitle],
        y = "Abundance"
      )
    if (x == color) {
      p <- p + guides(color = "none")
    }
    p
  }

  object <- drop_flagged(object, all_features)
  if (save) {
    save_feature_plots(
      object, file_path, format,
      title, subtitle, text_base_size, boxplot_fun, ...
    )
    log_text(paste("Saved group boxplots to:", file_path))
  } else {
    return(create_feature_plot_list(object, boxplot_fun))
    log_text("Created a list of group boxplots")
  }
}



#' Save beeswarm plots of each feature by group
#'
#' Draws a beeswarm plot of feature abundances in each group.
#' A separate plot is drawn and saved for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#' @param file_path character, a file path for PDF or prefix added to the file paths for other formats
#' @param format character, format in which the plots should be saved
#' @param x character, name of the column to be used as x-axis
#' @param add_boxplots logical, should boxplots be added to the figure?
#' @param title,subtitle column names from fData to use as plot title/filename and subtitle.
#' Set to NULL for no title/subtitle, this creates running numbered filenames
#' @param color character, name of the column to be used for coloring
#' @param color_scale the color scale as returned by a ggplot function
#' @param text_base_size integer, base size for text in figures
#' @param cex numeric, scaling for adjusting point spacing
#' @param size numeric, size of points
#' @param title_line_length integer, maximum length of the title line in characters, passed to stringr::str_wrap
#' @param theme a ggplot theme to be added to the plot
#' @param ... arguments to code{\link[ggplot2]{ggsave}}
#'
#' @examples
#' \dontrun{
#' # Default beeswarms by group
#' save_beeswarm_plots(drop_qcs(merged_sample),
#'   file_path = "./beeswarm_plots.pdf",
#'   format = "pdf"
#' )
#' # x and color can be a different variable
#' save_beeswarm_plots(drop_qcs(merged_sample)[1:10],
#'   file_path = "./beeswarm_plots/",
#'   format = "png",
#'   x = "Time",
#'   color = "Group"
#' )
#' }
#' # Plot one feature
#' save_beeswarm_plots(drop_qcs(merged_sample)[5, ], save = FALSE)
#' @export
save_beeswarm_plots <- function(object,
                                all_features = FALSE,
                                save = TRUE,
                                file_path = NULL,
                                format = "emf",
                                x = group_col(object),
                                add_boxplots = FALSE,
                                title = "Feature_ID",
                                subtitle = NULL,
                                color = group_col(object),
                                color_scale = getOption("notame.color_scale_dis"),
                                text_base_size = 14,
                                cex = 2,
                                size = 2,
                                title_line_length = 40,
                                theme = theme_bw(base_size = text_base_size),
                                ...) {
  beeswarm_fun <- function(object, fname) {
    data <- combined_data(object)
    p <- ggplot(data, aes(x = .data[[x]], y = .data[[fname]], color = .data[[color]]))

    if (add_boxplots) {
      p <- p +
        geom_boxplot(position = position_dodge(0.6), width = 0.5, lwd = .3) +
        stat_boxplot(geom = "errorbar", width = 0.5, lwd = .3)
    }
    p <- p +
      ggbeeswarm::geom_beeswarm(cex = cex, size = size) +
      color_scale +
      theme +
      labs(
        title = stringr::str_wrap(fData(object)[fname, title], title_line_length),
        subtitle = fData(object)[fname, subtitle],
        y = "Abundance"
      )
    if (x == color) {
      p <- p + guides(color = "none")
    }
    p
  }

  object <- drop_flagged(object, all_features)
  if (save) {
    save_feature_plots(
      object, file_path, format,
      title, subtitle, text_base_size, beeswarm_fun, ...
    )

    log_text(paste("Saved beeswarm plots to:", file_path))
  } else {
    return(create_feature_plot_list(object, beeswarm_fun))
    log_text("Created a list of beeswarm plots")
  }
}

#' Save scatter plots of each feature against a set variable
#'
#' Draws a scatterplots with a feature on y-axis and another variable on x-axis.
#' A separate plot is drawn and saved for each feature.
#'
#' @param object a MetaboSet object
#' @param x character, name of the column to be used as x-axis
#' @param file_path character, a file path for PDF or prefix added to the file paths for other formats
#' @param format character, format in which the plots should be saved
#' @param all_features logical, should all features be used? If FALSE
#' (the default), flagged features are removed before visualization.
#' @param color character, name of the column to be used for coloring
#' @param color_scale the color scale as returned by a ggplot function.
#' Set to NA to choose the appropriate scale based on the class of the coloring variable.
#' @param shape character, name of the column used for shape
#' @param title,subtitle column names from fData to use as plot title/filename and subtitle.
#' Set to NULL for no title/subtitle, this creates running numbered filenames
#' @param shape_scale the shape scale as returned by a ggplot function
#' @param text_base_size integer, base size for text in figures
#' @param point_size numeric, size of the points
#' @param title_line_length integer, maximum length of the title line in characters, passed to stringr::str_wrap
#' @param theme a ggplot theme to be added to the plot
#' @param ... arguments to code{\link[ggplot2]{ggsave}}
#'
#' @examples
#' \dontrun{
#' # Against injection order, colored by group
#' save_scatter_plots(
#'   object = merged_sample[1:10],
#'   x = "Injection_order",
#'   color = "Group",
#'   file_path = "./scatter_plots.pdf",
#'   format = "pdf"
#' )
#' }
#' # Plot one feature
#' save_scatter_plots(merged_sample[5, ], save = FALSE, color = "Group")
#' @export
save_scatter_plots <- function(object,
                               x = "Injection_order",
                               save = TRUE,
                               file_path = NULL,
                               format = "emf",
                               all_features = FALSE,
                               color = NULL,
                               color_scale = NA,
                               shape = NULL,
                               title = "Feature_ID",
                               subtitle = NULL,
                               shape_scale = getOption("notame.shape_scale"),
                               text_base_size = 14,
                               point_size = 2,
                               title_line_length = 40,
                               theme = theme_bw(base_size = text_base_size),
                               ...) {
  scatter_fun <- function(object, fname) {
    data <- combined_data(object)
    p <- scatter_plot(
      data = data,
      x = x,
      y = fname,
      color = color,
      color_scale = color_scale,
      shape = shape,
      shape_scale = shape_scale,
      point_size = point_size,
      fixed = FALSE,
      apply_theme_bw = FALSE
    ) +
      theme +
      labs(
        title = stringr::str_wrap(fData(object)[fname, title], title_line_length),
        subtitle = fData(object)[fname, subtitle],
        y = "Abundance"
      )
    p
  }
  object <- drop_flagged(object, all_features)
  if (save) {
    save_feature_plots(
      object, file_path, format,
      title, subtitle, text_base_size, scatter_fun, ...
    )
    log_text(paste("Saved scatter plots to:", file_path))
  } else {
    return(create_feature_plot_list(object, scatter_fun))
    log_text("Created a list of scatter plots")
  }
}


#' Save line plots with errorbars by group
#'
#' Plots the change in the feature abundances as a function of e.g. time.
#' A line is drawn for each group and error bars are added.
#' A separate plot is drawn for each feature.
#'
#' @param object a MetaboSet object
#' @param all_features logical, should all features be used? If FALSE (the default),
#' flagged features are removed before visualization.
#' @param file_path character, a file path for PDF or prefix added to the file paths for other formats
#' @param format character, format in which the plots should be saved
#' @param x character, name of the column to be used as x-axis
#' @param group character, name of the column containing group information, used for coloring
#' @param title,subtitle column names from fData to use as plot title/filename and subtitle.
#' Set to NULL for no title/subtitle, this creates running numbered filenames
#' @param fun.data passed to ggplot2::stat_summary and used for errorbars,
#' "A function that is given the complete data and should return a data frame with variables ymin, y, and ymax."
#' @param fun.min,fun,fun.max Alternative to fun.data, passed to ggplot2::stat_summary,
#' "supply three individual functions that are each passed a vector of x's and should return a single number"
#' @param position_dodge_amount numeric: how much the group mean points should dodge away from each other
#' @param color_scale the color scale as returned by a ggplot function
#' @param text_base_size integer, base size for text in figures
#' @param line_width numeric, width of the lines
#' @param point_size numeric, size of the points
#' @param title_line_length integer, maximum length of the title line in characters, passed to stringr::str_wrap
#' @param theme a ggplot theme to be added to the plot
#' @param ... arguments to code{\link[ggplot2]{ggsave}}
#'
#' @examples
#' \dontrun{
#' save_group_lineplots(drop_qcs(merged_sample),
#'   file_path = "./group_line_plots.pdf",
#'   format = "pdf"
#' )
#' save_group_lineplots(drop_qcs(merged_sample)[1:10],
#'   file_path = "./group_line_plots/",
#'   format = "png"
#' )
#' }
#' # Plot one feature
#' save_group_lineplots(drop_qcs(merged_sample[5, ]), save = FALSE)
#' @export
save_group_lineplots <- function(object,
                                 all_features = FALSE,
                                 save = TRUE,
                                 file_path = NULL,
                                 format = "emf",
                                 x = time_col(object),
                                 group = group_col(object),
                                 title = "Feature_ID",
                                 subtitle = NULL,
                                 fun.data = "mean_cl_boot", # nolint: object_name_linter.
                                 fun = NULL,
                                 fun.min = NULL, # nolint: object_name_linter.
                                 fun.max = NULL, # nolint: object_name_linter.
                                 position_dodge_amount = 0.2,
                                 color_scale = getOption("notame.color_scale_dis"),
                                 text_base_size = 14,
                                 line_width = 0.5,
                                 point_size = 4,
                                 title_line_length = 40,
                                 theme = theme_bw(base_size = text_base_size),
                                 ...) {
  if (is.na(group)) {
    stop("The group column is missing")
  }
  if (is.na(x)) {
    stop("The time column is missing")
  }

  line_fun <- function(object, fname) {
    data <- combined_data(object)
    p <- ggplot(
      data,
      aes(x = .data[[x]], y = .data[[fname]], group = .data[[group]], color = .data[[group]])
    ) +
      # Errorbars with solid lines
      stat_summary(
        fun.data = fun.data,
        geom = "errorbar", width = line_width,
        fun = fun,
        fun.min = fun.min,
        fun.max = fun.max,
        position = position_dodge(position_dodge_amount)
      ) +
      # Plot point to mean
      stat_summary(
        fun.data = fun.data,
        geom = "point",
        fun = fun,
        fun.min = fun.min,
        fun.max = fun.max,
        position = position_dodge(position_dodge_amount),
        size = point_size
      ) +
      # Line from mean to mean between for example timepoints
      stat_summary(
        fun.data = fun.data,
        geom = "line",
        position = position_dodge(position_dodge_amount), size = line_width,
        fun = fun,
        fun.min = fun.min,
        fun.max = fun.max
      ) +
      color_scale +
      theme +
      labs(
        title = stringr::str_wrap(fData(object)[fname, title], title_line_length),
        subtitle = fData(object)[fname, subtitle],
        y = "Abundance"
      )
    if (x == group) {
      p <- p + guides(color = "none")
    }
    p
  }

  object <- drop_flagged(object, all_features)
  if (save) {
    save_feature_plots(
      object, file_path, format,
      title, subtitle, text_base_size, line_fun, ...
    )

    log_text(paste("Saved line plots with mean line to:", file_path))
  } else {
    return(create_feature_plot_list(object, line_fun))
    log_text("Created a list of line plots with mean line")
  }
}
