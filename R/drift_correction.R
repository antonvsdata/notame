# Used to combine and return multiple objects from a foreach loop
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

#' Fit a cubic spline to correct drift
#'
#' Corrects the drift in the features by applying smoothed cubic spline regression
#' to each feature separately.
#'
#' @param object a MetaboSet object
#' @param log_transform logical, should drift correction be done on log-transformed values? See Details
#' @param spar smoothing parameter
#' @param spar_lower,spar_upper lower and upper limits for the smoothing parameter
#'
#' @return list with object = MetaboSet object as the one supplied, with drift corrected features
#' and predicted = matrix of the predicted values by the cubic spline (used in visualization)
#'
#' @details If \code{log_transform = TRUE}, the correction will be done on log-transformed values.
#' The correction formula depends on whether the correction is run on original values or log-transformed values.
#' In log-space: \eqn{corrected = original + mean of QCs - prediction by cubic spline}.
#' In original space: \eqn{corrected = original * prediction for first QC / prediction for current point}.
#' We recommend doing the correction in the log-space since the log-transfomred data better follows the
#' assumptions of cubic spline regression. The drift correction in the original space also sometimes results
#' in negative values, and results in rejection of the drift corrrection procedure.
#'
#' If \code{spar} is set to \code{NULL} (the default), the smoothing parameter will
#' be separately chosen for each feature from the range [\code{spar_lower, spar_upper}]
#' using cross validation.
#'
#' @examples
#' dc <- dc_cubic_spline(merged_sample)
#' corrected <- dc$object
#'
#' @seealso  \code{\link[stats]{smooth.spline}} for details about the regression,
#' \code{\link{inspect_dc}} for analysing the drift correction results,
#' \code{\link{save_dc_plots}} for plotting the drift correction process for each feature
#'
#' @importFrom Biobase exprs exprs<-
#'
#' @export
dc_cubic_spline <- function(object, log_transform = TRUE, spar = NULL, spar_lower = 0.5, spar_upper = 1.5) {

  # Start log
  log_text(paste("\nStarting drift correction at", Sys.time()))

  # Zero values do not behave correctly
  if (sum(exprs(object) == 0, na.rm = TRUE)) {
    log_text("Zero values in feature abundances detected. Zeroes will be replaced with 1.1")
    exprs(object)[exprs(object) == 0] <- 1.1
  }
  # Extract data and injection order for QC samples and the full dataset
  qc <- object[, object$QC == "QC"]
  qc_order <- qc$Injection_order
  qc_data <- exprs(qc)

  full_order <- object$Injection_order
  full_data <- exprs(object)

  # log-transform before fiting the cubic spline
  if (log_transform) {
    if (sum(exprs(object) == 1, na.rm = TRUE)) {
      log_text("Values of 1 in feature abundances detected. 1s will be replaced with 1.1.")
      exprs(object)[exprs(object) == 1] <- 1.1
    }
    qc_data <- log(qc_data)
    full_data <- log(full_data)
  }
  # comb needs a matrix with the right amount of columns
  n <- ncol(full_data)
  # Return both predicted values (for plotting) and drift corrected values for each feature
  dc_data <- foreach::foreach(i = seq_len(nrow(object)), .combine = comb) %dopar% {

    dnames <- list(rownames(full_data)[i], colnames(full_data))
    # Spline cannot be fitted if there are les than 4 QC values
    qc_detected <- !is.na(qc_data[i, ])
    if (sum(qc_detected) < 4) {
      return(list(corrected = matrix(NA_real_, nrow = 1, ncol = n, dimnames = dnames),
                  predicted = matrix(NA_real_, nrow = 1, ncol = n, dimnames = dnames)))
    }

    # Spline regression
    fit <- smooth.spline(x = qc_order[qc_detected], y = qc_data[i, qc_detected], all.knots = TRUE,
                         spar = spar, control.spar = list("low" = spar_lower, "high" = spar_upper))
    predicted <- predict(fit, full_order)$y
    # Correction
    # Substraction in log space, division in original space
    if (log_transform) {
      corrected <- full_data[i, ] + mean(qc_data[i, qc_detected]) - predicted
    } else {
      corr_factors <- predicted[1]/predicted
      corrected <- full_data[i, ] * corr_factors
    }

    # Each iteration of the loop returns one row to corrected and one row to predicted
    list(corrected = matrix(corrected, ncol = n, dimnames = dnames),
         predicted = matrix(predicted, ncol = n, dimnames = dnames))

  }
  corrected <- dc_data$corrected
  # Inverse the initial log transformation
  if (log_transform) {
    corrected <- exp(corrected)
  }
  exprs(object) <- corrected

  # Recompute quality metrics
  object <- assess_quality(object)

  log_text(paste("Drift correction performed at", Sys.time()))
  return(list(object = object, predicted = dc_data$predicted))
}


#' Flag the results of drift correction
#'
#' Determines whether the drift correction worked.
#' The primary reason is to search for features where there were too many missing values in the QCs,
#' so it was not possible to run drift correction. If the drift correction is run on the original
#' values (not log-transformed), then there is also a need to check that the correction did not result
#' in any negative values. This can sometimes happen if the prediction curve takes an extreme shape.
#' If quality is monitored,
#' a quality condition is checked for each feature. If the condition is fulfilled, the drift corrected feature is retained,
#' otherwise the original feature is retained and the drift corrected feature is discarded.
#' The result of this operation is recorded in the feature data.
#'
#' @param orig a MetaboSet object, before drift correction
#' @param dc a MetaboSet object, after drift correction
#' @param check_quality logical, whether quality should be monitored.
#' @param condition a character specifying the condition, see Details
#'
#' @return MeatboSet object
#'
#' @details The \code{condition} parameter should be a character giving a condition compatible
#' with dplyr::filter. The condition is applied on the \strong{changes} in the quality metrics
#' RSD, RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r < 0 and D_ratio_r < 0",
#' meaning that both RSD_r and D_ratio_r need to decrease in the drift correction, otherwise the
#' drift corrected feature is discarded and the original is retained.
#'
#' @seealso \code{\link{correct_drift}}, \code{\link{save_dc_plots}}
#'
#' @examples
#' dc <- dc_cubic_spline(merged_sample)
#' corrected <- dc$object
#' inspected <- inspect_dc(orig = merged_sample, dc = corrected,
#'                         check_quality = TRUE)
#'
#' @export
inspect_dc <- function(orig, dc, check_quality, condition = "RSD_r < 0 & D_ratio_r < 0") {

  if (is.null(quality(orig))) {
    orig <- assess_quality(orig)
  }
  if (is.null(quality(dc))) {
    dc <- assess_quality(dc)
  }

  orig_data <- exprs(orig)
  dc_data <- exprs(dc)
  fnames <- featureNames(orig)
  qdiff <- quality(dc)[2:5] - quality(orig)[2:5]

  log_text(paste("Inspecting drift correction results", Sys.time()))

  inspected <- foreach::foreach(i = seq_len(nrow(orig_data)), .combine = comb,
                                .export = c("%>%", "qdiff")) %dopar% {

                                  data <- orig_data[i, ]
                                  if (all(is.na(dc_data[i, ]))) {
                                    dc_note <- "Missing_QCS"
                                  } else if (any(dc_data[i, ] < 0, na.rm = TRUE)){
                                    dc_note <- "Negative_DC"
                                  } else if (check_quality) {
                                    pass <- paste0("qdiff[i, ] %>% dplyr::filter(", condition, ") %>% nrow() %>% as.logical()") %>%
                                      parse(text = .) %>% eval()
                                    if (!pass) {
                                      dc_note <- "Low_quality"
                                    } else {
                                      data <- dc_data[i, ]
                                      dc_note <- "Drift_corrected"
                                    }
                                  } else {
                                    data <- dc_data[i, ]
                                    dc_note <- "Drift_corrected"
                                  }

                                  list(data = matrix(data, nrow = 1, dimnames = list(fnames[i], names(data))),
                                       dc_notes = data.frame(Feature_ID = fnames[i],
                                                             DC_note = dc_note,
                                                             stringsAsFactors = FALSE))

                                }

  exprs(dc) <- inspected$data
  dc <- assess_quality(dc)
  dc <- join_fData(dc, inspected$dc_notes)

  log_text(paste("Drift correction results inspected at", Sys.time()))

  # Log information
  dc_note <- inspected$dc_notes$DC_note
  note_counts <- table(dc_note) %>% unname()
  note_percentage <- note_counts/sum(note_counts)
  note_percentage <- scales::percent(as.numeric(note_percentage))
  note_labels <- table(dc_note) %>% names()
  report <- paste(note_labels, note_percentage, sep = ": ", collapse = ",  ")
  log_text(paste0("\nDrift correction results inspected, report:\n", report))

  dc
}

#' Drift correction plots
#'
#' Plots the data before and after drift correction, with the regression line drawn with
#' the original data. If the drift correction was done on log-transformed data, then
#' plots of both the original and log-transformed data before and after correction are drawn.
#' The plot shows 2 standard deviation spread for both QC samples and regular samples.
#'
#' @param orig a MetaboSet object, before drift correction
#' @param dc a MetaboSet object, after drift correction as returned by correct_drift
#' @param predicted a matrix of predicted values, as returned by dc_cubic_spline
#' @param file path to the PDF file where the plots should be saved
#' @param log_transform logical, was the drift correction done on log-transformed data?
#' @param width,height width and height of the plots in inches
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param color_scale the color scale as returned by a ggplot function
#' @param shape_scale the shape scale as returned by a ggplot function
#'
#' @details If \code{shape} is set to \code{NULL} (the default), the column used for color
#' is also used for shape
#'
#' @seealso \code{\link{correct_drift}}, \code{\link{inspect_dc}}
#'
#' @examples
#' \dontrun{
#' dc <- dc_cubic_spline(merged_sample)
#' corrected <- dc$object
#' inspected <- inspect_dc(orig = merged_sample, dc = corrected,
#'                         check_quality = TRUE)
#' save_dc_plots(orig = merged_sample, dc = corrected, predicted = dc$predicted,
#'               file = "drift_plots.pdf")
#' }
#' @export
save_dc_plots <- function(orig, dc, predicted, file, log_transform = TRUE, width = 16, height = 8, color = "QC",
                          shape = color, color_scale = getOption("notame.color_scale_dis"),
                          shape_scale = scale_shape_manual(values = c(15, 16))) {

  if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
  }

  # Create a helper function for plotting
  dc_plot_helper <- function(data, fname, title = NULL) {
    p <- ggplot(mapping = aes_string(x = "Injection_order", y = fname)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      color_scale +
      shape_scale +
      labs(title = title)

    mean_qc <- finite_mean(data[data$QC == "QC", fname])
    sd_qc <- finite_sd(data[data$QC == "QC", fname])
    mean_sample <- finite_mean(data[data$QC != "QC", fname])
    sd_sample <- finite_sd(data[data$QC != "QC", fname])

    y_intercepts <- sort(c("-2 SD (Sample)" = mean_sample - 2 * sd_sample,
                           "-2 SD (QC)" = mean_qc - 2 * sd_qc,
                           "+2 SD (QC)" = mean_qc + 2 * sd_qc,
                           "+2 SD (Sample)" = mean_sample + 2 * sd_sample))

    for (yint in y_intercepts) {
      p <- p + geom_hline(yintercept = yint, color = "grey", linetype = "dashed")
    }
    p +
      scale_y_continuous(sec.axis = sec_axis(~ ., breaks = y_intercepts, labels = names(y_intercepts))) +
      geom_point(data = data, mapping = aes_string(color = color, shape = shape))
  }

  orig_data_log <- combined_data(log(orig))
  dc_data_log <- combined_data(log(dc))
  orig_data <- combined_data(orig)
  dc_data <- combined_data(dc)
  predictions <- as.data.frame(t(predicted))
  predictions$Injection_order <- orig_data$Injection_order

  pdf(file, width = width, height = height)

  for (fname in Biobase::featureNames(dc)) {

    p2 <- dc_plot_helper(data = dc_data, fname = fname,
                         title = "After")

    if (log_transform) {
      p1 <- dc_plot_helper(data = orig_data, fname = fname,
                           title = "Before")
      p3 <- dc_plot_helper(data = orig_data_log, fname = fname,
                           title = "Drift correction in log space") +
        geom_line(data = predictions, color = "grey")

      p4 <- dc_plot_helper(data = dc_data_log, fname = fname,
                           title = "Corrected data in log space")
      p <- cowplot::plot_grid(p1, p3, p2, p4, nrow = 2)
    } else {
      p1 <- dc_plot_helper(data = orig_data, fname = fname,
                           title = "Before (original values)") +
        geom_line(data = predictions, color = "grey")
      p <- cowplot::plot_grid(p1, p2, nrow = 2)
    }


    plot(p)
  }

  dev.off()

  log_text(paste("\nSaved drift correction plots to:", file))
}

#' Correct drift using cubic spline
#'
#' A wrapper function for applying cubic spline drift correction and saving
#' before and after plots
#'
#' @param object a MetaboSet object
#' @param log_transform logical, should drift correction be done on log-transformed values? See Details
#' @param spar smoothing parameter
#' @param spar_lower,spar_upper lower and upper limits for the smoothing parameter
#' @param check_quality logical, whether quality should be monitored.
#' @param condition a character specifying the condition used to decide whether drift correction
#' works adequately, see Details
#' @param plotting logical, whether plots should be drawn
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param color_scale,shape_scale the color and shape scales as returned by a ggplot function
#'
#' @return MetaboSet object as the one supplied, with drift corrected features
#'
#' @details If \code{log_transform = TRUE}, the correction will be done on log-transformed values.
#' The correction formula depends on whether the correction is run on original values or log-transformed values.
#' In log-space: \eqn{corrected = original + mean of QCs - prediction by cubic spline}.
#' In original space: \eqn{corrected = original * prediction for first QC / prediction for current point}.
#' We recommend doing the correction in the log-space since the log-transfomred data better follows the
#' assumptions of cubic spline regression. The drift correction in the original space also sometimes results
#' in negative values, and results in rejection of the drift corrrection procedure.
#' If \code{spar} is set to \code{NULL} (the default), the smoothing parameter will
#' be separately chosen for each feature from the range [\code{spar_lower, spar_upper}]
#' using cross validation. If  \code{check_quality = TRUE}, the \code{condition} parameter should be a character giving a condition compatible
#' with dplyr::filter. The condition is applied on the \strong{changes} in the quality metrics
#' RSD, RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r < 0 and D_ratio_r < 0",
#' meaning that both RSD_r and D_ratio_r need to decrease in the drift correction, otherwise the
#' drift corrected feature is discarded and the original is retained. If \code{shape} is set to \code{NULL} (the default), the column used for color
#' is also used for shape
#'
#' @examples
#' corrected <- correct_drift(merged_sample)
#'
#' @seealso  \code{\link{dc_cubic_spline}}, \code{\link[stats]{smooth.spline}} for details about the regression,
#' \code{\link{inspect_dc}} for analysing the drift correction results,
#' \code{\link{save_dc_plots}} for plotting the drift correction process for each feature
#'
#'
#' @export
correct_drift <- function(object, log_transform = TRUE, spar = NULL, spar_lower = 0.5, spar_upper = 1.5,
                          check_quality = FALSE, condition = "RSD_r < 0 & D_ratio_r < 0", plotting = FALSE,
                          file = NULL, width = 16, height = 8, color = "QC",
                          shape = NULL, color_scale = getOption("notame.color_scale_dis"),
                          shape_scale = scale_shape_manual(values = c(15, 16))) {
  # Fit cubic spline and correct
  corrected_list <- dc_cubic_spline(object, log_transform = log_transform, spar = spar,
                                    spar_lower = spar_lower, spar_upper = spar_upper)
  corrected <- corrected_list$object
  # Only keep corrected versions of features where drift correction increases quality
  inspected <- inspect_dc(orig = object, dc = corrected,
                          check_quality = check_quality, condition = condition)
  # Optionally save before and after plots
  if (plotting) {
    if (is.null(file)) {
      stop("File must be specified")
    }
    save_dc_plots(orig = object, dc = corrected, predicted = corrected_list$predicted,
                  file = file, log_transform = log_transform, width = width, height = height,
                  color = color, shape = shape, color_scale = color_scale, shape_scale = shape_scale)

  }
  # Return the final version
  inspected
}
