# Used to combine multiple obejcts returned from a foreach loop
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

#' Fit a cubic spline to correct drift
#'
#' Corrects the dirft in the features by applying smoothed cubic spline regression
#' to each feature separately.
#'
#' @param object a MetaboSet object
#' @param spar smoothing parameter
#' @param spar_lower,spar_upper lower and upper limits for the smoothing parameter
#'
#' @return MetaboSet object as the one supplied, with drift corrected fetures
#'
#' @details If \code{spar} is set to \code{NULL} (the default), the smoothing parameter will
#' be separately chosen for each feature from the range [\code{spar_lower, spar_upper}]
#' using cross validation.
#'
#' @seealso  \code{\link[stats]{smooth.spline}} for details about the regression,
#' \code{\link{inspect_dc}} for analysing the drift correction results,
#' \code{\link{plot_dc}} for plotting the drift correction process for each feature
#'
#' @importFrom Biobase exprs exprs<-
#' @seealso \code{\link{inspect_dc}}, \code{\link{save_dc_plots}}
#'
#' @export
dc_cubic_spline <- function(object, spar = NULL, spar_lower = 0.5, spar_upper = 1.5) {

  qc <- object[, object$QC == "QC"]
  qc_order <- qc$Injection_order
  qc_data <- exprs(qc)

  full_order <- object$Injection_order
  full_data <- exprs(object)

  n <- ncol(full_data)




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
    corr_factors <- predicted[1]/predicted
    corrected <- full_data[i, ] * corr_factors

    # Each iteration of the loop return one row to corrected and one row to predicted
    list(corrected = matrix(corrected, ncol = n, dimnames = dnames),
         predicted = matrix(predicted, ncol = n, dimnames = dnames))

  }

  exprs(object) <- dc_data$corrected
  predicted(object) <- dc_data$predicted
  object
}

#' Flag the results of drift correction
#'
#' Chooses whether the drift correction worked well enough by applying
#' the specified codition for each feature.
#' If the condition is fulfilled, the drift corrected feature is retained,
#' otherwise the original feature is retained and the drift corrrected feature is discarded.
#' The result of this operation is recorded in the feature data
#'
#' @param orig a MetaboSet object, before drift correction
#' @param dc a MetaboSet object, after drift correction
#' @param condition a character specifying the condition, see Details
#'
#' @return MeatboSet object
#'
#' @details The \code{condition} parameter should be a character giving a condition combatible
#' with dplyr::filter. The condition is applied on the \strong{changes} in the quality metrics
#' RSD, RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r < 0 and D_ratio_r < 0",
#' meaning that both RSD_r and D_ratio_r need to decrease in the drift correction, otherwise the
#' drift corrected feature is discarded and the original is retained.
#'
#' @seealso \code{\link{correct_drift}}, \code{\link{save_dc_plots}}
#'
#' @export
inspect_dc <- function(orig, dc, condition = "RSD_r < 0 & D_ratio_r < 0") {

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

  inspected <- foreach::foreach(i = seq_len(nrow(orig_data)), .combine = comb,
                                .export = c("%>%", "qdiff")) %dopar% {

    data = orig_data[i, ]
    if (all(is.na(dc_data[i, ]))) {
      dc_note <- "Missing_QCS"
    } else if (any(dc_data[i, ] < 0, na.rm = TRUE)){
      dc_note <- "Negative_DC"
    } else {
      pass <- paste0("qdiff[i, ] %>% dplyr::filter(", condition, ") %>% nrow() %>% as.logical()") %>%
        parse(text = .) %>% eval()
      if (!pass) {
        dc_note <- "Low_quality"
      } else {
        data <- dc_data[i, ]
        dc_note <- "Drift_corrected"
      }
    }

    list(data = matrix(data, nrow = 1, dimnames = list(fnames[i], names(data))),
         dc_notes = data.frame(Feature_ID = fnames[i],
                               DC_note = dc_note,
                               stringsAsFactors = FALSE))

  }

  exprs(dc) <- inspected$data
  dc <- assess_quality(dc)
  dc <- join_fdata(dc, inspected$dc_notes)
  dc

}

#' Drift correction plots
#'
#' Plots the data before and after drift correction, with the regression line drawn with
#' the original data.
#'
#' @param orig a MetaboSet object, before drift correction
#' @param dc a MetaboSet object, after drift correction as returned by correct_drift
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param color_scale the color scale as returned by a ggplot function
#'
#' @details If \code{shape} is set to \code{NULL} (the default), the column used for color
#' is also used for shape
#'
#' @seealso \code{\link{correct_drift}}, \code{\link{inspect_dc}}
#'
#' @export
save_dc_plots <- function(orig, dc, file, width = 8, height = 6, color = group_col(orig),
                    shape = NULL, color_scale = NULL) {
  # If color column not set, use QC column
  color <- color %||% "QC"
  shape <- shape %||% color
  color_scale <- color_scale %||% getOption("amp.color_scale_d")

  orig_data <- combined_data(orig)
  dc_data <- combined_data(dc)
  predictions <- as.data.frame(t(predicted(dc)))
  predictions$Injection_order <- orig_data$Injection_order

  pdf(file, width = width, height = height)

  for (fname in Biobase::featureNames(dc)) {
    p <- ggplot(mapping = aes_string(x = "Injection_order", y = fname)) +
      theme_bw() +
      color_scale


    p1 <- p +
      geom_point(data = orig_data, mapping = aes_string(color = color, shape = color)) +
      geom_line(data = predictions, color = "grey")

    p2 <- p +
      geom_point(data = dc_data, mapping = aes_string(color = color, shape = color))

    p <- gridExtra::arrangeGrob(p1, p2, nrow = 2)
    plot(p)
  }

  dev.off()

}

#' Correct drift using cubic spline
#'
#' A wrapper function for applying cubic spline drift correction and saving
#' before and after plots
#'
#' @param object a MetaboSet object
#' @param spar smoothing parameter
#' @param spar_lower,spar_upper lower and upper limits for the smoothing parameter
#' @param condition a character specifying the condition used to decide whether drift correction
#' works adequately, see Details
#' @param plotting logical, whether plots should be drawn
#' @param file path to the PDF file where the plots should be saved
#' @param width,height width and height of the plots in inches
#' @param color character, name of the column used for coloring the points
#' @param shape character, name of the column used for shape
#' @param color_scale the color scale as returned by a ggplot function
#'
#' @return MetaboSet object as the one supplied, with drift corrected fetures
#'
#' @details If \code{spar} is set to \code{NULL} (the default), the smoothing parameter will
#' be separately chosen for each feature from the range [\code{spar_lower, spar_upper}]
#' using cross validation. The \code{condition} parameter should be a character giving a condition combatible
#' with dplyr::filter. The condition is applied on the \strong{changes} in the quality metrics
#' RSD, RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r < 0 and D_ratio_r < 0",
#' meaning that both RSD_r and D_ratio_r need to decrease in the drift correction, otherwise the
#' drift corrected feature is discarded and the original is retained. If \code{shape} is set to \code{NULL} (the default), the column used for color
#' is also used for shape
#'
#' @seealso  \code{\link{correct_drift}}, \code{\link[stats]{smooth.spline}} for details about the regression,
#' \code{\link{inspect_dc}} for analysing the drift correction results,
#' \code{\link{plot_dc}} for plotting the drift correction process for each feature
#'
#'
#' @export
correct_drift <- function(object, spar = NULL, spar_lower = 0.5, spar_upper = 1.5,
                          condition = "RSD_r < 0 & D_ratio_r < 0", plotting = FALSE,
                          file = NULL, width = 8, height = 6, color = group_col(orig),
                          shape = NULL, color_scale = NULL) {

  # Fit cubic spline and correct
  corrected <- dc_cubic_spline(object, spar = spar,
                               spar_lower = spar_lower, spar_upper = spar_upper)
  # Only keep corrected versions of features where drift correction increases quality
  inspected <- inspect_dc(object, corrected, condition = condition)
  # Optionally save before and after plots
  if (plotting) {
    if (is.null(file)) {
      stop("File must be specified")
    }
    save_dc_plots(object, corrected, file = file, width = width, height = height,
                  color = color, shape = shape, color_scale = color_scale)

  }
  # Return the final version
  inspected
}
