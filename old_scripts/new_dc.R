save_dc_plots <- function(orig, dc, predicted, file, width = 8, height = 6, color = "QC",
                          shape = NULL, color_scale = NULL, shape_scale = NULL) {
  # If color column not set, use QC column
  color <- color %||% "QC"
  shape <- shape %||% color
  color_scale <- color_scale %||% getOption("amp.color_scale_dis")
  shape_scale <- shape_scale %||% scale_shape_manual(values = c(15, 16))

  # Create a helper function for plotting
  dc_plot_helper <- function(data, fname) {
    p <- ggplot(mapping = aes_string(x = "Injection_order", y = fname)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      color_scale +
      shape_scale

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
    p1 <- dc_plot_helper(data = orig_data_log, fname = fname) +
      geom_line(data = predictions, color = "grey")

    p2 <- dc_plot_helper(data = dc_data_log, fname = fname)
    p3 <- dc_plot_helper(data = orig_data, fname = fname)
    p4 <- dc_plot_helper(data = dc_data, fname = fname)

    p <- cowplot::plot_grid(p3, p1, p4, p2, nrow = 2)
    plot(p)
  }

  dev.off()

  log_text(paste("\nSaved drift correction plots to:", file))
}






dc_cubic_spline <- function(object, log_transform = TRUE, spar = NULL, spar_lower = 0.5, spar_upper = 1.5) {

  # Start log
  log_text(paste("\nStarting drift correction at", Sys.time()))

  # Zero values do not behave correctly
  if (sum(exprs(object) == 0, na.rm = TRUE)) {
    log_text("Zero values in feature abundances detected. Zeroes will be replaced with 1.1")
    exprs(object)[exprs(object) == 0] <- 1.1
  }

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
    #corr_factors <- predicted[1]/predicted
    #corrected <- full_data[i, ] * corr_factors
    corrected <- full_data[i, ] + mean(qc_data[i, qc_detected]) - predicted

    # Each iteration of the loop return one row to corrected and one row to predicted
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


ex <- mark_nas(merged_sample, 0)
exprs(ex)[exprs(ex) == 0] <- 1.1
ex <- ex[, 1:218]
dc <- dc_cubic_spline(ex)
predicted <- dc$predicted
dc <- dc$object

save_dc_plots(ex, dc, predicted, file = "xd.pdf", width = 18, height = 8)


