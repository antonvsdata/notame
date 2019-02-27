comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

#' @importFrom Biobase exprs exprs<-
correct_drift <- function(object, spar = NULL, spar_low = 0.5, spar_high = 1.5) {

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
                         spar = spar, control.spar = list("low" = spar_low, "high" = spar_high))
    predicted <- predict(fit, full_order)$y
    # Correction
    corr_factors <- predicted[1]/predicted
    corrected <- full_data[i, ] * corr_factors

    list(corrected = matrix(corrected, ncol = n, dimnames = dnames),
         predicted = matrix(predicted, ncol = n, dimnames = dnames))

  }

  exprs(object) <- dc_data$corrected
  predicted(object) <- dc_data$predicted
  object
}

inspect_dc <- function(orig, dc, condition) {

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


plot_dc <- function(orig, dc, file, width = 8, height = 6, group = group_col(orig),
                    color_scale = getOption("amp.color_scale_d")) {
  # If group column not set, use QC column
  group <- group %||% "QC"

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
      geom_point(data = orig_data, mapping = aes_string(color = group, shape = group)) +
      geom_line(data = predictions, color = "grey")

    p2 <- p +
      geom_point(data = dc_data, mapping = aes_string(color = group, shape = group))

    p <- gridExtra::arrangeGrob(p1, p2, nrow = 2)
    plot(p)
  }

  dev.off()

}



