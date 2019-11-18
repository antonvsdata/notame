devtools::load_all("~/amp")


hilic_neg_sample$Batch <- factor(c(rep(1, 110), rep(2, 111)))

orig <- corrected <- hilic_neg_sample

save_batch_plots <- function(orig, corrected, file, batch = "Batch", color = "Batch", shape = "QC",
                         color_scale = NULL, shape_scale = NULL) {

  color_scale <- color_scale %||% getOption("amp.color_scale_dis")
  shape_scale <- shape_scale %||% scale_shape_manual(values = c(15, 16))

  data_orig <- combined_data(orig)
  data_corr <- combined_data(corrected)

  batch_injections <- data_orig %>%
    dplyr::group_by(!! sym(batch)) %>%
    dplyr::summarise(start = min(Injection_order), end = max(Injection_order))

  batch_means_orig <- data_orig %>%
    dplyr::group_by(!! sym(batch)) %>%
    dplyr::summarise_at(featureNames(hilic_neg_sample), finite_mean) %>%
    dplyr::left_join(batch_injections, ., by = batch)

  batch_means_corr <- data_corr %>%
    dplyr::group_by(!! sym(batch)) %>%
    dplyr::summarise_at(featureNames(hilic_neg_sample), finite_mean) %>%
    dplyr::left_join(batch_injections, ., by = batch)


  batch_plot_helper <- function(data, fname, batch_means) {
    p <- ggplot() +
      geom_point(data = data, mapping = aes_string(x = "Injection_order", y = fname,
                                                   color = color, shape = shape)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      color_scale +
      shape_scale

    p <- p +
      geom_segment(data = batch_means, mapping = aes_string(x = "start", xend = "end",
                                                            y = fname, yend = fname))
    p
  }

  pdf(file)

  for (feature in featureNames(orig)) {
    p1 <- batch_plot_helper(data_orig, feature, batch_means_orig)

    p2 <- batch_plot_helper(data_corr, feature, batch_means_corr)

    p <- cowplot::plot_grid(p1, p2, nrow = 2)
    plot(p)
  }

  dev.off()
}


