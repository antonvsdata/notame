test_plot_saving_helper <- function(func, title) {
  prefix <- paste0(tempdir(), "\\test\\")
  func(drop_qcs(example_set), prefix = prefix, format = "pdf", title = title)
  if (is.null(title)) {
    expect_equal(all(list.files(path = prefix) %in% paste0(1:nrow(fData(example_set)), ".pdf")), TRUE)
  } else {
    expect_equal(all(list.files(path = prefix) %in% paste0(featureNames(example_set), ".pdf")), TRUE)
  }
  unlink(prefix, recursive = TRUE)
}


test_file_extension_helper <- function(fileext) {
  p <- ggplot()
  file <- tempfile(fileext = fileext)
  save_plot(p, file)
  expect_equal(file.exists(file), TRUE)
  unlink(file)
}
