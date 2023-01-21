context("Testing plotting and saving functions")

# Testing save_plot helper function ----

test_that("Recursive folder creation works", {
  p <- ggplot()
  folder <- paste0(tempdir(), "\\test\\recursive")
  expect_equal(file.exists(folder), FALSE)
  file <- tempfile(tmpdir = folder, fileext = ".pdf")
  save_plot(p, file)
  expect_equal(file.exists(folder), TRUE)
  unlink(folder, recursive = TRUE)
})

test_that("Creating pdf files work", {
  p <- ggplot()
  file <- tempfile(fileext = ".pdf")
  save_plot(p, file)
  expect_equal(file.exists(file), TRUE)
  unlink(file)
})

test_that("Creating emf files work", {
  test_file_extension_helper(".emf")
})

test_that("Creating svg files work", {
  test_file_extension_helper(".svg")
})

test_that("Creating png files work", {
  test_file_extension_helper(".png")
})

test_that("Creating tiff files work", {
  test_file_extension_helper(".tiff")
})

test_that("Giving invalid file format throws error", {
  p <- ggplot()
  file <- tempfile(fileext = ".jpeg")
  expect_error(save_plot(p, file), "is not valid")
  expect_equal(file.exists(file), FALSE)
  unlink(file)
})

# Testing save functions ----

test_that("Subject line plots are saved without title", {
  test_plot_saving_helper(save_subject_line_plots, title = NULL)
})

test_that("Subject line plot naming works", {
  test_plot_saving_helper(save_subject_line_plots, title = "Feature_ID")
})

test_that("Group boxplots are saved without title", {
  test_plot_saving_helper(save_group_boxplots, title = NULL)
})

test_that("Group boxplot naming works", {
  test_plot_saving_helper(save_group_boxplots, title = "Feature_ID")
})

test_that("Beeswarm plots are saved without title", {
  test_plot_saving_helper(save_beeswarm_plots, title = NULL)
})

test_that("Beeswarm plot naming works", {
  test_plot_saving_helper(save_beeswarm_plots, title = "Feature_ID")
})

test_that("Scatter plots are saved without title", {
  test_plot_saving_helper(save_scatter_plots, title = NULL)
})

test_that("Scatter plot naming works", {
  test_plot_saving_helper(save_scatter_plots, title = "Feature_ID")
})

test_that("Group lineplots are saved without title", {
  test_plot_saving_helper(save_group_lineplots, title = NULL)
})

test_that("Group lineplot naming works", {
  test_plot_saving_helper(save_group_lineplots, title = "Feature_ID")
})
