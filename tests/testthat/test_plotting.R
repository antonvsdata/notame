context("Testing plotting and saving functions")

# Testing save functions ----

test_that("Subject line plots are saved without title", {
  test_plot_saving_helper(save_subject_line_plots, title = NULL)
})

test_that("Subject line plot naming works", {
  test_plot_saving_helper(save_subject_line_plots, title = "Metabolite_name")
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
  test_plot_saving_helper(save_beeswarm_plots, title = "Metabolite_name")
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
