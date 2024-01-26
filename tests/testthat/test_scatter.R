context("Testing scatter.R plots")

library(notame)

test_that("mz_rt_plot uses correct data", {
  gg <- mz_rt_plot(merged_sample)
  expect_equal(fData(merged_sample), gg$data)
})

test_that("mz_rt_plot works with right objects", {
  expect_visible(mz_rt_plot(fData(merged_sample)))
  expect_visible(mz_rt_plot(merged_sample))
})

test_that("Low quality features are dropped when they should", {
  foreach::registerDoSEQ()
  flagged <- flag_quality(merged_sample)
  no_low <- drop_flagged(flagged)
  gg_flagged <- mz_rt_plot(flagged, all_features = FALSE)
  gg_all <- mz_rt_plot(flagged, all_features = TRUE)

  expect_equal(gg_flagged$data, fData(no_low))
  expect_equal(gg_all$data, fData(flagged))
})
