context("Testing statistics")

library(amp)

test_that("summary statistics work without grouping", {

  smry <- summary_statistics(mark_nas(example_set, 0),
                             grouping_cols = NA)
  ex <- exprs(mark_nas(example_set, 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(ex, 1, fun)), smry[, gsub("finite_", "", fun)])
  }
  expect_equal(unname(apply(ex, 1, finite_quantile, probs = 0.25)), smry$Q25)
  expect_equal(unname(apply(ex, 1, finite_quantile, probs = 0.75)), smry$Q75)
})


test_that("summary statistics work with grouping", {

  smry <- summary_statistics(mark_nas(example_set, 0))

  exa <- exprs(mark_nas(example_set[, example_set$Group == "A"], 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(exa, 1, fun)), smry[, gsub("finite", "Group_A", fun)])
  }
  expect_equal(unname(apply(exa, 1, finite_quantile, probs = 0.25)), smry$Group_A_Q25)
  expect_equal(unname(apply(exa, 1, finite_quantile, probs = 0.75)), smry$Group_A_Q75)

  exb <- exprs(mark_nas(example_set[, example_set$Group == "B"], 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(exb, 1, fun)), smry[, gsub("finite", "Group_B", fun)])
  }
  expect_equal(unname(apply(exb, 1, finite_quantile, probs = 0.25)), smry$Group_B_Q25)
  expect_equal(unname(apply(exb, 1, finite_quantile, probs = 0.75)), smry$Group_B_Q75)
})

test_that("summary statistics work with all NA features", {

  ex_set_na <- mark_nas(example_set, 0)
  exprs(ex_set_na)[1, ] <- NA

  smry <- summary_statistics(ex_set_na)

  expect_equal(nrow(smry), nrow(exprs(ex_set_na)))
  expect_equal(smry$Feature_ID, featureNames(ex_set_na))
  expect_true(all(is.na(smry[1, 2:ncol(smry)])))
})




