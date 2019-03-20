context("Testing transformations")

library(amp)

data("example_set")

test_that("Marking NAs works properly", {

  marked <- mark_nas(example_set, value = 0)
  zero_idx <- exprs(example_set) == 0
  na_idx <- is.na(exprs(marked))

  expect_equal(zero_idx, na_idx)
})


test_that("Imputation works as expected", {
  marked <- mark_nas(example_set, value = 0)
  imputed <- impute_rf(marked)

  # Check that all missing values are imputed
  expect_equal(sum(is.na(exprs(imputed))), 0)

  # Check that non-missing values are unchanged
  na_idx <- is.na(exprs(marked))

  non_na_marked <- exprs(marked)[!na_idx]
  non_na_imputed <- exprs(imputed)[!na_idx]
  expect_equal(non_na_imputed, non_na_marked)
})


test_that("Inverse normalization works as expected", {

  marked <- mark_nas(example_set, value = 0)
  imputed <- impute_rf(marked)
  normalized <- inverse_normalize(imputed)

  # Ranks should be identical
  orig_ranks <- exprs(imputed) %>%
    apply(1, rank) %>%
    t()
  norm_ranks <- exprs(normalized) %>%
    apply(1, rank) %>%
    t()
  expect_identical(norm_ranks, orig_ranks)

  # Check that all the rows follow standard normal distribution
  # Zero means up to 0.1 accuracy
  abs_means <- abs(rowMeans(exprs(normalized)))
  expect_true(all(abs_means < 0.1))
  # Unit variance
  sd_diffs <- abs(apply(exprs(normalized), 1, sd) - 1)
  expect_true(all(sd_diffs < 0.1))

  # Shapiro-Wilk normality tests are non-significant
  shapiro_p <- exprs(normalized) %>%
    apply(1, function(x) { shapiro.test(x)$p.value})
  expect_true(all(shapiro_p > 0.9))
})

