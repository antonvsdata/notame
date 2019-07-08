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

test_that("Flagging works as expected", {

  single_change <- matrix(1111, nrow = 1, ncol = 1,
                          dimnames = list("HILIC_neg_108_1065a2_6121",
                                          "Demo_1"))
  mrg <- merge_exprs(example_set, single_change)
  expect_equal(exprs(mrg)[2, 1], 1111)

  block_change <- matrix(runif(9), nrow = 3, ncol = 3,
                         dimnames = list(rownames(exprs(example_set))[4:6],
                                         colnames(exprs(example_set))[5:7]))
  mrg2 <- merge_exprs(example_set, block_change)
  expect_true(all(exprs(mrg2)[4:6, 5:7] == block_change))

  scattered_change <- matrix(runif(9), nrow = 3, ncol = 3,
                             dimnames = list(rownames(exprs(example_set))[c(1,5,7)],
                                             colnames(exprs(example_set))[c(2,5,8)]))
  mrg3 <- merge_exprs(example_set, scattered_change)
  expect_true(all(exprs(mrg3)[c(1,5,7), c(2,5,8)] == scattered_change))

  wrong_names <- matrix(runif(9), nrow = 3)
  expect_error(merge_exprs(example_set, wrong_names), "Column names")

  wrong_names2 <- matrix(runif(9), nrow = 3,
                         dimnames = list(letters[1:3],
                         colnames(exprs(example_set))[5:7]))
  expect_error(merge_exprs(example_set, wrong_names2), "Row names")

})


test_that("Flagged compounds are not imputed", {

  marked <- mark_nas(example_set, 0)
  results(marked)$Flag[c(1,4,6)] <- "Flagged"

  imputed <- impute_rf(marked)
  nas <- apply(exprs(imputed), 1, prop_na)

  expect_true(all(nas[c(1,4,6)] > 0))
  expect_true(all(nas[-c(1,4,6)] == 0))

  # All feature parameter test
  imputed <- impute_rf(marked, all_features = TRUE)
  nas <- apply(exprs(imputed), 1, prop_na)
  expect_true(all(nas == 0))

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

