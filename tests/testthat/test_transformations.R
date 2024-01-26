context("Testing transformations")

library(notame)

data("example_set")

test_that("Marking NAs works properly", {
  marked <- mark_nas(example_set, value = 0)
  zero_idx <- exprs(example_set) == 0
  na_idx <- is.na(exprs(marked))

  expect_equal(zero_idx, na_idx)
})


test_that("RF imputation works as expected", {
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
  single_change <- matrix(1111,
    nrow = 1, ncol = 1,
    dimnames = list(
      "HILIC_pos_108_1065a2_6121",
      "Demo_1"
    )
  )
  mrg <- merge_exprs(example_set, single_change)
  expect_equal(exprs(mrg)[2, 1], 1111)

  block_change <- matrix(runif(9),
    nrow = 3, ncol = 3,
    dimnames = list(
      rownames(exprs(example_set))[4:6],
      colnames(exprs(example_set))[5:7]
    )
  )
  mrg2 <- merge_exprs(example_set, block_change)
  expect_true(all(exprs(mrg2)[4:6, 5:7] == block_change))

  scattered_change <- matrix(runif(9),
    nrow = 3, ncol = 3,
    dimnames = list(
      rownames(exprs(example_set))[c(1, 5, 7)],
      colnames(exprs(example_set))[c(2, 5, 8)]
    )
  )
  mrg3 <- merge_exprs(example_set, scattered_change)
  expect_true(all(exprs(mrg3)[c(1, 5, 7), c(2, 5, 8)] == scattered_change))

  wrong_names <- matrix(runif(9), nrow = 3)
  expect_error(merge_exprs(example_set, wrong_names), "Column names")

  wrong_names2 <- matrix(runif(9),
    nrow = 3,
    dimnames = list(
      letters[1:3],
      colnames(exprs(example_set))[5:7]
    )
  )
  expect_error(merge_exprs(example_set, wrong_names2), "Row names")
})


test_that("Flagged compounds are not imputed", {
  marked <- mark_nas(example_set, 0)
  flag(marked)[c(1, 4, 6)] <- "Flagged"

  imputed <- impute_rf(marked)
  nas <- apply(exprs(imputed), 1, prop_na)

  expect_true(all(nas[c(1, 4, 6)] > 0))
  expect_true(all(nas[-c(1, 4, 6)] == 0))

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
    apply(1, function(x) {
      shapiro.test(x)$p.value
    })
  expect_true(all(shapiro_p > 0.9))
})

test_that("Simple imputation works as expected", {
  marked <- mark_nas(example_set, value = 0)
  imputed <- impute_simple(marked, value = 0)

  # Check that all missing values are imputed
  expect_equal(sum(is.na(exprs(imputed))), 0)

  # Check that non-missing values are unchanged
  na_idx <- is.na(exprs(marked))

  non_na_marked <- exprs(marked)[!na_idx]
  non_na_imputed <- exprs(imputed)[!na_idx]
  expect_equal(non_na_imputed, non_na_marked)
})

test_that("Simple imputation with only one feature works", {
  marked <- example_set
  exprs(marked)[2, 5] <- NA
  imputed <- impute_simple(marked, value = 0)
  #
  expect_equal(exprs(imputed)[2, 5], 0)
})

test_that("PQN normalization works correctly using median of QC samples as reference", {
  data <- exprs(example_set)
  # Calculate the median of QC samples
  reference <- apply(data[, example_set$QC == "QC"], 1, finite_median)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, exprs(pqn_normalization(example_set)))
})

test_that("PQN normalization works correctly using median of all samples as reference", {
  data <- exprs(example_set)
  # Calculate the median of all samples
  reference <- apply(data, 1, finite_median)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, exprs(pqn_normalization(example_set, ref = "all")))
})

test_that("PQN normalization works correctly using mean of QC samples as reference", {
  data <- exprs(example_set)
  # Calculate the mean of QC samples
  reference <- apply(data[, example_set$QC == "QC"], 1, finite_mean)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, exprs(pqn_normalization(example_set, method = "mean")))
})

test_that("PQN normalization works correctly using mean of all samples as reference", {
  data <- exprs(example_set)
  # Calculate the mean of all samples
  reference <- apply(data, 1, finite_mean)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, exprs(pqn_normalization(example_set, ref = "all", method = "mean")))
})

test_that("PQN normalization works with flagged features", {
  flagged_mset <- flag_quality(example_set)
  ref_data <- exprs(drop_flagged(flagged_mset))
  reference_spectrum <- apply(ref_data[, example_set$QC == "QC"], 1, finite_median)
  quotient <- ref_data / reference_spectrum
  quotient_md <- apply(quotient, 2, finite_median)

  data <- exprs(flagged_mset)
  pqn_data <- t(t(data) / quotient_md)

  pqn_mset <- pqn_normalization(flagged_mset)
  expect_equal(pqn_data, exprs(pqn_mset))
})
