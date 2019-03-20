context("Testing reading data")

library(amp)

test_that("Column conversion works", {

  df <- data.frame(Injection = 1:100,
                   Group = letters[1:4],
                   Time = 1:2,
                   Sample_ID = sample(letters, size = 100, replace = TRUE),
                   few_numbers = c(1.3, 2.5))
  converted <- best_classes(df)
  converted_classes <- unname(sapply(converted, class))
  expected_classes <- c("numeric", "factor", "factor", "character", "numeric")

  expect_equal(converted_classes, expected_classes)

})


test_that("Pheno data checking works", {

  df <- data.frame(Group = letters[1:2])
  expect_error(check_pheno_data(df), '"Injection_order" not found')

  df <- data.frame(Injection_order = c(1:5, 3:9))
  expect_error(check_pheno_data(df), "Injection_order is not unique")

  df <- data.frame(Injection_order = seq_len(10),
                   Sample_ID = c(letters[1:5], letters[1:5]))
  expect_error(check_pheno_data(df), "Sample_ID is not unique")

  df <- data.frame(Injection_order = seq_len(10),
                   Sample_ID = c(letters[1:5], rep("QC", 5)))
  checked <- check_pheno_data(df)
  expected <- data.frame(Sample_ID = c(letters[1:5], paste0("QC_", 1:5)),
                         Injection_order = seq_len(10),
                         stringsAsFactors = FALSE)
  rownames(expected) <- expected$Sample_ID
  expect_equal(checked, expected)

})
