context("Testing reading data")

library(amp)

test_that("Column conversion works", {

  set.seed(38)
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

  df <- data.frame(Injection_order = seq_len(5),
                   Sample_ID = c(letters[1:5]))
  expect_warning(check_pheno_data(df), "QC column not found")

  df <- data.frame(Injection_order = seq_len(10),
                   Sample_ID = c(letters[1:5], rep("QC", 5)),
                   QC = as.factor(rep(c("Sample", "QC"), each = 5)))
  checked <- check_pheno_data(df)
  expected <- data.frame(Sample_ID = c(letters[1:5], paste0("QC_", 1:5)),
                         Injection_order = seq_len(10),
                         stringsAsFactors = FALSE,
                         QC = as.factor(rep(c("Sample", "QC"), each = 5)))
  rownames(expected) <- expected$Sample_ID
  expect_equal(checked, expected)

})


test_that("Easy example data is read correctly", {
  # Pheno data
  pd <- data.frame(Sample_ID = paste0("TEST_", seq_len(12)),
                   Injection_order = seq_len(12),
                   Group = factor(rep(LETTERS[1:2], times = c(5,7))),
                   QC = as.factor("Sample"),
                   easy_Datafile = paste0("190102SR_RP_pos_0", 10:21),
                   stringsAsFactors = FALSE)
  rownames(pd) <- pd$Sample_ID

  # Feature data
  fd <- data.frame(Feature_ID = "",
                   Split = "easy",
                   Alignment = as.numeric(seq_len(10)),
                   Mass = 50 * seq_len(10),
                   RetentionTime = 0.5 *seq_len(10),
                   "MS/MS Spectrum" = c("(123.45; 678)", rep(NA, 9)),
                   stringsAsFactors = FALSE)
  fd <- name_features(fd)
  rownames(fd) <- fd$Feature_ID

  # Assay data
  ad <- matrix(0, nrow = 10, ncol = 12)
  for (i in seq_len(10)) {
    ad[i, ] <- 50000  -(i-1)*3000 + 1000 * (0:11)
  }
  dimnames(ad) <- list(rownames(fd), rownames(pd))

  # Read the file
  read <- read_from_excel(system.file("extdata", "easy_data.xlsx", package = "amp"),
                          sheet = 1,
                          corner_row = 4, corner_column = "D",
                          name = "easy", id_prefix = "TEST_")

  # Test that the parts are read as expected
  expect_equal(read$exprs, ad)
  expect_equal(read$pheno_data, pd)
  expect_equal(read$feature_data, fd)
})

test_that("Data is split correctly", {

  # Pheno data
  pd <- data.frame(Sample_ID = paste0("TEST_", seq_len(12)),
                   Injection_order = seq_len(12),
                   Group = rep(LETTERS[1:2], times = c(5,7)),
                   QC = "Sample",
                   Datafile = paste0("190102SR_RP_pos_0", 10:21),
                   stringsAsFactors = FALSE)
  rownames(pd) <- pd$Sample_ID
  qc_idx <- c(1, 6, 9, 12)
  pd$Group[qc_idx] <- "QC"
  pd$QC[qc_idx] <- "QC"
  pd$QC <- factor(pd$QC)



  # Feature data
  fd <- data.frame(Feature_ID = "",
                   Split = rep(c("RP_pos", "Hilic_pos", "Hilic_neg", "RP_neg"), each = 4),
                   Alignment = as.numeric(seq_len(16)),
                   Mass = 50 * seq_len(16),
                   RetentionTime = 0.5 *seq_len(16),
                   Column = rep(c("RP", "Hilic", "RP"), times = c(4, 8, 4)),
                   Mode = rep(c("pos", "neg"), each = 8),
                   "MS/MS Spectrum" = c("(123.45; 678)", rep(NA, 15)),
                   stringsAsFactors = FALSE)
  fd <- name_features(fd)
  rownames(fd) <- fd$Feature_ID

  # Assay data
  ad <- matrix(0, nrow = 16, ncol = 12)
  for (i in seq_len(16)) {
    ad[i, ] <- 50000  -(i-1)*3000 + 1000 * (0:11)
  }
  row_zeros <- c(1, 3, 3, 4, 5, 8, 8, 10, rep(11, 5))
  col_zeros <- c(3, 7, 12, 1, 5, 10, 12, 1, 1, 3, 5, 8, 12)
  for (j in seq_along(row_zeros)) {
    ad[row_zeros[j], col_zeros[j]] <- 0
  }
  dimnames(ad) <- list(rownames(fd), rownames(pd))


  # Read the file
  read <- read_from_excel(system.file("extdata", "split_data.xlsx", package = "amp"),
                          sheet = 1,
                          corner_row = 4, corner_column = "F",
                          split_by = c("Column", "Mode"), id_prefix = "TEST_")

  # Test that the parts are read as expected
  expect_equal(read$exprs, ad)
  expect_equal(read$pheno_data, pd)
  expect_equal(read$feature_data, fd)

})


