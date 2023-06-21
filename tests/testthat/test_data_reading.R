context("Testing reading data")

library(notame)

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
  # Check that error is thrown if injection order is not numeric
  df <- data.frame(Injection_order = c(1:2, "a", 4:9))
  expect_error(check_pheno_data(df), '"Injection_order" is not numeric')
  # Check if conversion to numeric works
  df <- data.frame(Injection_order = c(1:2, "3", 4:9))
  expect_warning(expect_visible(check_pheno_data(df, id_prefix = "ID_")),
    "QC not found"
    )
  # Check QC generator
  df <- data.frame(Injection_order = seq_len(10),
                   Sample_ID = c(letters[1:5], letters[1:5])
  )
  expect_warning(expect_error(check_pheno_data(df), "Sample_ID is not unique"),
                 "QC not found"
  )
  df <- data.frame(Injection_order = seq_len(5),
                   Sample_ID = c(letters[1:5])
  )
  expect_warning(check_pheno_data(df), "QC not found")

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


test_that("Feature data checking works", {
  df <- data.frame(Feature_ID = 1:5)
  expect_error(check_feature_data(df), "Numbers are not allowed as feature IDs")

  df <- data.frame(Feature_ID = c(letters[1:5], letters[3:5]))
  expect_error(check_feature_data(df), "Feature_ID values are not unique")

  df <- data.frame(Feature_ID = letters[1:9])
  df$Feature_ID[6] <- NA
  expect_error(check_feature_data(df), "Missing values in Feature IDs")

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
                   "MS_MS_Spectrum" = c("(123.45; 678)", rep(NA, 9)),
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
  read <- read_from_excel(system.file("extdata", "easy_data.xlsx", package = "notame"),
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
                   "MS_MS_Spectrum" = c("(123.45; 678)", rep(NA, 15)),
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
  read <- read_from_excel(system.file("extdata", "split_data.xlsx", package = "notame"),
                          sheet = 1,
                          corner_row = 4, corner_column = "F",
                          split_by = c("Column", "Mode"), id_prefix = "TEST_")

  # Test that the parts are read as expected
  expect_equal(read$exprs, ad)
  expect_equal(read$pheno_data, pd)
  expect_equal(read$feature_data, fd)

})

test_that("Splitting data works as expected", {
  split_by <- c("Ion mode", "gswregh") # Wrong column name
  expect_error(read_from_excel(system.file("extdata", "sample_data_whole.xlsx",
                                           package = "notame"),
                               corner_row = 4,
                               corner_column = "X",
                               split_by = split_by)
  )
})

test_that("Creating dummy injection order works as expected", {
  names <- list("HILIC_neg", "HILIC_pos", "RP_neg", "RP_pos")
  modes <- list()
  for (name in names) {
    file <- system.file("extdata", paste0(name, "_sample.xlsx"), package = "notame")
    mode <- read_from_excel(file, name = name)
    modes[name] <- construct_metabosets(mode$exprs, mode$pheno_data, mode$feature_data)
  }
  # Modify data
  modes$HILIC_neg$Injection_order <- modes$HILIC_neg$Injection_order + 1
  inj_ord_rn <- modes$RP_neg$Injection_order + 2
  modes$RP_neg$Injection_order <- inj_ord_rn
  inj_ord_rp <- modes$RP_pos$Injection_order[5:221] + 5
  modes$RP_pos$Injection_order[5:221] <- inj_ord_rp
  sampleNames(modes$HILIC_neg)[2] <- "ID_666"
  sampleNames(modes$RP_pos)[22] <- "ID_999"

  expect_warning(merged <- merge_metabosets(modes),
                 regexp = "Sample IDs are not identical|Unequal amount of samples"
  )
  # Dummy injection
  expect_equal(merged$Injection_order, -seq_along(merged$Sample_ID))
  # Original IOs
  expect_equal(sort(as.numeric(na.omit(merged$HILIC_neg_Injection_order))),
               modes$HILIC_neg$Injection_order
  )
  expect_equal(sort(as.numeric(na.omit(merged$HILIC_pos_Injection_order))),
               modes$HILIC_pos$Injection_order
  )
  expect_equal(sort(as.numeric(na.omit(merged$RP_neg_Injection_order))),
               modes$RP_neg$Injection_order
  )
  expect_equal(sort(as.numeric(na.omit(merged$RP_pos_Injection_order))),
               modes$RP_pos$Injection_order
  )
})


