context("Testing constructor class")

library(notame)

data("example_set")

test_that("Changing feature names updates all occurrences", {
  set <- example_set
  featureNames(set) <- letters[1:20]
  expect_equal(featureNames(set), featureNames(assayData(set)))
  expect_equal(featureNames(set), featureNames(featureData(set)))
  expect_equal(featureNames(set), fData(set)$Feature_ID)
})

test_that("Changing feature names only works if valid names are given", {
  set <- example_set
  expect_error(featureNames(set) <- NULL)
  # Numbers are not allowed
  expect_error(featureNames(set) <- 1:20)
  # Number of new names should equal number of rows
  expect_error(featureNames(set) <- letters[1:15])
  # Duplicates are not allowed
  names <- letters[1:20]
  names[2] <- "a"
  expect_warning(expect_error(featureNames(set) <- names),
                 "non-unique value when setting"
  )
  # Names are not allowed to start with numbers
  names[2] <- "2a"
  expect_error(featureNames(set) <- names)
  # NAs are not allowed
  names[2] <- NA
  expect_error(featureNames(set) <- names)
})

test_that("Changing sample names updates all occurrences", {
  set <- example_set
  sampleNames(set) <- tolower(sampleNames(set))
  expect_equal(sampleNames(set), sampleNames(assayData(set)))
  expect_equal(sampleNames(set), sampleNames(protocolData(set)))
  expect_equal(sampleNames(set), sampleNames(phenoData(set)))
  expect_equal(sampleNames(set), pData(set)$Sample_ID)
})

test_that("Changing sample names only works if valid names are given", {
  set <- example_set
  names <- sampleNames(set)
  expect_error(sampleNames(set) <- NULL)
  # Number of new names should equal number of columns
  expect_error(sampleNames(set) <- names[1:29])
  # Duplicates are not allowed
  names[2] <- names[1]
  expect_warning(expect_error(sampleNames(set) <- names),
                 "non-unique value"
  )
  # NAs are not allowed
  names[2] <- NA
  expect_error(sampleNames(set) <- names)
})

