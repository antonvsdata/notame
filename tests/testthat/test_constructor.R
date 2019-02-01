context("Testing class constructor")

library(amp)

# Toy data for testing
dada <- data.frame(Injection_order = 1:5,
                   Group = c("QC", "A", "A", "B", "B"),
                   c1 = rnorm(5),
                   c2 = rnorm(5))
response_info <- data.frame(Response = c("c1", "c2"), Mass = 1:2, RT = 3:4)
responses <- c("c1", "c2")

test_that("Works for toy example",{


  should <- list(data = dada, response_info = response_info, responses = responses,
                 response_name = "Response", order = "Injection_order",
                 group = "Group", sample_id = "Injection_order",
                 qc = "Group", time = NULL, subject_id = NULL)
  class(should) <- c("lcms_data", "list")

  is <- lcms_data(data = dada, response_info = response_info, responses = responses,
                  response_name_col = "Response", order_col = "Injection_order",
                  group_col = "Group")
  expect_identical(is, should)
})


test_that("Works with all columns defined", {
  dada_tmp <- dada
  dada_tmp$Visit <- 1:5
  dada_tmp$Subject_ID <- c(1,1,2,2,3)
  dada_tmp$Sample_ID <- 1:5
  dada_tmp$QC <- "QC"

  should <- list(data = dada_tmp, response_info = response_info, responses = responses,
                 response_name = "Response", order = "Injection_order",
                 group = "Group", sample_id = "Sample_ID",
                 qc = "QC", time = "Visit", subject_id = "Subject_ID")
  class(should) <- c("lcms_data", "list")

  is <- lcms_data(data = dada_tmp, response_info = response_info, responses = responses,
                  response_name_col = "Response", order_col = "Injection_order",
                  group_col = "Group", qc_col = "QC", time_col = "Visit", sample_id_col = "Sample_ID",
                  subject_id_col = "Subject_ID")
  expect_identical(is, should)
})

test_that("Correct errors for missing columns", {
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                  response_name_col = "Response", order_col = "Missing",
                  group_col = "Group"), "columns not found in data")
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Missing"), "columns not found in data")
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group", subject_id_col = "Missing"), "columns not found in data")
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection.order",
                         group_col = "Group", qc_col = "Missing"), "columns not found in data")
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group", time_col = "Missing"), "columns not found in data")
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group", sample_id_col = "Missing"), "columns not found in data")
  expect_error(lcms_data(data = dada, response_info = response_info, responses = responses,
                         response_name_col = "Mising", order_col = "Injection_order",
                         group_col = "Group"), "response_name_col")
})


test_that("Response info errors work", {
  expect_error(lcms_data(data = dada, response_info = response_info, responses = c("c1"),
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group"), "length of responses")

  response_info_tmp <- data.frame(Response = c("c1", "c2", "c3"), Mass = 1:3, RT = 3:5)

  expect_error(lcms_data(data = dada, response_info = response_info_tmp, responses = c("c1", "c2", "c4"),
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group"), "responses do not have matching")
  expect_error(lcms_data(data = dada, response_info = response_info_tmp, responses = c("c1", "c2", "c3"),
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group"), "responses not found in column names")

})

test_that("Problems with data detected", {
  dada_tmp <- dada
  dada_tmp$Injection_order <- 1
  expect_error(lcms_data(data = dada_tmp, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group"), "Injection order is not unique")

  dada_tmp <- dada
  dada_tmp$Group <- "A"
  expect_error(lcms_data(data = dada_tmp, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order",
                         group_col = "Group"), "not have any 'QC' labels")

  dada_tmp <- dada
  dada_tmp$Sample_ID <- c(1:4, 4)
  expect_error(lcms_data(data = dada_tmp, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order", sample_id_col = "Sample_ID",
                         group_col = "Group"), "Sample ID is not ")

  dada_tmp <- dada
  dada_tmp$QC <- 1:5
  expect_error(lcms_data(data = dada_tmp, response_info = response_info, responses = responses,
                         response_name_col = "Response", order_col = "Injection_order", qc_col = "QC",
                         group_col = "Group"), "not have any 'QC' labels")
})


