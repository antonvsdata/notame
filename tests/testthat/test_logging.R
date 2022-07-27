context("Testing global options")

library(notame)

test_that("Error handling changes when logging and is set back to default after logging", {
  expect_equal(options()$error, NULL)
  log_file <- tempfile()
  init_log(log_file)
  expect_equal(class(options()$error), "call")
  finish_log()
  expect_equal(options()$error, NULL)
})



