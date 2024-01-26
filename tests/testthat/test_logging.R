context("Testing logging")

library(notame)

test_that("Error/warning handling changes when logging and is set to NULL after logging", {
  # Save current error handling
  orig_error <- options("error")
  orig_warn <- options("warning.expression")
  log_file <- tempfile()
  init_log(log_file)
  # Initiating log should change error handling
  expect_failure(expect_identical(options("error")$error, orig_error$error))
  finish_log()
  # Finishing log should set error to "NULL"
  expect_null(options("error")$error)
  expect_null(options("warning.expression")$warning.expression)
})

test_that("Logging works as expected", {
  log_file <- tempfile()
  log_text("New line")
  # File doesn't exist yet
  expect_warning(expect_error(readLines(log_file), "cannot open the connection"))

  init_log(log_file)
  logs <- readLines(log_file)
  # Last line should match
  expect_identical(gsub("\\[.+\\] ", "", logs[length(logs)]), "INFO Starting logging")
  log_text("New line")
  logs <- readLines(log_file)
  # Last line should match
  expect_identical(gsub("\\[.+\\] ", "", logs[length(logs)]), "INFO New line")
  finish_log()
  logs <- readLines(log_file)
  log_len <- length(logs)
  # This shouldn't log to file anymore
  log_text("New line")
  logs <- readLines(log_file)
  expect_identical(log_len, length(logs))
})
