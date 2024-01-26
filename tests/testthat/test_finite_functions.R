context("Finite helper functions")

library(notame)

test_that("Work correctly for full vectors", {
  x <- 13:20
  expect_equal(finite_sd(x), sd(x))
  expect_equal(finite_mean(x), mean(x))
  expect_equal(finite_median(x), median(x))
  expect_equal(finite_mad(x), mad(x))
})

test_that("NAs and finites are excluded correctly", {
  x <- 1:10
  x[c(2, 7)] <- Inf
  x[c(4, 9)] <- NA
  y <- c(1, 3, 5, 6, 8, 10)
  expect_equal(finite_sd(x), sd(y))
  expect_equal(finite_mean(x), mean(y))
  expect_equal(finite_median(x), median(y))
  expect_equal(finite_mad(x), mad(y))
})
