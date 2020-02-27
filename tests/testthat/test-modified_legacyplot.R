library("dimension")

# Tests for modified_legacyplot default settting
# -----------------------------------------------

set.seed(seed = 1234)

test <- bcp(as.vector(c(rep(10, 10), 9.5, rep(0, 10))))

context("Argument input error")

test_that("Argument input error", {

  expect_error(
    modified_legacyplot(test, annotation = "1"), "numbers")
  expect_error(
    modified_legacyplot(test, annotation = 100), "less")
})
