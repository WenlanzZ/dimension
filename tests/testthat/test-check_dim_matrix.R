library("dimension")

# Tests for check_dim_matrix default settting
# --------------------------------------------
context("Valid check_dim_matrix input")

x1 <- x_sim(n = 150, p = 100, ncc = 10, var = 5)
x2 <- x_sim(n = 100, p = 150, ncc = 10, var = 5)

time_taken <- system.time({
  param_res1 <- check_dim_matrix(x1, rnk = 30)
  suppressWarnings(param_res2 <- check_dim_matrix(x2, rnk = 30))
})

test_that("invalid rnk input warnings returned", {
  expect_error(
    check_dim_matrix(x = NULL), regexp = "Invalid input")
  expect_error(
    check_dim_matrix(rnk = 300), regexp = "missing")
  expect_error(
    check_dim_matrix(x1, rnk = 300), regexp = "bounds")
  expect_error(
    check_dim_matrix(x1, rnk = -1), regexp = "positive")
  expect_error(
    check_dim_matrix(x1, rnk = 0), regexp = "positive")
  expect_message(
    check_dim_matrix(x1, verbose = TRUE),
    regexp = paste0("No component specified. ",
    "Calculating full singular value decomposition instead.\n"))
})


test_that("transpose messages", {

  expect_message(
    check_dim_matrix(x2, rnk = 30), regexp = "transpose")
})

context("Valid check_dim_matrix output")

test_that("output result", {

  expect_true(is.list(param_res1))
  expect_equal(param_res1$ndf, 150)
  expect_equal(param_res1$pdim, 100)
  expect_equal(param_res1$svr, 1.5)
  expect_equal(param_res1$rnk, 30)
  expect_equal(param_res1$transpose_flag, FALSE)

  expect_true(is.list(param_res2))
  expect_equal(param_res2$ndf, 150)
  expect_equal(param_res2$pdim, 100)
  expect_equal(param_res2$svr, 1.5)
  expect_equal(param_res2$rnk, 30)
  expect_equal(param_res2$transpose_flag, TRUE)
})
