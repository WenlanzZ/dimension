library("dimension")
context("Valid CheckDimMatrix input")

x <- x_sim(n = 150, p = 100, ncc = 10, var = 5)


test_that("invalid X input warnings returned", {

  expect_warning(
    check_dim_matrix(rnk = 30), regexp = "Invalid input X")

})


test_that("invalid rnk input warnings returned", {

  expect_warning(
    check_dim_matrix(x, rnk = 300), regexp = "Rnk out of bounds")

    expect_warning(
      check_dim_matrix(x, rnk = -1), regexp = "Rnk must be positive")

    expect_warning(
      check_dim_matrix(x, rnk = 0), regexp = "Rnk must be positive")

  	expect_warning(
  	  check_dim_matrix(x, rnk = 3.4), regexp = "Rnk must be positive")
})
