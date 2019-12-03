test_that("multiplication works", {

  X <- Xsim(n = 150, p = 100, ncc = 10, var = 5)
  params <- CheckDimMatrix(X, rnk = 30)


  expect_equal(2 * 2, 4)
})
