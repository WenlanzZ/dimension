library("dimension")
context("Valid CheckDimMatrix input")

X <- Xsim(n = 150, p = 100, ncc = 10, var = 5)


test_that("invalid X input warnings returned", {
  
  expect_warning(
  	CheckDimMatrix(rnk = 30),
  	regexp="Invalid input X")
 
})


test_that("invalid rnk input warnings returned", {
  
  expect_warning(
  	CheckDimMatrix(X, rnk = 300),
  	regexp="Rnk out of bounds")
 
    expect_warning(
  	CheckDimMatrix(X, rnk = -1),
  	regexp="Rnk must be positive")

    expect_warning(
  	CheckDimMatrix(X, rnk = 0),
  	regexp="Rnk must be positive")

  	expect_warning(
  	CheckDimMatrix(X, rnk = 3.4),
  	regexp="Rnk must be positive")
})

reticulate::py_install(packages = 'umap-learn')