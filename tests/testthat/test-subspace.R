library("dimension")

# Tests for check_dim_matrix default settting
# --------------------------------------------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 30, var = c(rep(10,5),rep(2,25)))

time_taken <- system.time({
	suppressWarnings(Subspace1 <- subspace(x1))
})

Subspace1_ref <- readRDS("reference_data/Subspace1.rds")

test_that("Default settings work as expected", {
	expect_is(Subspace1, "subspace")
  	expect_equal(Subspace1$ndf, 150)
  	expect_equal(Subspace1$pdim, 100)
  	expect_equal(Subspace1$components, 1:100)
  	expect_equivalent(Subspace1$var_correct, Subspace1_ref$var_correct, tolerance = 5e-2)
  	expect_equal(Subspace1$transpose_flag, TRUE)
  	expect_equivalent(Subspace1$irl$eigen, Subspace1_ref$irl$eigen, tolerance = 5e-2)
  	expect_equal(Subspace1$irl$dim, 1:100)
  	expect_equivalent(Subspace1$mp_irl$eigen, Subspace1_ref$mp_irl$eigen, tolerance = 5e-2)
  	expect_equal(Subspace1$mp_irl$dim, 1:100)
  	expect_equivalent(Subspace1$v, Subspace1_ref$v, tolerance = 5e-2)
  	expect_equivalent(Subspace1$u, Subspace1_ref$u, tolerance = 5e-2)
})


context("Argument input error")

test_that("Argument input error", {
expect_error(
		subspace(x1, "1"), "is.numeric")
	expect_error(
		subspace(x1, 1.1), "%%")
	expect_error(
		subspace(x1, 0), "larger")
	expect_error(
		subspace(x1, -1), "larger")
	expect_error(
		subspace(x1, 1000), "bounds")
	expect_error(
		subspace(x1, -1:20), "larger")
	expect_error(
		subspace(x1, 10:1000), "bounds")
	expect_error(
		subspace(x1, 0:1000), "larger")
})