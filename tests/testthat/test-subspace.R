library("dimension")

# Tests for subspace default settting
# -------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
x2 <- x_sim(n = 150, p = 100, ncc = 10, var = 5)

time_taken <- system.time({
  suppressWarnings(subspace1 <- subspace(x1))
})

subspace1_ref <- readRDS("reference_data/Subspace1.rds")

test_that("Default settings work as expected", {
  expect_is(subspace1, "subspace")
  expect_equal(subspace1$ndf, 150)
  expect_equal(subspace1$pdim, 100)
  expect_equal(subspace1$components, 1:100)
  expect_equivalent(subspace1$var_correct,
                    subspace1_ref$var_correct,
                    tolerance = 5e-2)
  expect_equal(subspace1$transpose_flag, TRUE)
  expect_equivalent(subspace1$irl$eigen,
                    subspace1_ref$irl$eigen,
                    tolerance = 5e-2)
  expect_equivalent(subspace1$irl$dim, 1:100)
  expect_equivalent(subspace1$mp_irl$eigen,
                    subspace1_ref$mp_irl$eigen,
                    tolerance = 5e-2)
  expect_equal(subspace1$mp_irl$dim, 1:100)
  expect_equivalent(subspace1$v, subspace1_ref$v, tolerance = 5e-2)
  expect_equivalent(subspace1$u, subspace1_ref$u, tolerance = 5e-2)
  expect_equal(print.subspace(subspace1), subspace1)
})


context("Argument input error")

test_that("class input error", {
 expect_message(subspace(x2, verbose = TRUE), "full")
 expect_error(correct_eigenvalues(x1), "type")
})


test_that("components input error", {
  expect_error(subspace(components = 1), "missing")
  expect_message(create_subspace(x2, verbose = TRUE), "specified")
  expect_message(subspace(x2, verbose = TRUE), "full")
  expect_message(subspace(x2, components = NULL, verbose = TRUE), "full")
  expect_message(subspace(x2, verbose = TRUE), "mp")
  expect_error(subspace(x1, components = "1"), "is.numeric")
  expect_error(subspace(x1, components = 1.1), "%%")
  expect_error(subspace(x1, components = 0), "larger")
  expect_error(subspace(x1, components = -1), "larger")
  expect_error(subspace(x1, components = 1000), "bounds")
  expect_error(subspace(x1, components = -1:20), "larger")
  expect_error(subspace(x1, components = 10:1000), "bounds")
  expect_error(subspace(x1, components = 0:1000), "larger")
  expect_message(subspace(x2, components = c(1, 3, 5)), "range")
})

test_that("num_est_samples input error", {
 expect_error(subspace(x2, num_est_samples = "1"), "is.numeric")
 expect_error(subspace(x2, num_est_samples = 1.1), "%%")
 expect_error(subspace(x2, num_est_samples = 1:5), "length")
 expect_error(subspace(x2, num_est_samples = -1), "positive")
 expect_error(subspace(x2, num_est_samples = 200), "smaller")
})

context("plot subspace input error")

test_that("plot subspace annotation error", {
 expect_true(inherits(plot(subspace1), "ggplot"))
 expect_true(inherits(plot(subspace1, changepoint = 10), "ggplot"))
 expect_error(plot(subspace1, annotation = "0"), "numbers")
 expect_error(plot(subspace1, annotation = -1), "positive")
 expect_error(plot(subspace1, annotation = 110), "less")
 expect_error(plot(subspace1, changepoint = "0"), "is.numeric")
 expect_error(plot(subspace1, changepoint = 1.5), "%%")
 expect_error(plot(subspace1, changepoint = c(0, 1)), "length")
 expect_error(plot(subspace1, changepoint = -1), "positive")
 expect_error(plot(subspace1, changepoint = 110), "less")
})



context("Memory allocate error")

# test_that("Memory allocate error", {
#   expect_message(subspace(x), "allocate")
# })
