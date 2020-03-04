library("dimension")

# Tests for subspace default settting
# -------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))

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
  expect_equal(subspace1$irl$dim, 1:100)
  expect_equivalent(subspace1$mp_irl$eigen,
                    subspace1_ref$mp_irl$eigen,
                    tolerance = 5e-2)
  expect_equal(subspace1$mp_irl$dim, 1:100)
  expect_equivalent(subspace1$v, subspace1_ref$v, tolerance = 5e-2)
  expect_equivalent(subspace1$u, subspace1_ref$u, tolerance = 5e-2)
})


context("Argument input error")

test_that("components input error", {
  expect_error(subspace(components = 1), "missing")
  expect_error(subspace(x1, components = "1"), "is.numeric")
  expect_error(subspace(x1, components = 1.1), "%%")
  expect_error(subspace(x1, components = 0), "larger")
  expect_error(subspace(x1, components = -1), "larger")
  expect_error(subspace(x1, components = 1000), "bounds")
  expect_error(subspace(x1, components = -1:20), "larger")
  expect_error(subspace(x1, components = 10:1000), "bounds")
  expect_error(subspace(x1, components = 0:1000), "larger")
})

test_that("times input error", {
 expect_error(subspace(x1, times = "1"), "is.numeric")
 expect_error(subspace(x1, times = 1.1), "%%")
 expect_error(subspace(x1, times = 1:5), "length")
 expect_error(subspace(x1, times = -1), "positive")
 expect_error(subspace(x1, times = 200), "less")
})

context("plot subspace input error")

test_that("plot subspace annotation error", {
 expect_true(inherits(plot(subspace1), "ggplot"))
 expect_error(plot(subspace1, annotation = "0"), "numbers")
 expect_error(plot(subspace1, annotation = -1), "positive")
 expect_error(plot(subspace1, annotation = 110), "less")
})
