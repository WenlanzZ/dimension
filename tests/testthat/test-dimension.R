library("dimension")

# Tests for dimension default settting
# ------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))

time_taken <- system.time({
  suppressWarnings(results1 <- dimension(x1))
})

subspace1_ref <- readRDS("reference_data/Subspace1.rds")

test_that("Default settings work as expected", {
  expect_is(results1$Subspace, "subspace")
  expect_equal(results1$Subspace$ndf, 150)
  expect_equal(results1$Subspace$pdim, 100)
  expect_equal(results1$Subspace$components, subspace1_ref$components)
  expect_equivalent(results1$Subspace$var_correct,
                    subspace1_ref$var_correct,
                    tolerance = 5e-2)
  expect_equal(results1$Subspace$transpose_flag, TRUE)
  expect_equivalent(results1$Subspace$irl$eigen,
                    subspace1_ref$irl$eigen,
                    tolerance = 5e-2)
  expect_equal(results1$Subspace$irl$dim, 1:100)
  expect_equivalent(results1$Subspace$mp_irl$eigen,
                    subspace1_ref$mp_irl$eigen,
                    tolerance = 5e-2)
  expect_equal(results1$Subspace$mp_irl$dim, 1:100)
  expect_equivalent(results1$Subspace$v, subspace1_ref$v, tolerance = 5e-2)
  expect_equivalent(results1$Subspace$u, subspace1_ref$u, tolerance = 5e-2)
  expect_equal(results1$dimension, 5)
})


context("Argument input error")

test_that("Argument input error", {

  expect_error(
    dimension(), "missing")
  expect_message(
    dimension(subspace_ = subspace1_ref), "already")

  subspace1_ref$mp_irl <- NULL
  expect_error(
    dimension(subspace_ = subspace1_ref), "missing")
})
