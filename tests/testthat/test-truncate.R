library("dimension")

# Tests for truncate default settting
# --------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
x2 <- x_sim(n = 150, p = 100, ncc = 10, var = 5)

subspace1_ref <- readRDS("reference_data/Subspace1.rds")

time_taken <- system.time({
 suppressWarnings(x_denoised1 <- truncate(x1,
                                          components = 20,
                                          method = "threshold",
                                          alpha = 0.9,
                                          zeroout = TRUE))
 suppressWarnings(x_denoised2 <- truncate(x1,
                                          components = 20,
                                          method = "hard",
                                          zeroout = FALSE))
 suppressWarnings(x_denoised3 <- truncate(x1,
                                          components = 100,
                                          method = "identity",
                                          location = c(1:15),
                                          zeroout = FALSE))
 suppressWarnings(x_denoised4 <- truncate(s = subspace1_ref,
                                          components = 20,
                                          method = "identity",
                                          location = c(1:15),
                                          zeroout = FALSE))
})

x_denoised1_ref <- readRDS("reference_data/x_denoised1.rds")

test_that("Threshold default settings work as expected", {
 expect_is(x_denoised1, "subspace_denoised")
 expect_equivalent(x_denoised1$xi_denoised,
                   x_denoised1_ref$xi_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised1$x_denoised,
                   x_denoised1_ref$x_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised1$v_denoised,
                   x_denoised1_ref$v_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised1$e_denoised,
                   x_denoised1_ref$e_denoised,
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised1$v), abs(x_denoised1_ref$v), tolerance = 5e-2)
 expect_equivalent(abs(x_denoised1$u), abs(x_denoised1_ref$u), tolerance = 5e-2)
})


x_denoised2_ref <- readRDS("reference_data/x_denoised2.rds")

test_that("Hard default settings work as expected", {
 expect_is(x_denoised2, "subspace_denoised")
 expect_equivalent(x_denoised2$xi_denoised,
                   x_denoised2_ref$xi_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised2$x_denoised,
                   x_denoised2_ref$x_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised2$v_denoised,
                   x_denoised2_ref$v_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised2$e_denoised,
                   x_denoised2_ref$e_denoised,
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised2$v),
                   abs(x_denoised2_ref$v),
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised2$u),
                   abs(x_denoised2_ref$u),
                   tolerance = 5e-2)
 # expect_message(truncate(x1, components = 20,
 #  method = "hard", zeroout = FALSE),
 #  c("denoised", "hard", "lambda_min", "truncated"))
})


x_denoised3_ref <- readRDS("reference_data/x_denoised3.rds")

test_that("Identity default settings work as expected", {
 expect_is(x_denoised3, "subspace_denoised")
 expect_equivalent(x_denoised3$xi_denoised,
                   x_denoised3_ref$xi_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$x_denoised,
                   x_denoised3_ref$x_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$v_denoised,
                   x_denoised3_ref$v_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$e_denoised,
                   x_denoised3_ref$e_denoised,
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised3$v),
                   abs(x_denoised3_ref$v),
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised3$u),
                   abs(x_denoised3_ref$u),
                   tolerance = 5e-2)

 expect_is(x_denoised4, "subspace_denoised")
 expect_equivalent(x_denoised3$xi_denoised,
                   x_denoised4$xi_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$x_denoised,
                   x_denoised4$x_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$v_denoised,
                   x_denoised4$v_denoised,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$e_denoised,
                   x_denoised4$e_denoised,
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised3$v),
                   abs(x_denoised4$v),
                   tolerance = 5e-2)
 expect_equivalent(abs(x_denoised3$u),
                   abs(x_denoised4$u),
                   tolerance = 5e-2)
 # expect_message(truncate(s = subspace1_ref, components = 20,
 #  method = "identity", location = c(1:15), zeroout = FALSE),
 #  c("denoised", "identity", "location", "truncated"))
})
context("Argument input error")

test_that("Argument input error", {
 expect_error(
   truncate(x = NULL), "Invalid input")
 expect_error(
   truncate(), "missing")
 expect_message(
   truncate(x = x2, method = "threshold", alpha = 0.9), "full")
 expect_error(
   truncate(x = x2, method = "threshold"), "specified")
 expect_error(
   truncate(x = x2, method = "threshold", alpha = -1), "positive")
 expect_error(
   truncate(x = x2, method = "threshold", alpha = 2), "less")
  expect_error(
   truncate(x = x2, method = "threshold", components = 1, alpha = 0.1), "preserved")
  expect_error(
   truncate(x = x2, method = "identity"), "location")
  expect_error(
   truncate(x = x2, method = "identity", location = "alpha"), "Invalid")
  expect_error(
   truncate(x = x2, method = "identity", location = 200), "smaller")
  expect_error(
   truncate(x = x2, method = "identity", location = -1), "bounds")
  expect_error(
   truncate(x = x2, method = "location"), "Invalid method input")
 expect_message(
   truncate(s = subspace1_ref, method = "hard"), "already")
})
