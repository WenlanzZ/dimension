library("dimension")

# Tests for truncate default settting
# --------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))

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
 suppressWarnings(x_denoised4 <- truncate(subspace_ = subspace1_ref,
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
 expect_equivalent(x_denoised1$v, x_denoised1_ref$v, tolerance = 5e-2)
 expect_equivalent(x_denoised1$u, x_denoised1_ref$u, tolerance = 5e-2)
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
 expect_equivalent(x_denoised2$v,
                   x_denoised2_ref$v,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised2$u,
                   x_denoised2_ref$u,
                   tolerance = 5e-2)
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
 expect_equivalent(x_denoised3$v,
                   x_denoised3_ref$v,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$u,
                   x_denoised3_ref$u,
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
 expect_equivalent(x_denoised3$v,
                   x_denoised4$v,
                   tolerance = 5e-2)
 expect_equivalent(x_denoised3$u,
                   x_denoised4$u,
                   tolerance = 5e-2)
})
context("Argument input error")

test_that("Argument input error", {
 expect_error(
   truncate(), "missing")
 expect_message(
   truncate(subspace_ = subspace1_ref, method = "hard"), "already")

 subspace1_ref$mp_irl <- NULL
 expect_error(
   dimension(subspace_ = subspace1_ref, method = "hard"), "missing")
})
