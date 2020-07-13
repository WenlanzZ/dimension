library("dimension")

# Tests for dimension default settting
# ------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x1 <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
x2 <- x_sim(n = 150, p = 100, ncc = 10, var = 5)

subspace1_ref <- readRDS("reference_data/Subspace1.rds")

time_taken <- system.time({
  suppressWarnings(results1 <- dimension(x1))
  suppressWarnings(results2 <- dimension(x1, method = "double_posterior"))
  suppressWarnings(results3 <- dimension(x1, method = "posterior"))
  suppressWarnings(results4 <- dimension(x1, method = "kmeans"))
  suppressWarnings(results5 <- dimension(x1, method = "ladle"))

  suppressWarnings(results2_ <- subspace1_ref %>% estimate_rank_double_posterior())
  suppressWarnings(results3_ <- subspace1_ref %>% estimate_rank_posterior())
  suppressWarnings(results4_ <- subspace1_ref %>% estimate_rank_kmeans())
  suppressWarnings(results5_ <- x1 %>% estimate_rank_ladle())
})


test_that("Default settings work as expected", {
  expect_is(results1$subspace, "subspace")
  expect_equal(results1$subspace$ndf, 150)
  expect_equal(results1$subspace$pdim, 100)
  expect_equal(results1$subspace$components, subspace1_ref$components)
  expect_equal(results1$subspace$transpose_flag, TRUE)
  expect_equivalent(results1$subspace$irl$eigen,
                    subspace1_ref$irl$eigen,
                    tolerance = 5e-2)
  expect_equal(results1$subspace$irl$dim, 1:100)
  expect_equivalent(results1$subspace$mp_irl$eigen,
                    subspace1_ref$mp_irl$eigen,
                    tolerance = 5e-2)
  expect_equal(results1$subspace$mp_irl$dim, 1:100)
  expect_equivalent(results1$subspace$v, subspace1_ref$v, tolerance = 5e-2)
  expect_equivalent(results1$subspace$u, subspace1_ref$u, tolerance = 5e-2)
  # expect_equal(results1$dimension, 10)
})

test_that("estimate_rank_double_posterior settings work as expected", {
  # expect_equal(results2$subspace, subspace1_ref)
  expect_equal(results2_$subspace, subspace1_ref)
  expect_equal(results2$dimension, results1$dimension)
  expect_equal(results2$dimension, results2_$dimension)
})

test_that("estimate_rank_posterior settings work as expected", {
  # expect_equal(results3$subspace, subspace1_ref)
  expect_equal(results3_$subspace, subspace1_ref)
  expect_equal(results3$dimension, results3_$dimension)
})

test_that("estimate_rank_kmeans settings work as expected", {
  # expect_equal(results4$subspace, subspace1_ref)
  expect_equal(results4_$subspace, subspace1_ref)
  expect_equal(results4$dimension, results4_$dimension)
  expect_true(inherits(km_plot(results4$within_var), "ggplot"))
    expect_error(km_plot(results4$within_var, annotation = "1"), "numbers")
  expect_error(km_plot(results4$within_var, annotation = 110), "less")
})

test_that("estimate_rank_ladle settings work as expected", {
  expect_equal(results5$d, results5_$d)
  expect_equal(results5$d, 5)
})

context("Argument input error")

test_that("Argument input error", {
  expect_error(
    dimension(x = NULL), "Invalid")
  expect_error(
    dimension(), "missing")
  expect_message(
    dimension(x2, verbose = TRUE), "full")
  expect_message(
    dimension(s = subspace1_ref, verbose = TRUE), "already")
  expect_error(
    dimension(s = x1, verbose = TRUE), "type")
})

context("modified legacyplot check")

test_that("modified legacyplot check", {
 expect_true(inherits(modified_legacyplot(results3$bcp_irl, annotation = 10), "gtable"))
 expect_error(modified_legacyplot(results3$bcp_irl, annotation = "0"), "numbers")
 expect_error(modified_legacyplot(results3$bcp_irl, annotation = 110), "less")
})


# test_that("print.dimension prints by default", {
#   expect_message(print(results1),"features")
# })
