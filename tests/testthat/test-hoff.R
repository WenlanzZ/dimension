library("dimension")

# Tests for hoff default settting
# ------------------------------------------
context("Default settings work as expected")

set.seed(seed = 1234)

x <- x_sim(n = 150, p = 100, ncc = 10, var = 6)

time_taken <- system.time({
  suppressWarnings(hoff_result <- hoff(y = x, NSCAN = 10))
})
