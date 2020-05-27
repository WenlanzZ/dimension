library(tibble)
library(dplyr)
library(tidyr)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
library(irlba)
library(ggplot2)
library(mvtnorm)
library(bcp)
library(cpm)
library(gridExtra)
library(ggrepel)
library(dr)
library(trelliscopejs)
library(purrr)
library(gtsummary)

mse <- function(x, x_hat) {
	mean((x - x_hat)^2)
}
mae <- function(x, x_hat) {
	mean(abs(x - x_hat))
}

setwd("/Users/wz262/Projects/dimension")
library(devtools)
document()

################--------------------> model 1 <--------------------------#################
bench_params <- 
  expand.grid(n = c(10, 100, 1000, 10000),
              p = c(10, 100, 1000),
              d = 3) %>% as_tibble()
#drop setting that ladle cannot run through
bench_params <- bench_params %>% slice(-c(1, 6, 7, 8, 10, 11, 12))

rank_ests <- tibble(
  rank_estimator_type = c("double_post", "double_post_cor", "kmeans", "kmeans_cor", 
                          "posterior", "posterior_cor", "ladle"),
  rank_estimator = list(estimate_rank_double_posterior,
                        estimate_rank_double_posterior_cor,
                        estimate_rank_kmeans,
                        estimate_rank_kmeans_cor,
                        estimate_rank_posterior,
                        estimate_rank_posterior_cor,
                        estimate_rank_ladle))
cnames <- c("dp_est", "dp_time", "dpc_est", "dpc_time","km_est", "km_time", "km_cor_est",
           "km_cor_time", "p_est", "p_time", "pc_est", "pc_time","lad_est", "lad_time")

bench_params$num_sim <- 1000
bench_params %>% print(n = Inf)

ncores <- detectCores()
registerDoParallel(ncores - 2)
registerDoRNG()

bench_params$runs <- foreach(i = seq_len(nrow(bench_params)), 
                            .options.RNG=1234) %dorng% {
  #run over all model setting
  foreach(s = seq_len(bench_params$num_sim[i]), 
          .combine = bind_rows,
          .errorhandling = "remove") %do% {
    #same model setting for d = 3
    M <- diag(c(2, 1, 1, rep(0, bench_params$p[i] - 3)))
    sigma <- diag(rep(0.54^2, bench_params$p[i])) + M
    cor_cols <- rmvnorm(bench_params$n[i], rep(0, bench_params$p[i]), 
                        sigma = sigma)
    #run all methods on the same matrix
    bench <- foreach(h = seq_len(nrow(rank_ests)-1), .combine = cbind) %do% {
      stime1 <- 
      suppressMessages(
        system.time(
          test1 <- cor_cols %>% create_subspace() %>% 
                    correct_eigenvalues() %>%
                    rank_ests$rank_estimator[[h]]()
      ))
      tibble(dim_est = test1$dimension, time = stime1[[3]])
    }
    #run ladle on the same matrix
    stime2 <- system.time(
          test2 <- cor_cols %>% rank_ests$rank_estimator[[nrow(rank_ests)]]()
    )
    bench <- bench %>% add_column(lad_est = test2$d, lad_time = stime2[[3]])
    colnames(bench) <- cnames
    cat("<-----------------------finish", s, "----------------------->\n")
    bench
  }
}

# Write out bench here. We can post process anything we on a smaller
# machine.

est_hist <- function(run) {
  ggplot(run, aes(x = dim_est)) + geom_histogram() + theme_bw()
}

# distribute results
bench<- foreach(i = seq_len(nrow(rank_ests)), .combine = bind_rows) %do% {
  bench <- bench_params
  bench$rank_estimator_type <- rank_ests$rank_estimator_type[i]
  bench$runs <- bench$runs %>% lapply(`[`, c(2 * i - 1, 2 * i)) %>% lapply(setNames, c("dim_est", "time"))
  bench
}

bench <- arrange(bench, -desc(n), -desc(p))
bench %>% print(n = Inf)
# calcualte metrics
bench <- bench %>%
  mutate(mse = map2_dbl(runs, d, 
                        ~ mean((.x$dim_est - rep(d, length(.x$dim_est)))^2)),
         mae = map2_dbl(runs, d, 
                        ~ mean(abs((.x$dim_est - rep(d, length(.x$dim_est)))))),
         bias = map2_dbl(runs, d, 
                             ~ mean(.x$dim_est - rep(d, length(.x$dim_est)))),
         mean_est = map_dbl(runs, ~ mean(.x$dim_est)),
         median_est = map_dbl(runs, ~ median(.x$dim_est)),
         plots = map(runs, est_hist))
bench %>% print(n = Inf)

bench_params_2_2_1<- bench_params
saveRDS(bench_params_2_2_1, "/Users/wz262/Projects/dimension/inst/benchmark/bench_params_2_2_1.rds")
bench_2_2_1<- bench
saveRDS(bench_2_2_1, "/Users/wz262/Projects/dimension/inst/benchmark/bench_2_2_1.rds")
filter(bench_params, n == 100, p == 100)$runs

bench %>% trelliscope(name = "Dimension Estimation", panel_col = "plots")

# Points from MK.
# Get this running on a single machine first.
# From the code, it's not clear that this is using the ladle estimator.
# Where does the 0.54^2 come from?

################--------------------> model 2 <--------------------------#################
bench_params <- 
  expand.grid(n = 10,
              p = c(10, 100, 1000),
              d = 3,
              sigma = c(2, 6)) %>% as_tibble()

bench_params <- 
expand.grid(n = c(100, 500, 5000, 50000),
            p = c(100, 500, 5000, 50000),
            d = 10,
            sigma = c(2, 6, 10)) %>% as_tibble()
#drop setting that ladle cannot run through
bench_params <- bench_params %>% slice(-c(1, 6, 7, 8, 10, 11, 12))


bench_params$num_sim <- 10
bench_params %>% print(n = Inf)

ncores <- detectCores()
registerDoParallel(ncores - 2)
registerDoRNG()

bench_params$runs <- foreach(i = seq_len(nrow(bench_params)), 
                            .options.RNG=1234) %dorng% {
  #run over all model setting
  foreach(s = seq_len(bench_params$num_sim[i]), 
          .combine = bind_rows,
          .errorhandling = "remove") %do% {
    #same model setting
    x <- x_sim(n = bench_params$n[i], p = bench_params$p[i], ncc = bench_params$d[i], var = bench_params$sigma[i])
    #run all methods on the same matrix
    bench <- foreach(h = seq_len(nrow(rank_ests)-1), .combine = cbind) %do% {
      stime1 <- 
      suppressMessages(
        system.time(
          test1 <- x %>% create_subspace() %>% 
                    correct_eigenvalues() %>%
                    rank_ests$rank_estimator[[h]]()
      ))
      tibble(dim_est = test1$dimension, time = stime1[[3]])
    }
    #run ladle on the same matrix
    stime2 <- system.time(
          test2 <- x %>% rank_ests$rank_estimator[[nrow(rank_ests)]]()
    )
    bench <- bench %>% add_column(lad_est = test2$d, lad_time = stime2[[3]])
    colnames(bench) <- cnames
    cat("<-----------------------finish", s, "----------------------->\n")
    bench
  }
}


# distribute results
bench<- foreach(i = seq_len(nrow(rank_ests)), .combine = bind_rows) %do% {
  bench <- bench_params
  bench$rank_estimator_type <- rank_ests$rank_estimator_type[i]
  bench$runs <- bench$runs %>% lapply(`[`, c(2 * i - 1, 2 * i)) %>% lapply(setNames, c("dim_est", "time"))
  bench
}

bench <- arrange(bench, -desc(n), -desc(p))
bench %>% print(n = Inf)
# calcualte metrics
bench <- bench %>%
  mutate(mse = map2_dbl(runs, d, 
                        ~ mean((.x$dim_est - rep(d, length(.x$dim_est)))^2)),
         mae = map2_dbl(runs, d, 
                        ~ mean(abs((.x$dim_est - rep(d, length(.x$dim_est)))))),
         bias = map2_dbl(runs, d, 
                             ~ mean(.x$dim_est - rep(d, length(.x$dim_est)))),
         mean_est = map_dbl(runs, ~ mean(.x$dim_est)),
         median_est = map_dbl(runs, ~ median(.x$dim_est)),
         plots = map(runs, est_hist))
bench %>% print(n = Inf)

bench_params_sigma<- bench_params
saveRDS(bench_params_sigma, "/Users/wz262/Projects/dimension/inst/benchmark/bench_params_sigma.rds")
bench_sigma<- bench
saveRDS(bench_sigma, "/Users/wz262/Projects/dimension/inst/benchmark/bench_sigma.rds")

bench %>% trelliscope(name = "Dimension Estimation", panel_col = "plots")
