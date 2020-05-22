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
library(trellisopejs)
library(purrr)
library(gtsummary)

mse <- function(x, x_hat) {
	mean((x-x_hat)^2)
}
mae <- function(x, x_hat) {
	mean(abs(x-x_hat))
}

library(devtools)
document()

two_one_one <- 
  expand.grid(n = c(10, 100, 1000, 10000),
              p = c(10, 100, 1000),
              d = c(3, 10),
              sigma = c(2, 6, 10)) %>% as_tibble()
two_one_one$num_sim <- 10

ncores <- detectCores()
registerDoParallel(ncores - 2)
registerDoRNG()

two_one_one$runs <- foreach(i = seq_len(nrow(two_one_one)), 
                            .options.RNG=1234) %dorng% {
  foreach(s = seq_len(two_one_one$num_sim[i]), 
          .combine = bind_rows,
          .errorhandling = "remove") %do% { 
    M <- diag(c(2, 1, 1, rep(0, two_one_one$p[i] - 3)))
    # What is 0.54^2?
    sigma <- diag(rep(0.54^2, two_one_one$p[i])) + M
    cor_cols <- rmvnorm(two_one_one$n[i], rep(0, two_one_one$p[i]), 
                        sigma = sigma)
    stime2 <- 
      suppressMessages(
        system.time(test2 <- dimension(cor_cols, verbose = FALSE)))
    tibble(dim_est = test2$dimension, time = stime2[[3]])
  }
  cat(i, "of", nrow(two_one_one), "\n")
}

# Write out two_one_one here. We can post process anything we on a smaller
# machine.

est_hist <- function(run) {
  ggplot(run, aes(x = dim_est)) + geom_histogram() + theme_bw()
}

two_one_one <- two_one_one %>%
  mutate(mse = map2_dbl(runs, d, 
                        ~ mean(.x$dim_est - rep(d, length(.x$dim_est))^2)),
         mae = map2_dbl(runs, d, 
                        ~ mean(abs((.x$dim_est - rep(d, length(.x$dim_est)))))),
         bias = map2_dbl(runs, d, 
                             ~ mean(.x$dim_est - rep(d, length(.x$dim_est)))),
         mean_est = map_dbl(runs, ~ mean(.x$dim_est)),
         median_est = map_dbl(runs, ~ median(.x$dim_est)),
         plots = map(runs, est_hist))


two_one_one %>% trelliscope(name = "Dimension Estimation", panel_col = "plots")

# Points from MK.
# Get this running on a single machine first.
# From the code, it's not clear that this is using the ladle estimator.
# Where does the 0.54^2 come from?


# #####------------------------> Setting in Ladle small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 100
# p <- 10
# d <- 3
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   M <- diag(c(2, 1, 1, rep(0, p-3)))
#   sigma <- diag(rep(0.54^2, p)) + M
#   cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
#   stime2 <- system.time(test2 <- dimension(cor_cols))
#   cat("dimension = ", test2$dimension, "time = ", stime2[[3]], "\n")
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_",n,"_p_",p,"_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)

# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_3.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()
# plot(test2$Subspace)
# ggsave(file.path(output_dir, "mkdim_scree_plot_est_est_3.pdf"))
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_3.RDS")))
within_var <- test2$within_var
plot(1:9, within_var[,1], "l")
ggplot() +
        theme_minimal() +
        geom_point(aes(x = 1:(rnk-1), y = within_var[,1])) +
        geom_line(aes(x = 1:(rnk-1), y = within_var[,1]), color = "black",cex = 0.5) +
        geom_text_repel(aes(x = 1:(rnk-1), y = within_var[,1],
                            label = 1:(rnk-1)), colour = "black", size = 5)
#ggsave(file.path(output_dir, paste0("mkdim_kmeans_customized_clust_est_", test2$dimension, ".pdf")))
##may 20 
#  3   4   5 
#988  10   2 
#####------------------------> Setting in Ladle small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 1000
# p <- 10
# d <- 3
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   M <- diag(c(2, 1, 1, rep(0, p-3)))
#   sigma <- diag(rep(0.54^2, p)) + M
#   cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
#   stime2 <- system.time(test2 <- dimension(cor_cols))

#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_",n,"_p_",p,"_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mae(accuracy_test_mkdim$mkdim_d, d)
# mse(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)

# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_3.RDS")))
#####------------------------> Setting in Ladle small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 10000
# p <- 10
# d <- 3
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   M <- diag(c(2, 1, 1, rep(0, p-3)))
#   sigma <- diag(rep(0.54^2, p)) + M
#   cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
#   stime2 <- system.time(test2 <- dimension(cor_cols))

#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_",n,"_p_",p,"_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_3.RDS")))
#####------------------------> Setting in Ladle large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 10
# p <- 100
# d <- 3
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   M <- diag(c(2, 1, 1, rep(0, p-3)))
#   sigma <- diag(rep(0.54^2, p)) + M
#   cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
#   stime2 <- system.time(test2 <- dimension(cor_cols))
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim, 
# 	file.path(output_dir, 
#  	paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_1000run_2_1_1.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_1.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_1.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()

# double posterior
  1   2   3   4   5   6   7 
234 240 210 177  91  36  12

## irlmax== postmax
  1   2   3   4   5   6   7   8   9 
195 127  67  65  79 111 163 137  56 

##prob_irl[post_max] > 0.9 * max(prob_irl)
  1   2   3   4   5   6   7   8   9 
210 104  74  67  75 139 142 153  36 
## prior posterior
  1   2   3   4   5   6   7   8 
 18 267 270 247 133  48  16   1 

   2   3   4   5   6   7   8   9 
378 257 153  92  68  44   6   2 
## kmeans
  1   2   3   4   5   6   7   8 
276 239 182 110 109  70   9   5 
## kmeans may 13
  1   2   3   4   5   6   7   8   9 
175 265 190 142  93  71  19   9   1 
## kmeans may 20
  1   2   3   4   5   6   7   8 
214 188 138 143 124  83  64  46 
###############
# wrong cases #
###############
plot(test2$Subspace, changepoint = test2$dimension, annotation = 10)
# ggsave(file.path(output_dir, "mkdim_scree_plot_est_3.pdf"))

# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_3.pdf"))
modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()

# pdf(file.path(output_dir, "mkdim_bcp_posterior_prob_plot_est_3.pdf"))
modified_legacyplot(test2$bcp_post, annotation = 10)
# dev.off()

#Bayesian Change Point
test2$Changepoint$bcp_irl
test2$Changepoint$bcp_post

#####------------------------> Setting in Ladle large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 10
# p <- 1000
# d <- 3
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   M <- diag(c(2, 1, 1, rep(0, p-3)))
#   sigma <- diag(rep(0.54^2, p)) + M
#   cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
#   stime2 <- system.time(test2 <- dimension(cor_cols))
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim, 
# 	file.path(output_dir, 
#  	paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_1000run_2_1_1.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_7.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_7.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()
## prior posterior
  1   2   3   4   5   6   7   8 
  3  51 138 286 286 175  51  10 

    2   3   4   5   6   7   8   9 
161 178 149 144 134 138  71  25

1   2   3   4   5   6   7   8   9 
 48 277 213 151 119 118  68   4   2
###kmeans
   1   2   3   4   5   6   7   8   9 
 75 160 153 134 124 167 113  55  19 
###kmeans may 20
   1   2   3   4   5   6   7   8 
 63 134 134 126 154 148 103 138
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 7.558
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 2.24
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.07237
#####------------------------> x_sim large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 10
# p <- 100
# d <- 3
# sigma <- 6
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
  
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_3.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_3.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()
  2   3   4   
408 283 127 

###kmeans
  1   2   3   4   5   6   7   8 
241 379 198  83  58  32   6   3
###kmeans may 20
  1   2   3   4   5   6   7   8 
165 231 209 155 110  61  46  23 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 3.346
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 1.418
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.070174
###kmeans may 20 sigma 10
  1   2   3   4   5   6   7   8 
170 314 272 145  54  24  15   6 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 1.961
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 1.069
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.073625

#####------------------------> x_sim large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 10
# p <- 1000
# d <- 3
# sigma <- 6
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
#   # plot(test2$Subspace, changepoint = test2$dimension, annotation = 10)
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_6.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_6.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()
###kmeans may 20
  1   2   3   4   5   6   7   8   9 
 91 144 121 144 136 125  95 143   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 7.452
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 2.218
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.072594
#####------------------------> x_sim large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 10
# p <- 10000
# d <- 3
# sigma <- 6
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
  
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_6.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_6.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 10)
# dev.off()
## kmeans
 1   2   3   4   5   6   7   8   9 
 29 200 232 152 110  86  84  71  35 
###kmeans may 20
   1   2   3   4   5   6   7   8   9 
 56 134 112 141 140 142 102 171   2 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 8.316
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 2.368
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.073928

#####------------------------> x_sim small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 500
# p <- 100
# d <- 10
# sigma <- 2
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
  
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)

# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_10.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_10.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()
# plot(test2$Subspace, changepoint = test2$dimension, annotation = 10)
# ggsave(file.path(output_dir, "mkdim_scree_plot_est_est_10.pdf"))


# 10  11  12  13  14  15  16  17  18  19  20  25 
#434  99 180 132  54  48  32   8   8   3   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 6.695
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 1.675
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.617829
#####------------------------> x_sim small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 5000
# p <- 100
# d <- 10
# sigma <- 2
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
  
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_10.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_10.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()
mkdim_d mkdim_user mkdim_system mkdim_elapsed
    <dbl>      <dbl>        <dbl>         <dbl>
1      10      0.499         0.09          9.56

10  11  12  13  14  15  16  17  18  20  22 
664  68 124  90  30  13   3   4   2   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 2.855
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 0.855
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.760581
#####------------------------> x_sim small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 50000
# p <- 100
# d <- 10
# sigma <- 2
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
  
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_10.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_10.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()
mkdim_d mkdim_user mkdim_system mkdim_elapsed
    <dbl>      <dbl>        <dbl>         <dbl>
1      10       1.81        0.218          51.4

 10  11  12  13  14  15 
835  26  64  50  22   3 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 1.159
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 0.407
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 4.036847
#####------------------------> ladle small p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 500
# p <- 100
# d <- 10
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   M <- diag(c(rep(2,10), rep(0, p-10)))
#   sigma <- diag(rep(0.54^2, p)) + M
#   cor_cols <- mvtnorm::rmvnorm(n, rep(0, p), sigma = sigma)
#   stime2 <- system.time(test2 <- dimension(cor_cols))

#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim, 
# 	file.path(output_dir, 
#  	paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_1000run_2_2_2.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_10.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_10.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 50)
# dev.off()
  mkdim_d mkdim_user mkdim_system mkdim_elapsed
    <dbl>      <dbl>        <dbl>         <dbl>
1      10      0.287        0.007         0.293

 10  11  12  13  14  15  16 
749  84 111  39  12   2   3 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 1.229
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 0.499
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.544973
#####------------------------> x_sim large p<------------------------#########
ncores <- detectCores()
registerDoParallel(ncores)
num_sim <- 1000
n <- 100
p <- 500
d <- 10
sigma <- 6
output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
if (file.exists(output_dir) == F) {
  dir.create(output_dir)    
}
accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
  x <- x_sim(n = n, p = p, ncc = d, var = sigma)
  stime2 <- system.time(test2 <- dimension(x))
  # tibble(post_prob = c(test2$bcp_irl$posterior.prob[-n], 0))
  tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
}
write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run_may_11_new.csv")))
table(accuracy_test_mkdim[,1])
mse(accuracy_test_mkdim$mkdim_d, d)
mae(accuracy_test_mkdim$mkdim_d, d)
mean(accuracy_test_mkdim$mkdim_elapsed)
### original
  2   3   4   5   6   7   8   9  10  11  12  13 
123 158 113 110 125  99 107  94  64   4   2   1

### threshold 0.90
2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  23  24 
6   6  30  48  59  85 163 209 255  53  40  18   8   2   6   3   5   2   1   1 

mse(accuracy_test_mkdim$mkdim_d, d)#7.291
mae(accuracy_test_mkdim$mkdim_d, d)#1.907

### kmenas cluster May 7
1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
  4  21  17  18  14  28  14  43 113 146 126 108 121  82  49  40  19   7  19   3 
 21  22  23  27 
  4   1   2   1 
mse(accuracy_test_mkdim$mkdim_d, d)#14.455
mae(accuracy_test_mkdim$mkdim_d, d)#2.861

### kmenas cluster May 11
 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
  1   2   5   5   3   8   9  35 120 207 122 150 136  78  48  39  12   6   6   6 
 21  22 
  1   1 
mse(accuracy_test_mkdim$mkdim_d, d)#9.299
mae(accuracy_test_mkdim$mkdim_d, d)#2.239
  ### kmenas cluster May 11 new
  1   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21 
  3   1   6   7   7   7  28 104 224 112 119 103  77  54  40  23  15  14  11  10 
 22  23  24  25  26  27  28  29  33  34  41  49  54  77 
  6   7   3   2   3   3   3   2   1   1   1   1   1   1 
mse(accuracy_test_mkdim$mkdim_d, d)#14.654
mae(accuracy_test_mkdim$mkdim_d, d)#2.672
  ### kmenas cluster May 13
  2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  22 
 11  16  27  30  51  62  95 171 224  95  64  50  31  27  20   7   8   4   2   2 
 24  26  32 
  1   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 10.43
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 2.26
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.365717
  ### kmenas cluster May 20
> table(accuracy_test_mkdim[,1])

  6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  22 
  8   6  31 123 266 170 134  93  65  38  22  27   7   7   2   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 7.436
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 1.874
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.617527


saveRDS(x, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_10.RDS")))
pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_10.pdf"))
modified_legacyplot(test2$bcp_irl, annotation = 50)
dev.off()

plot(test2$Subspace, changepoint = test2$dimension, annotation = 20)
ggsave(file.path(output_dir, "mkdim_scree_plot_est_est_1.pdf"))


var_sim <- foreach(1:1000, .combine = cbind) %dopar% {
  x <- x_sim(n = n, p = p, ncc = d, var = sigma)
  test2 <- dimension(x, threshold = 0.90)
  tibble(post_prob = c(test2$bcp_irl$posterior.prob[-n], 0))
  # tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
}
saveRDS(var_sim, file.path(output_dir,paste0("var_n_",n,"_p_",p,"_d_",d,"_est_all.RDS")))
# var_sim <- readRDS(file.path(output_dir, "var_n_100_p_500_d_10_est_9.RDS"))
apply(var_sim, 1, mean)
apply(var_sim, 1, sd)

# x change
# [1] 0.35879714 0.39209538 0.39031203 0.39576823 0.39548818 0.39790279
# [7] 0.39575521 0.38372846 0.38347799 0.36654979 0.28984743 0.26500669
# [13] 0.23426108 0.21824158 0.19580972 0.19054933 0.18608901 0.17212432
# [19] 0.16370118 0.15058173 0.14173178 0.14815341 0.14219765 0.13931335

# x fixed
# [1] 0.0012642721 0.0078454516 0.0298803743 0.0234185112 0.0003938928
# [6] 0.0056639922 0.0031538977 0.0075825424 0.0067569642 0.0426410908
# [11] 0.0607653676 0.0592031224 0.0502043461 0.0486981011 0.0285896358
# [16] 0.0516985200 0.0626299667 0.0400373563 0.0456028353 0.0365258785
# [21] 0.0550560559 0.0548410928 0.0486556692 0.0450780468 0.0309904580

prob_mean <-apply(var_sim, 1, mean)
sd_mean <- apply(var_sim, 1, sd)
error <- qnorm(0.975)*sd_mean/sqrt(n)
mean(error)
# [1] 0.02357066  # x change
#[1] 0.006932737	  # x fixed
left <- prob_mean-error
right <- prob_mean+error
left
right

plot(1:n, prob_mean, "l")
lines(1:n, left, col = "red")
lines(1:n, right, col = "red")

data <- tibble(x = 1:n, error = error, prob_mean = prob_mean, left = left, right = right)
ggplot() +
          theme_minimal() +
          geom_point(aes(x = data$x, y = data$error, colour = "red")) +
          geom_line(aes(x = data$x, y = data$error), color = "black",cex = 0.5) +
          geom_text_repel(aes(x = data$x, y = data$error,
                              label = data$x), colour = "black", size = 5)
ggsave(file.path(output_dir, "error_est_9.pdf"))

ggplot() +
          theme_minimal() +
          geom_point(aes(x = data$x, y = data$prob_mean, colour = "red")) +
          geom_line(aes(x = data$x, y = data$prob_mean), color = "black",cex = 0.5) +
          geom_line(aes(x = data$x, y = data$left), color = "green",cex = 0.5) +
          geom_line(aes(x = data$x, y = data$right), color = "green",cex = 0.5) +
          geom_text_repel(aes(x = data$x, y = data$prob_mean,
                              label = data$x), colour = "black", size = 5)
ggsave(file.path(output_dir, "confidence_interval_est_9.pdf"))

#####------------------------> x_sim large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 100
# p <- 5000
# d <- 10
# sigma <- 6
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))

#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run_MAY11.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_2.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_2.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()

### kmenas cluster May 7
  2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21 
235 168 142  84  64  53  42  24  12  32  23  11  10   6   9  11   5   7   7   2 
 22  23  24  25  26  27  28  29  30  31  32  33  34  36  37  40  41  43  44  45 
  4   7   3   2   2   1   1   3   2   4   1   4   1   2   2   1   1   2   1   1 
 47  48  56  60  93  95 
  1   2   1   2   1   1 
### kmenas cluster May 11
    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  21 
 35  94  96  97 110 114 107  77  71  54  33  34  21  23  14   6   5   2   3   2 
 22  24 
  1   1 
mse(accuracy_test_mkdim$mkdim_d, d)#25.453
mae(accuracy_test_mkdim$mkdim_d, d)#4.361

### kmenas cluster May 13
 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
 2  2 90 90 65 64 69 61 67 72 50 59 39 35 29 30 24 20 19 12 12 10  5  5  2  5 
26 27 28 29 30 32 33 34 36 37 38 39 40 41 42 45 46 47 52 61 63 69 71 79 83 87 
 3  6  2  4  4  3  3  2  1  2  1  1  1  2  1  1  1  1  1  2  1  1  1  1  1  1 
89 90 92 93 94 95 96 97 
 1  1  1  1  2  3  4  1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 181.803
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 6.881
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.913799
### kmenas cluster May 20
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
 7 62 77 77 81 73 56 59 64 38 41 35 40 21 32 18 19 15 10  3  8  9 11  7  8  7 
27 28 29 30 31 32 33 34 35 36 37 39 42 43 47 48 50 56 59 68 72 77 83 84 86 87 
 7  6  4  3  5  6  4  3  4  2  1  1  1  3  2  2  1  1  1  1  1  1  1  1  2  1 
88 90 91 92 93 94 95 96 97 98 
 1  2  1  9  5  7  9  8 12  3 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 514.544
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 11.004
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.853105

  #####------------------------> x_sim large p<------------------------#########
ncores <- detectCores()
registerDoParallel(ncores)
num_sim <- 1000
n <- 100
p <- 50000
d <- 10
sigma <- 6
output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
if (file.exists(output_dir) == F) {
  dir.create(output_dir)    
}
accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
  x <- x_sim(n = n, p = p, ncc = d, var = sigma)
  stime2 <- system.time(test2 <- dimension(x))
  
  tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
}
write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
table(accuracy_test_mkdim[,1])
mse(accuracy_test_mkdim$mkdim_d, d)
mae(accuracy_test_mkdim$mkdim_d, d)
mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_3.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_3.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()
# tmp <- tibble(bcp_irl_new = test2$Changepoint$bcp_irl$posterior.prob, bcp_irl = test2$Changepoint$prob_irl, bcp_prior = test2$Changepoint$prob_post)
 2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21 
129 190 130  88  46  60  25  40  21  28  15  23  18  12  16  15   8   2  12  10 
 22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41 
 10  11   4   4   3   5   6   3   1   7   4   3   3   2   2   1   1   3   2   1 
 42  45  46  51  55  56  57  64  65  66  67  68  78  84  88  93  94  95  97  98 
  1   2   1   1   2   3   2   1   1   1   1   2   1   2   1   2   1   3   5   3 


> table(accuracy_test_mkdim[,1])

  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
  2   4  77  91 108  83  54  59  71  41  43  43  45  32  26  22  19  14  14  17 
 20  21  22  23  24  25  26  27  28  29  30  31  33  34  35  37  38  40  41  43 
 11  12   6   4   6   7   4   2   2   4   6   5   2   1   4   2   2   3   1   1 
 44  45  47  48  52  54  56  57  59  60  63  65  67  69  73  77  82  83  85  86 
  1   1   1   1   1   1   1   2   1   1   3   1   1   1   2   1   1   1   3   2 
 87  89  91  93  94  95  96  97  98 
  1   1   2   4   5   1   5   3   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 302.428
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 8.74
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 6.177211
## may 20
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
13 54 47 40 50 35 30 22 30 20 23  9 12 15  9 11  8 12 13  8 10  7  6  6  6  4 
27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 48 49 50 51 52 53 
 2  5  4  4  1  5  2  4  1  4  6  3  3  1  4  2  1  6  1  2  1  1  1  1  1  4 
54 55 56 58 59 60 61 62 63 64 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
 2  1  2  1  2  1  6  4  2  3  3  3  2  3  6  1  2  3  5  2  1  5  3  4  5  9 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 
 3  7  5  9 11 13 13 14 20 21 23 30 38 40 39 61  2
 > mse(accuracy_test_mkdim$mkdim_d, d)
[1] 2767.448
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 38.452
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 4.263927
detect.spikes(tmp)
#####------------------------> x_sim large p<------------------------#########
ncores <- detectCores()
registerDoParallel(ncores)
num_sim <- 1000
n <- 100
p <- 500
d <- 10
sigma <- 10
output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
if (file.exists(output_dir) == F) {
  dir.create(output_dir)    
}
accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
  x <- x_sim(n = n, p = p, ncc = d, var = sigma)
  stime2 <- system.time(test2 <- dimension(x, p = 0.95))
  
  tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
}
write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
table(accuracy_test_mkdim[,1])
mse(accuracy_test_mkdim$mkdim_d, d)
mae(accuracy_test_mkdim$mkdim_d, d)
mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_10.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_10.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()
# modified_legacyplot(test2$Changepoint$bcp_irl, annotation = 10)
# plot(test2$Subspace,  annotation = 10)

# head(tibble(bcp_irl_new = test2$Changepoint$bcp_irl$posterior.prob, bcp_irl = test2$Changepoint$prob_irl, bcp_prior = test2$Changepoint$prob_post),15)

# var_sim <- accuracy_test_mkdim
# prob_mean <-apply(var_sim, 1, mean)
# sd_mean <- apply(var_sim, 1, sd)
# error <- qnorm(0.975)*sd_mean/sqrt(n)
# # > mean(error)
# # [1] 0.02357066
# left <- prob_mean-error
# right <- prob_mean+error
# left
# right

# plot(1:n, prob_mean, "l")
# lines(1:n, left, col = "red")
# lines(1:n, right, col = "red")
  3   4   5   6   7   8   9  10  11  12 
  3   1   5   7   9  24  87 859   2   3 
### may 13 kmeans overestimate a little bit
    8   9  10  11  12  13  14  15  16  17  18  19  20  21  25  30 
  3  31 660  43  78  68  38  20  19  15  10   6   4   3   1   1
  > mse(accuracy_test_mkdim$mkdim_d, d)
[1] 6.051
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 1.153
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.37478
## may 20
  7   9  10  11  12  13  14  15  16  17  18  19  23 
  1   3 501 125 173 109  36  29   7  11   3   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 4.344
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 1.258
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.567025

# #####------------------------> x_sim large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 100
# p <- 5000
# d <- 10
# sigma <- 10
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
  
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_3.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_3.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()
# plot(test2$Subspace, changepoint = test2$dimension, annotation = 10)
# ggsave(file.path(output_dir, "mkdim_scree_plot.pdf"))
 2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21 
295 233 141 109  89  43  34   9  10   8   4   3   3   3   1   3   1   2   1   1 
 27  28  31  37  47  51  93 
  1   1   1   1   1   1   1 
  ### may 13 kmeans 
1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
  2 160 119 100  90  71  77  75  41  51  39  24  27  24  19  15  15   4   6   8 
 21  22  23  24  25  28  31  32  34  35  38  41  53  86  87  89  92  95  96 
  3   1   2   3   4   1   1   3   2   1   1   1   4   1   1   1   1   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 83.066
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 5.674
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.80589
  ### may 20 kmeans 
  2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21 
  6  30  47  76 102 115  91  90  85  67  61  48  42  28  23  17  12  14   9  10 
 22  23  24  26  28  31  32  35  95  96  97 
  7   7   2   2   2   1   1   1   2   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 51.063
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 3.947
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 0.793513

#####------------------------> x_sim large p<------------------------#########
# ncores <- detectCores()
# registerDoParallel(ncores)
# num_sim <- 1000
# n <- 100
# p <- 50000
# d <- 10
# sigma <- 100
# output_dir <- file.path(output_path, paste0("n_", n, "_p_", p, "_d_", d))
# if (file.exists(output_dir) == F) {
#   dir.create(output_dir)    
# }
# accuracy_test_mkdim <- foreach(1:num_sim, .combine = rbind) %dopar% {
#   x <- x_sim(n = n, p = p, ncc = d, var = sigma)
#   stime2 <- system.time(test2 <- dimension(x))
#   # tibble(post_prob = c(test2$bcp_irl$posterior.prob[-n], 0))
#   # cat("dimension = ", test2$dimension, "time = ", stime2[[3]], "\n")
#   tibble(mkdim_d = test2$dimension, mkdim_user = stime2[[1]], mkdim_system = stime2[[2]], mkdim_elapsed = stime2[[3]])
# }
# write.csv(accuracy_test_mkdim,file.path(output_dir,paste0("accuracy_test_mkdim_n_", n, "_p_", p, "_var_", sigma, "_1000run.csv")))
# table(accuracy_test_mkdim[,1])
# mse(accuracy_test_mkdim$mkdim_d, d)
# mae(accuracy_test_mkdim$mkdim_d, d)
# mean(accuracy_test_mkdim$mkdim_elapsed)
# saveRDS(cor_cols, file.path(output_dir,paste0("X_n_",n,"_p_",p,"_d_",d,"_est_9.RDS")))
# pdf(file.path(output_dir, "mkdim_bcp_posterior_mean_plot_est_9.pdf"))
# modified_legacyplot(test2$bcp_irl, annotation = 20)
# dev.off()

modified_legacyplot(test2$bcp_irl, annotation = 10)
plot(test2$Subspace,  annotation = 10)

  2   3   4   5   6   7   8   9  10 
  3   2   5   4   7   7   9  52 911 

  0   9  10  11  12  13  14  15  16  20  96 
  1  18 957   6   3   7   2   3   1   1   1 

> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 7.838
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 0.186
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 6.974127
  ### may 20 kmeans 
7  10  11  12  13  14  15  16  17  18  19  20 
  3 650 136 117  52  19   8   9   3   1   1   1 
> mse(accuracy_test_mkdim$mkdim_d, d)
[1] 2.319
> mae(accuracy_test_mkdim$mkdim_d, d)
[1] 0.753
> mean(accuracy_test_mkdim$mkdim_elapsed)
[1] 4.621042


