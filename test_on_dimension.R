rm(list=lsf.str())
gc()

#auto generating documents after changing fxns
setwd("/Users/wz262/Projects/dimension")
library(roxygen2)
library(devtools)
document()

library(usethis)
use_build_ignore("./data-raw")

#library
setwd("/Users/wz262/Projects")
install("dimension")
3
library(dimension)

setwd("/Users/wz262/Projects/dimension")
library(devtools)
document()
load_all()
devtools::test()
devtools::check(vignettes = FALSE)
?x_sim
?check_dim_matrix
?subspace
?dimension
?truncate
?lung
#lintr
setwd("/Users/wz262/Projects/dimension")
library(devtools)
a <- lint(".", cache = FALSE)
library(dplyr)
a %>% as_tibble()
library(covr)
a <- covr::report()
a

subspace1_ref <- readRDS("./tests/testthat/reference_data/Subspace1.rds")
saveRDS(x_denoised4, "./tests/testthat/reference_data/x_denoised4.rds")
x_denoised1_ref <- readRDS("./tests/testthat/reference_data/x_denoised1.rds")
# test examples
x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
results <- x %>%
  create_subspace(components = 1:30) %>%
  correct_eigenvalues() %>%
  estimate_rank_kmeans()
km_plot(results$within_var)
str(results)


results <- dimension(x, components = 1:30)
Subspace <- subspace(x, components = 1:10)
results <- dimension(s = Subspace, method = "double_posterior")
results <- dimension(s = Subspace, method = "posterior")
results <- dimension(s = Subspace, method = "kmeans")
results <- dimension(x, method = "ladle")

str(results)
plot(results$subspace, changepoint = results$dimension,annotation = 10)
legacyplot(results$bcp_irl, annotation = 5)
legacyplot(results, annotation = 10)
plot(results, changepoint = results$dimension,annotation = 10)
x %>% dimension() %>% legacyplot()
x %>% dimension() %>% plot()
x %>% subspace(1:15) %>% truncate(location = 1:4)

km_plot(results$within_var)
# prob_irl   <- c(results$bcp_irl$posterior.prob[-10], 0)
# prob_irl_diff <- abs(diff(prob_irl))
# sign_irl_diff <- sign(diff(prob_irl))
# km_min <- which.min(results$within_var[,1])
# sign_irl_diff[km_min]


#Test on MKDim package
x <- x_sim(n = 150, p = 100, ncc = 30, var = c(rep(10,5),rep(3,25)))
t1 <- proc.time()
Subspace <- subspace(x, components = 1:30, times = 10)
print(proc.time() - t1)
gc()
Subspace$irl$eigen
results <- dimension(subspace_ = Subspace)
results <- dimension(x, components = 1:50, times = 10)

#Test on ladle function
stime1 <- system.time(test1 <- lx <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
                      adle(x = t(x), method = "pca"))
cat("d = ", test1$d, "; elapsed = ", stime1[[3]], "\n")
ggplot() + geom_line(aes(x = 1:length(test1$gn), y = test1$gn),colour="black") +
  geom_point(aes(x = 1:length(test1$gn), y = test1$gn), color = "red") +
  ggtitle(paste0(
    "The Ladle Plot", "\n", "n_", 100, "_p_", 500, "_d_", 10, "_est d_", test1$d)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(x = 1:length(test1$gn), y = test1$gn, label=1:length(test1$gn)),
                  colour = "black", size = 5)
# ggsave(file.path(output_dir,"ladle_obj_function_plot.pdf"))

#Test on hoff function
x <- x_sim(n = 100, p = 500, ncc = 10, var = 6)
hoff_result <- hoff(y = t(x), NSCAN = 10)
########################################################
#####test on check_dim_matrix#########
########################################################
#
# params <- check_dim_matrix(x, rnk = 50)
# params <- check_dim_matrix(x)
# str(params)
#
# #test on warning
# params <- check_dim_matrix(rnk = 30)
# params <- check_dim_matrix(x, rnk = -1)
# params <- check_dim_matrix(x, rnk = 2000)
#
#
# MarchenkoPasturPar(ndf = 150, pdim = 100, var = 1, svr = params$svr)

########################################################
#####test on subspace#########
########################################################
subspace1_ref <- readRDS("./tests/testthat/reference_data/Subspace1.rds")
Subspace <- subspace(x)
Subspace <- subspace(x, num_est_samples = 10)
Subspace <- subspace(x, components = 11:30)
Subspace <- subspace(x, components = c(2, 3, 6), num_est_samples = 10)
Subspace <- subspace(x, components = 1:20, num_est_samples = 10, mp= FALSE)
Subspace
str(Subspace)


#test on scree plot
plot(Subspace)
plot(Subspace, annotation = 0)
plot(Subspace, changepoint = 0, annotation = "0")
plot(Subspace, annotation = -1)
plot(Subspace, annotation = 110)
plot(Subspace, annotation = 1:10)

plot(Subspace, changepoint = "0")
plot(Subspace, changepoint = 1.1)
plot(Subspace, changepoint = 1:5)
plot(Subspace, changepoint = 110)

#test on warning
Subspace <- subspace(components = -1, times = 10)
Subspace <- subspace(x, times = -1)
Subspace <- subspace(x, components = -1, times = 10)


########################################################
#####test on dimension#########
########################################################

results <- dimension(x)
results <- dimension(x, components = 50, times = 10)
results <- dimension(x, Subspace)
results <- dimension(subspace_ = Subspace)
results <- dimension(subspace_ = subspace(x))
results <- dimension(x, components = 1:40, times = 10)
str(results)
plot.subspace(results$Subspace, changepoint = results$dimension, annotation=10)

#test on warning
results <- dimension(subspace_ = subspace(x, mp= FALSE))
results <- dimension(x, components=-1, times = -1)
results <- dimension(x, components=1:40, times = -1)
results <- dimension(x, times = 199)

#test on legacy plot
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 30)
modified_legacyplot(results$Changepoint$bcp_post, annotation = 30)
legacyplot(results$Changepoint$bcp_post)


#test on legacy plot warning
modified_legacyplot(results$Changepoint$bcp_irl)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 0)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = "0")
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 160)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = 1:2)
modified_legacyplot(results$Changepoint$bcp_irl, annotation = c(1:5, 55:60))
########################################################
#####test on clipped#########
########################################################

x_denoised <- truncate(x, components = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
x_denoised
str(x_denoised); x_denoised$xi_denoised
x_denoised<-truncate(x, components = 20, method = "hard", zeroout = TRUE)
str(x_denoised);x_denoised$xi_denoised
x_denoised
x_denoised<-truncate(x, components = 20, method = "hard", zeroout = FALSE)
str(x_denoised);x_denoised$xi_denoised
x_denoised
x_denoised<-truncate(x, components = 20, method = "identity", location = c(1:15), zeroout = TRUE,verbose = FALSE)
str(x_denoised);x_denoised$xi_denoised
x_denoised

Subspace <- subspace(x, components = 1:40, times = 10)
x_denoised<-truncate(x, Subspace,method="threshold",alpha = 0.9,zeroout = TRUE)
x_denoised<-truncate(subspace_ = Subspace,method = "threshold",alpha = 0.9,zeroout = TRUE)
x_denoised<-truncate(subspace_ = Subspace,method = "hard",zeroout = TRUE)
x_denoised<-truncate(subspace_ = Subspace,method = "identity",location = c(1:5),zeroout = TRUE)

#test on warning
x_denoised<-truncate(x,method = "threshold",alpha = 0,zeroout = TRUE)
x_denoised<-truncate(x,method = "threshold",alpha = -1,zeroout = TRUE)
x_denoised<-truncate(x,method = "threshold",alpha = 1.9,zeroout = TRUE)
x_denoised<-truncate(x,method = "threshold",zeroout = TRUE)


x_denoised<-truncate(x, components = 20,method = "hard",alpha = 0,zeroout = TRUE)
x_denoised<-truncate(x, components = 20,method = "hard",alpha = 1.9,zeroout = TRUE)
x_denoised<-truncate(x, components = 20,method = "hard",location = c(-1, 2),zeroout = FALSE)

x_denoised<-truncate(x, components = 20,method = "identity",location = c(0:5),zeroout = FALSE)
x_denoised<-truncate(x, components = 20,method = "identity",location = "zero",zeroout = FALSE)
x_denoised<-truncate(x, components = 20,method = "identity",zeroout = FALSE)


# #relation between cor and crossprod(x)
# xstd<-sweep(x, 2L, colMeans(x))
# dim(x)
# ec_cov<-cov(x)
# dim(ec_cov)
# params = check_dim_matrix(xstd)
# str(params);params$irl$eigen
# E_truncate = params$eigenvec%*%diag(params$irl$eigen*100)%*%t(params$eigenvec) / (params$ndf - 1L)
#
# dim(E_truncate)
# E_truncate[1:10,1:10]
# ec_cov[1:10,1:10]
#
# ## symmetric rescaling
# V_truncate<-E_truncate / tcrossprod(diag(E_truncate) ^ 0.5)
# V_truncate[1:10,1:10]
#
# ec_cor<-cor(x)
# ec_cor[1:10,1:10]
#
# x_denoised<-clipped(x,method="threshold",alpha=1,zeroout=TRUE)
# str(x_denoised)
# x_denoised$E_clipped[1:10,1:10]
#
# x_denoised<-clipped(x,rnk=20,method="threshold",alpha=1,zeroout=TRUE)
# str(x_denoised)
# x_denoised$E_clipped[1:10,1:10]
#
# #relation between cor and svd
# tmp<-svd(x)
# rnk=20
# d<-c(tmp$d[1:rnk],rep(0,length(tmp$d)-rnk))
# x_svd<-tmp$u%*%diag(d)%*%t(tmp$v)
# cor(x_svd)[1:10,1:10]

