rm(list=lsf.str())
gc()

#auto generating documents after changing fxns
setwd("/Users/wz262/Projects/dimension")
library(roxygen2)
library(devtools)
document()

#library
setwd("/Users/wz262/Projects")
install("dimension")
3
library(dimension)

setwd("/Users/wz262/Projects/dimension")
library(devtools)
load_all()
?x_sim
?check_dim_matrix
?subspace
?dimension
?clipped

#lintr
setwd("/Users/wz262/Projects/dimension")
library(devtools)
a <- lint(".", cache = FALSE)
library(dplyr)
a %>% as_tibble()



#Test on MKDim package
# x <- x_sim(n = 150, p = 100, ncc = 30, var = c(rep(10,5),rep(3,25)))
x <- x_sim(n = 100, p = 150, ncc = 10, var = 10)
t1 <- proc.time()
Subspace <- subspace(x, components = 1:30, times = 10)
print(proc.time() - t1)
gc()
Subspace$irl$eigen
########################################################
#####test on check_dim_matrix#########
########################################################

params <- check_dim_matrix(x, rnk = 50)
params <- check_dim_matrix(x)
str(params)

#test on warning
params <- check_dim_matrix(rnk = 30)
params <- check_dim_matrix(x, rnk = -1)
params <- check_dim_matrix(x, rnk = 2000)


MarchenkoPasturPar(ndf = 150, pdim = 100, var = 1, svr = params$svr)

########################################################
#####test on subspace#########
########################################################

Subspace <- subspace(x)
Subspace <- subspace(x, time = 10)
Subspace <- subspace(x, components = 11:30)
Subspace <- subspace(x, components = c(2, 3, 6), times = 10)
Subspace <- subspace(x, components = 1:20, times = 10, mp= FALSE)
Subspace
str(Subspace)

#test on scree plot
plot(Subspace)
plot(Subspace, changepoint = 0)
plot(Subspace, annotation = 0)
plot(Subspace, changepoint = 0, annotation = 5)

#test on warning
Subspace <- subspace(components = -1, times = 10)
Subspace <- subspace(x, time = -1)
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


########################################################
#####test on clipped#########
########################################################

x_clp <- clipped(x, components = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
x_clp
str(x_clp); x_clp$xi_clipped
x_clp<-clipped(x, components = 20, method = "hard", zeroout = TRUE)
str(x_clp);x_clp$xi_clipped
x_clp
x_clp<-clipped(x, components = 20, method = "hard", zeroout = FALSE)
str(x_clp);x_clp$xi_clipped
x_clp
x_clp<-clipped(x, components = 20, method = "identity", location = c(1:15), zeroout = TRUE,verbose = FALSE)
str(x_clp);x_clp$xi_clipped
x_clp

load_allSubspace <- subspace(x, components = 1:40, times = 10)
x_clp<-clipped(x,Subspace,method="threshold",alpha = 0.9,zeroout = TRUE)
x_clp<-clipped(subspace_ = Subspace,method = "threshold",alpha = 0.9,zeroout = TRUE)
x_clp<-clipped(subspace_ = Subspace,method = "hard",zeroout = TRUE)
x_clp<-clipped(subspace_ = Subspace,method = "identity",location = c(1:5),zeroout = TRUE)

#test on warning
x_clp<-clipped(x,method = "threshold",alpha = 0,zeroout = TRUE)
x_clp<-clipped(x,method = "threshold",alpha = -1,zeroout = TRUE)
x_clp<-clipped(x,method = "threshold",alpha = 1.9,zeroout = TRUE)
x_clp<-clipped(x,method = "threshold",zeroout = TRUE)


x_clp<-clipped(x,components = 20,method = "hard",alpha = 0,zeroout = TRUE)
x_clp<-clipped(x,components = 20,method = "hard",alpha = 1.9,zeroout = TRUE)
x_clp<-clipped(x,components = 20,method = "hard",location = c(-1, 2),zeroout = FALSE)

x_clp<-clipped(x,components = 20,method = "identity",location = c(0:5),zeroout = FALSE)
x_clp<-clipped(x,components = 20,method = "identity",location = "zero",zeroout = FALSE)
x_clp<-clipped(x,components = 20,method = "identity",zeroout = FALSE)




# #relation between cor and crossprod(x)
# xstd<-sweep(x, 2L, colMeans(x))
# dim(x)
# ec_cov<-cov(x)
# dim(ec_cov)
# params = check_dim_matrix(xstd)
# str(params);params$irl$eigen
# E_clipped = params$eigenvec%*%diag(params$irl$eigen*100)%*%t(params$eigenvec) / (params$ndf - 1L)
#
# dim(E_clipped)
# E_clipped[1:10,1:10]
# ec_cov[1:10,1:10]
#
# ## symmetric rescaling
# V_clipped<-E_clipped / tcrossprod(diag(E_clipped) ^ 0.5)
# V_clipped[1:10,1:10]
#
# ec_cor<-cor(x)
# ec_cor[1:10,1:10]
#
# x_clp<-clipped(x,method="threshold",alpha=1,zeroout=TRUE)
# str(x_clp)
# x_clp$E_clipped[1:10,1:10]
#
# x_clp<-clipped(x,rnk=20,method="threshold",alpha=1,zeroout=TRUE)
# str(x_clp)
# x_clp$E_clipped[1:10,1:10]
#
# #relation between cor and svd
# tmp<-svd(x)
# rnk=20
# d<-c(tmp$d[1:rnk],rep(0,length(tmp$d)-rnk))
# x_svd<-tmp$u%*%diag(d)%*%t(tmp$v)
# cor(x_svd)[1:10,1:10]
