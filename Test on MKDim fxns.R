rm(list=lsf.str())

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
?Xsim
?CheckDimMatrix
?subspace
?dimension
?clipped


#Test on MKDim package
X <- Xsim(n = 150, p = 100, ncc = 10, var = 2, fact = 1)

########################################################
#####test on CheckDimMatrix#########
########################################################

params = CheckDimMatrix(X, rnk = 30)
params = CheckDimMatrix(X)

#test on warning
params = CheckDimMatrix(rnk = 30)
params = CheckDimMatrix(X, rnk = -1)
params = CheckDimMatrix(X, rnk = 2000)
str(params)

MarchenkoPasturPar(ndf = 150, pdim = 100, var = 1, svr = params$svr)

########################################################
#####test on subspace#########
########################################################

subspace_ <- subspace(X)
subspace_ <- subspace(X, time = 10)
subspace_ <- subspace(X, rank = 1:10, basis = "eigen")
subspace_ <- subspace(X, rank = c(2, 3, 6), times = 10, basis = "eigen")
subspace_
str(subspace_)
plot(subspace_, Changepoint = 0)
plot(subspace_, annotation = 0)
plot(subspace_, Changepoint = 0, annotation = 5)
#test on warning
subspace_ <- subspace(X, time = -1)
subspace_ <- subspace(X, rank = -1, times = 10)
str(subspace_)

########################################################
#####test on OptimumDimension#########
########################################################

results = dimension(X)
results = dimension(X, rank = 30, times = 10)
results = dimension(X, subspace_)
results = dimension(subspace_ = subspace_)
results = dimension(X, rank=1:40, times = 10, basis="eigen")
#test on warning
results = dimension(X, rank=1:40, times = 10, basis="eigen")
results = dimension(X,times=99)
str(results)

modified_legacyplot(results$Changepoint$bcp_irl,annotation=10)
modified_legacyplot(results$Changepoint$bcp_post,annotation=10)
legacyplot(results$Changepoint$bcp_post)


########################################################
#####test on clipped#########
########################################################

X_clp<-clipped(X,rnk=20,method="threshold",alpha=0.9,zeroout=TRUE)
str(X_clp);X_clp$xi_clipped
X_clp<-clipped(X,rnk=20,method="hard",zeroout=TRUE)
str(X_clp);X_clp$xi_clipped
X_clp<-clipped(X,rnk=20,method="hard",zeroout=FALSE)
str(X_clp);X_clp$xi_clipped
X_clp<-clipped(X,rnk=20,method="identity",location=c(1:15),zeroout=FALSE)
str(X_clp);X_clp$xi_clipped

#test on warning
params <- CheckDimMatrix(X,rnk=20)
X_clp<-clipped(X,params,rnk=20,method="threshold",alpha=0.9,zeroout=TRUE)
X_clp<-clipped(X,params,rnk=40,method="threshold",alpha=0.9,zeroout=TRUE)
X_clp<-clipped(params=params,rnk=40,method="threshold",alpha=0.9,zeroout=TRUE)
X_clp<-clipped(params=params,method="threshold",alpha=0.9,zeroout=TRUE)

X_clp<-clipped(X,rnk=20,method="threshold",alpha=0,zeroout=TRUE)
X_clp<-clipped(X,rnk=20,method="threshold",alpha=-1,zeroout=TRUE)
X_clp<-clipped(X,rnk=20,method="threshold",alpha=1.9,zeroout=TRUE)
X_clp<-clipped(X,rnk=20,method="threshold",zeroout=TRUE)


X_clp<-clipped(X,rnk=20,method="hard",alpha=0,zeroout=TRUE)
X_clp<-clipped(X,rnk=20,method="hard",alpha=1.9,zeroout=TRUE)
X_clp<-clipped(X,rnk=20,method="hard",location=c(-1,2),zeroout=FALSE)

X_clp<-clipped(X,rnk=20,method="identity",location=c(0:5),zeroout=FALSE)
X_clp<-clipped(X,rnk=20,method="identity",location="zero",zeroout=FALSE)
X_clp<-clipped(X,rnk=20,method="identity",zeroout=FALSE)




#relation between cor and crossprod(X)
Xstd<-sweep(X, 2L, colMeans(X))
dim(X)
ec_cov<-cov(X)
dim(ec_cov)
params = CheckDimMatrix(Xstd)
str(params);params$irl$eigen
E_clipped = params$eigenvec%*%diag(params$irl$eigen*100)%*%t(params$eigenvec) / (params$ndf - 1L)

dim(E_clipped)
E_clipped[1:10,1:10]
ec_cov[1:10,1:10]

## symmetric rescaling
V_clipped<-E_clipped / tcrossprod(diag(E_clipped) ^ 0.5)
V_clipped[1:10,1:10]

ec_cor<-cor(X)
ec_cor[1:10,1:10]

X_clp<-clipped(X,method="threshold",alpha=1,zeroout=TRUE)
str(X_clp)
X_clp$E_clipped[1:10,1:10]

X_clp<-clipped(X,rnk=20,method="threshold",alpha=1,zeroout=TRUE)
str(X_clp)
X_clp$E_clipped[1:10,1:10]

#relation between cor and svd
tmp<-svd(X)
rnk=20
d<-c(tmp$d[1:rnk],rep(0,length(tmp$d)-rnk))
X_svd<-tmp$u%*%diag(d)%*%t(tmp$v)
cor(X_svd)[1:10,1:10]
