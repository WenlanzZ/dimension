library(devtools)
setwd("/Users/wz262/Projects")
install("MKDim")
library(MKDim)
?Xsim
?CheckDimMatrix
?MarcenkoPasturSample
?OptimumDimension
?clipped

#auto generating documents after changing fxns
setwd("/Users/wz262/Projects/MKDim")
library(roxygen2)
library(devtools)
document()

#Test on MKDim package
X <- Xsim(n=150,p=100,ncc=100,var=2,fact = 1,orthogonl = FALSE)

########################################################
#####test on CheckDimMatrix#########
########################################################

params = CheckDimMatrix(X,rnk=30)
params = CheckDimMatrix(X)

#test on warning
params = CheckDimMatrix(rnk=30)
params = CheckDimMatrix(X,rnk=-1)
params = CheckDimMatrix(X,rnk=2000)
str(params)

MarchenkoPasturPar(ndf=150,pdim=100,var=1,svr=params$svr)

########################################################
#####test on MarcenkoPasturSample#########
########################################################

MPSamples = MarcenkoPasturSample(X)
MPSamples = MarcenkoPasturSample(X,time=10)
MPSamples = MarcenkoPasturSample(X,rnk=10)
MPSamples = MarcenkoPasturSample(X,params,times=10)
MPSamples = MarcenkoPasturSample(params=params,times=10)
MPSamples = MarcenkoPasturSample(X,rnk=40,times=10)
#test on warning
MPSamples = MarcenkoPasturSample(X,time=-1)
MPSamples = MarcenkoPasturSample(X,rnk=-1,times=10)
MPSamples = MarcenkoPasturSample(X,params,rnk=40,times=10)
str(MPSamples)

########################################################
#####test on OptimumDimension#########
########################################################

results = OptimumDimension(X)
results = OptimumDimension(X,rnk=30)
results = OptimumDimension(X,MPSamples)
results = OptimumDimension(MPSamples=MPSamples)

#test on warning
results = OptimumDimension(X,rnk=30,times=100,p=0.95)
results = OptimumDimension(X,times=100)
str(results)

ScreePlot(results$MarcenkoPasturSample,Changepoint=results$Changepoint$changePoint,annotation=0)
modified_legacyplot(results$Changepoint$bcp_irl,annotation=10,medianfilter = FALSE)
modified_legacyplot(results$Changepoint$bcp_post,annotation=10,medianfilter = FALSE)
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
