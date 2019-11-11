library(devtools)
setwd("/Users/wz262/Projects")
install("MKDim")
library(MKDim)
?Xsim
?CheckDimMatrix
?MarcenkoPasturSample
?OptimumDimension

#auto generating documents after changing fxns
setwd("/Users/wz262/Projects/MKDim")
library(roxygen2)
library(devtools)
document()

#Test on MKDim package
X <- Xsim(n=1000,p=500,ncc=10,var=2,fact = 1,orthogonl = FALSE)

params = CheckDimMatrix(X,rnk=-1)
params = CheckDimMatrix(X)
params = CheckDimMatrix(rnk=30)
params = CheckDimMatrix(X,rnk=30)
str(params)

MPSamples = MarcenkoPasturSample(X)
MPSamples = MarcenkoPasturSample(X,time=-1)
MPSamples = MarcenkoPasturSample(X,time=10)
MPSamples = MarcenkoPasturSample(X,rnk=10)
MPSamples = MarcenkoPasturSample(X,rnk=-1,times=10)
MPSamples = MarcenkoPasturSample(X,params,times=10)
MPSamples = MarcenkoPasturSample(params=params,times=10)
MPSamples = MarcenkoPasturSample(X,rnk=40,times=10)
MPSamples = MarcenkoPasturSample(X,params,rnk=40,times=10)
str(MPSamples)

results = OptimumDimension(X,rnk=30,times=100)
results = OptimumDimension(X,rnk=30)
results = OptimumDimension(X,times=100)
results = OptimumDimension(X,MPSamples)
results = OptimumDimension(X)
results = OptimumDimension(MPSamples=MPSamples)
str(results)

ScreePlot(results$MarcenkoPasturSample,Changepoint=results$Changepoint$changePoint,annotation=0)
modified_legacyplot(results$Changepoint$bcp.irl,annotation=10,medianfilter = FALSE)
modified_legacyplot(results$Changepoint$bcp_post,annotation=10)
legacyplot(results$Changepoint$bcp_post)

#test on MKdim
n=200
p=40
d=3
M<-diag(c(2,1,1,rep(0,p-3)))
sigma <- diag(rep(0.54^2,p))+M
library(mvtnorm)
cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
rnk=10
system.time(test1<-OptimumDimension(cor_cols))
str(test1)

modified_legacyplot(test1$Changepoint$bcp.irl,annotation=10,medianfilter = FALSE)
modified_legacyplot(test1$Changepoint$bcp_post,annotation=10)
modified_legacyplot(test1$Changepoint$bcp_post,medianfilter = FALSE)
modified_legacyplot(test1$Changepoint$bcp_diff,medianfilter = FALSE)
modified_legacyplot(test1$Changepoint$bcp_diff,medianfilter = TRUE)
ScreePlot(test1$MarcenkoPasturSample,test1$Changepoint$changePoint,annotation=10)
