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
X <- Xsim(n=100,p=90,ncc=2,var=10,fact = 30,orthogonl=T)
params = CheckDimMatrix(X,rnk=40)
MPSamples = MarcenkoPasturSample(X,params,rnk=10,times=10)
str(MPSamples)
#test on MKdim
n=200
p=40
d=3
M<-diag(c(2,1,1,rep(0,p-3)))
sigma <- diag(rep(0.54^2,p))+M
cor_cols <- rmvnorm(n, rep(0, p), sigma = sigma)
rnk=10
system.time(test1<-OptimumDimension(cor_cols,rnk=10,times=1000))
str(test1)

modified_legacyplot(test1$Changepoint$bcp.irl,annotation=10)

bcp_post = bcp(as.vector(c(test1$Changepoint$bcp.irl$posterior.prob[-rnk],0)),p0=0.1)
modified_legacyplot(bcp_post,annotation=10)


ScreePlot(test1$MarcenkoPasturSample,Changepoint=1,annotation=0)
