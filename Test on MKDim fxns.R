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

params = CheckDimMatrix(X,rnk=-1)
params = CheckDimMatrix(X)
params = CheckDimMatrix(X,rnk=10)
str(params)

MPSamples = MarcenkoPasturSample(X,time=-1)
MPSamples = MarcenkoPasturSample(X,time=10)
MPSamples = MarcenkoPasturSample(X,rnk=10)
MPSamples = MarcenkoPasturSample(X,rnk=-1,times=10)
MPSamples = MarcenkoPasturSample(X,params,times=10)
MPSamples = MarcenkoPasturSample(X,rnk=40,times=10)
MPSamples = MarcenkoPasturSample(X,params,rnk=40,times=10)
str(MPSamples)

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

modified_legacyplot(test1$Changepoint$bcp.irl,annotation=10)
modified_legacyplot(test1$Changepoint$bcp_post,annotation=10)
modified_legacyplot(test1$Changepoint$bcp_diff)
ScreePlot(test1$MarcenkoPasturSample,test1$Changepoint$changePoint,annotation=10)


bcp.irl = bcp(as.vector(test1$MarcenkoPasturSample$irl$eigen-test1$MarcenkoPasturSample$MP_irl$eigen), p0 = 0.1)
modified_legacyplot(bcp.irl)

bcp_post = bcp(as.vector(c(bcp.irl$posterior.prob[-rnk],0)),p0=0.1)
modified_legacyplot(bcp_post,annotation=10)

bcp_diff = bcp(as.vector(-diff(test1$MarcenkoPasturSample$irl$eigen-test1$MarcenkoPasturSample$MP_irl$eigen)),p0=0.01)
bcp_diff
modified_legacyplot(bcp_diff)
legacyplot(bcp_diff)

#median filter
library(signal)
t <- seq(0, 1, len=100)                            # 1 second sample
x <- sin(2*pi*t*2.3) + 0.25*rlnorm(length(t), 0.5) # 2.3 Hz sinusoid+noise
plot(t, x, type = "l")
# 3-point filter
lines(t, medfilt1(x), col="red", lwd=2)
# 7-point filter
lines(t, filter(MedianFilter(7), x), col = "blue", lwd=2) # another way to call it

plot(1:10,c(bcp.irl$posterior.prob[-rnk],0),"l")
lines(1:10, medfilt1(c(bcp.irl$posterior.prob[-rnk],0)), col="red", lwd=2)
lines(1:10, filter(MedianFilter(7),c(bcp.irl$posterior.prob[-rnk],0)), col="red", lwd=2)
