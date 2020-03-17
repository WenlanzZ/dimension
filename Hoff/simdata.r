#source("svd.f.r")
##### generate simulated data
for(seed in 10 +10*c(1:9)) {
set.seed(seed)

m<-100
n<-100
K0<-5
phi0<-1
mu0<-sqrt(n+m+2*sqrt(n*m))

U0<-rnorm(m) ; U0<-U0/sqrt(sum(U0^2))
V0<-rnorm(n) ; V0<-V0/sqrt(sum(V0^2))
for(j in 1:(n-1)) {
u<-rnorm(m-j) ; u<-u/sqrt(sum(u^2)) ; U0<-cbind(U0,Null(U0)%*%u)
v<-rnorm(n-j) ; v<-v/sqrt(sum(v^2)) ; V0<-cbind(V0,Null(V0)%*%v)
                   }

D0<-diag ( c(runif(K0,mu0/2,3*mu0/2), rep(0,n-K0) )  )


M0<- U0%*%D0%*%t(V0)

Y0<-M0+matrix(rnorm(m*n),m,n)*sqrt(1/phi0)
Y<-Y0

######
ofname<-paste(m,".",n,"/OUT",seed,sep="")
cat(seed,"\n")
source("svd.r")
 }


