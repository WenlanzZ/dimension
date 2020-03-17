##### Precondition: Y= m*n  matrix with m>=n
#library(dimension)
#source("svd.f.r")

##### constants 
m<-dim(Y)[1]
n<-dim(Y)[2]
#####

##### hyperparameters

sY<-svd(Y)
s20.est<-var(c(Y))
t20.est<-0
mu0.est<-0
for(k in 1:n) { 
  ks<-seq(1,k,length=k)
  s20.est<-c(s20.est, 
             var(c((Y-sY$u[,ks]%*%diag(sY$d[ks],nrow=k)%*%t(sY$v[,ks]) )) ) )
    t20.est<-c(t20.est,var(sY$d[ks])  )
    mu0.est<-c(mu0.est,mean(sY$d[ks]) )
                   }
t20.est[2]<-0

## prior for phi
nu0<-2
s20<-mean(s20.est)


## prior for psi
eta0<-2
t20<-mean(t20.est)


## prior for mu
#kap0<-1
#mu0<-mean(mu0.est)
mu0<-mean(mu0.est)
premu0<-1/var(mu0.est)


##### starting values
K0<-0
U<-matrix(0,m,n) ; V<-matrix(0,n,n) 
U[,seq(1,K0,length=K0)]<-sY$u[,seq(1,K0,length=K0)] 
V[,seq(1,K0,length=K0)]<-sY$v[,seq(1,K0,length=K0)] 
D<-diag( c(sY$d[seq(1,K0,length=K0)],rep(0,n-K0))) 
phi<- 1/s20
psi<-1/t20
mu<-mu0
#####

##### MCMC
NSCAN<-10000
odens<-10
OUT<-matrix(nrow=NSCAN/odens,ncol=5)
MSE<-NULL
M.ps<-Y*0
D.ps<-rep(0,n)
for(ns in 1:NSCAN) {

gibbs.UVD.varrank(U,V,D,Y,phi,mu,psi,min(n,10))

### fixed rank update
gibbs.UVD.fixedrank(U,V,D,Y,phi,mu,psi)

### update phi
phi<-rgamma(1,  ( nu0+m*n)/2, (nu0*s20+sum( (Y-U%*%D%*%t(V))^2))/2 )

### update mu,psi
mu<-rnorm(1,(premu0*mu0+psi*sum(diag(D)))/
            (premu0+psi*sum(D!=0)),
	                1/sqrt(premu0+psi*sum(D!=0)) )
psi<-rgamma(1, (eta0+sum(D!=0) )/2,  (eta0*t20 + sum((D[D!=0]-mu)^2))/2 )




M.ps<-M.ps+U%*%D%*%t(V)
D.ps<-D.ps+ -sort(-diag(D))


# ### output
# if(ns %% odens==0) {
# out<-c(ns,sum(D!=0),mean((M0-M.ps/ns)^2),mean((M0-U%*%D%*%t(V))^2),
#        1/phi,mu,1/psi)
# cat(out,"\n")
# OUT[ns/odens,]<-out
#                     }

                    ### output
if(ns %% odens==0) {
out<-c(ns,sum(D!=0),1/phi,mu,1/psi)
cat(out,"\n")
OUT[ns/odens,]<-out
                    }
         }

##### end mcmc
# dput(OUT,ofname)




