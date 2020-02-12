#### do the whole sampling thing
gibbs.sample  <-  function(E, Nu, Nv, phi, mu, psi, lor, llb = 50, lub = 2000) {
  E2  <-  sum(E^2)
  n <- dim(E)[2]
  m <- dim(E)[1]
  E0 <-  t(Nu) %*% E %*% Nv
  n0 <- dim(E0)[2]
  m0 <- dim(E0)[1]
  svdE0 <- svd(E0)
  z <- svdE0$d^2
  E02 <- sum(z)
  theta <- z / E02
  del <- 0
  u <- rep(0, m0)
  v <- rep(0, n0)

  ### critical value for approximating bessel function
  amax <-  svdE0$d[1]^2 * phi^2 * E02 * max(theta)
  lcrit <-  .5 * (sqrt((n0 / 2 + 1)^2  + 4 * (amax-n0 / 2)) -(n0 / 2 + 1))
  lmax <-  1.25 * (lcrit)
  lmax <-  min(max(lmax, llb),  lub)
  lmax <- 2 * as.integer(lmax / 2)
  lseq <- 0:lmax

  ### sums
  la <-  lgamma(m0 / 2) - lgamma(m0 / 2 + lseq) - lgamma(lseq + 1) -lseq * log(4)  +
       lcr(theta, lmax)
  if(mu == 0) {
  lb <-   .5 * log(psi / (phi + psi))  + lseq * log(1 / (phi + psi)) +
       (lgamma(2 * lseq + 1) - lgamma(lseq + 1) -lseq * log(2))
  }
  if(mu != 0) {
    lb <-  ln2moment(mu * psi / (phi + psi), sqrt(1 / (psi + phi)), lmax) +
    .5 * log(psi / (phi + psi)) -
    .5 *  mu^2 * psi * phi / (phi + psi)
  }
  lc <-  lseq * (log(E02 * phi^2))
  lt <-  la + lb + lc
  lpe.r10 <- max(lt)  + log(sum(exp(lt - max(lt))))
  s <- rbinom(1, 1, 1 / (1 + exp( - (lor + lpe.r10))))
  if(s == 1) {
  ####sample del
    pdl <-  exp(lt - max(lt))
    l <- sample(x = lseq, size = 1, prob = pdl)
    if(l > 0 & mu != 0) {
      del <- rxl(mu * psi / (psi + phi),  sqrt(1 / (phi + psi)),  l)
    }
    if(l > 0 & mu == 0) {
      del <- del <- sqrt(rchisq(1, 2 * l + 1) / (phi + psi))
    }
    if(l == 0) {
      del <- rnorm(1, mu * psi / (psi + phi), sqrt(1 / (phi + psi)))
    }
  ####sample u, v
  tmp <- ruv.A(phi * del * E0)
  u <- tmp$u
  v <- tmp$v
  }
  list(del = del, u = u, v = v)
}

### computing the dirichlet average
# system("R CMD SHLIB -o svd.so svd.c")
dyn.load("svd.so")
lcr <- function(theta, L) {
.C("lcr", as.double(theta), as.integer(L), as.integer(length(theta)),
       lr = double(L + 1))$lr
}

rW <- function(kap, m) {
  .C("rW", kap = as.double(kap), m = as.integer(m), w = double(1))$w
}

ln2moment <- function(mu, sigma, lmax) {
  .C("ln2moment", as.double(abs(mu)), as.double(sigma^2), as.integer(lmax),
                 l2mom = double(lmax + 1))$l2mom
}

labc <- function(theta, lmax, m, E2, phi, mu, psi) {
  .C("labc", theta = as.double(theta), L = as.integer(lmax), m = as.integer(m),
            n = as.integer(length(theta)), E2 = as.double(E2), phi = as.double(phi),
      mu = as.double(abs(mu)), psi = as.double(psi), lt = double(lmax + 1),
         lmom = double(2 * lmax + 1), lr = double(lmax + 1))
}
lsumt <- function(theta, lmax, m, E2, phi, mu, psi) {
  .C("lsumt", theta = as.double(theta), L = as.integer(lmax), m = as.integer(m),
       n = as.integer(length(theta)), E2 = as.double(E2), phi = as.double(phi),
       mu = as.double(abs(mu)), psi = as.double(psi), lt = double(lmax + 1))$lt
}

rvmf <- function(kmu) {
  kap <- sqrt(sum(kmu^2))
  mu <- kmu / kap
  m <- length(mu)

  if(m == 1){
    u <-  (-1)^rbinom(1, 1, 1 / (1 + exp(2 * kap * mu)))
  }

  if(m > 1) {
    W <- rW(kap, m)
    V <- rnorm(m-1)
    V <- V / sqrt(sum(V^2))
    x <- c((1-W^2)^.5 * t(V), W)
    u <- cbind(Null(mu), mu) %*% x
  }
  u
}




###
gibbs.UVD.fixedrank <- function(U, V, D, Y, phi, mu, psi) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  for(j in resample((1:n)[diag(Dl) != 0])) {
  Nu <- Null(Ul[, -j])
  Nv <- Null(Vl[, -j])
  if(sum(Dl[-j, -j] != 0) == 0) {
    Nu <- diag(nrow = m)
    Nv <- diag(nrow = n)
    }
  E <-  Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
  Vl[, j] <-  Nv %*% rvmf(Dl[j, j] * t(Nv) %*% t(E) %*% Ul[, j] * phi)
  Ul[, j] <-  Nu %*% rvmf(Dl[j, j] * t(Nu) %*% E %*% Vl[, j] * phi)

  mn <- (t(Ul[, j]) %*% E %*% Vl[, j] * phi + mu * psi) / (phi + psi)
  se <- sqrt(1 / (phi + psi))
  d <- rnorm(1, mn, se)
  Dl[j, j] <- d
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}



gibbs.UVD.varrank <- function(U, V, D, Y, phi, mu, psi, nsamp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  for(j in resample(1:n, size = nsamp)) {
  Nu <- Null(Ul[, -j])
  Nv <- Null(Vl[, -j])
  if(sum(Dl[-j, -j] != 0) == 0) {
    Nu <- diag(nrow = m)
    Nv <- diag(nrow = n)
    }
  E <- Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
  lor <-  log((sum(diag(Dl)[-j] != 0)  + 1) / (n-sum(diag(Dl)[-j] != 0)))
  tmp <- gibbs.sample(E, Nu, Nv, phi, mu, psi, lor)
  Ul[, j] <- Nu %*% tmp$u
  Vl[, j] <- Nv %*% tmp$v
  Dl[j, j] <- tmp$del
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

################

gibbs.UVD.pvarrank <- function(U, V, D, Y, phi, mu, psi, nsamp, uperp, vperp) {
Ul <- U ; Vl <- V ; Dl <- D
m <- dim(Y)[1];  n <- dim(Y)[2]

for(j in resample(1:n, size = nsamp)) {
Nu <- Null(cbind(uperp, Ul[, -j])) ; Nv <- Null(cbind(vperp, Vl[, -j]))

E <- Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
lor <-  log((sum(diag(Dl)[-j] != 0)  + 1) / (n-sum(diag(Dl)[-j] != 0)))
tmp <- gibbs.sample(E, Nu, Nv, phi, mu, psi, lor)
Ul[, j] <- Nu %*% tmp$u ; Vl[, j] <- Nv %*% tmp$v ; Dl[j, j] <- tmp$del
                                    }
assign("U", Ul, env = parent.env(environment()))
assign("V", Vl, env = parent.env(environment()))
assign("D", Dl, env = parent.env(environment()))
                       }

gibbs.UVD.pfixedrank <- function(U, V, D, Y, phi, mu, psi, uperp, vperp) {
Ul <- U ; Vl <- V ; Dl <- D
m <- dim(Y)[1];  n <- dim(Y)[2]

for(j in resample((1:n)[diag(Dl) != 0])) {
Nu <- Null(cbind(uperp, Ul[, -j])) ; Nv <- Null(cbind(vperp, Vl[, -j]))

E <-  Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
Vl[, j] <-  Nv %*% rvmf(Dl[j, j] * t(Nv) %*% t(E) %*% Ul[, j] * phi)
Ul[, j] <-  Nu %*% rvmf(Dl[j, j] * t(Nu) %*% E %*% Vl[, j] * phi)

mn <-  (t(Ul[, j]) %*% E %*% Vl[, j] * phi + mu * psi) / (phi + psi)
se <- sqrt(1 / (phi + psi))
d <- rnorm(1, mn, se)
Dl[j, j] <- d
                                        }
assign("U", Ul, env = parent.env(environment()))
assign("V", Vl, env = parent.env(environment()))
assign("D", Dl, env = parent.env(environment()))

}

#############

gibbs.UVD.rc.varrank <- function(U, V, D, Y, phi, mu, psi, nsamp, up, vp){
Ul <- U ; Vl <- V ; Dl <- D
m <- dim(Y)[1];  n <- dim(Y)[2]
for(j in resample(3:n, size = nsamp)) {
Nu <- Null(cbind(up, Ul[, -j])) ; Nv <- Null(cbind(vp, Vl[, -j]))

E <- Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
lor <-  log((sum(diag(Dl)[-j] != 0)  + 1) / (n-sum(diag(Dl)[-j] != 0)))
tmp <- gibbs.sample(E, Nu, Nv, phi, mu, psi, lor)
Ul[, j] <- Nu %*% tmp$u ; Vl[, j] <- Nv %*% tmp$v ; Dl[j, j] <- tmp$del
}
assign("U", Ul, env = parent.env(environment()))
assign("V", Vl, env = parent.env(environment()))
assign("D", Dl, env = parent.env(environment()))
}

gibbs.UVD.rc.fixedrank <- function(U, V, D, Y, phi, mu, psi, up, vp) {
Ul <- U ; Vl <- V ; Dl <- D
m <- dim(Y)[1];  n <- dim(Y)[2]

for(j in resample((1:n)[diag(Dl) != 0])) {

E <-  Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])

if(j != 2) {
Nv <- Null(cbind(vp, Vl[, -j]))
Vl[, j] <-  Nv %*% rvmf(Dl[j, j] * t(Nv) %*% t(E) %*% Ul[, j] * phi)
}

if(j != 1) {
Nu <- Null(cbind(up, Ul[, -j]))
Ul[, j] <-  Nu %*% rvmf(Dl[j, j] * t(Nu) %*% E %*% Vl[, j] * phi)
}

mn <-  (t(Ul[, j]) %*% E %*% Vl[, j] * phi + mu * psi) / (phi + psi)
se <- sqrt(1 / (phi + psi))
d <- rnorm(1, mn, se)
Dl[j, j] <- d
}
assign("U", Ul, env = parent.env(environment()))
assign("V", Vl, env = parent.env(environment()))
assign("D", Dl, env = parent.env(environment()))

}



#############
gibbs.UVD.perp.varrank <- function(U, V, D, Y, phi, mu, psi, nsamp, up, vp){
Ul <- U ; Vl <- V ; Dl <- D
m <- dim(Y)[1];  n <- dim(Y)[2]
for(j in resample(1:n, size = nsamp)) {
Nu <- Null(cbind(up, Ul[, -j])) ; Nv <- Null(cbind(vp, Vl[, -j]))

E <- Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
lor <-  log((sum(diag(Dl)[-j] != 0)  + 1) / (n-sum(diag(Dl)[-j] != 0)))
tmp <- gibbs.sample(E, Nu, Nv, phi, mu, psi, lor)
Ul[, j] <- Nu %*% tmp$u ; Vl[, j] <- Nv %*% tmp$v ; Dl[j, j] <- tmp$del
}
assign("U", Ul, env = parent.env(environment()))
assign("V", Vl, env = parent.env(environment()))
assign("D", Dl, env = parent.env(environment()))
}

gibbs.UVD.perp.fixedrank <- function(U, V, D, Y, phi, mu, psi, up, vp) {
Ul <- U ; Vl <- V ; Dl <- D
m <- dim(Y)[1];  n <- dim(Y)[2]

for(j in resample((1:n)[diag(Dl) != 0])) {

E <-  Y-Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
Nv <- Null(cbind(vp, Vl[, -j]))
Vl[, j] <-  Nv %*% rvmf(Dl[j, j] * t(Nv) %*% t(E) %*% Ul[, j] * phi)
Nu <- Null(cbind(up, Ul[, -j]))
Ul[, j] <-  Nu %*% rvmf(Dl[j, j] * t(Nu) %*% E %*% Vl[, j] * phi)

mn <-  (t(Ul[, j]) %*% E %*% Vl[, j] * phi + mu * psi) / (phi + psi)
se <- sqrt(1 / (phi + psi))
d <- rnorm(1, mn, se)
Dl[j, j] <- d
}
assign("U", Ul, env = parent.env(environment()))
assign("V", Vl, env = parent.env(environment()))
assign("D", Dl, env = parent.env(environment()))

}


#############




### other functions
resample  <-  function(x, size, ...)
    if(length(x) < =  1) { if(!missing(size) && size == 0) x[FALSE] else x
    } else sample(x, size, ...)

Null <- function (M)
{
    tmp  <-  qr(M)
    set  <-  if (tmp$rank == 0)
        1:ncol(M)
    else -(1:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}
###


### link functions
lpY.pois <- function(Y, Theta) {  Y * Theta - exp(Theta) -lgamma(Y + 1) }

### link functions
lpY.binary <- function(Y, Theta) {  Y * Theta - log(1 + exp(Theta)) }


### link functions
lpY.tri <- function(Y, Theta) {
    Y * Theta + t(Y) * t(Theta) - log(1 + exp(Theta) + exp(t(Theta))) }


### reg stuff
XB <- function(X, beta) {
tmp <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
for(k in 1:length(beta)) { tmp <- tmp + beta[k] * X[, , k] }
tmp                   }


x.X <- function(X) {
x <- NULL
for(k in 1:dim(X)[3]) { x <- cbind(x, c(X[, , k])) }
x                 }

rmvnorm <- function(mu, Sig2){
R <- t(chol(Sig2))
R %*% (rnorm(length(mu), 0, 1))  + mu }

ldmvnorm <- function(y, mu, Sig){
       c(    -(length(y) / 2) * log(2 * pi) -.5 * log(det(Sig)) -.5 *
                   t(y-mu) %*% solve(Sig) %*% (y-mu) )   }



####




####
ruv.A <- function(A, nsamp = 25) {
m <- dim(A)[1]
n <- dim(A)[2]
tmp <- svd(A)

### sample a mode
j <- sample(x =  1:n, size = 1, prob = exp(tmp$d-max(tmp$d)))
v <-  (-1)^rbinom(1, 1, .5) * tmp$v[, j]

### sample around a mode
for(i in 1:nsamp) { u <- rvmf(A %*% v) ; v <- rvmf(t(A) %*% u) }

list(u = c(u), v = v)            }




####3


pK.X.bic <- function(X){
N <- dim(X)[2]
d <- dim(X)[1]
S <- t(X) %*% X / N
lam <- eigen(S)$val / N
lp <- NULL
for(k in 0:(d-1)) {
m <-  d * k - k * (k + 1) / 2
lp <- c(lp, -.5 * N * sum(log(lam[seq(1, k, length = k)])) - .5 * (m + k) * log(N) -
        .5 * N * (d-k) * log(mean(lam[seq(k + 1, d, length = d-k)])))
                         }
k <-  d
m <-  d * k - k * (k + 1) / 2
lp <- c(lp, -.5 * N * sum(log(lam[seq(1, k, length = k)])) - .5 * (m + k) * log(N))
pk <- exp(lp-max(lp))
pk / sum(pk)
}


pK.X.lap <- function(X){
N <- dim(X)[2]
d <- dim(X)[1]
S <- t(X) %*% X / N
lam <- eigen(S)$val / N
lp <-   NULL

for(k in 0:(d-1)) {

laz <- 0
for(i in seq(1, k, length = k)){
for(j in seq(i + 1, d, length = d-i)){
    laz <- laz + log(1 / lam[j]-1 / lam[i]) + log(lam[i]-lam[j]) + log(N) }}

lpu <- 0
if(k>0){
lpu <-  -k * log(2) + sum(lgamma((d-(1:k) + 1) / 2)  - .5 * (d-(1:k) + 1) * log(pi))
        }

m <-  d * k - k * (k + 1) / 2
lp <- c(lp, lpu -.5 * N * sum(log(lam[seq(1, k, length = k)])) -
               .5 * N * (d-k) * log(mean(lam[seq(k + 1, d, length = d-k)]))  +
               .5 * (m + k) * log(2 * pi) -.5 * laz -.5 * k * log(N) )
           }

k <- d
laz <- 0
for(i in seq(1, k, length = k)){
for(j in seq(i + 1, d, length = d-i)){
    laz <- laz + log(1 / lam[j]-1 / lam[i]) + log(lam[i]-lam[j]) + log(N) }}

lpu <- 0
if(k>0){
lpu <-  -k * log(2) + sum(lgamma((d-(1:k) + 1) / 2)  - .5 * (d-(1:k) + 1) * log(pi))
        }

m <-  d * k - k * (k + 1) / 2
lp <- c(lp, lpu -.5 * N * sum(log(lam[seq(1, k, length = k)]))  +
               .5 * (m + k) * log(2 * pi) -.5 * laz -.5 * k * log(N) )

pk <- exp(lp-max(lp))
pk / sum(pk)
}


##########
ln2moment <- function(mu, sigma, lmax) {
.C("ln2moment", as.double(abs(mu)), as.double(sigma^2), as.integer(lmax),
               l2mom = double(lmax + 1))$l2mom
                                  }



roots.cubic <- function(a2, a1, a0) {
Q <-  (3 * a1 -a2^2) / 9
R <-  (9 * a2 * a1-27 * a0-2 * a2^3) / 54
D <- Q^3 + R^2
S <- (R + sqrt(0i + D))^(1 / 3)
T <- (R-sqrt(0i + D))^(1 / 3)

c(-a2 / 3 + (S + T),
         -a2 / 3 -.5 * (S + T) + 1i * .5 * sqrt(3) * (S-T),
         -a2 / 3 -.5 * (S + T) - 1i * .5 * sqrt(3) * (S-T))
                                 }

roots.quartic <- function(a3, a2, a1, a0) {

z3 <- roots.cubic(-a2, a1 * a3-4 * a0, 4 * a2 * a0 -a1^2-a3^2 * a0)
y <-  Re(z3[ Im(z3) ==  0 ])[1]

R <-  sqrt(.25 * a3^2 -a2  + y)
D <-  sqrt(0i + .25 * 3 * a3^2 -R^2 -2 * a2 + .25 * (4 * a3 * a2-8 * a1-a3^3) / R)
E <-  sqrt(0i + .25 * 3 * a3^2 -R^2 -2 * a2-.25 * (4 * a3 * a2-8 * a1-a3^3) / R)

 c(-a3 / 4  + R / 2  + D / 2,
    -a3 / 4  + R / 2 -D / 2,
    -a3 / 4 -R / 2  + E / 2,
    -a3 / 4 -R / 2 -E / 2)
}

ldxl <- function(x, mu, sigma, l, ln2m =  ln2moment(mu, sigma, l)[l + 1]) {
l * log(x^2)  -log(sigma) -.5 * log(2 * pi) - .5 * ((x-mu) / sigma)^2 - ln2m }

rxl <- function(mu, sigma, l, nu = 1) {

theta <-  .5 * mu * (1 + sqrt(1 + 8 * l * sigma^2 / mu^2))
tau <-  1 / sqrt(1 / sigma^2 + 2 * l / (theta^2))

a <-  (-2 * theta-mu)
b <-  (theta^2  + 2 * mu * theta + nu * tau^2  + sigma^2 * (-nu-2 * l-1))
c <-  (-mu * theta^2 + (nu + 4 * l + 1) * sigma^2 * theta - mu * nu * tau^2)
d <-  (-2 * l * sigma^2  * theta^2 -2 * l * nu * sigma^2 * tau^2)

z4 <- roots.quartic(a, b, c, d)
xc <- Re(z4[Im(z4) == 0])

ln2m <- ln2moment(mu, sigma, l)[l + 1]

lM <-  max(ldxl(xc, mu, sigma, l, ln2m = ln2m)-(-log(tau) + log(dt((xc-theta) / tau, nu))))

samp <- T
while(samp) {
x <-  theta  + rt(1, nu) * tau
lrratio <-  ldxl(x, mu, sigma, l = l, ln2m = ln2m)-(-log(tau) + log(dt((x-theta) / tau, nu)))-lM
samp <-  (log(runif(1))>lrratio)
             }
x
}


