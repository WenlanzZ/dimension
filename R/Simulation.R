#' @title Gram-schmidt process
#'
#' @description This function uses Gram-schmidt process to generate orthorgonal basis for a vector space.
#' @param X a numeric matrix.
#' @examples
#' \donttest{
#' # 3 baskets, each with enrollement size 5
#' X <- Xsim(n=100,p=90,ncc=2,var=10)
#' gs1(X)
#' }
#' @export

gs1 = function(x){
  innovation.norms <- numeric(ncol(x))
  xx = x
  innovation.norms[1] <- sqrt(sum(x[,1]^2))
  xx[,1] = x[,1]/innovation.norms[1]
  r = ncol(x)
  for(j in 2:r){
    v = x[,j]
    for(k in 1:(j-1)){
      v = v - sum(x[,j]*xx[,k])*xx[,k]
    }
    innovation.norms[j] <- sqrt(sum(v^2))
    xx[,j] = v/innovation.norms[j]
  }
  return(list(ON.basis=xx, innovation.norms=innovation.norms))
}

#' @title simulate X matrix
#'
#' @description These functions provide multivariate normal distributed X matrix.
#' @param n number of rows for matrix X.
#' @param p number of columns for matrix X.
#' @param ncc number of correlated columns.
#' @param sigma covariance matrix.
#' @param orthogonl a logical value. Generate orthogonal noise if true.
#' @param seed the random number seed.
#' @examples
#' \donttest{
#' # 3 baskets, each with enrollement size 5
#' X <- Xsim(n=100,p=90,ncc=2,var=10)
#' }
#' @importFrom mvtnorm rmvnorm
#' @export

Xsim <- function(n = 100,
                        p = 90,
                        ncc = 2,
                        var = 10,
                        orthogonl = FALSE,
                        seed = 1000) {
  set.seed(seed)
  sigma <- matrix(rep(0,ncc*ncc), ncol = ncc)
  diag(sigma)<-var
  cor_cols <- rmvnorm(n, rep(0, ncc), sigma = sigma)
  ON<-gs1(cor_cols)$ON.basis
  xs <- cbind(ON, matrix(0, ncol=p-ncc, nrow = n))
  xn<-matrix(rnorm(n*p),n,p)
  if(orthogonl){
    xs <- cbind(ON, matrix(0, ncol=p-ncc, nrow = n))
    X <- xs + gs1(xn)$ON.basis
  }else{
    xs <- cbind(cor_cols, matrix(0, ncol=p-ncc, nrow = n))
    X <- xs + xn
  }
  return(X)
}

