#' @title Gram-schmidt process
#'
#' @description This function uses Gram-schmidt process to generate orthorgonal basis for a vector space.
#' @param X a numeric matrix.
#' @examples
#' \donttest{
#' # 3 baskets, each with enrollement size 5
#' X <- Xsim(n=100,p=90,ncc=2,var=10,fact = 30)
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
