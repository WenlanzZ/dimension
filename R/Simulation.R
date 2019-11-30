#' @title simulate X matrix
#'
#' @description These functions provide multivariate normal distributed X matrix.
#' @param n number of rows for matrix X.
#' @param p number of columns for matrix X.
#' @param ncc number of correlated columns.
#' @param var varaince.
#' @param fact scale on basis if orthogonal is true.
#' @param seed the random number seed.
#' @examples
#' \donttest{
#' # 3 baskets, each with enrollement size 5
#' X <- Xsim(n=100,p=90,ncc=2,var=10,fact = 30)
#' }
#' @importFrom mvtnorm rmvnorm
#' @export

Xsim <- function(n = 100,
                 p = 90,
                 ncc = 2,
                 var = 10,
                 fact = 30,
                 seed = 1000) {
  set.seed(seed)
  sigma <- matrix(rep(0, ncc * ncc), ncol = ncc)
  diag(sigma) <- var
  cor_cols <- rmvnorm(n, rep(0, ncc), sigma = sigma)
  xs <- cbind(cor_cols, matrix(0, ncol = p-ncc, nrow = n))
  xn <- matrix(rnorm(n * p), n, p)
  X <- xs + xn
  return(X)
}