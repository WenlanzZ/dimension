#' @title Simulate X matrix
#'
#' @description These functions provide multivariate normal distributed X matrix.
#' @param n The number of rows for matrix X.
#' @param p The number of columns for matrix X.
#' @param ncc The number of correlated columns.
#' @param var Varaince.
#' @param seed The random number seed.
#' @examples
#' \donttest{
#' # 3 baskets, each with enrollement size 5
#' X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
#' }
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @export

Xsim <- function(n = 100,
                 p = 90,
                 ncc = 2,
                 var = 10,
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