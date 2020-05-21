#' @title Simulate x matrix
#'
#' @description These functions provide multivariate normal
#'  distributed x matrix.
#' @param n The number of rows for matrix x.
#' @param p The number of columns for matrix x.
#' @param ncc The number of correlated columns.
#' @param var Varaince.
#' @examples
#' x <- x_sim(n = 150, p = 100, ncc = 10, var = 2)
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @export

x_sim <- function(n = 100,
                  p = 90,
                  ncc = 2,
                  var = 10) {
  sigma <- matrix(rep(0, ncc * ncc), ncol = ncc)
  diag(sigma) <- var
  cor_cols <- rmvnorm(n, rep(0, ncc), sigma = sigma)
  xs <- cbind(cor_cols, matrix(0, ncol = p - ncc, nrow = n))
  xn <- matrix(rnorm(n * p), n, p)
  x <- xs + xn
  return(x)
}
