#' @title ladle
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features (n > p).
#' @param nboot bootstrap sample size.

#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{fn0:}{ fn0 is the bootstrap eigenvector variability.}
#'   \item{fn:}{ fn is the renormalized fn0.}
#'   \item{lam:}{ lam is the sample eigenvalues.}
#'   \item{phin:}{ phin is the renormalized and shifted sample eigenvalues.}
#'   \item{gn:}{ gn is the objective function ("y" in ladle plot).}
#'   \item{d:}{ d is the ladle estimator.}
#' }
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' x %>% estimate_rank_ladle()
#' @importFrom stats var
#' @importFrom dr dr
#' @export
estimate_rank_ladle <- function(x) {
  p <- ncol(x)
  n <- nrow(x)
  if (p > 10) {
    r <- as.integer(p / log(p))
  } else {
      r <- p - 1
  }
  obs <- x
  mhat <- var(x)
  bhat <- eigen(mhat)$vectors[, 1:r]
  lam <- eigen(mhat)$values[1:(r + 1)]
  nboot <- as.integer(n / 2)

  fn0 <- rep(0, r + 1)
  for (j in 1:nboot) {
    u <- round(runif(n, min = -0.5, max = n + 0.5))
    bs <- obs[u, ]
    mstar <- var(bs)
    b_star <- eigen(mstar)$vectors[, 1:r]

    for (i in 1:r) {
      fn0[i + 1] <- fn0[i + 1] + 1 - abs(det(t(b_star[, 1:i]) %*% bhat[, 1:i]))
    }
  }
  fn0 <- fn0 / nboot
  res <- list(0)
  res$fn0 <- fn0
  res$fn <- fn0 / (1 + sum(fn0))
  res$lam <- lam
  res$phin <- lam / (1 + sum(lam))
  res$gn <- res$fn + res$phin
  res$d <- which.min(res$gn) - 1
  m3 <- paste0("dimension estimation = ", res$d, "\n")
  message(m3)
  return(res)
}