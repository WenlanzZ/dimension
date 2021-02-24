#' @title Signal subspace dimension estimation in high-dimensional matrix
#'  by double posterior.
#' @description Estimate the dimension of a signal-rich subspace in large,
#'  high-dimensional data.
#'
#' @param s a subspace class.
#' @param p threshold for selecting changepoint. Default is 0.90.
#' @param verbose output message
#' @param ... Extra parameters
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of x.}
#'   \item{pdim:}{ The number of dimensions of x.}
#'   \item{components:}{ A series of right singular vectors estimated.}
#'   \item{var_correct:}{ Corrected population variance
#'    for Marcenko-Pastur distribution.}
#'   \item{transpose_flag:}{ A logical value indicating
#'    whether the matrix x is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for
#'    specified components and corresponding dimensions.}
#'   \item{sigma_a:}{ A vector of corrected eigenvalues up to max(components).}
#'   \item{mp_irl:}{ A data frame of sampled expected eigenvalues from
#'    Marcenko-Pastur for specified components and corresponding dimensions.}
#'   \item{sigma_mp:}{ A vector of samped expected eigenvalues from
#'    Marcenko-Pastur up to max(components).}
#'   \item{v:}{ Right singular vectors of x matrix for specified components.}
#'   \item{u:}{ Left singular vectors of x matrix or specified components.}
#'   \item{dimension:}{ Estimated signal subspace dimension.}
#'   \item{bcp_irl:}{ Probability of change in mean and posterior means
#'    of eigenvalue difference between $x$ and $N$.}
#'   \item{bcp_post:}{ Probability of change in mean and posterior means
#'    of bcp_rirl.}
#' }
#' @section Details:
#'  We estimate the intrinsic dimension of a signal-rich subspace
#'  in large high-dimensional data by decomposing matrix into a
#'  signal-plus-noise space and approximate the signal-rich subspace
#'  with a rank K approximation \eqn{\hat{x}=\sum_{k=1}^{K}d_ku_k{v_k}^T}.
#'  To estimate rank K, we propose a simple procedure assuming that matrix
#'  x is composed of a low-rank signal matrix S and an average general noise
#'  random matrix \eqn{\bar{N}}. It has been shown that the average eigenvalues
#'  of random matrices N follows a universal Marchenko-Pastur (MP) distribution.
#'  We hypothesize that the deviation of eigenvalues of x from the MP
#'  distribution indicates the intrinsic dimension of signal-rich subspace.
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' results <- x %>%
#' create_subspace(components = 8:30) %>%
#' correct_eigenvalues() %>%
#' estimate_rank_double_posterior()
#'
#' str(results)
#' plot(results$subspace, changepoint = results$dimension,
#'      annotation = 10)
#' modified_legacyplot(results$bcp_irl, annotation = 10)
#' modified_legacyplot(results$bcp_post, annotation = 10)
#' @seealso [RMTstat] for details of Marchenko-Pastur distribution.
#' @seealso https://dracodoc.wordpress.com/2014/07/21/
#' a-simple-algorithm-to-detect-flat-segments-in-noisy-signals/ for detection
#' of flat and spike in noisy signals
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @export
estimate_rank_double_posterior <- function(s, p, verbose, ...) {
  UseMethod("estimate_rank_double_posterior", s)
}

#' @export
estimate_rank_double_posterior.default <- function(s, p, verbose, ...) {
  stop("Don't know how to estimate the rank for an object of type ",
       paste(class(s), collapse = " "), ".")
}

#' @export
estimate_rank_double_posterior.subspace <- function(s, p = 0.90,
                                                    verbose = FALSE, ...) {

  # -----------------------
  # Basic parameter set up
  # -----------------------
  rnk             <- max(s$components)
  sigma_a         <- s$sigma_a

  #Bayesian Change Point
  bcp_irl     <- bcp(as.vector(sigma_a[seq_len(rnk)]), p0 = 0.1)
  #Bayesian Posterior Prob Change Point
  bcp_post    <- bcp(as.vector(c(bcp_irl$posterior.prob[-rnk], 0)),
                     p0 = 0.1)

  prob_irl    <- c(bcp_irl$posterior.prob[-rnk], 0)
  prob_post   <- c(bcp_post$posterior.prob[-rnk], 0)

  #multiple changepoints in bcp_post
  #1. find all max changepoinst in bcp_post
  post_max    <- which(prob_post == max(prob_post))
  #2. Any of them in prob_irl >0.90 * irl_max
  irl_max     <- rnk + 1 - which.max(rev(prob_irl))
  threshold   <- prob_irl[post_max] > p * max(prob_irl)
  #3. If none in 2. then choose irl_max
  changepoint <- switch(2 - any(threshold),
                       max(post_max[threshold]),
                       irl_max)
  m3 <- paste0("dimension estimation = ", changepoint, ".\n")
  message(m3)
  ret <- list(subspace    = s,
              dimension   = changepoint,
              bcp_irl     = bcp_irl,
              bcp_post    = bcp_post,
              message     = list(m3))
  attr(ret, "class") <- "dimension"
  ret
}
