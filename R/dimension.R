#' @title Signal subspace dimension estimation in high-dimensional matrix
#'
#' @description Estimate the dimension of a signal-rich subspace in large,
#'  high-dimensional data.
#'
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p>n, a warning message is generated and
#'  the transpose of x is used.
#' @param subspace_ A subspace class.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x),ncol(x)).
#' @param times Split data into times-fold for parallel computation.
#' @param p Threshold for selecting changepoint. Default is 0.90.
#' @param verbose output message
#' @param ... Extra parameters
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of x.}
#'   \item{pdim:}{ The number of dimensions of x.}
#'   \item{components:}{ A series of right singular
#'    vectors estimated.}
#'   \item{var_correct:}{ Corrected population variance
#'    for Marcenko-Pastur distribution.}
#'   \item{transpose_flag:}{ A logical value indicating
#'    whether the matrix x is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for
#'    specified rank and corresponding dimensions.}
#'   \item{mp_irl:}{ A data frame of sampled expected eigenvalues from
#'    Marcenko-Pastur for specified rank and corresponding dimensions.}
#'   \item{v:}{ Right singular vectors of x matrix for specified rank.}
#'   \item{u:}{ Left singular vectors of x matrix or specified rank.}
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
#'  of random matrices N follows a universal Marcenko-Pastur (MP) distribution.
#'  We hypothesize that the deviation of eigenvalues of x from the MP
#'  distribution indicates the intrinsic dimension of signal-rich subspace.
#' @examples
#' \donttest{
#' x <- x_sim(n = 150, p = 100, ncc = 10, var = 2)
#' results <- dimension(x, components = 1:40, times = 10)
#'
#' #equivelantly, if subsapce is calcualted
#' Subspace <- subspace(x, components = 1:40, times = 10)
#' results <- dimension(subspace_ = Subspace)
#'
#' str(results)
#' plot(results$Subspace,
#'      Changepoint = results$dimension,
#'      annotation = 10)
#' modified_legacyplot(results$Changepoint$bcp_irl, annotation = 10)
#' modified_legacyplot(results$Changepoint$bcp_post, annotation = 10)
#' }
#' @seealso [RMTstat] for details of Marcenko-Pastur distribution.
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @export

dimension <- function(x,
                      subspace_ = NULL,
                      components = NA,
                      times = NA,
                      p = 0.90,
                      verbose = TRUE,
                      ...) {
# -----------------------
# Check input parameters
# -----------------------
  if (!missing(subspace_) & !is.null(subspace_$mp_irl)) {
    if (verbose) {
        cat("Subspace have already been calculated.\n")
    }
  } else {
    if (is.null(x)) {
      stop("Invalid input x")
    } else {
      # Checking for rank input
      if (missing(components)) {
        components <- 1:min(nrow(x), ncol(x))
        if (verbose) {
          cat("No rank specified.
              Calculating full singular value decomposition instead.\n")
        }
      }
      # Checking for times input
      if (missing(times)) {
      times <- 0
      }
      # Calcualte subspace
      subspace_ <- subspace(x, components = components, times = times)
      cat("Finish calculating subpsace.\n")
    }
  }
# -----------------------
# Basic parameter set up
# -----------------------
  rnk             <- max(subspace_$components)
  sigma_a         <- subspace_$sigma_a
  sigma_mp        <- subspace_$sigma_mp
# --------------------------
# Rank Estimation procedure
# --------------------------
  #Bayesian Change Point
  bcp_irl   <- bcp(as.vector(sigma_a - sigma_mp), p0 = 0.1)
  #Bayesian Posterior Prob Change Point
  bcp_post  <- bcp(as.vector(c(bcp_irl$posterior.prob[-rnk], 0)), p0 = 0.1)

  prob_post <- c(bcp_post$posterior.prob[-rnk], 0)
  prob_irl  <- c(bcp_irl$posterior.prob[-rnk], 0)
  post_max  <- which(prob_post == max(prob_post))
  num_max   <- length(post_max)

  if (num_max == 1) {
    #unimodal changepoint in bcp_post
    changepoint <- rnk + 1 - which.max(rev(prob_post))
  } else {
    #multiple changepoints in bcp_post
    #1. find all max changepoinst in bcp_post
    #2. Any of them in prob_irl >0.90*max(prob_irl)
    irl_max     <- rnk + 1 - which.max(rev(prob_irl))
    threshold   <- prob_irl[post_max] > p * max(prob_irl)
    #3. If none in 2. then choose irl_max
    changepoint <- switch(2 - any(threshold), max(post_max[threshold]), irl_max)
  }

  return(list(Subspace    = subspace_,
              dimension  = changepoint,
              Changepoint = list(bcp_irl    = bcp_irl,
                                 bcp_post   = bcp_post)))
}
