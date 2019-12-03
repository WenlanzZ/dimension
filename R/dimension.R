#' @title Signal subspace dimension estimation in high-dimensional matrix
#'
#' Estimate the dimension of a signal-rich subspace in large high-dimensional data.
#'
#' @param X A numeric real-valued matrix with n number of samples and p number of features. 
#'   If p>n, a warning message is generated and the transpose of X is used.
#' @param subspace_ A subspace class.
#' @param rank A series of right singular vectors to estimate. rank must be smaller or equal to min(nrow(X),ncol(X)).
#' @param basis Choose eigenvalue decomposition or singular value decomposition.
#' @param times Split data into X times for parallel computation.
#' @param p Threshold for selecting changepoint. default is 0.90.
#' @param verbose output message
#' @param ... Extra parameters
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of X.}
#'   \item{pdim:}{ The number of dimensions of X.}
#'   \item{rank:}{ A series of right singular vectors estimated.}
#'   \item{var_correct:}{ Corrected population variance for Marcenko-Pastur distribution.}
#'   \item{transpose_flag:}{ A logical value indicating whether the matrix X is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for specified rank and corresponding dimensions.}
#'   \item{MP_irl:}{ A data frame of samped expected eigenvalues from Marcenko-Pastur for specified rank and corresponding dimensions.}
#'   \item{v:}{ Right singular vectors of X matrix for specified rank.}
#'   \item{u:}{ Left singular vectors of X matrix or specified rank.}
#'   \item{bcp_irl:}{ Probability of change in mean and posterior means of eigenvalue difference between $X$ and $N$.}
#'   \item{bcp_post:}{ Probability of change in mean and posterior means of bcp_rirl.}
#'   \item{changePoint:}{ Estimated signal subspace dimension.}
#' }
#' @section Details:
#'   We estimate the intrinsic dimension of a signal-rich subspace in large high-dimensional data by decomposing matrix
#'   into a signal-plus-noise space and approximate the signal-rich subspace with a rank K approximation
#'   \eqn{\hat{X}=\sum_{k=1}^{K}d_ku_k{v_k}^T}. To estimate rank K, we propose a simple procedure assuming that matrix X is composed
#'   of a low-rank signal matrix S and an average general noise random matrix \eqn{\bar{N}}. It has been shown that
#'   the average eigenvalues of random matrices N follows a universal Marcenko-Pastur (MP) distribution.
#'   We hypothesize that the deviation of eigenvalues of X from the MP distribution indicates the intrinsic dimension of signal-rich subspace.
#' @examples
#' \donttest{
#' X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
#' results <- dimension(X, rank = 1:40, times = 10, basis="eigen")
#' 
#' #equivelantly, if subsapce is calcualted
#' Subspace <- subspace(X, rank = 1:40, times = 10,  basis = "eigen")
#' results <- dimension(subspace_ = Subspace)
#' 
#' str(results)
#' plot(results$Subspace, Changepoint = results$Changepoint$dimension, annotation = 10)
#' modified_legacyplot(results$Changepoint$bcp_irl, annotation = 10)
#' modified_legacyplot(results$Changepoint$bcp_post, annotation = 10)
#' }
#' @seealso [RMTstat] for details of Marcenko-Pastur distribution.
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @export

dimension <- function(X,                          # data matrix
                      subspace_ = NULL,               # subspace class
                      rank = NA,                      # number of singular vectors to estimate
                      basis = c("eigen","singular"),  # 
                      times = NA,                     # split data into X times for parallel computation.
                      p = 0.90,                       # threshold for selecting changepoint
                      verbose = TRUE,                 # output message
                      ...)
{
# ---------------------------------------------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------------------------------------------
  if (!missing(subspace_) & !is.null(subspace_$MP_irl)) {
    if (verbose) {
        cat("Subspace have already been calculated.\n")
    }
  } else {
    if (is.null(X)) {
      stop("Invalid input X")
    } else {
      # Checking for rank input
      if (missing(rank)) {
        rank <- 1:min(nrow(X), ncol(X))
        if (verbose) {
          cat("No rank specified. Calculating full singular value decomposition instead.\n")
        }
      } 
      # Checking for times input
      if (missing(times)) {
      times <- 0
      }
      # Calcualte subspace
      subspace_ <- subspace(X, rank = rank, times = times, basis = basis)
      cat("Finish calculating subpsace.\n")
    }
  }
# ---------------------------------------------------------------------------------------------------------
# Basic parameter set up
# ---------------------------------------------------------------------------------------------------------  
  ndf             <- subspace_$ndf
  pdim            <- subspace_$pdim
  rank            <- subspace_$rank
  rnk             <- max(rank)
  var_correct     <- subspace_$var_correct
  sigma_a         <- subspace_$sigma_a
  sigma_MP        <- subspace_$sigma_MP
# ---------------------------------------------------------------------------------------------------------
# Rank Estimation procedure
# ---------------------------------------------------------------------------------------------------------  
  #Bayesian Change Point
  bcp_irl   <- bcp(as.vector(sigma_a-sigma_MP), p0 = 0.1)
  #Bayesian Posterior Prob Change Point
  bcp_post  <- bcp(as.vector(c(bcp_irl$posterior.prob[-rnk], 0)), p0 = 0.1)

  prob_post <- c(bcp_post$posterior.prob[-rnk], 0)
  prob_irl  <- c(bcp_irl$posterior.prob[-rnk], 0)
  post_max  <- which(prob_post == max(prob_post))
  num_max   <- length(post_max)

  if (num_max == 1) {
    #unimodal changepoint in bcp_post
    changePoint <- rnk + 1 - which.max(rev(prob_post))
  } else {
    #multiple changepoints in bcp_post
    #1. find all max changepoinst in bcp_post 
    #2. Any of them in prob_irl >0.90*max(prob_irl)
    irl_max   <- rnk + 1 - which.max(rev(prob_irl))
    threshold <- prob_irl[post_max] > p*max(prob_irl)
    #3. If none in 2. then choose irl_max
    changePoint <- switch(2 - any(threshold), max(post_max[threshold]), irl_max)
  }

  return(list(Subspace    = subspace_,
              Changepoint = list(bcp_irl  = bcp_irl,
                                 bcp_post = bcp_post,
                                 dimension = changePoint)))
}