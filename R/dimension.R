#' @title Rank Estimation in High-Dimensional matrix X.
#'
#' @description Estimate a low rank approximation of a signal-rich subspace in large high-dimensional data.
#' @param X A numeric real- or real-valued sparse matrix with n number of samples and p number of features. If p > n, a warning message is generated and the transpose of X is used.
#' @param MPSamples A list. Samples from Marc\u{e}nko-Pastur (MP) distribution and eigenvalues of X. Output from MarcenkoPasturSample function.
#' @param rnk number of right singular vectors to estimate. rnk must be smaller or equal to max(nrow(X),ncol(X)).
#' @param times split data into X times for parallel computation.
#' @param p threshold for selecting changepoint.default is 0.90.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ number of degrees of freedom of X and the white Wishart matrix.}
#'   \item{pdim:}{ number of dimensions of X and the white Wishart matrix.}
#'   \item{rnk:}{ number of right singular vectors to estimate.}
#'   \item{var_correct:}{ population variance for Marcenko-Pastur distribution.}
#'   \item{transpose_flag:}{ whether the matrix X is transposed.}
#'   \item{irl:}{ a data frame of scaled eigenvalues and corresponding dimensions.}
#'   \item{MP_irl:}{ a data frame of samped expected eigenvalues from Marcenko-Pastur and corresponding dimensions.}
#'   \item{v:}{ right singular vectors of X matrix with truncation up to dimension rnk.}
#'   \item{u:}{ left singular vectors of X matrix with truncation up to dimension rnk.}
#'   \item{bcp_irl:}{ probability of change in mean and posterior means of eigenvalue difference between $X$ and $N$.}
#'   \item{bcp_post:}{ probability of change in mean and posterior means of bcp_rirl.}
#'   \item{changePoint:}{ estimated changepoint position by from bcp_post.}
#' }
#' @section Details:
#' We estimate a low rank approximation of a signal-rich subspace in large high-dimensional data by decomposing matrix
#' into a signal-plus-noise space and approximate the signal-rich subspace with a rank K approximation
#' \eqn{\hat{X}=\sum_{k=1}^{K}d_ku_k{v_k}^T}. To estimate rank K, we propose a simple procedure assuming that matrix X is composed
#' of a low-rank signal matrix S and an average general noise random matrix \eqn{\bar{N}}. It has been shown that
#' the average eigenvalues of random matrices N follows a universal Marc\u{e}nko-Pastur (MP) distribution.
#' We hypothesize that the deviation of eigenvalues of X from the MP distribution indicates the intrinsic dimension of signal-rich subspace.
#' @examples
#' \donttest{
#' X <- Xsim(n=1000,p=500,ncc=10,var=2,fact = 1,orthogonl = FALSE)
#' results<-dimension(X,rnk=10,times=100)
#' 
#' #equivelantly, if MarcenkoPasturSample is calcualted
#' MPSamples<-MarcenkoPasturSample(X,rnk=40,times=100)
#' results = dimension(MPSamples=MPSamples)
#' str(results)
#' ScreePlot(results$MarcenkoPasturSample,Changepoint=results$Changepoint$changePoint,annotation=10)
#' modified_legacyplot(results$Changepoint$bcp_post,annotation=10)
#' }
#' should import RMTstat after fix bug
#' @seealso \code{\link[RMTstat]} for details of Marcenko-Pastur distribution.
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @export

dimension <- function(X = X,                          # data matrix
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
  if (!missing(subspace_)) {
    cat("Subspace have already been calculated.\n")
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

  return(list(subspace    = subspace_,
              Changepoint = list(bcp_irl  = bcp_irl,
                                 bcp_post = bcp_post,
                                 changePoint = changePoint)))
}