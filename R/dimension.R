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
#'    for Marchenko-Pastur distribution.}
#'   \item{transpose_flag:}{ A logical value indicating
#'    whether the matrix x is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for
#'    specified rank and corresponding dimensions.}
#'   \item{mp_irl:}{ A data frame of sampled expected eigenvalues from
#'    Marchenko-Pastur for specified rank and corresponding dimensions.}
#'   \item{v:}{ Right singular vectors of x matrix for specified rank.}
#'   \item{u:}{ Left singular vectors of x matrix or specified rank.}
#'   \item{dimension:}{ Estimated signal subspace dimension.}
#'   \item{bcp_irl:}{ Probability of change in mean and posterior means
#'    of eigenvalue difference between $x$ and $N$.}
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
#' results <- dimension(x, components = 1:50)
#'
#' #equivelantly, if subsapce is calcualted
#' Subspace <- subspace(x, components = 1:50)
#' results <- dimension(subspace_ = Subspace)
#'
#' str(results)
#' plot(results$Subspace, changepoint = results$dimension,
#'      annotation = 10)
#' modified_legacyplot(results$bcp_irl, annotation = 10)
#' @seealso [RMTstat] for details of Marchenko-Pastur distribution.
#' @seealso https://dracodoc.wordpress.com/2014/07/21/
#' a-simple-algorithm-to-detect-flat-segments-in-noisy-signals/ for detection
#' of flat and spike in noisy signals
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @importFrom stringr str_locate
#' @export
dimension <- function(x,
                      subspace_ = NULL,
                      components = NA,
                      times = NA,
                      verbose = FALSE,
                      ...) {
# -----------------------
# Check input parameters
# -----------------------
  if (!missing(subspace_) && !is.null(subspace_$mp_irl)) {
    if (verbose) {
        message("Subspace has already been calculated.\n")
    }
  } else {
    if (is.null(x)) {
      stop("Invalid input x")
    } else {
      # Checking for rank input
      if (missing(components)) {
        components <- seq_len(min(nrow(x), ncol(x)))
        if (verbose) {
          message(paste0("No component specified. ",
              "Calculating full singular value decomposition instead.\n"))
        }
      }
      # Checking for times input
      if (missing(times)) {
        times <- 0
      }
      # Calcualte subspace
      subspace_ <- subspace(x, components = components, times = times)
    }
  }
# -----------------------
# Basic parameter set up
# -----------------------
  rnk             <- max(subspace_$components)
  sigma_a         <- subspace_$sigma_a
  sigma_mp        <- subspace_$sigma_mp
  sigma_a_adj     <- sigma_a - sigma_mp
# --------------------------
# Rank Estimation procedure
# --------------------------
  #Bayesian Change Point
  prob_prior <- (seq(0.9, 0, length.out = rnk))
  suppressWarnings(
    bcp_irl  <- bcp(as.vector(sigma_a_adj[1:rnk]),
                     p0 = prob_prior[1:rnk]))
  prob_irl   <- c(bcp_irl$posterior.prob[-rnk], 0)
  prob_irl_diff <- abs(diff(prob_irl))
  irl_max        <- rnk + 1 - which.max(rev(prob_irl[1:rnk]))
  if (verbose) {
    cat("irl_max = ", irl_max, "\n")
  }
  data <- tibble(diff = prob_irl_diff,
                 prob = prob_irl[-rnk])
  within_var <- km(data)
  changepoint <- which.min(within_var[,1])
  m3 <- paste0("Dimension estimation = ", changepoint, "\n")
  if (verbose) {
    message(m3)
  }

  ret <- list(Subspace    = subspace_,
              dimension   = changepoint,
              bcp_irl     = bcp_irl,
              data        = data,
              within_var  = within_var,
              message     = list(m3))
  attr(ret, "class") <- "dimension"
  ret
}

#' @title Print dimension
#'
#' A generic function.
#'
#' @param x a dimension class.
#' @param ... Extra parameters
#' @export
print.dimension <- function(x, ...) {
  cat("An object of class dimension estimated for",
      ifelse(x$transpose_flag, "transposed", ""),
      "X matrix with",
      x$Subspace$ndf,
      "samples and",
      x$Subspace$pdim,
      "features.\n")
  cat(x$message[[1]])
}
