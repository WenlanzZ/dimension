#' @title Signal subspace dimension estimation in high-dimensional matrix
#'
#' @description Estimate the dimension of a signal-rich subspace in large,
#'  high-dimensional data.
#'
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p > n, a warning message is generated and
#'  the transpose of x is used.
#' @param s A subspace class.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x),ncol(x)).
#' @param decomposition The method to be used; method = "svd"
#'  returns results from singular value decomposition; method = "eigen"
#'  returns results from eigenvalue decomposition.
#' @param method The method to be used; method = "double_posterior"
#'  returns results from function estimate_rank_double_posterior;
#'   method = "posterior" returns results from function estimate_rank_posterior;
#'   method = "kmeans" returns results from function estimate_rank_kmeans;
#'   method = "ladle" returns results from function estimate_rank_ladle.
#'   Default uses estimate_rank_double_posterior.
#' @param num_est_samples Split data into num_est_samples-fold
#'  for parallel computation.
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
#' results <- dimension(s = Subspace, method = "double_posterior")
#'
#' str(results)
#' plot(results$subspace, changepoint = results$dimension,
#'      annotation = 10)
#' modified_legacyplot(results$bcp_irl, annotation = 10)
#' @seealso [RMTstat] for details of Marchenko-Pastur distribution.
#' @seealso https://dracodoc.wordpress.com/2014/07/21/
#' a-simple-algorithm-to-detect-flat-segments-in-noisy-signals/ for detection
#' of flat and spike in noisy signals
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @importFrom stringr str_locate
#' @importFrom dplyr %>%
#' @export
dimension <- function(x,
                      s = NULL,
                      components = NA,
                      decomposition = c("svd", "eigen"), 
                      method = c("double_posterior", "posterior",
                                 "kmeans", "ladle"),
                      num_est_samples = NA,
                      verbose = FALSE,
                      ...) {
# -----------------------
# Check input parameters
# -----------------------
  if (!missing(s)) {
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
      if (missing(num_est_samples)) {
        num_est_samples <- 0
      }
      if (missing(decomposition)) {
        decomposition <- "svd"
      }
      # Calcualte subspace
      s <- subspace(x, components = components, decomposition = decomposition, 
                    num_est_samples = num_est_samples, mp = TRUE, verbose)
    }
  }
  if (missing(method)) {
    method <- "double_posterior"
    }
  switch(method,
        default = {
          s %>% estimate_rank_double_posterior()
        },
        double_posterior = {
          s %>% estimate_rank_double_posterior()
        },
        posterior = {
          s %>% estimate_rank_posterior()
        },
        kmeans = {
          s %>% estimate_rank_kmeans()
        },
        ladle = {
          x %>% estimate_rank_ladle()
        },
        stop("Invalid method input")
        )
}

#' @title Print dimension
#'
#' A generic function.
#'
#' @param x a dimension class.
#' @param ... Extra parameters
#' @export
print.dimension <- function(x, ...) {
  message("An object of class dimension estimated for",
      ifelse(x$transpose_flag, "transposed", ""),
      "X matrix with",
      x$Subspace$ndf,
      "samples and",
      x$Subspace$pdim,
      "features.\n")
  message(x$message[[1]])
  invisible(x)
}
