#' @title Signal subspace dimension estimation in high-dimensional matrix by double posterior.
#'
#' @description Estimate the dimension of a signal-rich subspace in large,
#'  high-dimensional data.
#'
#' @param s a subspace class.
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
#' estimate_rank_double_posterior_cor()
#'
#' str(results)
#' plot(results$Subspace, changepoint = results$dimension,
#'      annotation = 10)
#' modified_legacyplot(results$Changepoint$bcp_irl, annotation = 10)
#' modified_legacyplot(results$Changepoint$bcp_post, annotation = 10)
#' @seealso [RMTstat] for details of Marchenko-Pastur distribution.
#' @seealso https://dracodoc.wordpress.com/2014/07/21/
#' a-simple-algorithm-to-detect-flat-segments-in-noisy-signals/ for detection
#' of flat and spike in noisy signals
#' @importFrom bcp bcp
#' @importFrom  tibble tibble
#' @importFrom stringr str_locate
#' @export
estimate_rank_double_posterior_cor <- function(s, p, verbose, ...) {
  UseMethod("estimate_rank_double_posterior_cor", s)
}

#' @export
estimate_rank_double_posterior_cor.default <- function(s, p, verbose, ...) {
  stop("Don't know how to estimate the rank for an object of type ",
       paste(class(s), collapse = " "), ".")
}
#' @export
estimate_rank_double_posterior_cor.subspace <- function(s, p = 0.90, verbose = TRUE, ...) {  
  # -----------------------
  # Basic parameter set up
  # -----------------------
  rnk             <- max(s$components)
  sigma_a         <- s$sigma_a
  # --------------------------
  # Rank Estimation procedure
  # --------------------------
  if (rnk > 20) {
    #Trim unnecessary components when components large
    sigma_a_diff    <- diff(sigma_a / sigma_a[1])
    cutoff          <- min(min(sigma_a / sigma_a[1]), 0.005)
    m <- paste0("Cutoff value = ", cutoff, "\n")
    if (verbose) {
      message(m)
    }
    sigma_a_diff[abs(sigma_a_diff) < cutoff] <- 0
    sigma_a_diff[abs(sigma_a_diff) > cutoff] <- 1
    sigma_a_str     <- toString(sigma_a_diff)
    #Detect extreme cases with strange spikes
    flatstart       <- str_locate(sigma_a_str,
                                  paste0("1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ",
                                        "0, 0, 0, 0, 0, 0, 0, 0, 0, 0"))
    flatend         <- str_locate(sigma_a_str,
                                  paste0("0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ",
                                        "0, 0, 0, 0, 0, 0, 0, 0, 0, 1"))
    if (!is.na(sum(flatstart))) {
      cond_num      <- (flatstart[1] - 1) / 3 + 1
      cor_rnk   <- (flatstart[2] - 1) / 3 + 1
      m1 <- paste0("Detecting flat pattern from ",
            cond_num, " to ", cor_rnk,
            " and trim out components after ", cor_rnk, "\n")
      if (verbose) {
        message(m1)
      }
    } else if (!is.na(sum(flatend))) {
      spike_num     <- (flatend[2] - 1) / 3 + 1
      cor_rnk   <- (flatend[1] - 1) / 3 + 1
      m2 <- paste0("Detecting spike pattern from ",
            cor_rnk, " to ", spike_num,
            " and trim out components after ", cor_rnk, "\n")
      if (verbose) {
        message(m2)
      }
    } else {
      cor_rnk   <- rnk
    }
  } else {
      cor_rnk   <- rnk
      flatstart <- NA
  }

  #Bayesian Change Point
  bcp_irl     <- bcp(as.vector(sigma_a[seq_len(cor_rnk)]),
                     p0 = 0.1)
  #Bayesian Posterior Prob Change Point
  bcp_post    <- bcp(as.vector(c(bcp_irl$posterior.prob[-cor_rnk], 0)),
                     p0 = 0.1)

  prob_irl    <- c(bcp_irl$posterior.prob[-cor_rnk], 0)
  prob_post   <- c(bcp_post$posterior.prob[-cor_rnk], 0)

  #multiple changepoints in bcp_post
  #1. find all max changepoinst in bcp_post
  post_max    <- which(prob_post == max(prob_post))
  #2. Any of them in prob_irl >0.90 * irl_max
  irl_max     <- cor_rnk + 1 - which.max(rev(prob_irl))
  threshold   <- prob_irl[post_max] > p * max(prob_irl)
  #3. If none in 2. then choose irl_max
  changepoint <- ifelse(!is.na(sum(flatstart)) & irl_max == 1,
                        cond_num, switch(2 - any(threshold),
                                         max(post_max[threshold]),
                                         irl_max))
  m3 <- paste0("estimate_rank_double_posterior_cor estimation = ", changepoint, "\n")
  message(m3)
  ret <- list(subspace    = s,
              dimension   = changepoint,
              Changepoint = list(bcp_irl    = bcp_irl,
                                 bcp_post   = bcp_post),
              message    = list(m3))
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
      x$subspace$ndf,
      "samples and",
      x$subspace$pdim,
      "features.\n")
  cat(x$message[[1]])
}
