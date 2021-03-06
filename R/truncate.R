#' @title Eigenvalue truncation
#'
#' @description This function truncates the scaled eigenvalues
#'  of x in order to provide a denoised estimator e_denoised
#'  of the underlying correlation matrix. There are three methods
#' to chose from. The threshold method returns $(1-alpha)*rnk$
#' proportion of eigenvalues above threshold; the hard method returns
#' all the empirical eigenvalues greater than the upper limit of the
#' support to the Marchenko-Pastur spectrum; the identity method
#' returns eigenvalues specified in a location vector. While we
#' keep a proportion of eigenvalues, we can either shrink the
#' remaining ones by a trace-preserving constant (i.e.
#' $Tr(E\_denoised) = Tr(E)$) or set them all to zero.
#' This function is adapted from "Python for Random Matrix Theory"
#' credit to J.-P. Bouchaud and M. Potters.
#' @param x A subspace class or a numeric real-valued matrix with n number
#'  of samples and p number of features. If p>n,
#'  a warning message is generated and the transpose of
#'  x is used.
#' @param components A series of right singular vectors
#' to estimate. Components must be smaller or equal
#' to min(nrow(x),ncol(x)).
#' @param method The method to be used; method = "threshold"
#'  returns (1-alpha)*rnk proportion of eigenvalues above threshold;
#'   method = "hard" returns all the empirical eigenvalues greater
#'    than the upper limit of the support to the Marcenko-Pastur spectrum;
#'   method = "identity" returns eigenvalues specified in location vector.
#' @param alpha Determine the fraction to keep of the top
#'  eigenvalues in threshold method.
#' @param location Indicate the location of eigenvalues
#'  to keep in identity method.
#' @param zeroout A logical value to zero out eigenvalues
#'  when truncating. default is to set to FALSE.
#' @param verbose output message
#' @param ... Extra parameters
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{eigen_denoised:}{ A denoised estimator of eigenvalues
#'    through a simple eigenvalue truncation procedure (cf. reference below).}
#'   \item{x_denoised:}{ A denoised estimator of the true sigmal
#'    matrix underlying a noisy high dimensional matrix x.}
#'   \item{v_denoised:}{ A denoised estimator of the true covariance
#'    matrix underlying a noisy, in-sample estimate empirical
#'    correaltion matrix estimated from x (cf. reference below).}
#'   \item{e_denoised:}{ A denoised estimator of the true correlation
#'    matrix underlying a noisy, in-sample estimate empirical
#'    correaltion matrix estimated from x (cf. reference below).}
#'   \item{v:}{ Right singular vectors of x matrix for specified rank.}
#'   \item{u:}{ Left singular vectors of x matrix or specified rank.}
#' }
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' x_denoised <- truncate(x,
#'                        components = 20,
#'                        method = "threshold",
#'                        alpha = 0.9,
#'                        zeroout = TRUE)
#' x_denoised <- truncate(x,
#'                        components = 20,
#'                        method = "hard",
#'                        zeroout = FALSE)
#' x_denoised <- truncate(x,
#'                        components = 20,
#'                        method = "identity",
#'                        location = 1:15,
#'                        zeroout = FALSE)
#'
#' # equivalently, if subspace is calculated
#' Subspace  <- subspace(x,
#'                       components = 1:20)
#' x_denoised <- truncate(Subspace,
#'                        location = 1:15,
#'                        zeroout = FALSE)
#' @importFrom stats na.omit
#' @seealso
#' * [MarchenkoPasturPar()] calculates upper and lower
#' limits of Marcenko-Pastur distribution from RMTstat package.
#' @export

truncate <- function(x, components, method, alpha,
                     location, zeroout, verbose, ...) {
  UseMethod("truncate", x)
}

truncate.default <- function(x, components, method, alpha,
                             location, zeroout, verbose, ...) {
  stop("Don't know how to truncate from a class of type: ",
       class(x))
}

truncate.matrix <- function(x,
                            components = NA,
                            method = c("threshold", "hard", "identity"),
                            alpha = NA, location = NA, zeroout = FALSE,
                            verbose = TRUE, ...) {
  truncate(x = x, components = components, method = method, alpha = alpha,
           location = location, zeroout = zeroout, verbose = verbose, ... = ...)
}

truncate.Matrix <- function(x,
                            components = NA,
                            method = c("threshold", "hard", "identity"),
                            alpha = NA, location = NA, zeroout = FALSE,
                            verbose = TRUE, ...) {
  truncate(x = x, components = components, method = method, alpha = alpha,
           location = location, zeroout = zeroout, verbose = verbose, ... = ...)
}

truncate.subspace <- function(x,
                            components = NA,
                            method = c("threshold", "hard", "identity"),
                            alpha = NA, location = NA, zeroout = FALSE,
                            verbose = TRUE, ...) {
  truncate(x = x, components = components, method = method, alpha = alpha,
           location = location, zeroout = zeroout, verbose = verbose, ... = ...)
}

truncate <- function(x,
                    components = NA,
                    method = c("threshold", "hard", "identity"),
                    alpha = NA,
                    location = NA,
                    zeroout = FALSE,
                    verbose = TRUE,
                    ...) {
# ------------------------
# Check input parameters
# ------------------------

  if (any(class(x) == "subspace")) {
    if (verbose) {
        message("Subspace has already been calculated.\n")
    }
    s <- x
  } else {
    if (is.null(x)) {
      stop("Invalid input x")
    } else {
      if (missing(components)) {
        components <- 1:min(nrow(x), ncol(x))
        if (verbose) {
          message("No component specified.
              Calculating full singular value decomposition instead.\n")
        }
      }
      s <- subspace(x, components = components, mp = FALSE)
    }
  }
# -----------------------
# Basic parameter set up
# -----------------------
    ndf            <- s$ndf
    pdim           <- s$pdim
    svr            <- ndf / pdim
    components     <- s$components
    rnk            <- max(components)
    irl            <- s$irl
    v              <- s$v
    u              <- s$u
# ----------------------------------
# Determine eigenvalues to preserve
# ----------------------------------
    if (missing(method)) {
      method <- "identity"
    }
    switch(method,
        hard = {
          lambda_min <- marcenko_pastur_par(ndf, pdim, var = 1, svr = svr)$lower
          lambda_max <- marcenko_pastur_par(ndf, pdim, var = 1, svr = svr)$upper
          m <- paste0("Use method ", method, "\n",
                      "lambda_min = ", round(lambda_min, 4), "\n",
                      "lambda_max = ", round(lambda_max, 4), "\n")
          if (verbose) cat(m)
          eigen_denoised <- ifelse(irl$eigen >= lambda_max, irl$eigen, NA)
        },
        threshold = {
          if (missing(alpha)) {
            stop("Alpha must be specified")
          } else if (alpha <= 0) {
            stop("Alpha must be positive")
          } else if (alpha > 1) {
            stop("Alpha must be less or equal to 1")
          }
          eigen_denoised <- rep(NA, length(components))
          threshold <- ceiling(alpha * length(components))
          m <- paste0("Use method ", method, "\n",
                      "threshold = ", threshold, "\n")
          if (verbose) cat(m)
          if (threshold > 1) {
            eigen_denoised[1:threshold] <- irl$eigen[1:threshold]
            } else (stop("No eigenvalue preserved"))
        },
        identity = {
          if (!is.numeric(location)) {
            stop("Invalid location input")
          }
          if (missing(location)) {
            stop("Invalid location input")
          }
          if (max(location) > rnk) {
            stop("Location must be smaller than rnk")
          } else if (min(location) <= 0) {
            stop("Location out of bounds")
          }
          eigen_denoised <- rep(NA, length(components))
          eigen_denoised[location] <- irl$eigen[location]
          m <- paste0("Use method ", method, "\n",
                      "location = ", toString(location), "\n")
          if (verbose) cat(m)
        },
        stop("Invalid method input")
        )
# ----------------------------------------
# Zero out or average truncated eigenvalues
# ----------------------------------------
    numerator <- sum(irl$eigen) - sum(na.omit(eigen_denoised))
    denominator <- sum(is.na(eigen_denoised))
    gamma <- ifelse(zeroout,
                    0,
                    round(numerator / denominator, 4))
    mg <- paste0("The averaged truncated eigenvalues = ", gamma, "\n")
    if (verbose) cat(mg)
    eigen_denoised <- ifelse(is.na(eigen_denoised), gamma, eigen_denoised)
# ---------------------------------
# Calculate estimated X, COV, CORR
# ---------------------------------
    #calculate estimated X
    x_denoised <- u %*% diag(eigen_denoised * pdim) %*% t(v)
    #calculate empirical covariance
    v_denoised <- v %*% diag(eigen_denoised * pdim) %*% t(v) / (ndf - 1L)
    #symmetric rescaling to correlation matrix
    e_denoised <- v_denoised / tcrossprod(diag(v_denoised)^0.5)

    ret <- list(eigen_denoised = eigen_denoised,
                x_denoised  = x_denoised,
                v_denoised  = v_denoised,
                e_denoised  = e_denoised,
                v          = v,
                u          = u,
                message    = list(m, mg))
    attr(ret, "class") <- "subspace_denoised"
    ret
}

print.subspace_denoised <- function(x, ...) {
  cat(paste0("A denoised estimator of X, ",
              "correlation matrix and truncated eigenvalues.\n"))
  cat(x$message[[1]])
  cat(x$message[[2]])
}
