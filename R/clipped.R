#' @title Eigenvalue clipping Procedure
#'
#' @description This function clips the scaled eigenvalues of X in order to provide a cleaned estimator E_clipped 
#'   of the underlying correlation matrix. Proceeds by keeping the [N * alpha] top eigenvalues and shrinking 
#'   the remaining ones by a trace-preserving constant (i.e. Tr(E_clipped) = Tr(E)) or zero out remaining ones. 
#'   This function is adpated from "Python for Random Matrix Theory" credit to J.-P. Bouchaud and M. Potters.
#' @param X A numeric real-valued matrix with n number of samples and p number of features. 
#'   If p>n, a warning message is generated and the transpose of X is used.
#' @param subspace_ A subspace class.
#' @inheritParams rank A series of right singular vectors to estimate. rank must be smaller or equal to min(nrow(X),ncol(X)).
#' @param method The method to be used; method = "threshold" returns (1-alpha)*rnk proportion of eigenvalues above threshold; 
#'   method = "hard" returns all the empirical eigenvalues greater than the upper limit of the support to the Marcenko-Pastur spectrum;
#'   method = "identity" returns eigenvalued specified in location vector.
#' @param alpha Determine the fraction to keep of the top eigenvalues in threshold method. 
#' @param location Indicate the location of eigenvalues to keep in identity method.
#' @param zeroout A logical value to zero out eigenvalues when clipping. default is to set to FALSE.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{xi_clipped:}{ A cleaned estimator of eigenvalues through a simple eigenvalue clipping procedure (cf. reference below).}
#'   \item{X_clipped:}{ A cleaned estimator of the true sigmal matrix underlying a noisy high dimensional matrix X.}
#'   \item{E_clipped:}{ A cleaned estimator of the true correlation matrix underlying a noisy, in-sample estimate E (empirical correaltion matrix estimated from X (cf. reference below)).}
#'   \item{v:}{ Right singular vectors of X matrix for specified rank.}
#'   \item{u:}{ Left singular vectors of X matrix or specified rank.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
#' X_clp <- clipped(X, rank = 20, method = "threshold", alpha = 0.9, zeroout = TRUE)
#' x_clp<-clipped(x, rnk = 20, method = "hard", zeroout = FALSE)
#' x_clp<-clipped(x, rnk = 20, method = "identity", location = c(1:15), zeroout = FALSE)
#' 
#' #equivelantly, if Subspace is calcualted
#' Subspace <- subspace(X, rank = 1:40, times = 10, basis = "eigen")
#' X_clp <- clipped(subspace_ = Subspace, method = "identity", location = c(1:5), zeroout = TRUE)
#' }
#' Reference
#' ---------
#' "Financial Applications of Random Matrix Theory: a short review",
#' J.-P. Bouchaud and M. Potters
#' arXiv: 0910.1205 [q-fin.ST]
#' @seealso 
#' * [MarchenkoPasturPar()] calculates upper and lower limits of Marcenko-Pastur distribution from RMTstat package.
#' @export

clipped <- function(X,                                          # data matrix
                    subspace_ = NULL,                               # a list of ouput from function CheckDimMatrix
                    rank = NA,                                      # number of singular vectors to estimate
                    method = c("threshold", "hard", "identity"),    # choose method from c("threshold","hard")
                    alpha = NA,                                     # determining the fraction to keep of the top eigenvalues of an empirical correlation matrix. 
                    location = NA,                                  # preserve eigenvalues at these locations
                    zeroout = FALSE,                                # zero out eigenvalues when clipping. default is to set it to average.
                    verbose = TRUE,                                 # output message
                    ...)                                            # optional arguments
{
# ----------------------------------------------------------------------------------------------------------------------------
# Check input parameters
# ----------------------------------------------------------------------------------------------------------------------------

  if (!missing(subspace_)) {
    cat("Subspace have already been calculated.\n")
  } else {
    if (is.null(X)) {
      stop("Invalid input X")
    } else {
      if (missing(rank)) {
        rank <- 1:min(nrow(X), ncol(X))
        if (verbose) {
          cat("No rank specified. Calculating full singular value decomposition instead.\n")
        }
      }
      subspace_ <- subspace(X, rank = rank, MP = FALSE, basis = "eigen")
      if (verbose) {
      cat("Finish checking dimension of X and calculating eigenvalues of X.\n")
      }
    }
  }
# ---------------------------------------------------------------------------------------------------------
# Basic parameter set up
# ---------------------------------------------------------------------------------------------------------
    ndf  <- subspace_$ndf
    pdim <- subspace_$pdim
    svr  <- ndf/pdim
    rank <- subspace_$rank
    rnk  <- max(rank)
    irl  <- subspace_$irl
    v    <- subspace_$v
    u    <- subspace_$u
    transpose_flag <- subspace_$transpose_flag
# ----------------------------------------------------------------------------------------------------------
# Determine eigenvalues to preserve
# ----------------------------------------------------------------------------------------------------------  
    switch(method,
        hard = {
          lambda_min <- MarchenkoPasturPar(ndf, pdim, var = 1, svr = svr)$lower
          if (verbose) cat("lambda_min = ", lambda_min, "\n")
          lambda_max <- MarchenkoPasturPar(ndf, pdim, var = 1, svr = svr)$upper
          if (verbose) cat("lambda_max = ", lambda_max, "\n")
          xi_clipped <- ifelse(irl$eigen >= lambda_max, irl$eigen, NA)
        },
        threshold = {
          if (missing(alpha)) {
            stop("Alpha must be specified")
          } else if (alpha < 0) {
            stop("Alpha must be positive")
          } else if (alpha > 1) {
            stop("Alpha must be less or equal to 1")
          }
          xi_clipped <- rep(NA, length(rank))
          threshold <- ceiling(alpha * length(rank))
          if (verbose) cat("threshold = ", threshold, "\n")
          if (threshold > 0) xi_clipped[1:threshold] <- irl$eigen[1:threshold]
          else(stop("No eigenvalue preserved"))
        },
        identity = {
          if (!is.numeric(location)) {
            stop("Invalid location input")
          }
          if (missing(location)) {
            stop("Location must be specified")
          }
          if (max(location) > rnk) {
            stop ("Location must be smaller than rnk")
          } else if (min(location) <= 0) {
            stop ("Location out of bounds")
          }
          xi_clipped <- rep(NA, length(rank))
          xi_clipped[location] <- irl$eigen[location]
          },
        stop("Invalid method input")
    )
# ----------------------------------------------------------------------------------------------------------
# Zero out or average clipped eigenvalues
# ----------------------------------------------------------------------------------------------------------  
    gamma <- ifelse(zeroout, 0, (sum(irl$eigen) - sum(na.omit(xi_clipped))) / sum(is.na(xi_clipped))) 
    if (verbose) cat("gamma = ", gamma, "\n")   
    xi_clipped <- ifelse(is.na(xi_clipped), gamma, xi_clipped)
# ----------------------------------------------------------------------------------------------------------
# Calculate estimated X, COV, CORR
# ----------------------------------------------------------------------------------------------------------  
    #calculate estimated X
    X_clipped <- u %*% diag(xi_clipped * pdim) %*% t(v)
    #calculate empirical covariance
    V_clipped <- v %*% diag(xi_clipped * pdim) %*% t(v) / (ndf - 1L)
    ## symmetric rescaling to correlation matrix
    E_clipped <- V_clipped / tcrossprod(diag(V_clipped)^0.5)

    return(list(xi_clipped = xi_clipped, 
                X_clipped  = X_clipped, 
                E_clipped  = E_clipped, 
                v          = v, 
                u          = u))
}


