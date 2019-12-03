#' @title Return parameter setting for subspace function
#'
#' @param X A numeric real-valued matrix with n number of samples and p number of features. If p>n, a warning message is generated and the transpose of X is used.
#' @param rnk The number of right singular vectors to estimate. rnk must be smaller or equal to max(nrow(X),ncol(X)).
#' @param verbose output message
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of X.}
#'   \item{pdim:}{ The number of dimensions of X.}
#'   \item{svr:}{ The samples to variance ratio: equals to ndf divided by pdim.}
#'   \item{rnk:}{ The number of right singular vectors estimated.}
#'   \item{transpose_flag:}{ A logical value indicating whether the matrix X is transposed.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
#' params <- CheckDimMatrix(X, rnk = 40)
#' }
#' 
#' @seealso [checkDesignMatrix()] from Random Matrix Theory pacakge credit to Gregory Giecold and Lionel Ouaknin.
#' @importFrom  tibble tibble
#' @export

CheckDimMatrix <- function(X,                        # A data matrix
                           rnk = NA,                 # The number of singular vectors to estimate
                           verbose = TRUE)           # output message
{
# ----------------------------------------------------------------------------------------------------------
# Check input parameters
# ----------------------------------------------------------------------------------------------------------
    if (is.null(X)) {
        stop("Invalid input X")
    }

    ndf  <- nrow(X)
    pdim <- ncol(X)
    
    if (missing(rnk)) {
        rnk <- min(ndf, pdim)
        if (verbose) {
            cat("No rnk specified. Calculating full singular value decomposition instead.\n")
        }
    }

    if (rnk <= 0) {
        stop("Rnk must be positive")
    } else if (rnk > min(ndf, pdim)) {
        stop("Rnk out of bounds")
    }
# ----------------------------------------------------------------------------------------------------------
# Transpose matrix when p > n
# ----------------------------------------------------------------------------------------------------------
    transpose_flag <- FALSE
    if (nrow(X) < ncol(X)) {
        warning("The number of samples of X is smaller than the number of features of X. 
                A transpose of X is used instead.\n")
        X <- t(X)
        transpose_flag <-TRUE
        ndf <- nrow(X)
        pdim <- ncol(X)
    }

    svr <- ndf/pdim
    
    return(list(ndf  = ndf,
                pdim = pdim, 
                svr  = svr, 
                rnk  = rnk, 
                transpose_flag = transpose_flag))
}
