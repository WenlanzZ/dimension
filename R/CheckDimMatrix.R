#' @title Return parameter setting for MarcenkoPasturSample function.
#'
#' @param X A numeric real- or real-valued sparse matrix with n number of samples and p number of features. If p>n, a warning message is generated and the transpose of X is used.
#' @param rnk number of right singular vectors to estimate. rnk must be smaller or equal to max(nrow(X),ncol(X)).
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ number of degrees of freedom of X and the white Wishart matrix.}
#'   \item{pdim:}{ number of dimensions of X and the white Wishart matrix.}
#'   \item{svr:}{ samples to variance ratio; equals to ndf divided by pdim.}
#'   \item{rnk:}{ number of right singular vectors to estimate.}
#'   \item{transpose_flag:}{ whether the matrix X is transposed.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n=1000,p=500,ncc=10,var=2,fact = 1)
#' params <- CheckDimMatrix(X,rnk=40)
#' }
#' #should import RMTstat after fix bug
#' @seealso Random Matrix Theory pacakge credit to Gregory Giecold and Lionel Ouaknin
#' @importFrom  tibble tibble
#' @export

CheckDimMatrix <- function(X = X,                    # data matrix
                           rnk = NA,                 # number of singular vectors to estimate
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
