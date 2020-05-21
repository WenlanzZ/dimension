#' @title Return parameter setting for subspace function
#'
#' @param x A numeric real-valued matrix with n number of samples and
#' p number of features.
#' If p>n, a warning message is generated and the transpose of x is used.
#' @param rnk The number of right singular vectors to estimate.
#' rnk must be smaller or equal to max(nrow(x),ncol(x)).
#' @param verbose output message
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of x.}
#'   \item{pdim:}{ The number of dimensions of x.}
#'   \item{svr:}{ The samples to variance ratio: equals to ndf divided by pdim.}
#'   \item{rnk:}{ The number of right singular vectors estimated.}
#'   \item{transpose_flag:}{ A logical value indicating whether the
#'   matrix x is transposed.}
#' }
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' params <- check_dim_matrix(x, rnk = 40)
#' @seealso [checkDesignMatrix()] from Random Matrix Theory pacakge
#' credit to Gregory Giecold and Lionel Ouaknin.
#' @importFrom  tibble tibble
#' @export
check_dim_matrix <- function(x,
                             rnk = NA,
                             verbose = TRUE) {
# --------------------------------
# Check input parameters
# --------------------------------
    if (is.null(x)) {
        stop("Invalid input x.\n")
    }

    ndf  <- nrow(x)
    pdim <- ncol(x)

    if (min(ndf, pdim) <= 5) {
        stop("x matrix should have at least 5 rows or columns.\n")
    }

    if (missing(rnk)) {
        rnk <- min(ndf, pdim)
        if (verbose) {
            message(paste0("No component specified. ",
                "Calculating full singular value decomposition instead.\n"))
        }
    }

    if (rnk <= 0) {
        stop("Rnk must be positive.\n")
    } else if (rnk > min(ndf, pdim)) {
        stop("Rnk out of bounds.\n")
    }
# --------------------------------
# Transpose matrix when p > n
# --------------------------------
    transpose_flag <- FALSE

    if (nrow(x) < ncol(x)) {
        warning("nrow(x) < ncol(x). A transpose of x is used instead.\n")
        transpose_flag <- TRUE
        ndf <- ncol(x)
        pdim <- nrow(x)
    }

    svr <- ndf / pdim

    return(list(ndf  = ndf,
                pdim = pdim,
                svr  = svr,
                rnk  = rnk,
                transpose_flag = transpose_flag))
}
