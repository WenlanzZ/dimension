#' @title A constructor function for the subspace class
#'
#' @description This function calculates scaled eigenvalues
#'  and eigenvectors of x matrix.
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p > n, a warning message is generated and
#'  the transpose of x is used.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x), ncol(x)).
#' @param decomposition The method to be used; method = "svd"
#'  returns results from singular value decomposition; method = "eigen"
#'  returns results from eigenvalue decomposition.
#' @param verbose output message
#' @param ... Extra parameters
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of x.}
#'   \item{pdim:}{ The number of dimensions of x.}
#'   \item{components:}{ A series of right singular vectors estimated.}
#'   \item{transpose_flag:}{ A logical value indicating
#'    whether the matrix x is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for
#'    specified components and corresponding dimensions.}
#'   \item{sigma_a:}{ A vector of scaled eigenvalues up to max(components).}
#'   \item{v:}{ Right singular vectors of x matrix for specified components.}
#'   \item{u:}{ Left singular vectors of x matrix or specified components.}
#' }
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' x %>% create_subspace(components = 8:30) %>% str()
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import Matrix
#' @importFrom dplyr %>%
#' @export
create_subspace <- function(x, components = NULL, decomposition = c("svd", "eigen"), verbose = TRUE, ...) {
  # Checking for components input
  if (is.null(components)) {
    components <- seq_len(min(nrow(x), ncol(x)))
    if (verbose) {
      message(paste0("No component specified. ",
          "Calculating full singular value decomposition instead.\n"))
    }
  } else {
    components <- check_comp_input(components, nrow(x), ncol(x), verbose)
  }
  
  if (missing(decomposition)) {
    decomposition <- "svd"
  }
  switch(decomposition,
        svd = {
          # Check x matrix dimension
          params <- check_dim_matrix(x, rnk = max(components))
          ndf             <- params$ndf
          pdim            <- params$pdim
          rnk             <- params$rnk
          transpose_flag  <- params$transpose_flag
          # Col Mean center Matrix
          if (transpose_flag) {
            x <- t(x)
          }

          if (rnk >= pdim / 2) {
            tryCatch({
                x_std <- sweep(x, 2L, colMeans(x))
              },
              warning = function(w) {
                message(paste0("Cannot allocate matrix in memory, ",
                               "try transforming matrix to dgCMatrix",
                               "or a smaller proportion of eigenvalues.\n"))
              },
              error = function(e) {
                message("Caught an error in scaling matrix!\n")
              }
            )
            tmp <- svd(x_std)
          } else {
            tmp <- irlba(x, center = TRUE, nv = rnk)
          }

          irl     <- tibble(eigen = tmp$d[components]^2 / pdim, dim = components)
          sigma_a <- tmp$d[1:rnk]^2 / pdim
          v       <- tmp$v[, components]
          u       <- tmp$u[, components]

        },
        eigen = {
          ndf             <- nrow(x)
          pdim            <- ncol(x)
          rnk             <- max(components)
          transpose_flag  <- FALSE
          if (rnk >= pdim / 2) {
            tmp <- eigen(x)
          } else {
            tmp <- partial_eigen(x, n = rnk, symmetric = TRUE)
          }

          irl     <- tibble(eigen = tmp$values[components], dim = components)
          sigma_a <- tmp$values
          v       <- NULL
          u       <- tmp$vectors[, components]
          },
          stop("Invalid method input")
          )
  
  value <- (list(ndf  = ndf,
                 pdim = pdim,
                 components = components,
                 transpose_flag = transpose_flag,
                 irl = irl,
                 sigma_a = sigma_a,
                 v = v,
                 u = u))
  class(value) <- "subspace"
  value
}
