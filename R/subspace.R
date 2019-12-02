#' @title A constructor function for the "subspace" class
#'
#' @description This function calculates scaled eigenvalues and eigenvectors of X matrix, 
#'   as well as sampled eigenvalues from a random noise matrix N of the same dimension as X, which follows a 
#'   Marcenko-Pastur distribution with package "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param X A numeric real-valued matrix with n number of samples and p number of features. 
#'   If p>n, a warning message is generated and the transpose of X is used.
#' @param rank A series of right singular vectors to estimate. rank must be smaller or equal to min(nrow(X),ncol(X)).
#' @param basis Choose eigenvalue decomposition or singular value decomposition.
#' @param MP A logical value. If true, sample eigenvlaues from random noise matrix with MP distribution.
#' @param times Split data into X times for parallel computation.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of X.}
#'   \item{pdim:}{ The number of dimensions of X.}
#'   \item{rank:}{ A series of right singular vectors estimated.}
#'   \item{var_correct:}{ Corrected population variance for Marcenko-Pastur distribution.}
#'   \item{transpose_flag:}{ A logical value indicating whether the matrix X is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for specified rank and corresponding dimensions.}
#'   \item{sigma_a:}{ A vector of scaled eigenvalues up to max(rank).}
#'   \item{MP_irl:}{ A data frame of samped expected eigenvalues from Marcenko-Pastur for specified rank and corresponding dimensions.}
#'   \item{sigma_MP:}{ A vector of samped expected eigenvalues from Marcenko-Pastur up to max(rank).}
#'   \item{v:}{ Right singular vectors of X matrix for specified rank.}
#'   \item{u:}{ Left singular vectors of X matrix or specified rank.}
#'   \item{basis:}{ Basis for decomposition.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n = 150, p = 100, ncc = 10, var = 2)
#' Subspace <- subspace(X)
#' Subspace <- subspace(X, rank = 8:30, basis = "eigen")
#' Subspace <- subspace(X, rank = c(2, 3, 6, 16), times = 10, basis = "eigen")
#' Subspace <- subspace(X, rank = 1:20, MP= FALSE)
#' Subspace
#' plot(Subspace, Changepoint = 0, annotation = 15)
#' }
#' @seealso 
#' * [MarchenkoPasturPar()] calculates upper and lower limits of Marcenko-Pastur distribution from RMTstat package.
#' 
#' * [rmp()] sample scaled eigenvalues of random noise matrix from RMTstat package.
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export

subspace <- function(X,                              # data matrix
                     rank = NA,                      # a series of singular vectors to estimate
                     basis = c("eigen","singular"),  # chopse basis for decomposition
                     MP = TRUE,                      # Perform MP simulation for noise
                     times = NA,                     # split data into X times for parallel computation.
                     verbose = TRUE,                 # output message
                     ...)
{
# ---------------------------------------------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------------------------------------------
  if (is.null(X)) {
    stop("Invalid input X")
  } else {
    # Checking for rank input
    if (missing(rank)) {
      rank <- 1:min(nrow(X), ncol(X))
      if (verbose) {
        cat("No rank specified. Calculating full singular value decomposition instead.\n")
      }
    } else {
      rank <- CheckRankInput(rank, nrow(X), ncol(X), verbose)
    }
    # Checking for times input
    if (MP) {
      if (missing(times)) {
        times <- 0
      } else {
        CheckTimesInput(times, nrow(X), ncol(X), verbose)
      }
    }
    # Check X matrix dimension
    params <- CheckDimMatrix(X, rnk = max(rank))
  }
# ---------------------------------------------------------------------------------------------------------
# Basic parameter set up
# ---------------------------------------------------------------------------------------------------------
  ndf             <- params$ndf
  pdim            <- params$pdim
  svr             <- params$svr
  rnk             <- params$rnk
  transpose_flag  <- params$transpose_flag
# ----------------------------------------------------------------------------------------------------------
# Singular Value Decomposition of X matrix
# ----------------------------------------------------------------------------------------------------------  
  #Col Mean center Matrix
  Xstd <- sweep(X, 2L, colMeans(X))

  if (rnk > pdim/2) {
      tmp <- svd(Xstd)
  } else {
      tmp <- irlba(Xstd, nv = rnk)
  }

  irl <- tibble(eigen = tmp$d[rank]^2/(pdim), dim = rank)
  sigma_a <- tmp$d[1:rnk]^2/(pdim)
  v   <- tmp$v[, rank]
  u   <- tmp$u[, rank]
# ----------------------------------------------------------------------------------------------------------
# Simulate noise eigenvalues from rmp function
# ----------------------------------------------------------------------------------------------------------
  if (MP) {
    var_correct = min(irl$eigen)/MarchenkoPasturPar(ndf, pdim, var = 1, svr = svr)$upper
    if (verbose) {
      cat("The corrected variance of MP distribution is ",var_correct,".\n",sep = "")
    }
    if (times == 0) {
      system.time({sim = rmp(pdim, ndf = ndf, pdim = pdim, var = var_correct, svr = ndf/pdim)})
    } else {
      ncores <- detectCores()
      registerDoParallel(ncores) 
      system.time({
        sim <- foreach(1:times, .combine = c) %dopar% {
          tmp <- rmp(pdim/times, ndf = ndf, pdim = pdim, var = var_correct, svr = ndf/pdim)
        }})
    }

    MP_irl <- tibble(eigen = sim[order(sim, decreasing = T)][rank], dim = rank)
    sigma_MP <- sim[order(sim,decreasing = T)][1:rnk]
  }
  value <- (list(ndf  = ndf, 
                 pdim = pdim, 
                 rank = rank,
                 var_correct = ifelse(MP, var_correct, NA),
                 transpose_flag = transpose_flag, 
                 irl = irl, 
                 sigma_a = sigma_a,
                 MP_irl = switch(2 - MP, MP_irl, NULL),
                 sigma_MP = switch(2 - MP, sigma_MP, NULL),
                 v = v, 
                 u = u,
                 basis = basis))
  attr(value, "class") <- "subspace"
  value
}

#' @export
print.subspace <- function(obj) {
  cat("An object of class subspace within",ifelse(obj$transpose_flag,"transposed",""), "X matrix with", obj$ndf, "samples and", obj$pdim, "features.\n")
  if (all(diff(obj$rank)==1)) {
    cat("Estimated rank range from ",min(obj$rank), " to ", max(obj$rank),"\n")
  } else {
      cat("Estimated rank range", obj$rank, "\n")
  }
  
}

# #x_clipped <- x %–% subspace(x, 1:9, MP = FALSE, basis = “eigen”)
# #' @export
# `%-%` <- function(X, subspace_) {
# # ---------------------------------------------------------------------------------------------------------
# # Basic parameter set up
# # ---------------------------------------------------------------------------------------------------------  
#   AmbientSpace <- subspace(X, rnk = min(nrow(X),ncol(X)), MP = FALSE, basis = "eigen")

#   rank            <- subspace_$rank
#   var_correct     <- subspace_$var_correct
#   sigma_a         <- subspace_$sigma_a
#   sigma_MP        <- subspace_$sigma_MP
  
#   xi_clipped = AmbientSpace$sigma_a
#   xi_clipped[rank] = irl$eigen[rank]
#   xi_clipped = ifelse(is.na(xi_clipped),0,xi_clipped)
#   #calculate estimated X
#   X_clipped = u%*%diag(xi_clipped*pdim)%*%t(v)
#   #calculate empirical covariance
#   V_clipped = v%*%diag(xi_clipped*pdim)%*%t(v)/ (ndf - 1L)
#   ## symmetric rescaling to correlation matrix
#   E_clipped<-V_clipped / tcrossprod(diag(V_clipped) ^ 0.5)
#   return(list(xi_clipped=xi_clipped,X_clipped=X_clipped,E_clipped=E_clipped,v=v,u=u))
# }

# #' clipped by location
# #' @export
# setMethod('%-%',
#   signature(X = "matrix", subspace_="subspace"),
#   function(X, subspace_) return(SetCols.bm(x, j, value)))


#' @export
CheckRankInput <- function(rank, ndf, pdim, verbose=TRUE) {
  stopifnot(is.numeric(rank))
  stopifnot(rank%%1==0)
  if (length(rank) == 1) {
    if (rank < 1) {
      stop("Rank must be larger than 1.\n")
    }
    if (rank > min(ndf, pdim)) {
      stop("Rank out of bounds.\n")
    }
    if (verbose) {
      cat("Calculating rank from 1 to", rank, ".\n")
    }
    return (1:rank)
  } else if (length(rank) > 1) {
    if (min(rank) < 1) {
      stop("Rank must be larger than 1.\n")
    }
    if (max(rank) > min(ndf, pdim)) {
      stop("Rank out of bounds.\n")
    }
    if (verbose) {
      if (all(diff(rank)==1)) {
        cat("Calculating rank from", min(rank), "to", max(rank), ".\n")
      } else {
        cat("Calculating rank range", rank, "\n")
      }
    }
    return (rank[order(rank)])
  } else {
    stop ("Invalid rank input.\n")
  }
}

#' @export
CheckTimesInput <- function(times, ndf, pdim, verbose=TRUE) {
  stopifnot(is.numeric(times))
  stopifnot(times%%1==0)
  stopifnot(length(times) == 1)
  if (times < 0) {
    stop("Times must be positive")
  }
  else if (times > min(ndf - 1, pdim - 1)) {
    stop("Times must be strictly less than min(nrow(X), ncol(X))")
  }
}
