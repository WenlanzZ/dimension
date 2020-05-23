#' @title correct eigenvalues from a subspace object and return a subspace object
#'
#' @description This function correct eigenvalues from a subspace by substracting
#' sampling eigenvalues from a random noise matrix N of the same dimension as x,
#'  which follows a Marcenko-Pastur distribution with package
#'  "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param x A subspace object.
#' @param num_est_samples the number of resamples to take from the 
#' Marcenko-Pastur distribution to estimate the eigenvalues.
#' @param verbose output message
#' @param ... Extra parameters
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ The number of degrees of freedom of x.}
#'   \item{pdim:}{ The number of dimensions of x.}
#'   \item{components:}{ A series of right singular vectors estimated.}
#'   \item{var_correct:}{ Corrected population variance
#'    for Marcenko-Pastur distribution.}
#'   \item{transpose_flag:}{ A logical value indicating
#'    whether the matrix x is transposed.}
#'   \item{irl:}{ A data frame of scaled eigenvalues for
#'    specified components and corresponding dimensions.}
#'   \item{sigma_a:}{ A vector of scaled eigenvalues up to max(components).}
#'   \item{mp_irl:}{ A data frame of sampled expected eigenvalues from
#'    Marcenko-Pastur for specified components and corresponding dimensions.}
#'   \item{sigma_mp:}{ A vector of samped expected eigenvalues from
#'    Marcenko-Pastur up to max(components).}
#'   \item{v:}{ Right singular vectors of x matrix for specified components.}
#'   \item{u:}{ Left singular vectors of x matrix or specified components.}
#' }
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' x %>% create_subspace(components = 8:30) %>% correct_eigenvalues() %>% plot()
#' @seealso
#' * [MarchenkoPasturPar()] calculates upper and lower limits
#'  of Marcenko-Pastur distribution from RMTstat package.
#'
#' * [rmp()] sample scaled eigenvalues of random noise matrix
#'  from RMTstat package.
#' @export
correct_eigenvalues <- function(subspace, num_est_samples, verbose, ...) {
  UseMethod("correct_eigenvalues", subspace)
}

#' @export
correct_eigenvalues.default <- function(subspace, num_est_samples, verbose,
                                        ...) {
  stop("Don't know how to correct eigenvalues for an object of type ",
       paste(class(subspace), collapse = " "),
       ". Did you mean to call subspace()?")
}

#' @importFrom tibble tibble
#' @import foreach
#' @export
correct_eigenvalues.subspace <- 
  function(subspace, num_est_samples = min(c(subspace$ndf, subspace$pdim)-1), 
           verbose = FALSE, ...) {
  #check if it is a subspace object
  check_times_input(num_est_samples, subspace$ndf, subspace$pdim, 
                    verbose = verbose)

  # ----------------------
  # Basic parameter set up
  # ----------------------

  ndf             <- subspace$ndf
  pdim            <- subspace$pdim
  components      <- subspace$components
  transpose_flag  <- subspace$transpose_flag
  
  denominator <- marcenko_pastur_par(ndf, pdim, var = 1, svr = ndf / pdim)$upper
  var_correct <- min(subspace$irl$eigen) / denominator
  if (verbose) {
  	cat("The corrected variance of mp distribution is ",
      var_correct,
      ".\n",
      sep = "")
  }

  sim <- foreach(seq_len(num_est_samples), .combine = c) %dorng% {
    rmp(pdim / num_est_samples, ndf = ndf, pdim = pdim, var = var_correct,
        svr = ndf / pdim)
  }

  mp_irl <- tibble(eigen = sim[order(sim, decreasing = TRUE)][components],
                 dim = components)
  sigma_mp <- sim[order(sim, decreasing = T)][seq_len(max(components))]

  ret <- list(ndf  = ndf,
	            pdim = pdim,
	            components = components,
	            var_correct = var_correct,
	            transpose_flag = transpose_flag,
	            irl = subspace$irl,
	            sigma_a = subspace$sigma_a - sigma_mp,
	            mp_irl = mp_irl,
	            sigma_mp = sigma_mp,
	            v = subspace$v,
	            u = subspace$u)
  class(ret) <- c("eigen_corrected_subspace", "subspace")
  ret
}
