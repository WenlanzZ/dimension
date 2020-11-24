#' @title A constructor function for the subspace class
#'
#' @description This function calculates scaled eigenvalues
#'  and eigenvectors of x matrix, as well as sampled eigenvalues
#'  from a random noise matrix N of the same dimension as x,
#'  which follows a Marcenko-Pastur distribution with package
#'  "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p > n, a warning message is generated and
#'  the transpose of x is used.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x), ncol(x)).
#' @param decomposition The method to be used; method = "svd"
#'  returns results from singular value decomposition; method = "eigen"
#'  returns results from eigenvalue decomposition.
#' @param mp A logical value. If true, sample eigenvlaues from random noise
#'  matrix with mp distribution.
#' @param num_est_samples Split data into num_est_samples-fold for
#'  parallel computation.
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
#'   \item{sigma_a:}{ A vector of corrected eigenvalues up to max(components)
#'    if mp is true.}
#'   \item{mp_irl:}{ A data frame of sampled expected eigenvalues from
#'    Marcenko-Pastur for specified components and corresponding dimensions.}
#'   \item{sigma_mp:}{ A vector of samped expected eigenvalues from
#'    Marcenko-Pastur up to max(components).}
#'   \item{v:}{ Right singular vectors of x matrix for specified components.}
#'   \item{u:}{ Left singular vectors of x matrix or specified components.}
#' }
#' @examples
#' x <- x_sim(n = 100, p = 150, ncc = 10, var = c(rep(10, 5), rep(1, 5)))
#' Subspace <- subspace(x)
#' Subspace <- subspace(x, components = 8:30)
#' Subspace <- subspace(x, components = c(2, 3, 6, 16))
#' Subspace <- subspace(x, components = 1:20, mp = FALSE)
#' Subspace
#' plot(Subspace, changepoint = 0, annotation = 1:10)
#' @seealso
#' * [MarchenkoPasturPar()] calculates upper and lower limits
#'  of Marcenko-Pastur distribution from RMTstat package.
#'
#' * [rmp()] sample scaled eigenvalues of random noise matrix
#'  from RMTstat package.
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import doRNG
#' @import Matrix
#' @export
subspace <- function(x, components, decomposition, mp, num_est_samples, verbose, ...) {
  UseMethod("subspace", x)
}

#' @title A constructor function for the subspace class
#'
#' @description This function calculates scaled eigenvalues
#'  and eigenvectors of x matrix, as well as sampled eigenvalues
#'  from a random noise matrix N of the same dimension as x,
#'  which follows a Marcenko-Pastur distribution with package
#'  "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p > n, a warning message is generated and
#'  the transpose of x is used.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x), ncol(x)).
#' @param decomposition The method to be used; method = "svd"
#'  returns results from singular value decomposition; method = "eigen"
#'  returns results from eigenvalue decomposition.
#' @param mp A logical value. If true, sample eigenvlaues from random noise
#'  matrix with mp distribution.
#' @param num_est_samples Split data into num_est_samples-fold for
#'  parallel computation.
#' @param verbose output message
#' @param ... Extra parameters
#' @seealso
#' * [MarchenkoPasturPar()] calculates upper and lower limits
#'  of Marcenko-Pastur distribution from RMTstat package.
#'
#' * [rmp()] sample scaled eigenvalues of random noise matrix
#'  from RMTstat package.
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import doRNG
#' @import Matrix
#' @export
subspace.default <- function(x, components, decomposition, mp, num_est_samples, verbose, ...) {
  stop("Don't know how to create a subspace object from a class of type: ",
       class(x))
}

#' @title A constructor function for the subspace class
#'
#' @description This function calculates scaled eigenvalues
#'  and eigenvectors of x matrix, as well as sampled eigenvalues
#'  from a random noise matrix N of the same dimension as x,
#'  which follows a Marcenko-Pastur distribution with package
#'  "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p > n, a warning message is generated and
#'  the transpose of x is used.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x), ncol(x)).
#' @param decomposition The method to be used; method = "svd"
#'  returns results from singular value decomposition; method = "eigen"
#'  returns results from eigenvalue decomposition.
#' @param mp A logical value. If true, sample eigenvlaues from random noise
#'  matrix with mp distribution.
#' @param num_est_samples Split data into num_est_samples-fold for
#'  parallel computation.
#' @param verbose output message
#' @param ... Extra parameters
#' @seealso
#' * [MarchenkoPasturPar()] calculates upper and lower limits
#'  of Marcenko-Pastur distribution from RMTstat package.
#'
#' * [rmp()] sample scaled eigenvalues of random noise matrix
#'  from RMTstat package.
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import doRNG
#' @import Matrix
#' @export
subspace.matrix <- function(x, components = NULL, decomposition = c("svd", "eigen"), mp = TRUE,
                            num_est_samples = NA, verbose = FALSE, ...) {
  subspace(x = x, components = components, decomposition = c("svd", "eigen"), mp = mp,
           num_est_samples = num_est_samples, verbose = verbose, ... = ...)
}

#' @title A constructor function for the subspace class
#'
#' @description This function calculates scaled eigenvalues
#'  and eigenvectors of x matrix, as well as sampled eigenvalues
#'  from a random noise matrix N of the same dimension as x,
#'  which follows a Marcenko-Pastur distribution with package
#'  "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features. If p > n, a warning message is generated and
#'  the transpose of x is used.
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x), ncol(x)).
#' @param decomposition The method to be used; method = "svd"
#'  returns results from singular value decomposition; method = "eigen"
#'  returns results from eigenvalue decomposition.
#' @param mp A logical value. If true, sample eigenvlaues from random noise
#'  matrix with mp distribution.
#' @param num_est_samples Split data into num_est_samples-fold for
#'  parallel computation.
#' @param verbose output message
#' @param ... Extra parameters
#' @seealso
#' * [MarchenkoPasturPar()] calculates upper and lower limits
#'  of Marcenko-Pastur distribution from RMTstat package.
#'
#' * [rmp()] sample scaled eigenvalues of random noise matrix
#'  from RMTstat package.
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import doRNG
#' @import Matrix
#' @export
subspace.Matrix <- function(x, components = NULL, decomposition = c("svd", "eigen"), mp = TRUE,
                            num_est_samples = NA, verbose = FALSE, ...) {
  subspace(x, components = components, decomposition = c("svd", "eigen"), mp = mp,
           num_est_samples = num_est_samples, verbose = verbose, ... = ...)
}

subspace <- function(x, components = NULL, decomposition = c("svd", "eigen"), mp = TRUE,
                     num_est_samples = NA, verbose = FALSE, ...) {

  # ----------------------
  # Check input parameters
  # ----------------------

  # Checking for components input
  if (is.null(components)) {
    components <- seq_len(min(nrow(x), ncol(x)))
    if (verbose) {
      message(paste0("No component specified. ",
          "Calculating full singular value decomposition instead.\n"))
    }
  } else {
    components <- check_comp_input(components, nrow(x), ncol(x), verbose = TRUE)
  }
  if (missing(decomposition)) {
    decomposition <- "svd"
  }
  s <- create_subspace(x, components = components, decomposition = decomposition, verbose)

  if (mp & decomposition == "svd") {
    if (missing(num_est_samples)) {
      num_est_samples <- 0
    } else {
      check_num_est_samples_input(num_est_samples, s$ndf, s$pdim,
                      verbose = TRUE)
    }
    s <- correct_eigenvalues(s, num_est_samples = num_est_samples, verbose)
  }
  s
}


#' @title Print subsapce
#'
#' A generic function.
#'
#' @param x a subspace class.
#' @param ... Extra parameters
#' @export
print.subspace <- function(x, ...) {
  message("An object of class subspace within",
      ifelse(x$transpose_flag, "transposed", ""),
      "X matrix with ",
      x$ndf,
      " samples and ",
      x$pdim,
      " features.\n")
  if (all(diff(x$components) == 1)) {
    message("Estimated components range from ",
        min(x$components),
        " to ",
        max(x$components),
        ".\n")
  } else {
      message("Estimated components range", x$components, ".\n")
  }
  invisible(x)
}

#' @title Scree plot of scaled eigenvalues of x and random noise matrix N
#'
#' @description This is a generic function for supspace class to plot scree
#'  plot of scaled eigenvalues of X and sampled scaled expected eigenvalues
#'   from Marcenko-Pastur distribution.
#' @param x A subsapce class.
#' @param changepoint A number. Estimated changepoint in dimension function.
#' @param annotation A number. Choose to label points up to annotation number.
#'  Set to 0 with no annotation.
#' @param verbose output message
#' @param ... Extra parameters
#' @examples
#' \donttest{
#' plot(Subspace, changepoint = 0, annotation = 15)
#' }
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom ggrepel geom_text_repel
#' @export

plot.subspace <- function(x,
                          changepoint = NULL,
                          annotation = NULL,
                          verbose = TRUE,
                          ...) {
# -----------------------
# Basic parameter set up
# -----------------------
  ndf             <- x$ndf
  pdim            <- x$pdim
  components      <- x$components
  rnk             <- max(components)
  var_correct     <- x$var_correct
  transpose_flag  <- x$transpose_flag
  irl             <- x$irl
  mp_irl          <- x$mp_irl
# ----------------------
# Check input parameters
# ----------------------
  if (missing(annotation)) {
    annotation <- rnk
    if (verbose) {
      cat("Anotating for all points. Set annotation = 0 to stop annotation.\n")
    }
  }
  if (!is.numeric(annotation)) {
    stop("Anotation must be numbers.\n")
  }
  if (!is.null(annotation) & min(annotation) < 0) {
    stop("Annotation number must be positive numbers.\n")
  }
  if (!is.null(annotation) & max(annotation) > rnk) {
    stop(paste0("Annotation number must be strictly less or equal",
            " to than maximum components.\n"))
  }
  if (length(annotation) > 1) {
    mark <- rep("", rnk)
    mark[annotation] <- annotation
  } else {
    mark <- ifelse(annotation >= components, components, "")
    if (annotation == 0) {
      mark <- rep("", rnk)
    }
  }
# -----------------------------------------------------------------------
# Scree plot for both eigenvalues of X and simulated eigenvalues of noise
# -----------------------------------------------------------------------
  scree <- ggplot() +
            geom_line(aes(x = dim, y = eigen), irl, colour = "black") +
            geom_point(aes(x = dim, y = eigen), irl, color = "red") +
            theme_minimal() +
            xlab("Dimension") +
            ylab("Eigenvalue Scaled") +
            theme(plot.title = element_text(hjust = 0.5)) +
            geom_text_repel(aes(x = irl$dim, y = irl$eigen, label = mark),
                            colour = "black",
                            size = 5) +
            ggtitle(
              paste0("Scree Plot\n",
                     "N = ",
                     ifelse(transpose_flag, pdim, ndf),
                     ", P = ",
                     ifelse(transpose_flag, ndf, pdim),
                     ", Var = ",
                     ifelse(is.null(var_correct), NA, round(var_correct, 2)))
              )
  if (!is.null(mp_irl)) {
    scree <- scree +
              geom_line(aes(x = dim, y = eigen), mp_irl, colour = "black") +
              geom_point(aes(x = dim, y = eigen), mp_irl, color = "blue")
  }

  if (!missing(changepoint)) {
    check_changepoint_input(changepoint, ndf, pdim, verbose)
    scree <- scree +
              ggtitle(
                paste0("Scree Plot\n",
                       "N = ",
                       ifelse(transpose_flag, pdim, ndf),
                       ", P = ",
                       ifelse(transpose_flag, ndf, pdim),
                       ", Var = ",
                       ifelse(is.null(var_correct), NA, round(var_correct, 2)),
                       ", Changepoint est = ",
                       changepoint))
  }

  return(scree)
}



#' @title Check components Input
#'
#' A generic function.
#'
#' @param components A series of right singular vectors to estimate.
#'  Components must be smaller or equal to min(nrow(x),ncol(x)).
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param verbose output message
#' @export
check_comp_input <- function(components, ndf, pdim, verbose = TRUE) {
  stopifnot(is.numeric(components))
  stopifnot(components %% 1 == 0)
  if (length(components) == 1) {
    if (components < 1) {
      stop("Components must be larger than 1.\n")
    }
    if (components > min(ndf, pdim)) {
      stop("Components out of bounds.\n")
    }
    if (verbose) {
      message("Calculating components from 1 to ", components, ".\n")
    }
    return(1:components)
  } else if (length(components) > 1) {
    if (min(components) < 1) {
      stop("Components must be larger than 1.\n")
    }
    if (max(components) > min(ndf, pdim)) {
      stop("Components out of bounds.\n")
    }
    if (verbose) {
      if (all(diff(components) == 1)) {
        message("Calculating components from ",
            min(components),
            " to ",
            max(components),
            ".\n")
      } else {
        message("Calculating components range", components, ".\n")
      }
    }
    return(components[order(components)])
  } else {
    stop("Invalid components input.\n")
  }
}

#' @title Check num_est_samples Input
#'
#' A generic function.
#'
#' @param num_est_samples Split data into num_est_samples-fold
#'  for parallel computation.
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param verbose output message
#' @export
check_num_est_samples_input <- function(num_est_samples,
                                        ndf, pdim, verbose = TRUE) {
  stopifnot(is.numeric(num_est_samples))
  stopifnot(num_est_samples %% 1 == 0)
  stopifnot(length(num_est_samples) == 1)
  if (num_est_samples < 0) {
    stop("num_est_samples must be positive")
  }
  else if (num_est_samples > min(ndf, pdim)) {
    stop("num_est_samples must be smaller or euqal to min(nrow(x), ncol(x))")
  }
}

#' @title Check changepoint Input
#'
#' A generic function.
#'
#' @param changepoint estimated changepoint in dimension function
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param verbose output message
#' @export
check_changepoint_input <- function(changepoint, ndf, pdim, verbose = TRUE) {
  stopifnot(is.numeric(changepoint))
  stopifnot(changepoint %% 1 == 0)
  stopifnot(length(changepoint) == 1)
  if (changepoint < 0) {
    stop("Changepoint must be positive")
  }
  else if (changepoint > min(ndf - 1, pdim - 1)) {
    stop("Changepoint must be strictly less than min(nrow(x), ncol(x))")
  }
}
