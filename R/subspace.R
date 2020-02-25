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
#' @param mp A logical value. If true, sample eigenvlaues from random noise
#'  matrix with mp distribution.
#' @param times Split data into times-fold for parallel computation.
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
#' \donttest{
#' x <- x_sim(n = 150, p = 100, ncc = 10, var = 2)
#' Subspace <- subspace(x)
#' Subspace <- subspace(x, components = 8:30)
#' Subspace <- subspace(x, components = c(2, 3, 6, 16), times = 10)
#' Subspace <- subspace(x, components = 1:20, mp= FALSE)
#' Subspace
#' plot(Subspace, changepoint = 0, annotation = 15)
#' }
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
#' @export

subspace <- function(x,
                     components = NA,
                     mp = TRUE,
                     times = NA,
                     verbose = TRUE,
                     ...) {
# ----------------------
# Check input parameters
# ----------------------
  if (is.null(x)) {
    stop("Invalid input x")
  } else {
    # Checking for components input
    if (missing(components)) {
      components <- 1:min(nrow(x), ncol(x))
      if (verbose) {
        message(paste0("No component specified. ",
            "Calculating full singular value decomposition instead.\n"))
      }
    } else {
      components <- check_comp_input(components, nrow(x), ncol(x), verbose)
    }
    # Checking for times input
    if (mp) {
      if (missing(times)) {
        times <- 0
      } else {
        check_times_input(times, nrow(x), ncol(x), verbose)
      }
    }
    # Check X matrix dimension
    params <- check_dim_matrix(x, rnk = max(components))
  }
# ----------------------
# Basic parameter set up
# ----------------------
  ndf             <- params$ndf
  pdim            <- params$pdim
  svr             <- params$svr
  rnk             <- params$rnk
  transpose_flag  <- params$transpose_flag
# ----------------------------------------
# Singular Value Decomposition of X matrix
# ----------------------------------------
  #Col Mean center Matrix
  if (transpose_flag) x <- t(x)

  if (rnk > pdim / 2) {
    tryCatch({
      x_std <- sweep(x, 2L, colMeans(x))
    }, warning = function(w) {
      message(paste0("Cann't allocate matrix in memory, try transforming matrix",
              " or a smaller proportion of eigenvalues.\n"))
    }, error = function(e) {
      message("Caught an error!\n")
    }
    )
      tmp <- svd(x_std)
  } else {
      tmp <- irlba(x, center = colMeans(x), nv = rnk)
  }

  irl     <- tibble(eigen = tmp$d[components]^2 / pdim, dim = components)
  sigma_a <- tmp$d[1:rnk]^2 / pdim
  v       <- tmp$v[, components]
  u       <- tmp$u[, components]
# --------------------------------------------
# Simulate noise eigenvalues from rmp function
# --------------------------------------------
  if (mp) {
    denominator <- marcenko_pastur_par(ndf, pdim, var = 1, svr = svr)$upper
    var_correct <- min(irl$eigen) / denominator
    if (verbose) {
      cat("The corrected variance of mp distribution is ",
          var_correct,
          ".\n",
          sep = "")
    }
    if (times == 0) {
      sim <- rmp(pdim,
                ndf = ndf,
                pdim = pdim,
                var = var_correct,
                svr = ndf / pdim)
    } else {
      ncores <- detectCores()
      registerDoParallel(ncores)
      sim <- foreach(1:times, .combine = c) %dopar% {
        tmp <- rmp(pdim / times,
                   ndf = ndf,
                   pdim = pdim,
                   var = var_correct,
                   svr = ndf / pdim)
      }
    }

    mp_irl <- tibble(eigen = sim[order(sim, decreasing = T)][components],
                     dim = components)
    sigma_mp <- sim[order(sim, decreasing = T)][1:rnk]
  }
  value <- (list(ndf  = ndf,
                 pdim = pdim,
                 components = components,
                 var_correct = ifelse(mp, var_correct, NA),
                 transpose_flag = transpose_flag,
                 irl = irl,
                 sigma_a = sigma_a,
                 mp_irl = switch(2 - mp, mp_irl, NULL),
                 sigma_mp = switch(2 - mp, sigma_mp, NULL),
                 v = v,
                 u = u))
  attr(value, "class") <- "subspace"
  value
}


#' @title Print subsapce
#'
#' A generic function.
#'
#' @param x a subspace class.
#' @param ... Extra parameters
#' @export
print.subspace <- function(x, ...) {
  cat("An object of class subspace within",
      ifelse(x$transpose_flag, "transposed", ""),
      "X matrix with",
      x$ndf,
      "samples and",
      x$pdim,
      "features.\n")
  if (all(diff(x$components) == 1)) {
    cat("Estimated components range from ",
        min(x$components),
        " to ",
        max(x$components),
        "\n")
  } else {
      cat("Estimated components range", x$components, "\n")
  }
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
                     round(var_correct, 2)))

  if (!is.null(mp_irl)) {
    scree <- scree +
              geom_line(aes(x = dim, y = eigen), mp_irl, colour = "black") +
              geom_point(aes(x = dim, y = eigen), mp_irl, color = "blue")
  }

  if (!missing(changepoint)) {
    scree <- scree +
              ggtitle(
                paste0("Scree Plot\n",
                       "N = ",
                       ifelse(transpose_flag, pdim, ndf),
                       ", P = ",
                       ifelse(transpose_flag, ndf, pdim),
                       ", Var = ",
                       round(var_correct, 2),
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
      cat("Calculating components from 1 to", components, ".\n")
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
        cat("Calculating components from",
            min(components),
            "to",
            max(components),
            ".\n")
      } else {
        cat("Calculating components range", components, "\n")
      }
    }
    return(components[order(components)])
  } else {
    stop("Invalid components input.\n")
  }
}

#' @title Check Times Input
#'
#' A generic function.
#'
#' @param times Split data into times-fold for parallel computation.
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param verbose output message
#' @export
check_times_input <- function(times, ndf, pdim, verbose = TRUE) {
  stopifnot(is.numeric(times))
  stopifnot(times %% 1 == 0)
  stopifnot(length(times) == 1)
  if (times < 0) {
    stop("Times must be positive")
  }
  else if (times > min(ndf - 1, pdim - 1)) {
    stop("Times must be strictly less than min(nrow(x), ncol(x))")
  }
}
