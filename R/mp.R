#' @title Marcenko-Pastur distribution within package "RMTsata".
#' @description All credit to Iain M. Johnstone, Zongming Ma,
#'  Patrick O. Perry and Morteza Shahram
#'  (https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param var Population variance.
#' @param svr ndfxpdim.
#' @export
MarchenkoPasturPar <- function(ndf = NA, pdim = NA, var = 1, svr = ndf / pdim) {
    gamma          <- svr

    inv_gamma_sqrt <- sqrt(gamma)
    a              <- var * (1 - inv_gamma_sqrt)^2
    b              <- var * (1 + inv_gamma_sqrt)^2

    list(lower = a, upper = b)
}

#' @title dmp within package "RMTsata".
#'
#' @description All credit to Iain M. Johnstone, Zongming Ma,
#'  Patrick O. Perry and Morteza Shahram
#'  (https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param x Vector of quantiles.
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param var Population variance.
#' @param svr ndf/pdim.
#' @param log A logical value.
#' @export
dmp <- function(x, ndf = NA, pdim = NA, var = 1, svr = ndf / pdim,
                log = FALSE) {
    gamma  <- svr

    params <- MarchenkoPasturPar(ndf, pdim, var, svr)
    a      <- params$lower
    b      <- params$upper

    if (!log) {
        # we have to handle + / - zero carefully when gamma = 1
        d <- ifelse(gamma == 1 & x == 0 & 1 / x > 0, Inf,
                 ifelse(x <= a | x >= b, 0,
                     suppressWarnings(
                         gamma / (2 * pi * var * x) * sqrt((x - a) * (b - x))
                        )))
    } else {
        d <- ifelse(gamma == 1 & x == 0 & 1 / x > 0, Inf,
                 ifelse(x <= a | x >= b, -Inf,
                     suppressWarnings(
                         log(gamma) -
                         (log(2) + log(pi) + log(var) + log(x))
                         +
                         0.5 * log(x - a)
                         +
                         0.5 * log(b - x)
                        )))
    }
    d
}

#' @title pmp within package "RMTsata".
#'
#' @description All credit to Iain M. Johnstone, Zongming Ma,
#'  Patrick O. Perry and Morteza Shahram
#'  (https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param q Vector of quantiles.
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param var Population variance.
#' @param svr ndfxpdim.
#' @param lower_tail A logical value.
#' @param log.p A logical value.
#' @export
pmp <- function(q, ndf = NA, pdim = NA, var = 1, svr = ndf / pdim,
                 lower_tail = TRUE, log.p = FALSE) {
    params <- MarchenkoPasturPar(ndf, pdim, var, svr)
    a <- params$lower
    b <- params$upper

    f <- function(x) {
             dmp(x, ndf, pdim, var, svr)
         }

    if (lower_tail) {
        p <- ifelse(q <= a, 0,
                 ifelse(q >= b, 1,
                     integrate(f, a, q)$value
                )
            )
        p <- ifelse(svr < 1 && q >= 0, p + (1 - svr), p)
    } else {
        p <- ifelse(q >= b, 0,
                 ifelse(q <= a, min(1, svr),
                     integrate(f, q, b)$value
                )
            )
        p <- ifelse(svr < 1 && q <= 0, p + (1 - svr), p)
    }

    res <- ifelse(log.p, log(p), p)
    res
}
pmp <- Vectorize(pmp)

#' @title qmp within package "RMTsata".
#'
#' @description All credit to Iain M. Johnstone,
#'  Zongming Ma, Patrick O. Perry and Morteza Shahram
#'  (https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param p Vector of probabilities.
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param var Population variance.
#' @param svr ndfxpdim.
#' @param lower_tail A logical value.
#' @param log.p A logical value.
#' @export
qmp <- function(p, ndf = NA, pdim = NA, var = 1, svr = ndf / pdim,
                 lower_tail = TRUE, log.p = FALSE) {
    p      <- ifelse(log.p, exp(p), p)
    p      <- ifelse(lower_tail, p, 1 - p)
    params <- MarchenkoPasturPar(ndf, pdim, var, svr)

    q      <- NULL
    if (p <= 0) {
        q <- ifelse(svr <= 1, -0, params$lower)
    } else if (p >= 1) {
        q <- params$upper
    }

    if (svr < 1) {
        if (p < 1 - svr) {
            q <- -0
        } else if (p == 1 - svr) {
            q <- 0
        }
    }

    if (is.null(q)) {
        f <- function(x) {
                 pmp(x, ndf, pdim, var, svr) - p
             }
        q <- uniroot(f, interval = c(params$lower, params$upper))$root
    }

    q
}
qmp <- Vectorize(qmp)

#' @title rmp within package "RMTsata".
#'
#' @description All credit to Iain M. Johnstone, Zongming Ma,
#'  Patrick O. Perry and Morteza Shahram
#'   (https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param n Number of observations.
#' @param ndf The number of degrees of freedom of x.
#' @param pdim The number of dimensions of x.
#' @param var Population variance.
#' @param svr ndfxpdim.
#' @importFrom stats runif
#' @export
rmp <- function(n, ndf = NA, pdim = NA, var = 1, svr = ndf / pdim) {
    u <- runif(n)
    qmp(u, ndf, pdim, var, svr)
}
