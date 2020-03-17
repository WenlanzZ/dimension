#' @title stand
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x numbers.
stand <- function(x) {
  n   <- nrow(x)
  p   <- ncol(x)
  xb  <- apply(x, 2, mean)
  xb  <- t(matrix(xb, p, n))
  x1  <- x - xb
  sigma <- t(x1) %*% (x1) / (n - 1)
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1 / sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}

#' @title matpower
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param a an original matrix
#' @param alpha power
matpower <- function(a, alpha) {
  a   <- (a + t(a)) / 2
  tmp <- eigen(a)
  return(tmp$vectors %*% diag((tmp$values)^alpha) %*%
           t(tmp$vectors))
}

#' @title mhat_cca
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x numbers.
#' @param y numbers.
#' @importFrom stats cov
mhat_cca <- function(x, y) {
  n   <- nrow(x)
  p   <- ncol(x)
  q   <- ncol(y)
  nsx <- matpower(cov(x), -0.5)
  nvy <- matpower(cov(y), -1)
  cxy <- cov(x, y)
  m   <- nsx %*% cxy %*% nvy %*% t(cxy) %*% nsx
  return(m)
}

#' @title mhat_ica
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x numbers.
mhat_ica <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  z <- stand(x)
  w <- apply(z^2, 1, sum)
  m <- t(z) %*% diag(w) %*% z / n - (p + 2) * diag(p)
  return(m %*% m)
}

#' @title slicing
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param y numbers.
#' @param h numbers.
#' @param n numbers.
#' @importFrom stats quantile
#' @importFrom stats quantile
slicing <- function(y, h, n) {
  if (length(levels(as.factor(y))) > h) {
    y_tilde <- rep(0, h + 1)
    y_tilde[1] <- min(y)
    for (h in 1:(h - 1)) {
      y_tilde[h + 1] <- quantile(y, h / h)
    }
  }
  if (length(levels(as.factor(y))) <= h) {
    h <- length(levels(as.factor(y)))
    y_tilde <- rep(0, h + 1)
    y_tilde[1] <- min(y)
    for (h in 1:(h - 1)) {
      y_tilde[h + 1] <- min(y[y > y_tilde[h]])
    }
  }
  y_tilde[h + 1] <- max(y) + 1
  prop <- rep(1, h)
  for (i in 1:h) {
    prop[i] <- sum((y >= y_tilde[i]) & (y < y_tilde[i + 1])) / n
  }
  res <- list()
  res$h <- h
  res$y_tilde <- y_tilde
  res$prop <- prop
  return(res)
}

#' @title mhat_dr
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x numbers.
#' @param y numbers.
#' @param h numbers.
mhat_dr <- function(x, y, h) {
  n <- nrow(x)
  p <- ncol(x)
  z <- stand(x)
  dy <- slicing(y, h)
  h <- dy$h
  y_tilde <- dy$y_tilde
  prop <- dy$prop
  ind <- matrix(0, n, h)
  zbar <- matrix(0, p, h)
  for (j in 1:(h - 1)) {
    ind[, j] <- ((y >= y_tilde[j]) & (y < y_tilde[j + 1]))
    zbar[, j] <- (t(z) %*% (ind[, j])) / sum(ind[, j])
  }
  ind[, h] <- (y >= y_tilde[h])
  zbar[, h] <- (t(z) %*% (ind[, h])) / sum(ind[, h])
  a <- matrix(0, p, p)
  b <- matrix(0, p, p)
  c <- 0
  for (q in 1:h) {
    z <- (t(z))[, ind[, q] == 1] - zbar[, q]
    a <- a + prop[q] * ((z %*% t(z) / (sum(ind[, q]) - 1)
          + zbar[, q] %*% t(zbar[, q])) %*%
          (z %*% t(z) / (sum(ind[, q]) - 1)
           + zbar[, q] %*% t(zbar[, q])) - diag(1, p))
    b <- b + sqrt(prop[j]) * (zbar[, q] %*% t(zbar[, q]))
    c <- c + sqrt(prop[j]) * (t(zbar[, q]) %*% zbar[, q])
  }
  c <- as.vector(c)
  m <- 2 * a + 2 * (b %*% b) + 2 * b * c
  return(m)
}

#' @title mhat_ksir
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x numbers.
#' @param phi numbers.
#' @param y numbers.
#' @param nslices nslices is needed for sufficient dimension reduction.
#' @importFrom stats var
mhat_ksir <- function(x = x, phi = phi, y = y, nslices = nslices) {
  n <- length(y)
  h <- nslices
  if ((missing(phi))) {
    p <- ncol(x)
    phi <- cbind(1, x)
    for (i in 1:p) {
      phi <- cbind(phi, x[, i] * x[, 1:i])
      }
  }
  oy <- order(y)
  ophi <- phi[oy, ]
  bslice <- h
  wslice <- n / bslice
  ephiy <- numeric()
  for (i in 1:bslice) {
    ephiy <- rbind(ephiy,
                   apply(ophi[((i - 1) * wslice + 1):(i * wslice), ],
                         2, mean))
  }
  ephi <- apply(phi, 2, mean)
  m <- var(t(t(ephiy) - ephi))
  vd <- svd(t(t(ephiy) - ephi))
  v <- vd$v[, order((vd$d)^2, decreasing = TRUE)]
  vc <- svd(diag(1, ncol(phi)) - v %*% t(v))$u
  lam <- (vd$d)^2 / (h - 1)
  res <- list()
  res$mhat <- m
  res$nx <- phi
  res$evectors <- cbind(v, vc)
  res$evalues <- lam
  return(res)
}

#' @title ladle
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param x A numeric real-valued matrix with n number of samples and
#'  p number of features (n > p).
#' @param y A vector.
#' @param nslices nslices is needed for sufficient dimension reduction.
#' @param nboot bootstrap sample size.
#' @param method options for method include: "pca", "cca", "ica", "sir",
#      "save", "dr" and "kernel sir", in which "pca" includes
#     (polynomial) kernel principal component analysis.
#' @param order the order of the polynomial in polynomial
#      principal component analysis.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{fn0:}{ fn0 is the bootstrap eigenvector variability.}
#'   \item{fn:}{ fn is the renormalized fn0.}
#'   \item{lam:}{ lam is the sample eigenvalues.}
#'   \item{phin:}{ phin is the renormalized and shifted sample eigenvalues.}
#'   \item{gn:}{ gn is the objective function ("y" in ladle plot).}
#'   \item{d:}{ d is the ladle estimator.}
#' }
#' @importFrom stats var
#' @importFrom dr dr
#' @export
ladle <- function(x = x, y = y, nslices = nslices, nboot = nboot,
                method = method, order = order) {
  p <- ncol(x)
  n <- nrow(x)
  z <- stand(x)
  if (p > 10) {
    r <- as.integer(p / log(p))
  } else {
      r <- p - 1
  }
  if (missing(y)) {
    obs <- x
  } else {
      obs <- cbind(y, x)
  }
  if (!(missing(nslices))) h <- nslices
  if (missing(nboot)) nboot <- as.integer(n / 2)
  if (method == "pca") {
    if (missing(order)) {
      nx <- x
      } else {
      nx <- x
      if (order > 1) {
        for (j in 2:order) nx <- cbind(nx, x^j)
        }
      }
    mhat <- var(nx)
    obs <- nx
    p <- ncol(nx)
  }
  if (method == "cca") {
    if (is.vector(y) == TRUE) y <- matrix(y, , 1)
    q <- ncol(y)
    if (q > 10) {
      r <- min(r, as.integer(q / log(q)))
    }
    mhat <- mhat_cca(x, y)
  }
  if (method == "ica") mhat <- mhat_ica(x)
  if (method == "sir") {
    obs <- cbind(y, z)
    sdr <- dr(y ~ z, method = method, nslices = h, numdir = r)
    mhat <- sdr$evectors %*% diag(sdr$evalues) %*% t(sdr$evectors)
    r <- min(r, h - 1)
  }
  if (method == "save") {
    obs <- cbind(y, z)
    sdr <- dr(y ~ z, method = method, nslices = h, numdir = r)
    mhat <- sdr$evectors %*% diag(sdr$evalues) %*% t(sdr$evectors)
  }
  if (method == "dr") {
    mhat <- mhat_dr(z, y, h)
    }
  if (method != "kernel sir") {
    bhat <- eigen(mhat)$vectors[, 1:r]
    lam <- eigen(mhat)$values[1:(r + 1)]
  } else {
    m <- mhat_ksir(x = x, y = y, nslices = h)
    r <- min(r, h - 2)
    #only consider nonzero sample eigenvalues
    lam <- m$evalues[1:(r + 1)]
    nx <- (m$nx) %*% m$evectors[, 1:p]
    obs <- cbind(y, nx)
    bhat <- diag(1, p)
  }
  fn0 <- rep(0, r + 1)
  for (j in 1:nboot) {
    u <- round(runif(n, min = -0.5, max = n + 0.5))
    bs <- obs[u, ]
    if (method == "pca") {
      mstar <- var(bs)
      }
    if (method == "cca") {
      mstar <- mhat_cca(bs[, - (1:q)], bs[, 1:q])
      }
    if (method == "ica") {
      mstar <- mhat_ica(bs)
      }
    if ((method == "sir") | (method == "save")) {
      bsdr <- dr(bs[, 1] ~ bs[, -1], method = method, nslices = h)
      vsdr <- qr.Q(qr(bsdr$evectors))
      mstar <- vsdr %*% diag(bsdr$evalues) %*% t(vsdr)
    }
    if (method == "dr") {
      mstar <- mhat_dr(stand(bs[, -1]), bs[, 1], h)
      }
    if (method != "kernel sir") {
      b_star <- eigen(mstar)$vectors[, 1:r]
    } else {
      b_star <- mhat_ksir(phi = bs[, -1],
                          y = bs[, 1],
                          nslices = h)$evectors[, 1:r]
      }
    for (i in 1:r) {
      fn0[i + 1] <- fn0[i + 1] + 1 - abs(det(t(b_star[, 1:i]) %*% bhat[, 1:i]))
    }
  }
  fn0 <- fn0 / nboot
  res <- list(0)
  res$fn0 <- fn0
  res$fn <- fn0 / (1 + sum(fn0))
  res$lam <- lam
  res$phin <- lam / (1 + sum(lam))
  res$gn <- res$fn + res$phin
  res$d <- which.min(res$gn) - 1
  return(res)
}
