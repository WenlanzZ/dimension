#' @title hoff_gibbs_sample
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param e An \tt{R} object.
#' @param nu An \tt{R} object.
#' @param nv An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param lor An \tt{R} object.
#' @param llb An \tt{R} object.
#' @param lub An \tt{R} object.
#' @export
hoff_gibbs_sample <- function(e,
                              nu,
                              nv,
                              phi,
                              mu,
                              psi,
                              lor,
                              llb = 50,
                              lub = 2000) {
  e2 <- sum(e^2)
  n <- dim(e)[2]
  m <- dim(e)[1]
  e0 <- t(nu) %*% e %*% nv
  n0 <- dim(e0)[2]
  m0 <- dim(e0)[1]
  svd_e0 <- svd(e0)
  z <- svd_e0$d^2
  e02 <- sum(z)
  theta <- z / e02
  del <- 0
  u <- rep(0, m0)
  v <- rep(0, n0)

  ### critical value for approximating bessel function
  amax <- svd_e0$d[1]^2 * phi^2 * e02 * max(theta)
  lcrit <- .5 * (sqrt((n0 / 2 + 1)^2 + 4 * (amax - n0 / 2)) - (n0 / 2 + 1))
  lmax <- 1.25 * (lcrit)
  lmax <- min(max(lmax, llb), lub)
  lmax <- 2 * as.integer(lmax / 2)
  lseq <- 0:lmax

  ### sums
  la <- lgamma(m0 / 2) - lgamma(m0 / 2 + lseq) - lgamma(lseq + 1) -
    lseq * log(4) + lcr(theta, lmax)
  if (mu == 0) {
    lb <- .5 * log(psi / (phi + psi)) + lseq * log(1 / (phi + psi)) +
      (lgamma(2 * lseq + 1) - lgamma(lseq + 1) - lseq * log(2))
  }
  if (mu != 0) {
    lb <- ln2moment(mu * psi / (phi + psi), sqrt(1 / (psi + phi)), lmax) +
      .5 * log(psi / (phi + psi)) -
      .5 * mu^2 * psi * phi / (phi + psi)
  }
  lc <- lseq * (log(e02 * phi^2))
  lt <- la + lb + lc
  lpe_r10 <- max(lt) + log(sum(exp(lt - max(lt))))
  s <- rbinom(1, 1, 1 / (1 + exp( - (lor + lpe_r10))))
  if (s == 1) {
    #### sample del
    pdl <- exp(lt - max(lt))
    l <- sample(x = lseq, size = 1, prob = pdl)
    if (l > 0 & mu != 0) {
      del <- rxl(mu * psi / (psi + phi), sqrt(1 / (phi + psi)), l)
    }
    if (l > 0 & mu == 0) {
      del <- del <- sqrt(rchisq(1, 2 * l + 1) / (phi + psi))
    }
    if (l == 0) {
      del <- rnorm(1, mu * psi / (psi + phi), sqrt(1 / (phi + psi)))
    }
    #### sample u, v
    tmp <- ruv_A(phi * del * e0)
    u <- tmp$u
    v <- tmp$v
  }
  list(del = del, u = u, v = v)
}

#' @title lcr
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param theta An \tt{R} object.
#' @param L An \tt{R} object.
lcr <- function(theta, L) {
  .C("lcr", as.double(theta), as.integer(L), as.integer(length(theta)),
    lr = double(L + 1)
  )$lr
}

#' @title rW
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param kap An \tt{R} object.
#' @param m An \tt{R} object.
rW <- function(kap, m) {
  .C("rW", kap = as.double(kap), m = as.integer(m), w = double(1))$w
}

#' @title ln2moment
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param mu An \tt{R} object.
#' @param sigma An \tt{R} object.
#' @param lmax An \tt{R} object.
ln2moment <- function(mu, sigma, lmax) {
  .C("ln2moment", as.double(abs(mu)), as.double(sigma^2), as.integer(lmax),
    l2mom = double(lmax + 1)
  )$l2mom
}

#' @title labc
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param theta An \tt{R} object.
#' @param lmax An \tt{R} object.
#' @param m An \tt{R} object.
#' @param e2 An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
labc <- function(theta, lmax, m, e2, phi, mu, psi) {
  .C("labc",
    theta = as.double(theta), L = as.integer(lmax), m = as.integer(m),
    n = as.integer(length(theta)), e2 = as.double(e2), phi = as.double(phi),
    mu = as.double(abs(mu)), psi = as.double(psi), lt = double(lmax + 1),
    lmom = double(2 * lmax + 1), lr = double(lmax + 1)
  )
}

#' @title lsumt
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param theta An \tt{R} object.
#' @param lmax An \tt{R} object.
#' @param m An \tt{R} object.
#' @param e2 An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
lsumt <- function(theta, lmax, m, e2, phi, mu, psi) {
  .C("lsumt",
    theta = as.double(theta), L = as.integer(lmax), m = as.integer(m),
    n = as.integer(length(theta)), e2 = as.double(e2), phi = as.double(phi),
    mu = as.double(abs(mu)), psi = as.double(psi), lt = double(lmax + 1)
  )$lt
}

#' @title rvmf
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param kmu An \tt{R} object.
rvmf <- function(kmu) {
  kap <- sqrt(sum(kmu^2))
  mu <- kmu / kap
  m <- length(mu)

  if (m == 1) {
    u <- (-1)^rbinom(1, 1, 1 / (1 + exp(2 * kap * mu)))
  }

  if (m > 1) {
    W <- rW(kap, m)
    V <- rnorm(m - 1)
    V <- V / sqrt(sum(V^2))
    x <- c((1 - W^2)^.5 * t(V), W)
    u <- cbind(Null(mu), mu) %*% x
  }
  u
}

#' @title hoff_gibbs_UVD_fixedrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @export
hoff_gibbs_UVD_fixedrank <- function(U, V, D, Y, phi, mu, psi) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  for (j in resample((1:n)[diag(Dl) != 0])) {
    nu <- Null(Ul[, -j])
    nv <- Null(Vl[, -j])
    if (sum(Dl[-j, -j] != 0) == 0) {
      nu <- diag(nrow = m)
      nv <- diag(nrow = n)
    }
    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    Vl[, j] <- nv %*% rvmf(Dl[j, j] * t(nv) %*% t(e) %*% Ul[, j] * phi)
    Ul[, j] <- nu %*% rvmf(Dl[j, j] * t(nu) %*% e %*% Vl[, j] * phi)

    mn <- (t(Ul[, j]) %*% e %*% Vl[, j] * phi + mu * psi) / (phi + psi)
    se <- sqrt(1 / (phi + psi))
    d <- rnorm(1, mn, se)
    Dl[j, j] <- d
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_varrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param nsamp An \tt{R} object.
#' @export
hoff_gibbs_UVD_varrank <- function(U, V, D, Y, phi, mu, psi, nsamp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  for (j in resample(1:n, size = nsamp)) {
    nu <- Null(Ul[, -j])
    nv <- Null(Vl[, -j])
    if (sum(Dl[-j, -j] != 0) == 0) {
      nu <- diag(nrow = m)
      nv <- diag(nrow = n)
    }
    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    lor <- log((sum(diag(Dl)[-j] != 0) + 1) / (n - sum(diag(Dl)[-j] != 0)))
    tmp <- hoff_gibbs_sample(e, nu, nv, phi, mu, psi, lor)
    Ul[, j] <- nu %*% tmp$u
    Vl[, j] <- nv %*% tmp$v
    Dl[j, j] <- tmp$del
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_pvarrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param nsamp An \tt{R} object.
#' @param uperp An \tt{R} object.
#' @param vperp An \tt{R} object.
#' @export
hoff_gibbs_UVD_pvarrank <- function(U,
                                    V,
                                    D,
                                    Y,
                                    phi,
                                    mu,
                                    psi,
                                    nsamp,
                                    uperp,
                                    vperp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]

  for (j in resample(1:n, size = nsamp)) {
    nu <- Null(cbind(uperp, Ul[, -j]))
    nv <- Null(cbind(vperp, Vl[, -j]))

    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    lor <- log((sum(diag(Dl)[-j] != 0) + 1) / (n - sum(diag(Dl)[-j] != 0)))
    tmp <- hoff_gibbs_sample(e, nu, nv, phi, mu, psi, lor)
    Ul[, j] <- nu %*% tmp$u
    Vl[, j] <- nv %*% tmp$v
    Dl[j, j] <- tmp$del
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_pfixedrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param uperp An \tt{R} object.
#' @param vperp An \tt{R} object.
#' @export
hoff_gibbs_UVD_pfixedrank <- function(U,
                                      V,
                                      D,
                                      Y,
                                      phi,
                                      mu,
                                      psi,
                                      uperp,
                                      vperp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]

  for (j in resample((1:n)[diag(Dl) != 0])) {
    nu <- Null(cbind(uperp, Ul[, -j]))
    nv <- Null(cbind(vperp, Vl[, -j]))

    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    Vl[, j] <- nv %*% rvmf(Dl[j, j] * t(nv) %*% t(e) %*% Ul[, j] * phi)
    Ul[, j] <- nu %*% rvmf(Dl[j, j] * t(nu) %*% e %*% Vl[, j] * phi)

    mn <- (t(Ul[, j]) %*% e %*% Vl[, j] * phi + mu * psi) / (phi + psi)
    se <- sqrt(1 / (phi + psi))
    d <- rnorm(1, mn, se)
    Dl[j, j] <- d
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_rc_varrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param nsamp An \tt{R} object.
#' @param up An \tt{R} object.
#' @param vp An \tt{R} object.
#' @export
hoff_gibbs_UVD_rc_varrank <- function(U,
                                      V,
                                      D,
                                      Y,
                                      phi,
                                      mu,
                                      psi,
                                      nsamp,
                                      up,
                                      vp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  for (j in resample(3:n, size = nsamp)) {
    nu <- Null(cbind(up, Ul[, -j]))
    nv <- Null(cbind(vp, Vl[, -j]))

    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    lor <- log((sum(diag(Dl)[-j] != 0) + 1) / (n - sum(diag(Dl)[-j] != 0)))
    tmp <- hoff_gibbs_sample(e, nu, nv, phi, mu, psi, lor)
    Ul[, j] <- nu %*% tmp$u
    Vl[, j] <- nv %*% tmp$v
    Dl[j, j] <- tmp$del
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_rc_fixedrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param up An \tt{R} object.
#' @param vp An \tt{R} object.
#' @export
hoff_gibbs_UVD_rc_fixedrank <- function(U,
                                        V,
                                        D,
                                        Y,
                                        phi,
                                        mu,
                                        psi,
                                        up,
                                        vp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]

  for (j in resample((1:n)[diag(Dl) != 0])) {
    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])

    if (j != 2) {
      nv <- Null(cbind(vp, Vl[, -j]))
      Vl[, j] <- nv %*% rvmf(Dl[j, j] * t(nv) %*% t(e) %*% Ul[, j] * phi)
    }

    if (j != 1) {
      nu <- Null(cbind(up, Ul[, -j]))
      Ul[, j] <- nu %*% rvmf(Dl[j, j] * t(nu) %*% e %*% Vl[, j] * phi)
    }

    mn <- (t(Ul[, j]) %*% e %*% Vl[, j] * phi + mu * psi) / (phi + psi)
    se <- sqrt(1 / (phi + psi))
    d <- rnorm(1, mn, se)
    Dl[j, j] <- d
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_perp_varrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param nsamp An \tt{R} object.
#' @param up An \tt{R} object.
#' @param vp An \tt{R} object.
#' @export
hoff_gibbs_UVD_perp_varrank <- function(U,
                                        V,
                                        D,
                                        Y,
                                        phi,
                                        mu,
                                        psi,
                                        nsamp,
                                        up,
                                        vp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  for (j in resample(1:n, size = nsamp)) {
    nu <- Null(cbind(up, Ul[, -j]))
    nv <- Null(cbind(vp, Vl[, -j]))

    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    lor <- log((sum(diag(Dl)[-j] != 0) + 1) / (n - sum(diag(Dl)[-j] != 0)))
    tmp <- hoff_gibbs_sample(e, nu, nv, phi, mu, psi, lor)
    Ul[, j] <- nu %*% tmp$u
    Vl[, j] <- nv %*% tmp$v
    Dl[j, j] <- tmp$del
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title hoff_gibbs_UVD_perp_fixedrank
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param U An \tt{R} object.
#' @param V An \tt{R} object.
#' @param D An \tt{R} object.
#' @param Y An \tt{R} object.
#' @param phi An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param psi An \tt{R} object.
#' @param up An \tt{R} object.
#' @param vp An \tt{R} object.
#' @export
hoff_gibbs_UVD_perp_fixedrank <- function(U,
                                          V,
                                          D,
                                          Y,
                                          phi,
                                          mu,
                                          psi,
                                          up,
                                          vp) {
  Ul <- U
  Vl <- V
  Dl <- D
  m <- dim(Y)[1]
  n <- dim(Y)[2]

  for (j in resample((1:n)[diag(Dl) != 0])) {
    e <- Y - Ul[, -j] %*% Dl[-j, -j] %*% t(Vl[, -j])
    nv <- Null(cbind(vp, Vl[, -j]))
    Vl[, j] <- nv %*% rvmf(Dl[j, j] * t(nv) %*% t(e) %*% Ul[, j] * phi)
    nu <- Null(cbind(up, Ul[, -j]))
    Ul[, j] <- nu %*% rvmf(Dl[j, j] * t(nu) %*% e %*% Vl[, j] * phi)

    mn <- (t(Ul[, j]) %*% e %*% Vl[, j] * phi + mu * psi) / (phi + psi)
    se <- sqrt(1 / (phi + psi))
    d <- rnorm(1, mn, se)
    Dl[j, j] <- d
  }
  assign("U", Ul, env = parent.env(environment()))
  assign("V", Vl, env = parent.env(environment()))
  assign("D", Dl, env = parent.env(environment()))
}

#' @title resample
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param x An \tt{R} object.
#' @param size An \tt{R} object.
resample <- function(x, size, ...) {
  if (length(x) <= 1) {
    if (!missing(size) && size == 0) x[FALSE] else x
  } else {
    sample(x, size, ...)
  }
}

#' @title Null
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param M An \tt{R} object.
#' @export
Null <- function(M) {
  tmp <- qr(M)
  set <- if (tmp$rank == 0) {
    1:ncol(M)
  } else {
    -(1:tmp$rank)
  }
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}

#' @title lpY_pois
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param Y An \tt{R} object.
#' @param Theta An \tt{R} object.
lpY_pois <- function(Y, Theta) {
  Y * Theta - exp(Theta) - lgamma(Y + 1)
}

#' @title lpY_binary
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param Y An \tt{R} object.
#' @param Theta An \tt{R} object.
lpY_binary <- function(Y, Theta) {
  Y * Theta - log(1 + exp(Theta))
}

#' @title lpY_tri
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param Y An \tt{R} object.
#' @param Theta An \tt{R} object.
lpY_tri <- function(Y, Theta) {
  Y * Theta + t(Y) * t(Theta) - log(1 + exp(Theta) + exp(t(Theta)))
}

#' @title XB
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param X An \tt{R} object.
#' @param beta An \tt{R} object.
XB <- function(X, beta) {
  tmp <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
  for (k in 1:length(beta)) {
    tmp <- tmp + beta[k] * X[, , k]
  }
  tmp
}

#' @title x.X
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param X An \tt{R} object.
x.X <- function(X) {
  x <- NULL
  for (k in 1:dim(X)[3]) {
    x <- cbind(x, c(X[, , k]))
  }
  x
}

#' @title hoff_rmvnorm
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param mu An \tt{R} object.
#' @param Sig2 An \tt{R} object.
hoff_rmvnorm <- function(mu, Sig2) {
  R <- t(chol(Sig2))
  R %*% (rnorm(length(mu), 0, 1)) + mu
}

#' @title ldmvnorm
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param y An \tt{R} object.
#' @param mu An \tt{R} object.
#' @param Sig An \tt{R} object.
ldmvnorm <- function(y, mu, Sig) {
  c(-(length(y) / 2) * log(2 * pi) - .5 * log(det(Sig)) - .5 *
    t(y - mu) %*% solve(Sig) %*% (y - mu))
}

#' @title ruv_A
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param A An \tt{R} object.
#' @param nsamp An \tt{R} object.
ruv_A <- function(A, nsamp = 25) {
  m <- dim(A)[1]
  n <- dim(A)[2]
  tmp <- svd(A)

  ### sample a mode
  j <- sample(x = 1:n, size = 1, prob = exp(tmp$d - max(tmp$d)))
  v <- (-1)^rbinom(1, 1, .5) * tmp$v[, j]

  ### sample around a mode
  for (i in 1:nsamp) {
    u <- rvmf(A %*% v)
    v <- rvmf(t(A) %*% u)
  }

  list(u = c(u), v = v)
}

#' @title pK_X_bic
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param X An \tt{R} object.
pK_X_bic <- function(X) {
  N <- dim(X)[2]
  d <- dim(X)[1]
  S <- t(X) %*% X / N
  lam <- eigen(S)$val / N
  lp <- NULL
  for (k in 0:(d - 1)) {
    m <- d * k - k * (k + 1) / 2
    lp <- c(lp, -.5 * N * sum(log(lam[seq(1, k, length = k)])) - .5 * (m + k) * log(N) -
      .5 * N * (d - k) * log(mean(lam[seq(k + 1, d, length = d - k)])))
  }
  k <- d
  m <- d * k - k * (k + 1) / 2
  lp <- c(lp, -.5 * N * sum(log(lam[seq(1, k, length = k)])) - .5 * (m + k) * log(N))
  pk <- exp(lp - max(lp))
  pk / sum(pk)
}

#' @title pK_X_lap
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @param X An \tt{R} object.
pK_X_lap <- function(X) {
  N <- dim(X)[2]
  d <- dim(X)[1]
  S <- t(X) %*% X / N
  lam <- eigen(S)$val / N
  lp <- NULL

  for (k in 0:(d - 1)) {
    laz <- 0
    for (i in seq(1, k, length = k)) {
      for (j in seq(i + 1, d, length = d - i)) {
        laz <- laz + log(1 / lam[j] - 1 / lam[i]) + log(lam[i] - lam[j]) + log(N)
      }
    }

    lpu <- 0
    if (k > 0) {
      lpu <- -k * log(2) + sum(lgamma((d - (1:k) + 1) / 2) - .5 * (d - (1:k) + 1) * log(pi))
    }

    m <- d * k - k * (k + 1) / 2
    lp <- c(lp, lpu - .5 * N * sum(log(lam[seq(1, k, length = k)])) -
      .5 * N * (d - k) * log(mean(lam[seq(k + 1, d, length = d - k)])) +
      .5 * (m + k) * log(2 * pi) - .5 * laz - .5 * k * log(N))
  }

  k <- d
  laz <- 0
  for (i in seq(1, k, length = k)) {
    for (j in seq(i + 1, d, length = d - i)) {
      laz <- laz + log(1 / lam[j] - 1 / lam[i]) + log(lam[i] - lam[j]) + log(N)
    }
  }

  lpu <- 0
  if (k > 0) {
    lpu <- -k * log(2) + sum(lgamma((d - (1:k) + 1) / 2) - .5 * (d - (1:k) + 1) * log(pi))
  }

  m <- d * k - k * (k + 1) / 2
  lp <- c(lp, lpu - .5 * N * sum(log(lam[seq(1, k, length = k)])) +
    .5 * (m + k) * log(2 * pi) - .5 * laz - .5 * k * log(N))

  pk <- exp(lp - max(lp))
  pk / sum(pk)
}

#' @title ln2moment
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param mu An \tt{R} object.
#' @param sigma An \tt{R} object.
#' @param lmax An \tt{R} object.
ln2moment <- function(mu, sigma, lmax) {
  .C("ln2moment", as.double(abs(mu)), as.double(sigma^2), as.integer(lmax),
    l2mom = double(lmax + 1)
  )$l2mom
}

#' @title roots_cubic
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param a2 An \tt{R} object.
#' @param a1 An \tt{R} object.
#' @param a0 An \tt{R} object.
roots_cubic <- function(a2, a1, a0) {
  Q <- (3 * a1 - a2^2) / 9
  R <- (9 * a2 * a1 - 27 * a0 - 2 * a2^3) / 54
  D <- Q^3 + R^2
  S <- (R + sqrt(0i + D))^(1 / 3)
  T <- (R - sqrt(0i + D))^(1 / 3)

  c(
    -a2 / 3 + (S + T),
    -a2 / 3 - .5 * (S + T) + 1i * .5 * sqrt(3) * (S - T),
    -a2 / 3 - .5 * (S + T) - 1i * .5 * sqrt(3) * (S - T)
  )
}

#' @title roots_quartic
#' @description All credit to Hoff, Peter D. "Model averaging
#' and dimension selection for the singular value decomposition."
#' Journal of the American Statistical Association 102.478 (2007): 674-685.
#' @useDynLib dimension
#' @param a3 An \tt{R} object.
#' @param a2 An \tt{R} object.
#' @param a1 An \tt{R} object.
#' @param a0 An \tt{R} object.
roots_quartic <- function(a3, a2, a1, a0) {
  z3 <- roots_cubic(-a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1^2 - a3^2 * a0)
  y <- Re(z3[Im(z3) == 0])[1]

  R <- sqrt(.25 * a3^2 - a2 + y)
  D <- sqrt(0i + .25 * 3 * a3^2 - R^2 - 2 * a2 + .25 * (4 * a3 * a2 - 8 * a1 - a3^3) / R)
  e <- sqrt(0i + .25 * 3 * a3^2 - R^2 - 2 * a2 - .25 * (4 * a3 * a2 - 8 * a1 - a3^3) / R)

  c(
    -a3 / 4 + R / 2 + D / 2,
    -a3 / 4 + R / 2 - D / 2,
    -a3 / 4 - R / 2 + e / 2,
    -a3 / 4 - R / 2 - e / 2
  )
}

ldxl <- function(x, mu, sigma, l, ln2m = ln2moment(mu, sigma, l)[l + 1]) {
  l * log(x^2) - log(sigma) - .5 * log(2 * pi) - .5 * ((x - mu) / sigma)^2 - ln2m
}

rxl <- function(mu, sigma, l, nu = 1) {
  theta <- .5 * mu * (1 + sqrt(1 + 8 * l * sigma^2 / mu^2))
  tau <- 1 / sqrt(1 / sigma^2 + 2 * l / (theta^2))

  a <- (-2 * theta - mu)
  b <- (theta^2 + 2 * mu * theta + nu * tau^2 + sigma^2 * (-nu - 2 * l - 1))
  c <- (-mu * theta^2 + (nu + 4 * l + 1) * sigma^2 * theta - mu * nu * tau^2)
  d <- (-2 * l * sigma^2 * theta^2 - 2 * l * nu * sigma^2 * tau^2)

  z4 <- roots_quartic(a, b, c, d)
  xc <- Re(z4[Im(z4) == 0])

  ln2m <- ln2moment(mu, sigma, l)[l + 1]

  lM <- max(ldxl(xc, mu, sigma, l, ln2m = ln2m) - (-log(tau) + log(dt((xc - theta) / tau, nu))))

  samp <- T
  while (samp) {
    x <- theta + rt(1, nu) * tau
    lrratio <- ldxl(x, mu, sigma, l = l, ln2m = ln2m) - (-log(tau) + log(dt((x - theta) / tau, nu))) - lM
    samp <- (log(runif(1)) > lrratio)
  }
  x
}

#' @title hoff
#' @description All credit to Luo, W. and Li, B. (2016),
#' Combining Eigenvalues and Variation of Eigenvectors for Order
#' Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
#' @param y A numeric real-valued matrix with n number of samples and
#'  p number of features (n > p).
#' @param NSCAN A numberic value decide MCMC sample size (size = NSCAN / 10)
#' @export
hoff <- function(y = y, NSCAN = 10000) {
  # source("svd.r")
  ##### constants
  m <- dim(y)[1]
  n <- dim(y)[2]
  #####

  ##### hyperparameters
  sY <- svd(y)
  s20.est <- var(c(y))
  t20.est <- 0
  mu0.est <- 0
  for (k in 1:n) {
    ks <- seq(1,  k,  length = k)
    s20.est <- c(s20.est,
                 var(c((y - sY$u[, ks] %*%
                          diag(sY$d[ks], nrow = k) %*%
                          t(sY$v[, ks])))))
    t20.est <- c(t20.est, var(sY$d[ks]))
    mu0.est <- c(mu0.est, mean(sY$d[ks]))
  }
  t20.est[2] <- 0
  ## prior for phi
  nu0 <- 2
  s20 <- mean(s20.est)
  ## prior for psi
  eta0 <- 2
  t20 <- mean(t20.est)
  ## prior for mu
  # kap0 <- 1
  # mu0 <- mean(mu0.est)
  mu0 <- mean(mu0.est)
  premu0 <- 1 / var(mu0.est)
  ##### starting values
  K0 <- 0
  U <- matrix(0, m, n)
  V <- matrix(0, n, n)
  U[, seq(1, K0, length = K0)] <- sY$u[, seq(1, K0, length = K0)]
  V[, seq(1, K0, length = K0)] <- sY$v[, seq(1, K0, length = K0)]
  D <- diag(c(sY$d[seq(1, K0, length = K0)], rep(0, n - K0)))
  phi <-  1 / s20
  psi <- 1 / t20
  mu <- mu0
  #####

  ##### MCMC
  # NSCAN <- 10000
  odens <- 10
  OUT <- matrix(nrow = NSCAN / odens, ncol = 5)
  colnames(OUT) <- c("sample id", "dim", "inv_phi", "mu", "inv_psi")
  MSE <- NULL
  M.ps <- y * 0
  D.ps <- rep(0, n)
  for (ns in 1:NSCAN) {
    hoff_gibbs_UVD_varrank(U, V, D, y, phi, mu, psi, min(n, 10))
    ### fixed rank update
    hoff_gibbs_UVD_fixedrank(U, V, D, y, phi, mu, psi)
    ### update phi
    phi <- rgamma(1, (nu0 + m * n) / 2,
                  (nu0 * s20 + sum((y - U %*% D %*% t(V))^2)) / 2)
    ### update mu, psi
    mu <- rnorm(1,(premu0 * mu0 + psi * sum(diag(D))) /
                  (premu0 + psi * sum(D != 0)),
                1 / sqrt(premu0 + psi * sum(D != 0)))
    psi <- rgamma(1,(eta0 + sum(D != 0)) / 2,
                  (eta0 * t20 + sum((D[D != 0] - mu)^2)) / 2)
    M.ps <- M.ps + U %*% D %*% t(V)
    D.ps <- D.ps + -sort(-diag(D))
    # ### output
    # if(ns %% odens==0) {
    # out <- c(ns, sum(D!=0), mean((M0-M.ps/ns)^2), mean((M0-U%*%D%*%t(V))^2),
    #        1/phi, mu, 1/psi)
    # cat(out, "\n")
    # OUT[ns/odens, ] <- out
    #                     }

    ### output
    if(ns %% odens==0) {
      # out <- tibble(ns = ns, dim = sum(D != 0), inv_phi = 1 / phi, mu = mu, inv_psi = 1 / psi)
      out <- c(ns, sum(D!=0), 1/phi, mu, 1/psi)
      cat(out, "\n")
      OUT[ns / odens, ] <- out
    }
  }
  return(OUT)
}
