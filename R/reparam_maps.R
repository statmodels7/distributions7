#' @include distrib.R generics.R
NULL

# The explicit map derivatives of the second parametrizations. Every partial
# is a written formula: the univariate chains are Faa di Bruno to fourth
# order spelled out, and the one genuinely bivariate map (the generalized
# gamma's) goes through the written-out two-variable composition template
# below -- the same fixed algebra the PIG kernel carries in C++. No jets
# anywhere.

#' Univariate Composition to Fourth Order
#'
#' @description
#' The derivatives of \eqn{h(u(x))} for scalar chains, Faa di Bruno written
#' out: \eqn{(h \circ u)'' = h''u_1^2 + h'u_2} and so on to order four.
#'
#' @param h A list of the outer derivatives \code{h1} to \code{h4} at the
#'   inner value.
#' @param u A list of the inner derivatives \code{u1} to \code{u4}.
#' @return A list with the four derivatives of the composition.
#' @keywords internal
fdb1 <- function(h, u) {
  list(
    d1 = h[[1]] * u[[1]],
    d2 = h[[2]] * u[[1]]^2 + h[[1]] * u[[2]],
    d3 = h[[3]] * u[[1]]^3 + 3 * h[[2]] * u[[1]] * u[[2]] + h[[1]] * u[[3]],
    d4 = h[[4]] * u[[1]]^4 + 6 * h[[3]] * u[[1]]^2 * u[[2]] +
      h[[2]] * (3 * u[[2]]^2 + 4 * u[[1]] * u[[3]]) + h[[1]] * u[[4]]
  )
}

#' Bivariate Composition to Fourth Order
#'
#' @description
#' The fourteen partials of \eqn{h(u(x, z))} for a scalar outer \eqn{h} and
#' a bivariate inner \eqn{u}, Faa di Bruno written out component by
#' component. The inner partials arrive as a named list with entries
#' \code{x}, \code{z}, \code{xx}, \code{xz}, \code{zz}, \code{xxx},
#' \code{xxz}, \code{xzz}, \code{zzz}, \code{xxxx}, \code{xxxz},
#' \code{xxzz}, \code{xzzz}, \code{zzzz}; missing entries count as zero.
#'
#' @param h A list of the outer derivatives \code{h1} to \code{h4} at the
#'   inner value.
#' @param u The named list of inner partials.
#' @return A named list of the fourteen partials of the composition, keyed
#'   like the input.
#' @keywords internal
fdb2 <- function(h, u) {
  g <- function(k) if (is.null(u[[k]])) 0 else u[[k]]
  ux <- g("x"); uz <- g("z")
  uxx <- g("xx"); uxz <- g("xz"); uzz <- g("zz")
  uxxx <- g("xxx"); uxxz <- g("xxz"); uxzz <- g("xzz"); uzzz <- g("zzz")
  uxxxx <- g("xxxx"); uxxxz <- g("xxxz"); uxxzz <- g("xxzz")
  uxzzz <- g("xzzz"); uzzzz <- g("zzzz")
  h1 <- h[[1]]; h2 <- h[[2]]; h3 <- h[[3]]; h4 <- h[[4]]
  list(
    x = h1 * ux,
    z = h1 * uz,
    xx = h2 * ux^2 + h1 * uxx,
    xz = h2 * ux * uz + h1 * uxz,
    zz = h2 * uz^2 + h1 * uzz,
    xxx = h3 * ux^3 + 3 * h2 * uxx * ux + h1 * uxxx,
    xxz = h3 * ux^2 * uz + h2 * (uxx * uz + 2 * uxz * ux) + h1 * uxxz,
    xzz = h3 * ux * uz^2 + h2 * (uzz * ux + 2 * uxz * uz) + h1 * uxzz,
    zzz = h3 * uz^3 + 3 * h2 * uzz * uz + h1 * uzzz,
    xxxx = h4 * ux^4 + 6 * h3 * uxx * ux^2 +
      h2 * (3 * uxx^2 + 4 * uxxx * ux) + h1 * uxxxx,
    xxxz = h4 * ux^3 * uz + h3 * (3 * uxz * ux^2 + 3 * uxx * ux * uz) +
      h2 * (3 * uxx * uxz + 3 * uxxz * ux + uxxx * uz) + h1 * uxxxz,
    xxzz = h4 * ux^2 * uz^2 +
      h3 * (uxx * uz^2 + uzz * ux^2 + 4 * uxz * ux * uz) +
      h2 * (uxx * uzz + 2 * uxz^2 + 2 * uxzz * ux + 2 * uxxz * uz) +
      h1 * uxxzz,
    xzzz = h4 * ux * uz^3 + h3 * (3 * uxz * uz^2 + 3 * uzz * ux * uz) +
      h2 * (3 * uzz * uxz + 3 * uxzz * uz + uzzz * ux) + h1 * uxzzz,
    zzzz = h4 * uz^4 + 6 * h3 * uzz * uz^2 +
      h2 * (3 * uzz^2 + 4 * uzzz * uz) + h1 * uzzzz
  )
}

# a partial table entry: parent index i, psi-position tuple, value vector
.md_entry <- function(store, i, tup, val) {
  key <- paste(sort(tup), collapse = ",")
  store[[i]][[key]] <- val
  store
}

#' Explicit Map Derivatives of the Second Parametrizations
#'
#' @description
#' Each function returns, for its family's map \eqn{\theta = h(\psi)}, the
#' non-zero partials \eqn{\partial^B \theta_i / \partial \psi_B} to fourth
#' order: a list over parent parameters, each a list keyed by the sorted
#' tuple of \eqn{\psi} positions. A missing key is an exact zero. Every
#' formula is derived by hand and validated against one numerical pass per
#' order in the tests.
#'
#' @param psi The aligned list of the new parameters.
#' @return A list over parent parameters of keyed partial tables.
#' @name reparam_map_derivs
#' @keywords internal
NULL

#' @rdname reparam_map_derivs
#' @keywords internal
md_lognormal2 <- function(psi) {
  m <- psi[[1]]; w <- psi[[2]]
  # t = 1 + w/m^2, s2 = log t; mu = log m - s2/2, sigma2 = s2
  t <- 1 + w / m^2
  u <- list(
    x = -2 * w / m^3, z = 1 / m^2,
    xx = 6 * w / m^4, xz = -2 / m^3,
    xxx = -24 * w / m^5, xxz = 6 / m^4,
    xxxx = 120 * w / m^6, xxxz = -24 / m^5
  )
  h <- list(1 / t, -1 / t^2, 2 / t^3, -6 / t^4)
  s2 <- fdb2(h, u)
  out <- list(list(), list())
  # sigma2 = s2
  for (k in names(s2)) {
    tup <- chartr("xz", "12", strsplit(k, "")[[1]])
    out <- .md_entry(out, 2L, tup, s2[[k]])
  }
  # mu = log m - s2/2: pure-m log derivatives plus -s2/2 everywhere
  for (k in names(s2)) {
    tup <- chartr("xz", "12", strsplit(k, "")[[1]])
    out <- .md_entry(out, 1L, tup, -s2[[k]] / 2)
  }
  lm <- list(1 / m, -1 / m^2, 2 / m^3, -6 / m^4)
  for (r in 1:4) {
    key <- paste(rep("1", r), collapse = ",")
    out[[1L]][[key]] <- out[[1L]][[key]] + lm[[r]]
  }
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_weibull3 <- function(psi) {
  m <- psi[[1]]; sg <- psi[[2]]
  # mu = m * G(sigma), G = exp(-lgamma(1 + 1/sigma)); sigma passes through
  z <- 1 + 1 / sg
  zd <- list(-1 / sg^2, 2 / sg^3, -6 / sg^4, 24 / sg^5)
  wd <- fdb1(list(digamma(z), trigamma(z), psigamma(z, 2), psigamma(z, 3)), zd)
  G <- exp(-lgamma(z))
  w1 <- wd$d1; w2 <- wd$d2; w3 <- wd$d3; w4 <- wd$d4
  G1 <- -w1 * G
  G2 <- (w1^2 - w2) * G
  G3 <- (-w1^3 + 3 * w1 * w2 - w3) * G
  G4 <- (w1^4 - 6 * w1^2 * w2 + 3 * w2^2 + 4 * w1 * w3 - w4) * G
  out <- list(list(), list())
  out <- .md_entry(out, 1L, "1", G)
  out <- .md_entry(out, 1L, "2", m * G1)
  out <- .md_entry(out, 1L, c("1", "2"), G1)
  out <- .md_entry(out, 1L, c("2", "2"), m * G2)
  out <- .md_entry(out, 1L, c("1", "2", "2"), G2)
  out <- .md_entry(out, 1L, c("2", "2", "2"), m * G3)
  out <- .md_entry(out, 1L, c("1", "2", "2", "2"), G3)
  out <- .md_entry(out, 1L, c("2", "2", "2", "2"), m * G4)
  out <- .md_entry(out, 2L, "2", rep_len(1, length(m)))
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_student_t2 <- function(psi) {
  sg <- psi[[2]]; nu <- psi[[3]]
  # sigma_scale = sigma * h(nu), h = sqrt(1 - 2/nu); mu and nu pass through
  u <- 1 - 2 / nu
  ud <- list(2 / nu^2, -4 / nu^3, 12 / nu^4, -48 / nu^5)
  s <- sqrt(u)
  hd <- fdb1(list(0.5 / s, -0.25 / s^3, 0.375 / s^5, -0.9375 / s^7), ud)
  one <- rep_len(1, length(sg))
  out <- list(list(), list(), list())
  out <- .md_entry(out, 1L, "1", one)
  out <- .md_entry(out, 3L, "3", one)
  out <- .md_entry(out, 2L, "2", s * one)
  out <- .md_entry(out, 2L, "3", sg * hd$d1)
  out <- .md_entry(out, 2L, c("2", "3"), hd$d1 * one)
  out <- .md_entry(out, 2L, c("3", "3"), sg * hd$d2)
  out <- .md_entry(out, 2L, c("2", "3", "3"), hd$d2 * one)
  out <- .md_entry(out, 2L, c("3", "3", "3"), sg * hd$d3)
  out <- .md_entry(out, 2L, c("2", "3", "3", "3"), hd$d3 * one)
  out <- .md_entry(out, 2L, c("3", "3", "3", "3"), sg * hd$d4)
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_gengamma2 <- function(psi) {
  m <- psi[[1]]; d <- psi[[2]]; p <- psi[[3]]
  # a = m * H(d, p), H = exp(W), W = lgamma(d/p) - lgamma((d+1)/p);
  # d and p pass through. x stands for d and z for p in fdb2.
  q_part <- function(dd) {
    q <- dd / p
    list(value = q, u = list(
      x = 1 / p, z = -dd / p^2, xz = -1 / p^2, zz = 2 * dd / p^3,
      xzz = 2 / p^3, zzz = -6 * dd / p^4, xzzz = -6 / p^4,
      zzzz = 24 * dd / p^5
    ))
  }
  q1 <- q_part(d); q2 <- q_part(d + 1)
  lg <- function(q) list(digamma(q), trigamma(q), psigamma(q, 2), psigamma(q, 3))
  W1 <- fdb2(lg(q1$value), q1$u)
  W2 <- fdb2(lg(q2$value), q2$u)
  W <- Map(`-`, W1, W2)
  # H = exp(W): compose exp with the W partials through fdb2, whose inner
  # is W itself (value drops out of the derivative formulas)
  H0 <- exp(lgamma(q1$value) - lgamma(q2$value))
  H <- fdb2(list(H0, H0, H0, H0), W)
  out <- list(list(), list(), list())
  one <- rep_len(1, length(m))
  out <- .md_entry(out, 2L, "2", one)
  out <- .md_entry(out, 3L, "3", one)
  # a = m * H(d, p): the Leibniz split by whether position 1 appears
  key2 <- function(k) chartr("xz", "23", strsplit(k, "")[[1]])
  out <- .md_entry(out, 1L, "1", H0 * one)
  for (k in names(H)) {
    tup <- key2(k)
    out <- .md_entry(out, 1L, tup, m * H[[k]])
    if (length(tup) < 4L) out <- .md_entry(out, 1L, c("1", tup), H[[k]])
  }
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_invgauss2 <- function(psi) {
  lam <- psi[[2]]
  one <- rep_len(1, length(lam))
  out <- list(list(), list())
  out <- .md_entry(out, 1L, "1", one)
  out <- .md_entry(out, 2L, "2", -1 / lam^2)
  out <- .md_entry(out, 2L, c("2", "2"), 2 / lam^3)
  out <- .md_entry(out, 2L, c("2", "2", "2"), -6 / lam^4)
  out <- .md_entry(out, 2L, c("2", "2", "2", "2"), 24 / lam^5)
  out
}


#' @rdname reparam_map_derivs
#' @keywords internal
md_skewnormal2 <- function(psi) {
  sg <- psi[[2]]; g <- psi[[3]]
  b <- sn_b()
  # r = (2 gamma1 / (4 - pi))^(1/3), the real cube root; the whole map is
  # xi = mu - sigma r, omega = sigma sqrt(1 + r^2),
  # alpha = r / sqrt(b^2 + (b^2 - 1) r^2). The power rule written through
  # r/gamma keeps every derivative real on both signs of the skewness.
  r <- sign(g) * (2 * abs(g) / (4 - pi))^(1 / 3)
  r1 <- r / (3 * g)
  r2 <- -2 * r / (9 * g^2)
  r3 <- 10 * r / (27 * g^3)
  r4 <- -80 * r / (81 * g^4)

  u <- 1 + r^2
  ud <- list(2 * r * r1, 2 * (r1^2 + r * r2), 2 * (3 * r1 * r2 + r * r3),
             2 * (3 * r2^2 + 4 * r1 * r3 + r * r4))
  sq <- sqrt(u)
  qd <- fdb1(list(0.5 / sq, -0.25 / sq^3, 0.375 / sq^5, -0.9375 / sq^7), ud)

  cb <- b^2 - 1
  D <- b^2 + cb * r^2
  ad_r <- list(b^2 * D^-1.5,
               -3 * b^2 * cb * r * D^-2.5,
               -3 * b^2 * cb * (D - 5 * cb * r^2) * D^-3.5,
               15 * b^2 * cb^2 * r * (3 * D - 7 * cb * r^2) * D^-4.5)
  al <- fdb1(ad_r, list(r1, r2, r3, r4))

  one <- rep_len(1, length(g))
  out <- list(list(), list(), list())
  # xi = mu - sigma r
  out <- .md_entry(out, 1L, "1", one)
  out <- .md_entry(out, 1L, "2", -r)
  out <- .md_entry(out, 1L, "3", -sg * r1)
  out <- .md_entry(out, 1L, c("2", "3"), -r1)
  out <- .md_entry(out, 1L, c("3", "3"), -sg * r2)
  out <- .md_entry(out, 1L, c("2", "3", "3"), -r2)
  out <- .md_entry(out, 1L, c("3", "3", "3"), -sg * r3)
  out <- .md_entry(out, 1L, c("2", "3", "3", "3"), -r3)
  out <- .md_entry(out, 1L, c("3", "3", "3", "3"), -sg * r4)
  # omega = sigma sqrt(1 + r^2)
  out <- .md_entry(out, 2L, "2", sq * one)
  out <- .md_entry(out, 2L, "3", sg * qd$d1)
  out <- .md_entry(out, 2L, c("2", "3"), qd$d1 * one)
  out <- .md_entry(out, 2L, c("3", "3"), sg * qd$d2)
  out <- .md_entry(out, 2L, c("2", "3", "3"), qd$d2 * one)
  out <- .md_entry(out, 2L, c("3", "3", "3"), sg * qd$d3)
  out <- .md_entry(out, 2L, c("2", "3", "3", "3"), qd$d3 * one)
  out <- .md_entry(out, 2L, c("3", "3", "3", "3"), sg * qd$d4)
  # alpha, a function of gamma1 alone
  out <- .md_entry(out, 3L, "3", al$d1)
  out <- .md_entry(out, 3L, c("3", "3"), al$d2)
  out <- .md_entry(out, 3L, c("3", "3", "3"), al$d3)
  out <- .md_entry(out, 3L, c("3", "3", "3", "3"), al$d4)
  out
}


#' @rdname reparam_map_derivs
#' @keywords internal
md_gaussian2 <- function(psi) {
  v <- psi[[2]]
  # mu passes through; sigma = sqrt(sigma2)
  list(
    list("1" = rep_len(1, length(v))),
    list("2" = 0.5 / sqrt(v), "2,2" = -0.25 / v^1.5,
         "2,2,2" = 0.375 / v^2.5, "2,2,2,2" = -0.9375 / v^3.5)
  )
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_gaussian3 <- function(psi) {
  tau <- psi[[2]]
  # mu passes through; sigma = tau^(-1/2)
  list(
    list("1" = rep_len(1, length(tau))),
    list("2" = -0.5 / tau^1.5, "2,2" = 0.75 / tau^2.5,
         "2,2,2" = -1.875 / tau^3.5, "2,2,2,2" = 6.5625 / tau^4.5)
  )
}
