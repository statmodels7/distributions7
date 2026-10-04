#' @include distrib.R generics.R
NULL

# The explicit map derivatives of the second parametrizations. Every partial
# is a written formula: the univariate chains are Faa di Bruno to fourth
# order spelled out, and the one genuinely bivariate map (the generalized
# gamma's) goes through the written-out two-variable composition template
# below -- the same fixed algebra the PIG kernel carries in C++. No jets
# anywhere.

#' Univariate Composition, One Order
#'
#' @description
#' The derivative of order `order` of \eqn{h(u(x))} for scalar chains, Faa di
#' Bruno written out: \eqn{(h \circ u)'' = h''u_1^2 + h'u_2} and so on to
#' order four. Each order is its own formula and forms nothing above it.
#'
#' @param h A list of the outer derivatives `h1` to at least `h_order` at the
#'   inner value.
#' @param u A list of the inner derivatives `u1` to at least `u_order`.
#' @param order An integer from 1 to 4.
#' @return A numeric vector, the derivative of order `order` of the
#'   composition.
#'
#' @seealso [fdb2()] for a bivariate inner function;
#'   [reparam_map_derivs()], which consumes it.
#' @keywords internal
fdb1 <- function(h, u, order) {
  switch(order,
    h[[1]] * u[[1]],
    h[[2]] * u[[1]]^2 + h[[1]] * u[[2]],
    h[[3]] * u[[1]]^3 + 3 * h[[2]] * u[[1]] * u[[2]] + h[[1]] * u[[3]],
    h[[4]] * u[[1]]^4 + 6 * h[[3]] * u[[1]]^2 * u[[2]] +
      h[[2]] * (3 * u[[2]]^2 + 4 * u[[1]] * u[[3]]) + h[[1]] * u[[4]]
  )
}

#' Bivariate Composition, One Order
#'
#' @description
#' The partials of order `order` of \eqn{h(u(x, z))} for a scalar outer
#' \eqn{h} and a bivariate inner \eqn{u}, Faa di Bruno written out component
#' by component: two partials at order one, three at order two, four at order
#' three and five at order four. The inner partials arrive as a named list
#' with entries `x`, `z`, `xx`, `xz`, `zz`, `xxx`, `xxz`, `xzz`, `zzz`,
#' `xxxx`, `xxxz`, `xxzz`, `xzzz`, `zzzz`; missing entries count as zero, and
#' only the entries of order up to `order` are read. A caller that needs
#' several orders calls it once per order.
#'
#' @param h A list of the outer derivatives `h1` to at least `h_order` at the
#'   inner value.
#' @param u The named list of inner partials.
#' @param order An integer from 1 to 4.
#' @return A named list of the partials of order `order` of the composition,
#'   keyed as `u` is.
#'
#' @seealso [fdb1()] for a univariate inner function;
#'   [gengamma_components()], [gpd_components()] and [reparam_map_derivs()],
#'   which consume it.
#' @keywords internal
fdb2 <- function(h, u, order) {
  g <- function(k) if (is.null(u[[k]])) 0 else u[[k]]
  ux <- g("x"); uz <- g("z")
  h1 <- h[[1]]
  if (order == 1L) return(list(x = h1 * ux, z = h1 * uz))
  uxx <- g("xx"); uxz <- g("xz"); uzz <- g("zz"); h2 <- h[[2]]
  if (order == 2L) {
    return(list(
      xx = h2 * ux^2 + h1 * uxx,
      xz = h2 * ux * uz + h1 * uxz,
      zz = h2 * uz^2 + h1 * uzz
    ))
  }
  uxxx <- g("xxx"); uxxz <- g("xxz"); uxzz <- g("xzz"); uzzz <- g("zzz")
  h3 <- h[[3]]
  if (order == 3L) {
    return(list(
      xxx = h3 * ux^3 + 3 * h2 * uxx * ux + h1 * uxxx,
      xxz = h3 * ux^2 * uz + h2 * (uxx * uz + 2 * uxz * ux) + h1 * uxxz,
      xzz = h3 * ux * uz^2 + h2 * (uzz * ux + 2 * uxz * uz) + h1 * uxzz,
      zzz = h3 * uz^3 + 3 * h2 * uzz * uz + h1 * uzzz
    ))
  }
  uxxxx <- g("xxxx"); uxxxz <- g("xxxz"); uxxzz <- g("xxzz")
  uxzzz <- g("xzzz"); uzzzz <- g("zzzz")
  h4 <- h[[4]]
  list(
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
#' non-zero partials \eqn{\partial^B \theta_i / \partial \psi_B} up to order
#' `order`: a list over parent parameters, each a list keyed by the sorted
#' tuple of \eqn{\psi} positions. A missing key is an exact zero. A partial of
#' higher order than `order` is not formed. Every formula is derived by hand
#' and validated against one numerical pass per order in the tests.
#'
#' @param psi The aligned list of the new parameters.
#' @param order The highest order of partial to form, an integer from 1 to 4.
#' @return A list over parent parameters of keyed partial tables.
#' @name reparam_map_derivs
#' @keywords internal
NULL

#' @rdname reparam_map_derivs
#' @keywords internal
md_lognormal2 <- function(psi, order) {
  m <- psi[[1]]; w <- psi[[2]]
  # t = 1 + w/m^2, s2 = log t; mu = log m - s2/2, sigma2 = s2
  t <- 1 + w / m^2
  u <- list(x = -2 * w / m^3, z = 1 / m^2)
  if (order >= 2L) { u$xx <- 6 * w / m^4; u$xz <- -2 / m^3 }
  if (order >= 3L) { u$xxx <- -24 * w / m^5; u$xxz <- 6 / m^4 }
  if (order >= 4L) { u$xxxx <- 120 * w / m^6; u$xxxz <- -24 / m^5 }
  h <- list(function() 1 / t, function() -1 / t^2, function() 2 / t^3,
            function() -6 / t^4)
  h <- lapply(h[seq_len(order)], function(f) f())
  s2 <- do.call(c, lapply(seq_len(order), function(o) fdb2(h, u, o)))
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
  for (r in seq_len(order)) {
    key <- paste(rep("1", r), collapse = ",")
    out[[1L]][[key]] <- out[[1L]][[key]] +
      (-1)^(r - 1L) * factorial(r - 1L) / m^r
  }
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_weibull3 <- function(psi, order) {
  m <- psi[[1]]; sg <- psi[[2]]
  # mu = m * G(sigma), G = exp(-lgamma(1 + 1/sigma)); sigma passes through
  z <- 1 + 1 / sg
  zd <- lapply(seq_len(order), function(j) (-1)^j * factorial(j) / sg^(j + 1L))
  hz <- lapply(seq_len(order), function(j) psigamma(z, j - 1L))
  wd <- lapply(seq_len(order), function(o) fdb1(hz, zd, o))
  G <- exp(-lgamma(z))
  # G = exp(-w): its derivatives in sigma, one order at a time
  Gd <- list(function() -wd[[1]] * G,
             function() (wd[[1]]^2 - wd[[2]]) * G,
             function() (-wd[[1]]^3 + 3 * wd[[1]] * wd[[2]] - wd[[3]]) * G,
             function() (wd[[1]]^4 - 6 * wd[[1]]^2 * wd[[2]] + 3 * wd[[2]]^2 +
                           4 * wd[[1]] * wd[[3]] - wd[[4]]) * G)
  Gd <- lapply(Gd[seq_len(order)], function(f) f())
  out <- list(list(), list())
  out <- .md_entry(out, 1L, "1", G)
  out <- .md_entry(out, 2L, "2", rep_len(1, length(m)))
  for (k in seq_len(order)) {
    out <- .md_entry(out, 1L, rep("2", k), m * Gd[[k]])
    if (k < order) out <- .md_entry(out, 1L, c("1", rep("2", k)), Gd[[k]])
  }
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_student_t2 <- function(psi, order) {
  sg <- psi[[2]]; nu <- psi[[3]]
  # sigma_scale = sigma * h(nu), h = sqrt(1 - 2/nu); mu and nu pass through
  u <- 1 - 2 / nu
  ud <- list(function() 2 / nu^2, function() -4 / nu^3,
             function() 12 / nu^4, function() -48 / nu^5)
  ud <- lapply(ud[seq_len(order)], function(f) f())
  s <- sqrt(u)
  hs <- list(function() 0.5 / s, function() -0.25 / s^3,
             function() 0.375 / s^5, function() -0.9375 / s^7)
  hs <- lapply(hs[seq_len(order)], function(f) f())
  hd <- lapply(seq_len(order), function(o) fdb1(hs, ud, o))
  one <- rep_len(1, length(sg))
  out <- list(list(), list(), list())
  out <- .md_entry(out, 1L, "1", one)
  out <- .md_entry(out, 3L, "3", one)
  out <- .md_entry(out, 2L, "2", s * one)
  for (k in seq_len(order)) {
    out <- .md_entry(out, 2L, rep("3", k), sg * hd[[k]])
    if (k < order) out <- .md_entry(out, 2L, c("2", rep("3", k)), hd[[k]] * one)
  }
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_gengamma2 <- function(psi, order) {
  m <- psi[[1]]; d <- psi[[2]]; p <- psi[[3]]
  # a = m * H(d, p), H = exp(W), W = lgamma(d/p) - lgamma((d+1)/p);
  # d and p pass through. x stands for d and z for p in fdb2.
  q_part <- function(dd) {
    u <- list(x = 1 / p, z = -dd / p^2)
    if (order >= 2L) { u$xz <- -1 / p^2; u$zz <- 2 * dd / p^3 }
    if (order >= 3L) { u$xzz <- 2 / p^3; u$zzz <- -6 * dd / p^4 }
    if (order >= 4L) { u$xzzz <- -6 / p^4; u$zzzz <- 24 * dd / p^5 }
    list(value = dd / p, u = u)
  }
  q1 <- q_part(d); q2 <- q_part(d + 1)
  lg <- function(q) lapply(seq_len(order), function(j) psigamma(q, j - 1L))
  l1 <- lg(q1$value); l2 <- lg(q2$value)
  W <- do.call(c, lapply(seq_len(order), function(o) {
    Map(`-`, fdb2(l1, q1$u, o), fdb2(l2, q2$u, o))
  }))
  # H = exp(W): compose exp with the W partials through fdb2, whose inner
  # is W itself (value drops out of the derivative formulas)
  H0 <- exp(lgamma(q1$value) - lgamma(q2$value))
  eh <- rep(list(H0), order)
  H <- do.call(c, lapply(seq_len(order), function(o) fdb2(eh, W, o)))
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
    if (length(tup) < order) out <- .md_entry(out, 1L, c("1", tup), H[[k]])
  }
  out
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_invgauss2 <- function(psi, order) {
  lam <- psi[[2]]
  out <- list(list(), list())
  out <- .md_entry(out, 1L, "1", rep_len(1, length(lam)))
  # phi = 1/lambda, whose j-th derivative is (-1)^j j! lambda^(-1-j)
  for (j in seq_len(order)) {
    out <- .md_entry(out, 2L, rep("2", j), (-1)^j * factorial(j) / lam^(j + 1L))
  }
  out
}


#' @rdname reparam_map_derivs
#' @keywords internal
md_skewnormal2 <- function(psi, order) {
  sg <- psi[[2]]; g <- psi[[3]]
  b <- sn_b()
  # r = (2 gamma1 / (4 - pi))^(1/3), the real cube root; the whole map is
  # xi = mu - sigma r, omega = sigma sqrt(1 + r^2),
  # alpha = r / sqrt(b^2 + (b^2 - 1) r^2). The power rule written through
  # r/gamma keeps every derivative real on both signs of the skewness. Only
  # the partials up to `order` are formed.
  r <- sign(g) * (2 * abs(g) / (4 - pi))^(1 / 3)
  rd <- list(function() r / (3 * g), function() -2 * r / (9 * g^2),
             function() 10 * r / (27 * g^3), function() -80 * r / (81 * g^4))
  rd <- lapply(rd[seq_len(order)], function(f) f())
  u <- 1 + r^2
  ud <- list(function() 2 * r * rd[[1]],
             function() 2 * (rd[[1]]^2 + r * rd[[2]]),
             function() 2 * (3 * rd[[1]] * rd[[2]] + r * rd[[3]]),
             function() 2 * (3 * rd[[2]]^2 + 4 * rd[[1]] * rd[[3]] + r * rd[[4]]))
  ud <- lapply(ud[seq_len(order)], function(f) f())
  sq <- sqrt(u)
  hq <- list(function() 0.5 / sq, function() -0.25 / sq^3,
             function() 0.375 / sq^5, function() -0.9375 / sq^7)
  hq <- lapply(hq[seq_len(order)], function(f) f())
  qd <- lapply(seq_len(order), function(o) fdb1(hq, ud, o))

  cb <- b^2 - 1
  # b^2 + (b^2-1) r^2 is b^2 (1 - delta^2) exactly, and the second form
  # is the one that survives the top of the range: the first reaches
  # -1.11e-16 at the largest skewness the link can produce, where the
  # quantity is 9.42e-17.
  D <- b^2 * sn_one_minus_delta2(g, sign(g))
  ad_r <- list(function() b^2 * D^-1.5,
               function() -3 * b^2 * cb * r * D^-2.5,
               function() -3 * b^2 * cb * (D - 5 * cb * r^2) * D^-3.5,
               function() 15 * b^2 * cb^2 * r * (3 * D - 7 * cb * r^2) * D^-4.5)
  ad_r <- lapply(ad_r[seq_len(order)], function(f) f())
  al <- lapply(seq_len(order), function(o) fdb1(ad_r, rd, o))

  one <- rep_len(1, length(g))
  out <- list(list(), list(), list())
  put <- function(out, i, tup, val) {
    if (length(tup) <= order) .md_entry(out, i, tup, val) else out
  }
  # xi = mu - sigma r
  out <- put(out, 1L, "1", one)
  out <- put(out, 1L, "2", -r)
  out <- put(out, 1L, "3", -sg * rd[[1]])
  for (k in seq_len(order)[-1]) {
    out <- put(out, 1L, c("2", rep("3", k - 1)), -rd[[k - 1]])
    out <- put(out, 1L, rep("3", k), -sg * rd[[k]])
  }
  if (order >= 4L) out <- put(out, 1L, c("2", "3", "3", "3"), -rd[[3]])
  # omega = sigma sqrt(1 + r^2)
  out <- put(out, 2L, "2", sq * one)
  for (k in seq_len(order)) {
    out <- put(out, 2L, rep("3", k), sg * qd[[k]])
    out <- put(out, 2L, c("2", rep("3", k)), qd[[k]] * one)
  }
  # alpha, a function of gamma1 alone
  for (k in seq_len(order)) out <- put(out, 3L, rep("3", k), al[[k]])
  out
}


#' @rdname reparam_map_derivs
#' @keywords internal
md_laplace2 <- function(psi, order) {
  lam <- psi[[2]]
  # mu passes through; sigma = 1/lambda, whose j-th derivative is
  # (-1)^j j! lambda^(-1-j)
  s <- lapply(seq_len(order), function(j) (-1)^j * factorial(j) / lam^(j + 1L))
  names(s) <- vapply(seq_len(order), function(j)
    paste(rep("2", j), collapse = ","), "")
  list(list("1" = rep_len(1, length(lam))), s)
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_gaussian2 <- function(psi, order) {
  v <- psi[[2]]
  # mu passes through; sigma = sqrt(sigma2), whose j-th derivative is
  # c_j v^(1/2 - j) with c = 1/2, -1/4, 3/8, -15/16
  cf <- c(0.5, -0.25, 0.375, -0.9375)
  s <- lapply(seq_len(order), function(j) cf[j] * v^(0.5 - j))
  names(s) <- vapply(seq_len(order), function(j)
    paste(rep("2", j), collapse = ","), "")
  list(list("1" = rep_len(1, length(v))), s)
}

#' @rdname reparam_map_derivs
#' @keywords internal
md_gaussian3 <- function(psi, order) {
  tau <- psi[[2]]
  # mu passes through; sigma = tau^(-1/2), whose j-th derivative is
  # c_j tau^(-1/2 - j) with c = -1/2, 3/4, -15/8, 105/16
  cf <- c(-0.5, 0.75, -1.875, 6.5625)
  s <- lapply(seq_len(order), function(j) cf[j] * tau^(-0.5 - j))
  names(s) <- vapply(seq_len(order), function(j)
    paste(rep("2", j), collapse = ","), "")
  list(list("1" = rep_len(1, length(tau))), s)
}
