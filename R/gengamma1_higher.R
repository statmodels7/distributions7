#' @include gengamma1_distrib.R reparam_maps.R
NULL

# The generalized gamma's third and fourth derivatives, written out.
#
# The log-density splits into five terms, and each one is either elementary or
# a composition of a univariate function with a two-variable inner map, so
# the written-out template of fdb2() covers it:
#
#   l = log p - d log a - lgamma(d/p) + (d - 1) log y - exp(p L),  L = log(y/a)
#
# The two compositions do not share a variable pair: -lgamma(d/p) involves
# (d, p) and -exp(pL) involves (a, p), so every component of the three-variable
# derivative is one term of one of them plus the elementary pieces. The
# assembly is checked at order two against the compiled Hessian, which was
# written independently.

#' Derivative Components of the Generalized Gamma
#'
#' @description
#' The components of \eqn{\partial^{\alpha+\beta+\gamma}\ell /
#' \partial a^\alpha \partial d^\beta \partial p^\gamma} at any order from one
#' to four, assembled from the five terms of the log-density.
#'
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of component vectors, keyed as
#'   \code{\link{deriv_names}}.
#'
#' @seealso \code{\link{gengamma1_distrib}}, \code{\link{fdb2}}
#' @keywords internal
gengamma_components <- function(y, theta, order) {
  a <- theta[[1]]
  d <- theta[[2]]
  p <- theta[[3]]
  n <- max(length(y), lengths(theta[1:3]))
  y <- rep_len(y, n)
  a <- rep_len(a, n); d <- rep_len(d, n); p <- rep_len(p, n)
  one <- rep_len(1, n)

  L <- base::log(y) - base::log(a)
  w <- exp(p * L)
  k <- d / p

  # -lgamma(d/p): outer derivatives at k, inner k(d, p) with x = d, z = p
  gam <- fdb2(
    list(-digamma(k), -trigamma(k), -psigamma(k, 2L), -psigamma(k, 3L)),
    list(x = one / p, z = -d / p^2,
         xz = -one / p^2, zz = 2 * d / p^3,
         xzz = 2 * one / p^3, zzz = -6 * d / p^4,
         xzzz = -6 * one / p^4, zzzz = 24 * d / p^5)
  )
  # -exp(v) with v = p L(a): outer derivatives all equal exp(v) = w, inner
  # v(a, p) with x = a, z = p
  ex <- fdb2(
    list(w, w, w, w),
    list(x = -p / a, z = L,
         xx = p / a^2, xz = -one / a,
         xxx = -2 * p / a^3, xxz = one / a^2,
         xxxx = 6 * p / a^4, xxxz = -2 * one / a^3)
  )

  # d^m log(a) / da^m
  dlog_a <- function(m) (-1)^(m - 1L) * factorial(m - 1L) / a^m
  key <- function(xc, zc) {
    paste0(strrep("x", xc), strrep("z", zc))
  }

  comp <- function(al, be, ga) {
    out <- numeric(n)
    # log p
    if (al == 0L && be == 0L && ga >= 1L) {
      out <- out + (-1)^(ga - 1L) * factorial(ga - 1L) / p^ga
    }
    # -d log a
    if (ga == 0L) {
      if (be == 0L && al >= 1L) out <- out - d * dlog_a(al)
      if (be == 1L && al >= 1L) out <- out - dlog_a(al)
      if (be == 1L && al == 0L) out <- out - base::log(a)
    }
    # (d - 1) log y
    if (al == 0L && be == 1L && ga == 0L) out <- out + base::log(y)
    # -lgamma(d/p), which does not involve a
    if (al == 0L && be + ga >= 1L) {
      v <- gam[[key(be, ga)]]
      if (!is.null(v)) out <- out + v
    }
    # -exp(pL), which does not involve d
    if (be == 0L && al + ga >= 1L) {
      v <- ex[[key(al, ga)]]
      if (!is.null(v)) out <- out - v
    }
    out
  }

  nms <- deriv_names(c("a", "d", "p"), order)
  stats::setNames(lapply(nms, function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    comp(sum(parts == "a"), sum(parts == "d"), sum(parts == "p"))
  }), nms)
}


#' @title Generalized Gamma Third and Fourth Derivatives
#' @name distrib_deriv3.GenGamma1Distrib
#' @description
#' Closed form at both orders, from \code{\link{gengamma_components}}: the
#' log-density is elementary apart from \eqn{\lgamma(d/p)} and
#' \eqn{\exp(p\log(y/a))}, and each of those is a univariate function of a
#' two-variable map, so the written-out composition covers every component.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx The approximation used when \code{expected} is \code{TRUE}.
#' @param nsim Monte Carlo draws when \code{approx = "mc"}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_deriv3, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  gengamma_components(y, theta, 3L)
}

#' @rdname distrib_deriv3.GenGamma1Distrib
#' @name distrib_deriv4.GenGamma1Distrib
#' @return A named list of fourth-derivative components.
S7::method(distrib_deriv4, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  gengamma_components(y, theta, 4L)
}
