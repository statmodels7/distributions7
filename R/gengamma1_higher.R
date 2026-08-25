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
#' Returns the components of
#' \eqn{\partial^{\alpha+\beta+\gamma}\ell / \partial a^\alpha \partial d^\beta
#' \partial p^\gamma} at any order from one to four, assembled term by term
#' from the five pieces of the log-density.
#'
#' @details
#' # Why the assembly is short
#'
#' With \eqn{L = \log(y/a)} the log-density splits into
#' \deqn{\ell = \log p - d\log a - \log\Gamma(d/p) + (d-1)\log y - e^{pL},}
#' and each piece is either elementary or a univariate function composed with a
#' **two-variable** inner map, which the written-out template of [fdb2()]
#' covers.
#'
#' What keeps the sum from growing is that the two compositions do not share a
#' variable pair: \eqn{-\log\Gamma(d/p)} involves \eqn{(d, p)} and
#' \eqn{-e^{pL}} involves \eqn{(a, p)}. Every component of the three-variable
#' derivative is therefore one term of one composition plus the elementary
#' pieces, and no genuinely three-variable expansion is ever formed. A
#' component naming both \eqn{a} and \eqn{d} comes from the elementary
#' \eqn{-d\log a} alone.
#'
#' The assembly is checked at order two against the compiled Hessian, which was
#' written independently, so the orders that cannot be checked against a
#' hand-written form rest on the orders that can.
#'
#' @param y A numeric vector of positive observations.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`, all strictly positive. Shorter
#'   components are recycled to the common length.
#' @param order The derivative order, an integer from 1 to 4.
#'
#' @return A named list of component vectors, one per distinct multi-index of
#'   the given order and keyed as [deriv_names()] keys them: three at order 1,
#'   six at order 2, ten at order 3 and fifteen at order 4. Each has the
#'   recycled length of the inputs.
#'
#' @seealso [distrib_deriv3.GenGamma1Distrib()] and
#'   [distrib_deriv4.GenGamma1Distrib()], which call this;
#'   [fdb2()] for the two-variable composition template; and
#'   [gengamma1_distrib()] for the family.
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


#' @title Generalized Gamma Third-Order Derivatives
#' @name distrib_deriv3.GenGamma1Distrib
#' @description
#' Computes the ten distinct third derivatives of the generalized gamma
#' log-density in \eqn{a}, \eqn{d} and \eqn{p}, **in closed form**, through
#' [gengamma_components()]. The log-density is elementary apart from
#' \eqn{\log\Gamma(d/p)} and \eqn{\exp\{p\log(y/a)\}}, and each of those is a
#' univariate function of a two-variable map, so the written-out composition
#' [fdb2()] covers every component without forming a three-variable expansion.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead: the
#' expected third derivatives have no closed form. That is the one place on
#' this page where `approx` and `nsim` are read.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of positive observations. With `expected = TRUE`
#'   only its length is read.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`, all strictly positive. A
#'   component of length 1 is recycled.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`, the strategy [expected_derivative()] uses. Read only when
#'   `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of ten numeric vectors, `a_a_a` through `p_p_p`, each
#'   of length `max(length(y), lengths(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{a > 0} the scale,
#' \eqn{d > 0} and \eqn{p > 0} the two shapes, and \eqn{\Gamma} the gamma
#' function.
#'
#' @seealso [distrib_hessian.GenGamma1Distrib()] for the order below,
#'   [distrib_deriv4.GenGamma1Distrib()] for the order above,
#'   [gengamma_components()] for the assembly, and [distrib_deriv3()] for the
#'   generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' y <- c(0.6, 1.4, 3.1)
#' th <- list(a = 2, d = 1.5, p = 1.3)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-scale component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(a = 2 + eps, d = 1.5, p = 1.3))$a_a
#' dn <- distrib_hessian(d, y, list(a = 2 - eps, d = 1.5, p = 1.3))$a_a
#' all.equal((up - dn) / (2 * eps), d3$a_a_a, tolerance = 1e-6)
#'
#' # And the fully mixed component, which the two compositions never both
#' # contribute to.
#' up <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 + eps))$a_d
#' dn <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 - eps))$a_d
#' all.equal((up - dn) / (2 * eps), d3$a_d_p, tolerance = 1e-6)
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

#' @title Generalized Gamma Fourth-Order Derivatives
#' @name distrib_deriv4.GenGamma1Distrib
#' @description
#' Computes the fifteen distinct fourth derivatives of the generalized gamma
#' log-density in \eqn{a}, \eqn{d} and \eqn{p}, **in closed form**, by the
#' construction [distrib_deriv3.GenGamma1Distrib()] describes carried one order
#' further: the two two-variable compositions of [fdb2()] plus the elementary
#' pieces of the log-density.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead: the
#' expected fourth derivatives have no closed form. That is the one place on
#' this page where `approx` and `nsim` are read.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of positive observations. With `expected = TRUE`
#'   only its length is read.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`, all strictly positive. A
#'   component of length 1 is recycled.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of fifteen numeric vectors named for the multi-index
#'   they carry, from `a_a_a_a` to `p_p_p_p`, each of length
#'   `max(length(y), lengths(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{a > 0} the scale,
#' \eqn{d > 0} and \eqn{p > 0} the two shapes.
#'
#' @seealso [distrib_deriv3.GenGamma1Distrib()] for the order below and the
#'   construction, [gengamma_components()] for the assembly, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' y <- c(0.6, 1.4, 3.1)
#' th <- list(a = 2, d = 1.5, p = 1.3)
#' d4 <- distrib_deriv4(d, y, th)
#' length(d4)
#' names(d4)[1:4]
#'
#' # A central difference of the third order reproduces the pure-scale
#' # component.
#' eps <- 1e-4
#' up <- distrib_deriv3(d, y, list(a = 2 + eps, d = 1.5, p = 1.3))$a_a_a
#' dn <- distrib_deriv3(d, y, list(a = 2 - eps, d = 1.5, p = 1.3))$a_a_a
#' all.equal((up - dn) / (2 * eps), d4$a_a_a_a, tolerance = 1e-5)
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
