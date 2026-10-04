#' @include gengamma1_distrib.R reparam_maps.R y_higher.R cross2_derivatives.R cross_theta2_derivatives.R higher_derivatives.R link_scale.R
NULL

# The generalized gamma's derivatives of orders three to five, and those in
# the response, come from the compiled kernels of src/gengamma1.cpp, generated
# by stabilita/gen_gengamma1.py. gengamma_components() below is an independent
# assembly of orders one to four, kept as the reference the tests compare the
# kernels with.
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
#' @seealso [distrib_deriv3.GenGamma1Distrib()] for the compiled kernels
#'   this assembly is checked against;
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

  # Every component of order `order` reads the partials of the two
  # compositions below at that order and at no other, so fdb2() forms that
  # order alone, from the outer derivatives and inner partials up to it.
  # -lgamma(d/p): outer derivatives at k, inner k(d, p) with x = d, z = p
  uk <- list(x = one / p, z = -d / p^2)
  if (order >= 2L) { uk$xz <- -one / p^2; uk$zz <- 2 * d / p^3 }
  if (order >= 3L) { uk$xzz <- 2 * one / p^3; uk$zzz <- -6 * d / p^4 }
  if (order >= 4L) { uk$xzzz <- -6 * one / p^4; uk$zzzz <- 24 * d / p^5 }
  gam <- fdb2(lapply(seq_len(order), function(j) -psigamma(k, j - 1L)),
              uk, order)
  # -exp(v) with v = p L(a): outer derivatives all equal exp(v) = w, inner
  # v(a, p) with x = a, z = p
  uv <- list(x = -p / a, z = L)
  if (order >= 2L) { uv$xx <- p / a^2; uv$xz <- -one / a }
  if (order >= 3L) { uv$xxx <- -2 * p / a^3; uv$xxz <- one / a^2 }
  if (order >= 4L) { uv$xxxx <- 6 * p / a^4; uv$xxxz <- -2 * one / a^3 }
  ex <- fdb2(rep(list(w), order), uv, order)

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


#' @title Generalized Gamma Derivatives of Orders Three to Five
#' @name distrib_deriv3.GenGamma1Distrib
#' @aliases distrib_deriv4.GenGamma1Distrib distrib_deriv5.GenGamma1Distrib
#'
#' @description
#' Return the third, fourth and fifth derivatives of the log-density in
#' \eqn{(a, d, p)}, and for orders three and four their expectations, each
#' from its own compiled kernel.
#'
#' @details
#' With \eqn{k = d/p} and \eqn{U = p\log(y/a)}, the log-density is
#' \deqn{\ell = \log p - \log y + kU - e^{U} - \log\Gamma(k),}
#' and \eqn{e^{U} = (y/a)^p} is a gamma variable with shape \eqn{k} and unit
#' rate. Every component is a closed form derived offline and written out: a
#' polynomial in \eqn{U}, \eqn{e^{U}} and \eqn{a/y} whose coefficients depend
#' on \eqn{(d, p)} alone, multiplied by \eqn{a^{-r}} for \eqn{r} derivatives
#' in \eqn{a}. The coefficients contain the polygamma functions at
#' \eqn{k + 1}, through
#' \eqn{\psi^{(n)}(k) = \psi^{(n)}(k+1) + (-1)^{n+1} n!/k^{n+1}}, so that the
#' poles in \eqn{1/k} which the derivatives combine cancel in the closed form
#' and not in floating point. The expected derivatives are polynomials in the
#' central moments of \eqn{U}, whose cumulants are the polygamma functions at
#' \eqn{k}; `approx` and `nsim` are not read.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`, all strictly positive.
#' @param expected Logical of length 1; for orders three and four, whether the
#'   expectation under the model is returned. Defaults to `FALSE`.
#' @param scale `"parameter"` or `"link"`; the link scale is applied by the
#'   generic for orders three and four and by [deriv5_scale()] for order five.
#' @param approx,nsim Accepted for the generic's signature; the expectations
#'   are closed forms.
#' @param ... Unused.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list with one numeric vector per component, keyed as
#'   [deriv_names()] keys them: 10, 15 and 21 components at orders three, four
#'   and five.
#'
#' @seealso [distrib_hessian.GenGamma1Distrib()] for the order below,
#'   [gengamma_components()] for an independent assembly of orders one to four,
#'   and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' y <- c(0.6, 1.4, 3.1)
#' th <- list(a = 2, d = 1.5, p = 1.3)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces a mixed component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 + eps))$d_p
#' dn <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 - eps))$d_p
#' all.equal((up - dn) / (2 * eps), d3$d_p_p, tolerance = 1e-6)
#'
#' # The expected fourth derivatives, closed form.
#' distrib_deriv4(d, 1, th, expected = TRUE)$p_p_p_p
#'
#' # Order five on the link scale.
#' distrib_deriv5(d, y, th, scale = "link")$p_p_p_p_p
S7::method(distrib_deriv3, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gengamma1_deriv3_expected_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  } else {
    gengamma1_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

S7::method(distrib_deriv4, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gengamma1_deriv4_expected_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  } else {
    gengamma1_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

S7::method(distrib_deriv5, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ..., threads = 1L) {
  deriv5_scale(distrib, y, theta,
               gengamma1_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads),
               match.arg(scale))
}


#' @title Generalized Gamma Derivatives in the Response
#' @name distrib_grad_y.GenGamma1Distrib
#' @aliases distrib_hess_y.GenGamma1Distrib distrib_deriv3_y.GenGamma1Distrib
#'   distrib_deriv4_y.GenGamma1Distrib distrib_cross_y.GenGamma1Distrib
#'   distrib_cross2_y.GenGamma1Distrib distrib_grad_y_hess.GenGamma1Distrib
#'   distrib_hess_y_hess.GenGamma1Distrib
#'
#' @description
#' Return the derivatives of the log-density in the response to order four,
#' and the mixed derivatives of orders one and two in the response and one and
#' two in \eqn{(a, d, p)}, each a closed form from its own compiled kernel,
#' written as [distrib_deriv3.GenGamma1Distrib()] describes.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`, all strictly positive.
#' @param scale `"parameter"` or `"link"`, for the mixed derivatives; the link
#'   scale is applied by the generic.
#' @param ... Unused.
#'
#' @return A numeric vector for the derivatives in the response; a named list,
#'   keyed by parameter or by parameter pair in [hess_names()]'s order, for the
#'   mixed derivatives.
#'
#' @seealso [distrib_gradient.GenGamma1Distrib()] for the derivatives in the
#'   parameters.
#'
#' @examples
#' d <- gengamma1_distrib()
#' th <- list(a = 2, d = 3, p = 1.5)
#' y <- c(0.5, 1.5, 4)
#'
#' # The score in y written out: ((d - 1) - p (y/a)^p) / y.
#' all.equal(distrib_grad_y(d, y, th), (2 - 1.5 * (y / 2)^1.5) / y)
#'
#' distrib_hess_y_hess(d, y, th)$p_p
S7::method(distrib_grad_y, GenGamma1Distrib) <- function(distrib, y, theta, ...) {
  gengamma1_dy1_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}

S7::method(distrib_hess_y, GenGamma1Distrib) <- function(distrib, y, theta, ...) {
  gengamma1_dy2_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}

S7::method(distrib_deriv3_y, GenGamma1Distrib) <- function(distrib, y, theta, ...) {
  gengamma1_dy3_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}

S7::method(distrib_deriv4_y, GenGamma1Distrib) <- function(distrib, y, theta, ...) {
  gengamma1_dy4_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}

S7::method(distrib_cross_y, GenGamma1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma1_cross_y_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_cross2_y, GenGamma1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma1_cross2_y_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_grad_y_hess, GenGamma1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma1_grad_y_hess_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_hess_y_hess, GenGamma1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma1_hess_y_hess_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}
