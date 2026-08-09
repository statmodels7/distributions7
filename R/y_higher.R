#' @include y_derivatives.R generics.R
#' @include gaussian1_distrib.R gaussian2_distrib.R gaussian3_distrib.R
#' @include cauchy_distrib.R logistic_distrib.R laplace_distrib.R
#' @include laplace2_distrib.R enet_distrib.R pseudohuber_distrib.R
#' @include student_t1_distrib.R
#' @include skewnormal1_distrib.R skewnormal2_distrib.R skewt_distrib.R
#' @include gumbel_distrib.R
NULL

# ===========================================================================
# Third and fourth derivatives of the log-density with respect to the
# response.
#
# For a family in which the response enters only as y - mu there is nothing
# to derive: differentiating in y and in the location are the same operation
# up to a sign,
#
#   d^k l / dy^k = (-1)^k  d^k l / dmu^k,
#
# so every location family inherits these orders from the parameter
# derivatives it already has. What is left is the families on a half line,
# where the two are unrelated, and those take one stencil of the requested
# order applied to the log-density -- never a difference of a difference.
# ===========================================================================

#' Numerical Third and Fourth Response Derivatives
#'
#' @description
#' One central stencil of the requested order applied to
#' \code{distrib_pdf(..., log = TRUE)}.
#'
#' @details
#' The stencil reaches two steps either side, so the step is clamped to half
#' of what \code{\link{fd_steps_y}} allows: a Gamma observation near zero
#' would otherwise be differentiated at a point outside the support. The
#' relative step is \eqn{\varepsilon^{1/(k+2)}}, which balances the
#' \eqn{h^{2}} truncation against the \eqn{\varepsilon/h^{k}} rounding.
#'
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 3 or 4.
#' @param h_rel Relative finite-difference step.
#'
#' @return A numeric vector the length of \code{y}.
#'
#' @seealso \code{\link{numerical_hess_y}}
#'
#' @examples
#' numerical_deriv_y(gaussian1_distrib(), c(-1, 0, 1),
#'                   list(mu = 0, sigma = 1), order = 3)
#'
#' @export
numerical_deriv_y <- function(distrib, y, theta, order,
                              h_rel = .Machine$double.eps^(1 / (order + 2))) {
  h <- fd_steps_y(y, distrib@bounds, h_rel) / 2
  ld <- function(s) distrib_pdf(distrib, y + s * h, theta, log = TRUE)
  if (order == 3L) {
    (-0.5 * ld(-2) + ld(-1) - ld(1) + 0.5 * ld(2)) / h^3
  } else {
    (ld(-2) - 4 * ld(-1) + 6 * ld(0) - 4 * ld(1) + ld(2)) / h^4
  }
}

#' Third and Fourth Derivatives With Respect to the Response
#'
#' @description
#' \eqn{\partial^{3}\ell/\partial y^{3}} and
#' \eqn{\partial^{4}\ell/\partial y^{4}}, completing the sequence begun by
#' \code{\link{distrib_grad_y}} and \code{\link{distrib_hess_y}}.
#'
#' @details
#' A family whose response enters only as \eqn{y - \mu} gets these from the
#' derivatives it already has in the location, since
#' \eqn{\partial^{k}\ell/\partial y^{k} = (-1)^{k}\,
#' \partial^{k}\ell/\partial\mu^{k}}; the others take one stencil of the
#' requested order on the log-density.
#'
#' As with the orders below, a discrete family has no such derivative and the
#' generic rejects rather than returning a difference across the lattice.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param ... Passed to methods.
#'
#' @return A numeric vector the length of \code{y}.
#'
#' @seealso \code{\link{distrib_hess_y}}, \code{\link{numerical_deriv_y}}
#'
#' @examples
#' distrib_deriv3_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
distrib_deriv3_y <- S7::new_generic(
  "distrib_deriv3_y", "distrib",
  function(distrib, y, theta, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @rdname distrib_deriv3_y
#'
#' @examples
#' distrib_deriv4_y(logistic_distrib(), 0.5, list(mu = 0, sigma = 1))
#'
#' @export
distrib_deriv4_y <- S7::new_generic(
  "distrib_deriv4_y", "distrib",
  function(distrib, y, theta, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @title Default Third and Fourth Response Derivatives
#' @name distrib_deriv3_y.continuous_distrib
#' @description One stencil of the order asked for, on the log-density.
#' @param distrib A \code{continuous_distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_deriv3_y, continuous_distrib) <- function(distrib, y, theta,
                                                             ...) {
  numerical_deriv_y(distrib, y, theta, 3L)
}

#' @rdname distrib_deriv3_y.continuous_distrib
#' @name distrib_deriv4_y.continuous_distrib
S7::method(distrib_deriv4_y, continuous_distrib) <- function(distrib, y, theta,
                                                             ...) {
  numerical_deriv_y(distrib, y, theta, 4L)
}

#' The Response Derivative of a Location Family
#'
#' @description
#' Builds the order-\code{k} response derivative of a family whose response
#' enters only as \eqn{y - \mu}, as \eqn{(-1)^{k}} times the pure derivative
#' in the location.
#'
#' @details
#' The identity is exact and needs no formula of its own, which is the point:
#' the location derivative is already written, often as a compiled kernel, and
#' the response derivative is the same number with a sign. It also fixes the
#' two orders below, where \eqn{\partial\ell/\partial y = -\partial\ell/
#' \partial\mu} and \eqn{\partial^{2}\ell/\partial y^{2} =
#' \partial^{2}\ell/\partial\mu^{2}}, and the tests check the new orders
#' against exactly that.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @keywords internal
loc_deriv_y_k <- function(order) {
  key <- paste(rep("mu", order), collapse = "_")
  sgn <- (-1)^order
  function(distrib, y, theta, ...) {
    d <- if (order == 3L) distrib_deriv3(distrib, y, theta) else
                          distrib_deriv4(distrib, y, theta)
    rep_len(sgn * d[[key]], length(y))
  }
}

# The families whose response enters only as y - mu. A family on a half line
# is not among them: there the location and the response are unrelated
# directions and the fallback stencil is what applies.
for (.cls in list(Gaussian1Distrib, Gaussian2Distrib, Gaussian3Distrib,
                  CauchyDistrib, LogisticDistrib, LaplaceDistrib,
                  Laplace2Distrib, EnetDistrib, PseudoHuberDistrib,
                  StudentT1Distrib, SkewNormal1Distrib,
                  SkewNormal2Distrib, SkewTDistrib, GumbelDistrib)) {
  S7::method(distrib_deriv3_y, .cls) <- loc_deriv_y_k(3L)
  S7::method(distrib_deriv4_y, .cls) <- loc_deriv_y_k(4L)
}
rm(.cls)
