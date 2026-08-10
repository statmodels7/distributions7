#' @include cross_theta2_derivatives.R cross_derivatives_families.R y_higher.R
NULL

# ===========================================================================
# The location-scale identity, one order further.
#
# A family whose response enters only through z = (y - mu)/sigma has
#
#   l^(y)  = A(z)/sigma^2 * sigma = A(z)/sigma      with A = g'
#   l^(yy) = B(z)/sigma^2                           with B = g'' = A'
#
# and dz/dmu = -1/sigma, dz/dsigma = -z/sigma. Everything the location and
# the scale do to either quantity follows by the chain rule:
#
#   d2 l^(y)/dmu^2      =  A''/sigma^3
#   d2 l^(y)/dmu dsigma = (z A'' + 2 A')/sigma^3
#   d2 l^(y)/dsigma^2   = (z^2 A'' + 4 z A' + 2 A)/sigma^3
#
#   d2 l^(yy)/dmu^2      =  B''/sigma^4
#   d2 l^(yy)/dmu dsigma = (z B'' + 3 B')/sigma^4
#   d2 l^(yy)/dsigma^2   = (z^2 B'' + 6 z B' + 6 B)/sigma^4
#
# and A, A', A'', B, B', B'' are the family's own response derivatives times
# a power of sigma: A = sigma l^(y), A' = B = sigma^2 l^(yy), A'' = B' =
# sigma^3 l^(yyy), B'' = sigma^4 l^(yyyy). So no new algebra is needed at all
# -- the same bargain the first-order block already takes.
#
# A shape parameter is not covered and falls back, as it does there.
# ===========================================================================

#' The Location and Scale Block of a Second-Order Mixed Derivative
#'
#' @description
#' The three components in \eqn{(\mu, \sigma)} of
#' \eqn{\partial^2\ell^{(y)}/\partial\theta^2} or of
#' \eqn{\partial^2\ell^{(yy)}/\partial\theta^2}, for a family whose response
#' enters only through \eqn{z = (y-\mu)/\sigma}.
#'
#' @details
#' The quantities are the family's own response derivatives scaled by a power
#' of \eqn{\sigma}: writing \eqn{A = \sigma\ell^{(y)}}, its derivatives in
#' \eqn{z} are \eqn{\sigma^2\ell^{(yy)}} and \eqn{\sigma^3\ell^{(yyy)}}, and
#' for \eqn{B = \sigma^2\ell^{(yy)}} they are \eqn{\sigma^3\ell^{(yyy)}} and
#' \eqn{\sigma^4\ell^{(yyyy)}}. Nothing new is differentiated.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, location first and scale second.
#' @param order \code{1} for the block of \code{\link{distrib_grad_y}},
#'   \code{2} for that of \code{\link{distrib_hess_y}}.
#'
#' @return A named list of three components, keyed
#'   \code{mu_mu}, \code{sigma_sigma}, \code{mu_sigma} under the family's own
#'   parameter names.
#'
#' @seealso \code{\link{distrib_grad_y_hess}}
#' @keywords internal
loc_scale_theta2_block <- function(distrib, y, theta, order = 1L) {
  params <- distrib@params
  s <- theta[[2]]
  z <- (y - theta[[1]]) / s
  ly <- distrib_grad_y(distrib, y, theta) + 0 * y
  lyy <- distrib_hess_y(distrib, y, theta) + 0 * y
  l3 <- distrib_deriv3_y(distrib, y, theta) + 0 * y
  if (order == 1L) {
    a <- s * ly
    a1 <- s^2 * lyy
    a2 <- s^3 * l3
    out <- list(a2 / s^3,
                (z^2 * a2 + 4 * z * a1 + 2 * a) / s^3,
                (z * a2 + 2 * a1) / s^3)
  } else {
    l4 <- distrib_deriv4_y(distrib, y, theta) + 0 * y
    b <- s^2 * lyy
    b1 <- s^3 * l3
    b2 <- s^4 * l4
    out <- list(b2 / s^4,
                (z^2 * b2 + 6 * z * b1 + 6 * b) / s^4,
                (z * b2 + 3 * b1) / s^4)
  }
  stats::setNames(out, c(paste0(params[1], "_", params[1]),
                         paste0(params[2], "_", params[2]),
                         paste0(params[1], "_", params[2])))
}


#' The Location and Scale Components of the Response Curvature's Derivative
#'
#' @description
#' \eqn{\partial\ell^{(yy)}/\partial\mu = -B'/\sigma^3} and
#' \eqn{\partial\ell^{(yy)}/\partial\sigma = -(zB' + 2B)/\sigma^3}, with
#' \eqn{B = \sigma^2\ell^{(yy)}} and \eqn{B' = \sigma^3\ell^{(yyy)}}, for a
#' family whose response enters only through \eqn{z = (y-\mu)/\sigma}.
#'
#' @details
#' The same identity as \code{\link{loc_scale_cross_block}} one derivative
#' further in the response, and it reads the family's own third response
#' derivative rather than differencing its second.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, location first and scale second.
#'
#' @return A list of two component vectors, unnamed.
#'
#' @seealso \code{\link{distrib_cross2_y}}
#' @keywords internal
loc_scale_cross2_block <- function(distrib, y, theta) {
  s <- theta[[2]]
  z <- (y - theta[[1]]) / s
  lyy <- distrib_hess_y(distrib, y, theta) + 0 * y
  l3 <- distrib_deriv3_y(distrib, y, theta) + 0 * y
  b <- s^2 * lyy
  b1 <- s^3 * l3
  list(-b1 / s^3, -(z * b1 + 2 * b) / s^3)
}

#' Second-Response Mixed Derivatives of a Location-Scale Family
#'
#' @description
#' The \code{\link{distrib_cross2_y}} body shared by the families that are
#' location-scale in both their parameters, and the partial form for those
#' with a shape parameter beyond the two.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per parameter.
#'
#' @seealso \code{\link{loc_scale_cross2_block}}
#' @keywords internal
loc_scale_cross2_y <- function(distrib, y, theta,
                               scale = c("parameter", "link"), ...) {
  stats::setNames(loc_scale_cross2_block(distrib, y, theta), distrib@params)
}

#' @rdname loc_scale_cross2_y
#' @keywords internal
partial_loc_scale_cross2_y <- function(distrib, y, theta,
                                       scale = c("parameter", "link"), ...) {
  params <- distrib@params
  rest <- params[-(1:2)]
  out <- c(loc_scale_cross2_block(distrib, y, theta),
           numerical_cross2_y(distrib, y, theta, which = rest))
  stats::setNames(out, params)
}


#' Second-Order Mixed Derivatives of a Location-Scale Family
#'
#' @description
#' The \code{\link{distrib_grad_y_hess}} and
#' \code{\link{distrib_hess_y_hess}} bodies shared by the families that are
#' location-scale in both their parameters.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#'
#' @return A named list keyed by parameter pair.
#'
#' @seealso \code{\link{loc_scale_theta2_block}}
#' @keywords internal
loc_scale_grad_y_hess <- function(distrib, y, theta,
                                  scale = c("parameter", "link"), ...) {
  loc_scale_theta2_block(distrib, y, theta, 1L)[hess_names(distrib@params)]
}

#' @rdname loc_scale_grad_y_hess
#' @keywords internal
loc_scale_hess_y_hess <- function(distrib, y, theta,
                                  scale = c("parameter", "link"), ...) {
  loc_scale_theta2_block(distrib, y, theta, 2L)[hess_names(distrib@params)]
}


#' Second-Order Mixed Derivatives When Only Two Parameters Are Location-Scale
#'
#' @description
#' The location and scale pairs in closed form and every pair touching a shape
#' parameter by one central difference of the analytic first-order component.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#'
#' @return A named list keyed by parameter pair.
#'
#' @seealso \code{\link{loc_scale_theta2_block}}
#' @keywords internal
partial_loc_scale_grad_y_hess <- function(distrib, y, theta,
                                          scale = c("parameter", "link"),
                                          ...) {
  partial_theta2(distrib, y, theta, 1L)
}

#' @rdname partial_loc_scale_grad_y_hess
#' @keywords internal
partial_loc_scale_hess_y_hess <- function(distrib, y, theta,
                                          scale = c("parameter", "link"),
                                          ...) {
  partial_theta2(distrib, y, theta, 2L)
}


#' Splice the Closed Location-Scale Pairs Into the Fallback
#'
#' @description
#' Takes the differenced components and replaces the three that the
#' location-scale identity gives in closed form.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order \code{1} or \code{2}.
#'
#' @return A named list keyed by parameter pair.
#'
#' @keywords internal
partial_theta2 <- function(distrib, y, theta, order) {
  inner <- if (order == 1L) {
    function(th) distrib_cross_y(distrib, y, th)
  } else {
    function(th) distrib_cross2_y(distrib, y, th)
  }
  out <- numerical_theta2_y(distrib, y, theta, inner)
  closed <- loc_scale_theta2_block(distrib, y, theta, order)
  out[names(closed)] <- closed
  out
}


# --- the location-scale families -------------------------------------------
#
# The same list the first-order block is registered on, minus the ones whose
# second parameter is not a scale. A family with a shape parameter beyond the
# two takes the partial form, where the pairs touching the shape fall back.

#' @name distrib_grad_y_hess.LogisticDistrib
#' @rdname loc_scale_grad_y_hess
#' @keywords internal
S7::method(distrib_cross2_y, LogisticDistrib) <- loc_scale_cross2_y
S7::method(distrib_cross2_y, CauchyDistrib) <- loc_scale_cross2_y
S7::method(distrib_cross2_y, GumbelDistrib) <- loc_scale_cross2_y
S7::method(distrib_cross2_y, LaplaceDistrib) <- loc_scale_cross2_y
S7::method(distrib_cross2_y, PseudoHuberDistrib) <- partial_loc_scale_cross2_y
S7::method(distrib_cross2_y, SkewNormal1Distrib) <- partial_loc_scale_cross2_y
S7::method(distrib_cross2_y, SkewTDistrib) <- partial_loc_scale_cross2_y

S7::method(distrib_grad_y_hess, LogisticDistrib) <- loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, LogisticDistrib) <- loc_scale_hess_y_hess

S7::method(distrib_grad_y_hess, CauchyDistrib) <- loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, CauchyDistrib) <- loc_scale_hess_y_hess

S7::method(distrib_grad_y_hess, GumbelDistrib) <- loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, GumbelDistrib) <- loc_scale_hess_y_hess

S7::method(distrib_grad_y_hess, StudentT1Distrib) <-
  partial_loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, StudentT1Distrib) <-
  partial_loc_scale_hess_y_hess

S7::method(distrib_grad_y_hess, PseudoHuberDistrib) <-
  partial_loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, PseudoHuberDistrib) <-
  partial_loc_scale_hess_y_hess

S7::method(distrib_grad_y_hess, SkewNormal1Distrib) <-
  partial_loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, SkewNormal1Distrib) <-
  partial_loc_scale_hess_y_hess

S7::method(distrib_grad_y_hess, SkewTDistrib) <- partial_loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, SkewTDistrib) <- partial_loc_scale_hess_y_hess
