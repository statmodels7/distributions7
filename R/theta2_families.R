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

#' @title The Location and Scale Block of a Second-Order Mixed Derivative
#'
#' @description
#' Computes the three components in \eqn{(\mu, \sigma)} of
#' \eqn{\partial^2\ell^{(y)}/\partial\theta^2} or of
#' \eqn{\partial^2\ell^{(yy)}/\partial\theta^2}, for a family whose response
#' enters only through \eqn{z = (y-\mu)/\sigma}. It is the arithmetic behind
#' every location-scale method of [distrib_grad_y_hess()] and
#' [distrib_hess_y_hess()], and it derives nothing new.
#'
#' @details
#' # The identity
#'
#' With \eqn{A = \sigma\ell^{(y)}} and \eqn{B = \sigma^2\ell^{(yy)}} functions
#' of \eqn{z} alone, and \eqn{\partial z/\partial\mu = -1/\sigma},
#' \eqn{\partial z/\partial\sigma = -z/\sigma}, the chain rule gives
#' \deqn{\frac{\partial^2\ell^{(y)}}{\partial\mu^2} = \frac{A''}{\sigma^3},
#'   \qquad \frac{\partial^2\ell^{(y)}}{\partial\mu\,\partial\sigma}
#'     = \frac{zA'' + 2A'}{\sigma^3}, \qquad
#'   \frac{\partial^2\ell^{(y)}}{\partial\sigma^2}
#'     = \frac{z^2A'' + 4zA' + 2A}{\sigma^3},}
#' and at the next order in the response
#' \deqn{\frac{\partial^2\ell^{(yy)}}{\partial\mu^2} = \frac{B''}{\sigma^4},
#'   \qquad \frac{\partial^2\ell^{(yy)}}{\partial\mu\,\partial\sigma}
#'     = \frac{zB'' + 3B'}{\sigma^4}, \qquad
#'   \frac{\partial^2\ell^{(yy)}}{\partial\sigma^2}
#'     = \frac{z^2B'' + 6zB' + 6B}{\sigma^4}.}
#'
#' # Why nothing new is differentiated
#'
#' \eqn{A}, \eqn{A'}, \eqn{A''} and \eqn{B}, \eqn{B'}, \eqn{B''} are the
#' family's OWN response derivatives times a power of \eqn{\sigma}:
#' \eqn{A = \sigma\ell^{(y)}}, \eqn{A' = B = \sigma^2\ell^{(yy)}},
#' \eqn{A'' = B' = \sigma^3\ell^{(yyy)}} and \eqn{B'' = \sigma^4\ell^{(yyyy)}}.
#' A family that already carries [distrib_deriv3_y()] and [distrib_deriv4_y()]
#' therefore gets both orders in the parameters for free, which is the bargain
#' [loc_scale_cross_block()] already takes at first order.
#'
#' A shape parameter beyond the two is NOT covered and falls back; see
#' [partial_theta2()], which splices these three components into the
#' differenced ones.
#'
#' @param distrib An object inheriting from class `distrib`, whose first two
#'   parameters are a location and a scale in that order.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, location first and scale second.
#' @param order `1` for the block of [distrib_grad_y()], `2` for that of
#'   [distrib_hess_y()]. Default `1`.
#'
#' @return A named list of three numeric vectors, keyed `mu_mu`,
#'   `sigma_sigma` and `mu_sigma` under the family's own parameter names.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{z = (y-\mu)/\sigma} the standardized residual, and \eqn{A} and \eqn{B}
#' the two standardized quantities \eqn{\sigma\ell^{(y)}} and
#' \eqn{\sigma^2\ell^{(yy)}}, whose derivatives are taken in \eqn{z}.
#'
#' @seealso [loc_scale_grad_y_hess()] and [partial_loc_scale_grad_y_hess()],
#'   the two method bodies built on it, [loc_scale_cross2_block()] for the
#'   first order in \eqn{\theta}, and [distrib_grad_y_hess()] for the generic.
#'
#' @keywords internal
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' blk <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
#' blk$mu_mu
#'
#' # The identity written out: A'' / sigma^3, with A'' the family's own third
#' # response derivative scaled.
#' s <- 1.2
#' A2 <- s^3 * distrib_deriv3_y(d, y, theta)
#' A2 / s^3
#'
#' # And the scale pair, which carries all three of A, A' and A''.
#' z <- (y - 0.3) / s
#' A <- s * distrib_grad_y(d, y, theta)
#' A1 <- s^2 * distrib_hess_y(d, y, theta)
#' c(reported = blk$sigma_sigma[1],
#'   formula = ((z^2 * A2 + 4 * z * A1 + 2 * A) / s^3)[1])
#'
#' # Order 2 reads one derivative further and is the same shape.
#' distributions7:::loc_scale_theta2_block(d, y, theta, 2L)$mu_mu
#' (s^4 * distrib_deriv4_y(d, y, theta)) / s^4
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


#' @title The Location and Scale Parts of the Response Curvature's Derivative
#'
#' @description
#' Computes \eqn{\partial\ell^{(yy)}/\partial\mu = -B'/\sigma^3} and
#' \eqn{\partial\ell^{(yy)}/\partial\sigma = -(zB' + 2B)/\sigma^3}, with
#' \eqn{B = \sigma^2\ell^{(yy)}} and \eqn{B' = \sigma^3\ell^{(yyy)}}, for a
#' family whose response enters only through \eqn{z = (y-\mu)/\sigma}. It is
#' the arithmetic behind every location-scale method of [distrib_cross2_y()].
#'
#' @details
#' The same identity [loc_scale_cross_block()] uses, one derivative further in
#' the response: it reads the family's own [distrib_deriv3_y()] instead of
#' differencing its [distrib_hess_y()]. A shape parameter beyond the two is not
#' covered, and [partial_loc_scale_cross2_y()] differences those components
#' instead.
#'
#' @param distrib An object inheriting from class `distrib`, whose first two
#'   parameters are a location and a scale in that order.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, location first and scale second.
#'
#' @return A list of two numeric vectors, UNNAMED, in location-then-scale
#'   order. The callers name them.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{z = (y-\mu)/\sigma} the standardized residual, and \eqn{A} and \eqn{B}
#' the two standardized quantities \eqn{\sigma\ell^{(y)}} and
#' \eqn{\sigma^2\ell^{(yy)}}, whose derivatives are taken in \eqn{z}.
#'
#' @seealso [loc_scale_cross2_y()] and [partial_loc_scale_cross2_y()], the two
#'   method bodies built on it, [loc_scale_theta2_block()] for the next order
#'   in \eqn{\theta}, and [distrib_cross2_y()] for the generic.
#'
#' @keywords internal
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' blk <- distributions7:::loc_scale_cross2_block(d, y, theta)
#' vapply(blk, function(z) z[1], numeric(1))
#'
#' # The identity written out.
#' s <- 1.2
#' z <- (y - 0.3) / s
#' B <- s^2 * distrib_hess_y(d, y, theta)
#' B1 <- s^3 * distrib_deriv3_y(d, y, theta)
#' c(mu = (-B1 / s^3)[1], sigma = (-(z * B1 + 2 * B) / s^3)[1])
#'
#' # It is what the family's distrib_cross2_y method returns, named.
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
loc_scale_cross2_block <- function(distrib, y, theta) {
  s <- theta[[2]]
  z <- (y - theta[[1]]) / s
  lyy <- distrib_hess_y(distrib, y, theta) + 0 * y
  l3 <- distrib_deriv3_y(distrib, y, theta) + 0 * y
  b <- s^2 * lyy
  b1 <- s^3 * l3
  list(-b1 / s^3, -(z * b1 + 2 * b) / s^3)
}

#' @title Second-Response Mixed Derivatives of a Location-Scale Family
#'
#' @description
#' The [distrib_cross2_y()] method bodies shared by the families whose response
#' enters only through \eqn{z = (y-\mu)/\sigma}. `loc_scale_cross2_y()` serves
#' a family that is location-scale in BOTH its parameters, and
#' `partial_loc_scale_cross2_y()` one carrying shape parameters beyond the two:
#' the location and scale components come from [loc_scale_cross2_block()] in
#' closed form, and each shape component from one central difference of
#' [distrib_hess_y()] through [numerical_cross2_y()].
#'
#' @details
#' Registered on the logistic, Cauchy, Gumbel and Laplace for the first body,
#' and on the pseudo-Huber, skew normal and skew t for the second. The
#' gaussian and the Student t write their own out, and every remaining
#' continuous family takes [distrib_cross2_y.continuous_distrib()].
#'
#' @param distrib An object inheriting from class `distrib`, whose first two
#'   parameters are a location and a scale in that order.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{z = (y-\mu)/\sigma} the standardized residual, and \eqn{A} and \eqn{B}
#' the two standardized quantities \eqn{\sigma\ell^{(y)}} and
#' \eqn{\sigma^2\ell^{(yy)}}, whose derivatives are taken in \eqn{z}.
#'
#' @seealso [loc_scale_cross2_block()] for the closed components,
#'   [numerical_cross2_y()] for the differenced ones, and [distrib_cross2_y()]
#'   for the generic.
#'
#' @aliases distrib_cross2_y.LogisticDistrib distrib_cross2_y.CauchyDistrib
#'   distrib_cross2_y.GumbelDistrib distrib_cross2_y.LaplaceDistrib
#'   distrib_cross2_y.PseudoHuberDistrib distrib_cross2_y.SkewNormal1Distrib
#'   distrib_cross2_y.SkewTDistrib
#'
#' @keywords internal
#'
#' @examples
#' y <- c(-0.7, 0.3, 1.4)
#'
#' # Fully location-scale: both components closed.
#' d <- cauchy_distrib()
#' theta <- list(mu = 0.3, sigma = 1.2)
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::grad(f, c(0.3, 1.2))
#'
#' # With a shape parameter, the third component is differenced and the first
#' # two are not.
#' ds <- skewnormal1_distrib()
#' ts <- list(mu = 0.3, sigma = 1.2, alpha = 0.7)
#' vapply(distrib_cross2_y(ds, y, ts), function(z) z[1], numeric(1))
#' g <- function(v) distrib_hess_y(ds, y[1], list(mu = v[1], sigma = v[2],
#'                                                alpha = v[3]))
#' numDeriv::grad(g, c(0.3, 1.2, 0.7))
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


#' @title Second-Order Mixed Derivatives of a Location-Scale Family
#'
#' @description
#' The [distrib_grad_y_hess()] and [distrib_hess_y_hess()] method bodies shared
#' by the families that are location-scale in BOTH their parameters. Each one
#' calls [loc_scale_theta2_block()] at its own order and returns the three
#' components in [hess_names()]'s order, so neither derives anything of its
#' own.
#'
#' @details
#' Registered on the logistic, Cauchy, Gumbel and Laplace. A family carrying a
#' shape parameter beyond the two takes
#' [partial_loc_scale_grad_y_hess()] instead, where the pairs touching that
#' shape are differenced.
#'
#' Both bodies read the family's own third and fourth response derivatives, so
#' a family supplying [distrib_deriv3_y()] and [distrib_deriv4_y()] gets these
#' in closed form with no further algebra.
#'
#' @param distrib An object inheriting from class `distrib`, whose first two
#'   parameters are a location and a scale in that order.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{z = (y-\mu)/\sigma} the standardized residual, and \eqn{A} and \eqn{B}
#' the two standardized quantities \eqn{\sigma\ell^{(y)}} and
#' \eqn{\sigma^2\ell^{(yy)}}, whose derivatives are taken in \eqn{z}.
#'
#' @seealso [loc_scale_theta2_block()], which does the work,
#'   [partial_loc_scale_grad_y_hess()] for a family with a shape parameter, and
#'   [distrib_grad_y_hess()] and [distrib_hess_y_hess()] for the generics.
#'
#' @aliases distrib_grad_y_hess.LogisticDistrib
#'   distrib_grad_y_hess.CauchyDistrib distrib_grad_y_hess.GumbelDistrib
#'   distrib_grad_y_hess.LaplaceDistrib distrib_hess_y_hess.LogisticDistrib
#'   distrib_hess_y_hess.CauchyDistrib distrib_hess_y_hess.GumbelDistrib
#'   distrib_hess_y_hess.LaplaceDistrib
#'
#' @keywords internal
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' g3 <- distrib_grad_y_hess(d, y, theta)
#' vapply(g3, function(z) z[1], numeric(1))
#'
#' # Against a numerical Hessian of the response gradient.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(f, c(0.3, 1.2))
#'
#' # The fourth order takes the same body one order up.
#' g4 <- distrib_hess_y_hess(d, y, theta)
#' vapply(g4, function(z) z[1], numeric(1))
#' h <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(h, c(0.3, 1.2))
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


#' @title Second-Order Mixed Derivatives With a Shape Parameter
#'
#' @description
#' The [distrib_grad_y_hess()] and [distrib_hess_y_hess()] method bodies for a
#' family carrying shape parameters beyond its location and scale. The three
#' pairs in \eqn{(\mu, \sigma)} come from [loc_scale_theta2_block()] in closed
#' form, and every pair touching a shape parameter from one central difference
#' of the analytic first-order component. Both delegate to [partial_theta2()],
#' which does the splicing.
#'
#' @details
#' Registered on the Student t, the pseudo-Huber, the skew normal and the skew
#' t. Measured against a numerical Hessian of the response derivative, the
#' closed pairs agree to about \eqn{10^{-11}} and the differenced ones to
#' between \eqn{10^{-7}} and \eqn{10^{-6}}, which is one stencil's accuracy.
#'
#' @param distrib An object inheriting from class `distrib`, whose first two
#'   parameters are a location and a scale in that order and whose remaining
#'   parameters are shapes.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{z = (y-\mu)/\sigma} the standardized residual, and \eqn{A} and \eqn{B}
#' the two standardized quantities \eqn{\sigma\ell^{(y)}} and
#' \eqn{\sigma^2\ell^{(yy)}}, whose derivatives are taken in \eqn{z}.
#'
#' @seealso [partial_theta2()], which splices the two halves,
#'   [loc_scale_grad_y_hess()] for a family with no shape parameter, and
#'   [distrib_grad_y_hess()] and [distrib_hess_y_hess()] for the generics.
#'
#' @aliases distrib_grad_y_hess.StudentT1Distrib
#'   distrib_grad_y_hess.PseudoHuberDistrib
#'   distrib_grad_y_hess.SkewNormal1Distrib distrib_grad_y_hess.SkewTDistrib
#'   distrib_hess_y_hess.StudentT1Distrib
#'   distrib_hess_y_hess.PseudoHuberDistrib
#'   distrib_hess_y_hess.SkewNormal1Distrib distrib_hess_y_hess.SkewTDistrib
#'
#' @keywords internal
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma = 1.2, nu = 6)
#'
#' g3 <- distrib_grad_y_hess(d, y, theta)
#' names(g3)
#' vapply(g3, function(z) z[1], numeric(1))
#'
#' # Against a numerical Hessian of the response gradient. The three pairs in
#' # (mu, sigma) are closed and the three touching nu are differenced.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2],
#'                                               nu = v[3]))
#' numDeriv::hessian(f, c(0.3, 1.2, 6))
#'
#' # The closed three are exactly the location-scale block's.
#' cl <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
#' all(vapply(names(cl), function(k) identical(g3[[k]], cl[[k]]), logical(1)))
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


#' @title Splice the Closed Location-Scale Pairs Into the Fallback
#'
#' @description
#' Differences every component through [numerical_theta2_y()] and then replaces
#' the three that [loc_scale_theta2_block()] gives in closed form. It is what
#' [partial_loc_scale_grad_y_hess()] and [partial_loc_scale_hess_y_hess()] both
#' run, at order 1 and 2 respectively.
#'
#' @details
#' Differencing every component and overwriting three of them costs three
#' unnecessary differences, and is how the two halves are kept keyed the same
#' way without a second enumeration: [numerical_theta2_y()] returns
#' [hess_names()]'s full set, and the closed block's three names index into it
#' directly.
#'
#' @param distrib An object inheriting from class `distrib`, whose first two
#'   parameters are a location and a scale.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order `1` to splice into the derivative of [distrib_grad_y()], `2`
#'   into that of [distrib_hess_y()]. It selects both the first-order reference
#'   differenced, [distrib_cross_y()] or [distrib_cross2_y()], and the order of
#'   the closed block.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{z = (y-\mu)/\sigma} the standardized residual, and \eqn{A} and \eqn{B}
#' the two standardized quantities \eqn{\sigma\ell^{(y)}} and
#' \eqn{\sigma^2\ell^{(yy)}}, whose derivatives are taken in \eqn{z}.
#'
#' @seealso [loc_scale_theta2_block()] for the closed three,
#'   [numerical_theta2_y()] for the differenced rest, and
#'   [partial_loc_scale_grad_y_hess()], its caller.
#'
#' @keywords internal
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma = 1.2, alpha = 0.7, nu = 6)
#'
#' out <- distributions7:::partial_theta2(d, y, theta, 1L)
#' names(out)
#'
#' # The three location-scale pairs are spliced in unchanged.
#' cl <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
#' all(vapply(names(cl), function(k) identical(out[[k]], cl[[k]]), logical(1)))
#'
#' # And it is exactly what the family's method returns.
#' identical(out, distrib_grad_y_hess(d, y, theta))
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
