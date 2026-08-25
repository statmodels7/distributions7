#' @include theta2_mapped.R reparam_maps.R
NULL

# ===========================================================================
# What the census left, in the two shapes that need no new algebra.
#
# The FIRST is a family written out in its own parametrization whose map onto
# a parent is already tabulated in reparam_maps.R. gaussian2, gaussian3,
# laplace2 and invgauss2 have their own classes and their own kernels, so they
# are not ReparamContinuousDistrib and do not inherit its methods, but the map
# tables are there and the chain rule is the same one.
#
# The SECOND is a family with no location, where z = y/sigma. The derivation
# of the scale components never used dz/dmu, only that sigma is a scale, so
# every sigma formula holds unchanged. That covers the exponential, the
# Weibull's scale and the generalized Pareto's -- the same three the
# first-order block already reaches this way.
# ===========================================================================

#' @title The Mixed Grid of a Family Written Out Against a Tabulated Map
#'
#' @description
#' Builds the [distrib_cross2_y()], [distrib_grad_y_hess()] and
#' [distrib_hess_y_hess()] method bodies for a family that carries its own
#' class and kernels while its map onto a parent is tabulated. Each body is the
#' chain rule of [mapped_cross2_y()] or [mapped_theta2()] with the parent, the
#' map and the partial tables closed over, so no new algebra is involved.
#'
#' @details
#' The families it serves have their own classes and their own compiled
#' kernels, so they are NOT `ReparamContinuousDistrib` and do not inherit its
#' methods; what they share with it is the map, which `reparam_maps.R` already
#' tabulates. Registered on `gaussian2`, `gaussian3`, `laplace2` and
#' `invgauss2`.
#'
#' @param parent The parent distribution, built once and closed over.
#' @param th_par A function of the new parameters returning the parent's, as a
#'   named list.
#' @param tables A function of the new parameters returning the map's keyed
#'   partial tables, one of the `md_*` functions of `reparam_maps.R`.
#'
#' @return A list of three method bodies, named `cross2`, `grad2` and `hess2`,
#'   each with the generics' signature.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{\sigma} the scale and \eqn{z} the standardized response, \eqn{y/\sigma}
#' where there is no location and \eqn{(y-\mu)/\sigma} where there is.
#'
#' @seealso [mapped_theta2()] and [mapped_cross2_y()], which do the work,
#'   [reparam_map_derivs()] for the tables, and [scale_only_theta2_methods()]
#'   for the other shape this file covers.
#'
#' @aliases distrib_cross2_y.Gaussian2Distrib
#'   distrib_grad_y_hess.Gaussian2Distrib distrib_hess_y_hess.Gaussian2Distrib
#'   distrib_cross2_y.Gaussian3Distrib distrib_grad_y_hess.Gaussian3Distrib
#'   distrib_hess_y_hess.Gaussian3Distrib distrib_cross2_y.Laplace2Distrib
#'   distrib_grad_y_hess.Laplace2Distrib distrib_hess_y_hess.Laplace2Distrib
#'   distrib_cross2_y.InvGauss2Distrib distrib_grad_y_hess.InvGauss2Distrib
#'   distrib_hess_y_hess.InvGauss2Distrib
#'
#' @keywords internal
#'
#' @examples
#' # gaussian2 carries (mu, sigma2) against a parent in (mu, sigma).
#' d <- gaussian2_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma2 = 1.44)
#'
#' g3 <- distrib_grad_y_hess(d, y, theta)
#' vapply(g3, function(z) z[1], numeric(1))
#'
#' # Against a numerical Hessian of the response gradient in the family's own
#' # parameters, which shares none of the map's arithmetic.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(0.3, 1.44))
#'
#' # laplace2 carries the rate, so its map is 1 / lambda.
#' dl <- laplace2_distrib()
#' vapply(distrib_cross2_y(dl, c(-0.7, 1.1, 1.4), list(mu = 0.3, lambda = 0.8)),
#'        function(z) z[1], numeric(1))
mapped_theta2_methods <- function(parent, th_par, tables) {
  list(
    cross2 = function(distrib, y, theta, scale = c("parameter", "link"),
                      ...) {
      mapped_cross2_y(distrib, parent, th_par(theta), tables(theta), y)
    },
    grad2 = function(distrib, y, theta, scale = c("parameter", "link"),
                     ...) {
      th <- th_par(theta)
      mapped_theta2(distrib, parent, th, tables(theta), y,
                    distrib_cross_y(parent, y, th),
                    distrib_grad_y_hess(parent, y, th))
    },
    hess2 = function(distrib, y, theta, scale = c("parameter", "link"),
                     ...) {
      th <- th_par(theta)
      mapped_theta2(distrib, parent, th, tables(theta), y,
                    distrib_cross2_y(parent, y, th),
                    distrib_hess_y_hess(parent, y, th))
    })
}


#' @title The Scale Component of a Family With No Location
#'
#' @description
#' Computes \eqn{\partial^2\ell^{(y)}/\partial\sigma^2} or
#' \eqn{\partial^2\ell^{(yy)}/\partial\sigma^2} for a family whose response
#' enters only through \eqn{z = y/\sigma}. It is [loc_scale_theta2_block()]'s
#' scale component, and the same formula: the derivation never used the
#' location, only that \eqn{\sigma} is a scale.
#'
#' @details
#' With \eqn{A = \sigma\ell^{(y)}} and \eqn{B = \sigma^2\ell^{(yy)}},
#' \deqn{\frac{\partial^2\ell^{(y)}}{\partial\sigma^2}
#'   = \frac{z^2A'' + 4zA' + 2A}{\sigma^3}, \qquad
#'   \frac{\partial^2\ell^{(yy)}}{\partial\sigma^2}
#'   = \frac{z^2B'' + 6zB' + 6B}{\sigma^4},}
#' and the primed quantities are the family's own third and fourth response
#' derivatives times a power of \eqn{\sigma}, so nothing is differentiated
#' here. Any parameter beyond the scale is a SHAPE and is not covered;
#' [scale_only_theta2_methods()] differences those.
#'
#' @param distrib An object inheriting from class `distrib`, with a scale and
#'   no location.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order `1` for the block of [distrib_grad_y()], `2` for that of
#'   [distrib_hess_y()]. Default `1`.
#' @param at The index of the scale parameter in `distrib@params`. Default `1`.
#'
#' @return A numeric vector as long as `y`: the scale's own pair alone, not a
#'   list.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{\sigma} the scale and \eqn{z} the standardized response, \eqn{y/\sigma}
#' where there is no location and \eqn{(y-\mu)/\sigma} where there is.
#'
#' @seealso [loc_scale_theta2_block()], the location-scale form,
#'   [scale_only_theta2_methods()], which splices it into the fallback, and
#'   [scale_only_cross2_method()] for the first order in \eqn{\theta}.
#'
#' @keywords internal
#'
#' @examples
#' # The exponential has a scale and nothing else.
#' d <- exponential_distrib()
#' y <- c(0.4, 1.1, 2.3)
#' theta <- list(mu = 1.5)
#'
#' distributions7:::scale_only_theta2(d, y, theta, 1L, 1L)
#'
#' # The identity written out.
#' s <- 1.5
#' z <- y / s
#' A <- s * distrib_grad_y(d, y, theta)
#' A1 <- s^2 * distrib_hess_y(d, y, theta)
#' A2 <- s^3 * distrib_deriv3_y(d, y, theta)
#' (z^2 * A2 + 4 * z * A1 + 2 * A) / s^3
#'
#' # And what the family's own method reports.
#' distrib_grad_y_hess(d, y, theta)$mu_mu
scale_only_theta2 <- function(distrib, y, theta, order = 1L, at = 1L) {
  s <- theta[[at]]
  z <- y / s
  lyy <- distrib_hess_y(distrib, y, theta) + 0 * y
  l3 <- distrib_deriv3_y(distrib, y, theta) + 0 * y
  if (order == 1L) {
    ly <- distrib_grad_y(distrib, y, theta) + 0 * y
    a <- s * ly
    a1 <- s^2 * lyy
    a2 <- s^3 * l3
    (z^2 * a2 + 4 * z * a1 + 2 * a) / s^3
  } else {
    l4 <- distrib_deriv4_y(distrib, y, theta) + 0 * y
    b <- s^2 * lyy
    b1 <- s^3 * l3
    b2 <- s^4 * l4
    (z^2 * b2 + 6 * z * b1 + 6 * b) / s^4
  }
}


#' @title The Mixed Grid of a Family With No Location
#'
#' @description
#' Builds the [distrib_grad_y_hess()] and [distrib_hess_y_hess()] method bodies
#' for a family with a scale and no location: the scale's own pair comes from
#' [scale_only_theta2()] in closed form, and every pair touching a shape
#' parameter from one central difference of the analytic first-order component
#' through [numerical_theta2_y()].
#'
#' @details
#' A one-parameter family returns the closed component alone and never
#' differences anything. Registered on the exponential, which is that case, and
#' on the Weibull and the generalized Pareto, whose second parameter is a
#' shape.
#'
#' @param at The index of the scale parameter in `distrib@params`. Default `1`.
#'
#' @return A list of two method bodies, named `grad2` and `hess2`, each with
#'   the generics' signature.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{\sigma} the scale and \eqn{z} the standardized response, \eqn{y/\sigma}
#' where there is no location and \eqn{(y-\mu)/\sigma} where there is.
#'
#' @seealso [scale_only_theta2()] for the closed component,
#'   [numerical_theta2_y()] for the differenced ones, and [partial_theta2()],
#'   the location-scale counterpart.
#'
#' @aliases distrib_grad_y_hess.ExponentialDistrib
#'   distrib_hess_y_hess.ExponentialDistrib
#'   distrib_grad_y_hess.Weibull1Distrib distrib_hess_y_hess.Weibull1Distrib
#'   distrib_grad_y_hess.GPDDistrib distrib_hess_y_hess.GPDDistrib
#'
#' @keywords internal
#'
#' @examples
#' y <- c(0.4, 1.1, 2.3)
#'
#' # One parameter: the closed component alone, nothing differenced.
#' d <- exponential_distrib()
#' distrib_grad_y_hess(d, y, list(mu = 1.5))
#'
#' # Two: the scale pair closed, the two touching the shape differenced.
#' dw <- weibull1_distrib()
#' theta <- list(mu = 1.5, sigma = 1.3)
#' g <- distrib_grad_y_hess(dw, y, theta)
#' names(g)
#' vapply(g, function(z) z[1], numeric(1))
#'
#' # Against a numerical Hessian of the response gradient.
#' f <- function(v) distrib_grad_y(dw, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(f, c(1.5, 1.3))
#'
#' # The scale pair is exactly the closed component.
#' identical(g$mu_mu,
#'           distributions7:::scale_only_theta2(dw, y, theta, 1L, 1L))
scale_only_theta2_methods <- function(at = 1L) {
  body <- function(order) {
    function(distrib, y, theta, scale = c("parameter", "link"), ...) {
      p <- distrib@params
      nm <- paste(p[at], p[at], sep = "_")
      closed <- scale_only_theta2(distrib, y, theta, order, at)
      if (length(p) == 1L) return(stats::setNames(list(closed), nm))
      inner <- if (order == 1L) {
        function(th) distrib_cross_y(distrib, y, th)
      } else {
        function(th) distrib_cross2_y(distrib, y, th)
      }
      out <- numerical_theta2_y(distrib, y, theta, inner)
      out[[nm]] <- closed
      out
    }
  }
  list(grad2 = body(1L), hess2 = body(2L))
}


#' @title The Scale Component of the Response Curvature's Derivative
#'
#' @description
#' Builds the [distrib_cross2_y()] method body for a family with a scale and no
#' location. The scale's component is the identity
#' [loc_scale_cross2_block()] uses, read at \eqn{z = y/\sigma}:
#' \eqn{-(zB' + 2B)/\sigma^3} with \eqn{B = \sigma^2\ell^{(yy)}} and
#' \eqn{B' = \sigma^3\ell^{(yyy)}}. Any shape component is differenced through
#' [numerical_cross2_y()].
#'
#' @details
#' A one-parameter family returns the closed component alone. Registered on the
#' exponential, the Weibull and the generalized Pareto.
#'
#' @param at The index of the scale parameter in `distrib@params`. Default `1`.
#'
#' @return A method body with [distrib_cross2_y()]'s signature, returning a
#'   named list with one numeric vector per parameter.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} its first and second derivatives in the response,
#' \eqn{\sigma} the scale and \eqn{z} the standardized response, \eqn{y/\sigma}
#' where there is no location and \eqn{(y-\mu)/\sigma} where there is.
#'
#' @seealso [loc_scale_cross2_block()], the location-scale form,
#'   [scale_only_theta2_methods()] for the next order in \eqn{\theta}, and
#'   [distrib_cross2_y()] for the generic.
#'
#' @aliases distrib_cross2_y.ExponentialDistrib
#'   distrib_cross2_y.Weibull1Distrib distrib_cross2_y.GPDDistrib
#'
#' @keywords internal
#'
#' @examples
#' y <- c(0.4, 1.1, 2.3)
#'
#' # The generalized Pareto: a scale and a shape.
#' d <- gpd_distrib()
#' theta <- list(sigma = 1.5, xi = 0.3)
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(sigma = v[1], xi = v[2]))
#' numDeriv::grad(f, c(1.5, 0.3))
#'
#' # The exponential has one parameter, so nothing is differenced.
#' distrib_cross2_y(exponential_distrib(), y, list(mu = 1.5))
scale_only_cross2_method <- function(at = 1L) {
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    p <- distrib@params
    s <- theta[[at]]
    z <- y / s
    lyy <- distrib_hess_y(distrib, y, theta) + 0 * y
    l3 <- distrib_deriv3_y(distrib, y, theta) + 0 * y
    closed <- -(z * s^3 * l3 + 2 * s^2 * lyy) / s^3
    if (length(p) == 1L) return(stats::setNames(list(closed), p))
    out <- numerical_cross2_y(distrib, y, theta)
    out[[p[at]]] <- closed
    out
  }
}


# --- the four written-out reparametrizations --------------------------------

.m <- mapped_theta2_methods(
  gaussian1_distrib(),
  function(th) list(mu = th[[1]], sigma = sqrt(th[[2]])),
  md_gaussian2)
S7::method(distrib_cross2_y, Gaussian2Distrib) <- .m$cross2
S7::method(distrib_grad_y_hess, Gaussian2Distrib) <- .m$grad2
S7::method(distrib_hess_y_hess, Gaussian2Distrib) <- .m$hess2

.m <- mapped_theta2_methods(
  gaussian1_distrib(),
  function(th) list(mu = th[[1]], sigma = 1 / sqrt(th[[2]])),
  md_gaussian3)
S7::method(distrib_cross2_y, Gaussian3Distrib) <- .m$cross2
S7::method(distrib_grad_y_hess, Gaussian3Distrib) <- .m$grad2
S7::method(distrib_hess_y_hess, Gaussian3Distrib) <- .m$hess2

.m <- mapped_theta2_methods(
  laplace_distrib(),
  function(th) list(mu = th[[1]], sigma = 1 / th[[2]]),
  md_laplace2)
S7::method(distrib_cross2_y, Laplace2Distrib) <- .m$cross2
S7::method(distrib_grad_y_hess, Laplace2Distrib) <- .m$grad2
S7::method(distrib_hess_y_hess, Laplace2Distrib) <- .m$hess2

.m <- mapped_theta2_methods(
  invgauss1_distrib(),
  function(th) list(mu = th[[1]], phi = 1 / th[[2]]),
  md_invgauss2)
S7::method(distrib_cross2_y, InvGauss2Distrib) <- .m$cross2
S7::method(distrib_grad_y_hess, InvGauss2Distrib) <- .m$grad2
S7::method(distrib_hess_y_hess, InvGauss2Distrib) <- .m$hess2

rm(.m)


# --- the three families with no location ------------------------------------

.s <- scale_only_theta2_methods(1L)
S7::method(distrib_cross2_y, ExponentialDistrib) <- scale_only_cross2_method(1L)
S7::method(distrib_grad_y_hess, ExponentialDistrib) <- .s$grad2
S7::method(distrib_hess_y_hess, ExponentialDistrib) <- .s$hess2

S7::method(distrib_cross2_y, Weibull1Distrib) <- scale_only_cross2_method(1L)
S7::method(distrib_grad_y_hess, Weibull1Distrib) <- .s$grad2
S7::method(distrib_hess_y_hess, Weibull1Distrib) <- .s$hess2

S7::method(distrib_cross2_y, GPDDistrib) <- scale_only_cross2_method(1L)
S7::method(distrib_grad_y_hess, GPDDistrib) <- .s$grad2
S7::method(distrib_hess_y_hess, GPDDistrib) <- .s$hess2

rm(.s)


# --- the lognormal, which is a gaussian at log y ----------------------------
#
# With t = log y the density is that of a gaussian in t, times the Jacobian:
#
#   l(y) = l_N(t; mu, sigma2) - log y
#
# The transformation carries NO parameter, so t does not move with theta and
# every theta-derivative is the gaussian's own at t. What the response
# derivatives carry is the Jacobian:
#
#   l^(y)  = (g' - 1)/y
#   l^(yy) = (g'' - g' + 1)/y^2
#
# with g the gaussian's log-density in t. Differentiating in theta kills the
# constants and leaves
#
#   d l^(y)/dth   = (dg'/dth)/y
#   d l^(yy)/dth  = (dg''/dth - dg'/dth)/y^2
#
# and the same two lines at second order. The parent is gaussian2, which
# carries the same (mu, sigma2) the lognormal does, so theta passes through
# unchanged.

#' @title The Mixed Grid of the Lognormal
#'
#' @description
#' Computes the gaussian's components at \eqn{t = \log y}, carried by the
#' Jacobian of the transformation. It serves all three of the lognormal's mixed
#' methods, [distrib_cross2_y()], [distrib_grad_y_hess()] and
#' [distrib_hess_y_hess()], with `order` selecting the response derivative and
#' `second` the order in \eqn{\theta}.
#'
#' @details
#' # Why theta passes through untouched
#'
#' With \eqn{t = \log y} the log-density is the gaussian's in \eqn{t} plus the
#' log-Jacobian, \eqn{\ell(y) = \ell_N(t;\mu,\sigma^2) - \log y}. The
#' transformation carries NO parameter, so \eqn{t} does not move with
#' \eqn{\theta} and every \eqn{\theta}-derivative is the gaussian's own at
#' \eqn{t}. The parent is `gaussian2`, which carries the same
#' \eqn{(\mu, \sigma^2)} the lognormal does, so `theta` needs no map at all.
#'
#' # What the Jacobian does carry
#'
#' Writing \eqn{g} for the gaussian's log-density in \eqn{t},
#' \deqn{\ell^{(y)} = \frac{g' - 1}{y}, \qquad
#'   \ell^{(yy)} = \frac{g'' - g' + 1}{y^2},}
#' and differentiating in \eqn{\theta} kills the constants, leaving
#' \deqn{\frac{\partial\ell^{(y)}}{\partial\theta} = \frac{1}{y}
#'     \frac{\partial g'}{\partial\theta}, \qquad
#'   \frac{\partial\ell^{(yy)}}{\partial\theta} = \frac{1}{y^2}\left(
#'     \frac{\partial g''}{\partial\theta} -
#'     \frac{\partial g'}{\partial\theta}\right),}
#' with the same two lines at second order in \eqn{\theta}.
#'
#' @param y A numeric vector of observations, strictly positive.
#' @param theta A named list containing `mu` and `sigma2`, read by the gaussian
#'   unchanged.
#' @param order `1` for the derivative of [distrib_grad_y()], `2` for that of
#'   [distrib_hess_y()].
#' @param second `FALSE` for one derivative in \eqn{\theta}, keyed by
#'   parameter, and `TRUE` for two, keyed by parameter pair.
#'
#' @return A named list of numeric vectors, keyed by parameter when `second` is
#'   `FALSE` and by parameter pair when it is `TRUE`.
#'
#' @section Notation:
#' \eqn{t = \log y}, \eqn{g} is the gaussian's log-density in \eqn{t}, and a
#' prime denotes a derivative in \eqn{t}.
#'
#' @seealso [distrib_cross2_y.Lognormal1Distrib()], the three methods it
#'   serves, and [mapped_theta2()], the route a family whose transformation
#'   DOES carry a parameter takes instead.
#'
#' @keywords internal
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.4, 1.1, 2.3)
#' theta <- list(mu = 0.2, sigma2 = 0.6)
#'
#' distributions7:::lognormal_theta_chain(y, theta, 2L, FALSE)
#'
#' # Which is what the family's distrib_cross2_y reports.
#' distrib_cross2_y(d, y, theta)
#'
#' # And it is the gaussian's own component at t = log y, over y^2.
#' g2 <- gaussian2_distrib()
#' t <- log(y)
#' first <- distrib_cross_y(g2, t, theta)
#' hi <- distrib_cross2_y(g2, t, theta)
#' (hi$mu - first$mu) / y^2
lognormal_theta_chain <- function(y, theta, order, second) {
  t <- base::log(y)
  g2 <- gaussian2_distrib()
  first <- if (second) distrib_grad_y_hess(g2, t, theta) else
    distrib_cross_y(g2, t, theta)
  if (order == 1L) return(lapply(first, function(v) (v + 0 * t) / y))
  hi <- if (second) distrib_hess_y_hess(g2, t, theta) else
    distrib_cross2_y(g2, t, theta)
  stats::setNames(lapply(names(hi), function(nm)
    ((hi[[nm]] - first[[nm]]) + 0 * t) / y^2), names(hi))
}

#' @title Lognormal Second-Response and Second-Order Mixed Derivatives
#' @name distrib_cross2_y.Lognormal1Distrib
#'
#' @description
#' Reads the gaussian's own components at \eqn{t = \log y} and divides by a
#' power of \eqn{y}. The transformation carries no parameter, so \eqn{t} does
#' not move with \eqn{\theta} and only the Jacobian of the RESPONSE derivatives
#' enters; the parent is `gaussian2`, which carries the same
#' \eqn{(\mu, \sigma^2)}, so `theta` passes through unchanged. All three mixed
#' methods take this route, through [lognormal_theta_chain()].
#'
#' @param distrib A `Lognormal1Distrib` object.
#' @param y A numeric vector of observations, strictly positive.
#' @param theta A named list containing `mu` and `sigma2`.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, keyed by parameter for
#'   [distrib_cross2_y()] and by parameter pair for [distrib_grad_y_hess()] and
#'   [distrib_hess_y_hess()].
#'
#' @section Notation:
#' \eqn{t = \log y}, and \eqn{\ell^{(y)}} and \eqn{\ell^{(yy)}} are the first
#' and second derivatives of the log-density in the response.
#'
#' @seealso [lognormal_theta_chain()], which does the work,
#'   [distrib_cross2_y()] for the generic, and [transformation()], the general
#'   wrapper for a change of variable.
#'
#' @keywords internal
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.4, 1.1, 2.3)
#' theta <- list(mu = 0.2, sigma2 = 0.6)
#'
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::grad(f, c(0.2, 0.6))
#'
#' # The second order in theta, and its numerical Hessian.
#' g3 <- distrib_grad_y_hess(d, y, theta)
#' vapply(g3, function(z) z[1], numeric(1))
#' h <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(h, c(0.2, 0.6))
S7::method(distrib_cross2_y, Lognormal1Distrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    lognormal_theta_chain(y, theta, 2L, FALSE)
  }

#' @rdname distrib_cross2_y.Lognormal1Distrib
#' @name distrib_grad_y_hess.Lognormal1Distrib
#' @keywords internal
S7::method(distrib_grad_y_hess, Lognormal1Distrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    lognormal_theta_chain(y, theta, 1L, TRUE)
  }

#' @rdname distrib_cross2_y.Lognormal1Distrib
#' @name distrib_hess_y_hess.Lognormal1Distrib
#' @keywords internal
S7::method(distrib_hess_y_hess, Lognormal1Distrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    lognormal_theta_chain(y, theta, 2L, TRUE)
  }
