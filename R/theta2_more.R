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

#' The Mixed Grid of a Family Written Out Against a Tabulated Map
#'
#' @description
#' `distrib_cross2_y`, `distrib_grad_y_hess` and
#' `distrib_hess_y_hess` for a family that carries its own class and
#' kernels while its map onto a parent is tabulated.
#'
#' @param parent The parent distribution.
#' @param th_par A function of the new parameters returning the parent's.
#' @param tables A function of the new parameters returning the map's keyed
#'   partial tables.
#'
#' @return A list of three method bodies, named `cross2`, `grad2`
#'   and `hess2`.
#'
#' @seealso [mapped_theta2()], [reparam_map_derivs()]
#' @keywords internal
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


#' The Scale Component of a Family With No Location
#'
#' @description
#' \eqn{\partial^2\ell^{(y)}/\partial\sigma^2} and
#' \eqn{\partial^2\ell^{(yy)}/\partial\sigma^2} for a family whose response
#' enters only through \eqn{z = y/\sigma}.
#'
#' @details
#' The derivation of the scale components of
#' [loc_scale_theta2_block()] never used the location, only that
#' \eqn{\sigma} is a scale, so the formulas hold unchanged with \eqn{z =
#' y/\sigma}. Any other parameter is a shape and is not covered.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order `1` for the block of [distrib_grad_y()],
#'   `2` for that of [distrib_hess_y()].
#' @param at The index of the scale parameter.
#'
#' @return A numeric vector.
#'
#' @seealso [loc_scale_theta2_block()]
#' @keywords internal
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


#' The Mixed Grid of a Family With No Location
#'
#' @description
#' The scale's own pair in closed form and every pair touching a shape by one
#' central difference of the analytic first-order component.
#'
#' @param at The index of the scale parameter.
#'
#' @return A list of two method bodies, named `grad2` and `hess2`.
#'
#' @seealso [scale_only_theta2()]
#' @keywords internal
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


#' The Scale Component of the Response Curvature's Derivative
#'
#' @description
#' [distrib_cross2_y()] for a family with no location: the same
#' identity as [loc_scale_cross2_block()] with \eqn{z = y/\sigma}.
#'
#' @param at The index of the scale parameter.
#'
#' @return A method body.
#'
#' @seealso [loc_scale_cross2_block()]
#' @keywords internal
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

#' The Mixed Grid of the Lognormal
#'
#' @description
#' The gaussian's components at \eqn{t = \log y}, carried by the Jacobian of
#' the transformation.
#'
#' @param y A numeric vector of observations.
#' @param theta A named list containing `mu` and `sigma2`.
#' @param order `1` for the block of [distrib_grad_y()],
#'   `2` for that of [distrib_hess_y()].
#' @param second Whether the second-order theta components are wanted.
#'
#' @return A named list, keyed by parameter or by parameter pair.
#'
#' @seealso [distrib_grad_y_hess()]
#' @keywords internal
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
#' @description
#' The gaussian's own components at \eqn{t = \log y}: the transformation
#' carries no parameter, so \eqn{t} does not move with \eqn{\theta} and only
#' the Jacobian of the response derivatives enters.
#' @param distrib A `Lognormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `sigma2`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list, keyed by parameter or by parameter pair.
#' @keywords internal
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
