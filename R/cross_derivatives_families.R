#' @include cross_derivatives.R reparametrize.R reparam_maps.R
#' @include logistic_distrib.R cauchy_distrib.R laplace_distrib.R
#' @include gumbel_distrib.R pseudohuber_distrib.R skewnormal1_distrib.R
#' @include skewt_distrib.R exponential_distrib.R weibull1_distrib.R
#' @include gpd_distrib.R gengamma1_distrib.R invgauss2_distrib.R
#' @include vonmises2_distrib.R
NULL

# ===========================================================================
# Mixed response-parameter derivatives in closed form.
#
# One identity does most of the work. If the response enters the log-density
# only through z = (y - mu)/sigma, with mu a pure location and sigma a pure
# scale, then differentiating l_y = g'(z)/sigma gives
#
#   d^2 l / dy dmu    = -l_yy
#   d^2 l / dy dsigma = -z l_yy - l_y/sigma
#
# so both components come from the response derivatives the family already
# provides, and no new algebra is needed. The same sigma formula holds when
# there is no location at all and z = y/sigma, since it only uses that sigma
# is a scale: that covers the exponential, the Weibull's scale and the
# generalized Pareto's. A shape parameter is not covered and is written out
# per family, or differenced where it has no elementary form.
# ===========================================================================

#' The Location and Scale Components of a Mixed Derivative
#'
#' @description
#' \eqn{\partial^2\ell/\partial y\partial\mu = -\ell^{(yy)}} and
#' \eqn{\partial^2\ell/\partial y\partial\sigma = -z\ell^{(yy)} -
#' \ell^{(y)}/\sigma}, for a family whose response enters only through
#' \eqn{z = (y-\mu)/\sigma}.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, location first and scale second.
#'
#' @return A list of two component vectors, unnamed.
#'
#' @seealso \code{\link{distrib_cross_y}}
#' @keywords internal
loc_scale_cross_block <- function(distrib, y, theta) {
  s <- theta[[2]]
  z <- (y - theta[[1]]) / s
  ly <- distrib_grad_y(distrib, y, theta) + 0 * y
  lyy <- distrib_hess_y(distrib, y, theta) + 0 * y
  list(-lyy, -z * lyy - ly / s)
}

#' Mixed Derivatives of a Location-Scale Family
#'
#' @description
#' The \code{\link{distrib_cross_y}} body shared by the families that are
#' location-scale in both their parameters.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per parameter.
#'
#' @seealso \code{\link{loc_scale_cross_block}}
#' @keywords internal
loc_scale_cross_y <- function(distrib, y, theta,
                              scale = c("parameter", "link"), ...) {
  stats::setNames(loc_scale_cross_block(distrib, y, theta), distrib@params)
}

#' Mixed Derivatives When Only Two Parameters Are Location-Scale
#'
#' @description
#' The location and scale components in closed form and the remaining shape
#' components by one central difference of \code{\link{distrib_grad_y}}.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per parameter.
#'
#' @seealso \code{\link{loc_scale_cross_block}}
#' @keywords internal
partial_loc_scale_cross_y <- function(distrib, y, theta,
                                      scale = c("parameter", "link"), ...) {
  params <- distrib@params
  rest <- params[-(1:2)]
  out <- c(loc_scale_cross_block(distrib, y, theta),
           numerical_cross_y(distrib, y, theta, which = rest))
  stats::setNames(out, params)
}


# --- the location-scale families -------------------------------------------

#' @title Logistic Mixed Derivatives
#' @name distrib_cross_y.LogisticDistrib
#' @description
#' Closed form from the location-scale identity: the location component is
#' \eqn{-\ell^{(yy)}} and the scale one \eqn{-z\ell^{(yy)} - \ell^{(y)}/\sigma}.
#' @param distrib A \code{LogisticDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, LogisticDistrib) <- loc_scale_cross_y

#' @rdname distrib_cross_y.LogisticDistrib
#' @name distrib_cross_y.CauchyDistrib
#' @keywords internal
S7::method(distrib_cross_y, CauchyDistrib) <- loc_scale_cross_y

#' @rdname distrib_cross_y.LogisticDistrib
#' @name distrib_cross_y.GumbelDistrib
#' @keywords internal
S7::method(distrib_cross_y, GumbelDistrib) <- loc_scale_cross_y

#' @title Laplace Mixed Derivatives
#' @name distrib_cross_y.LaplaceDistrib
#' @description
#' Closed form from the location-scale identity, away from the kink at
#' \eqn{y = \mu}, where the density has no derivative in the response and the
#' quantity does not exist.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, LaplaceDistrib) <- loc_scale_cross_y

#' @title Pseudo-Huber Mixed Derivatives
#' @name distrib_cross_y.PseudoHuberDistrib
#' @description
#' Closed form in the location and the scale; the shape is differenced.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{sigma} and \code{nu}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, PseudoHuberDistrib) <- partial_loc_scale_cross_y

#' @rdname distrib_cross_y.PseudoHuberDistrib
#' @name distrib_cross_y.SkewNormal1Distrib
#' @keywords internal
S7::method(distrib_cross_y, SkewNormal1Distrib) <- partial_loc_scale_cross_y

#' @rdname distrib_cross_y.PseudoHuberDistrib
#' @name distrib_cross_y.SkewTDistrib
#' @keywords internal
S7::method(distrib_cross_y, SkewTDistrib) <- partial_loc_scale_cross_y


# --- families with a pure scale and, sometimes, a shape --------------------

#' @title Exponential Mixed Derivatives
#' @name distrib_cross_y.ExponentialDistrib
#' @description
#' Closed form: the mean is a pure scale, so the identity
#' \eqn{-z\ell^{(yy)} - \ell^{(y)}/\mu} applies with \eqn{z = y/\mu} and
#' returns \eqn{1/\mu^2}.
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector.
#' @keywords internal
S7::method(distrib_cross_y, ExponentialDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  list(mu = rep_len(1 / theta[[1]]^2, length(y)))
}

#' @title Weibull Mixed Derivatives
#' @name distrib_cross_y.Weibull1Distrib
#' @description
#' Closed form at both parameters. The scale follows the pure-scale identity
#' and reduces to \eqn{\sigma^2 w/(\mu y)} with \eqn{w = (y/\mu)^\sigma}; the
#' shape comes from differentiating
#' \eqn{\ell^{(y)} = ((\sigma-1) - \sigma w)/y}, giving
#' \eqn{(1 - w - \sigma w L)/y} with \eqn{L = \log(y/\mu)}.
#' @param distrib A \code{Weibull1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Weibull1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ...) {
  mu <- theta[[1]]
  sg <- theta[[2]]
  r <- y / mu
  w <- r^sg
  L <- base::log(r)
  list(mu = sg^2 * w / (mu * y),
       sigma = (1 - w - sg * w * L) / y)
}

#' @title Generalized Pareto Mixed Derivatives
#' @name distrib_cross_y.GPDDistrib
#' @description
#' Closed form at both parameters, from
#' \eqn{\ell^{(y)} = -(\xi+1)/(\sigma t)} with \eqn{t = 1 + \xi y/\sigma}:
#' the scale gives \eqn{(\xi+1)(t - \xi z)/(\sigma^2 t^2)} and the shape
#' \eqn{-1/(\sigma t) + (\xi+1)y/(\sigma^2 t^2)}. Neither carries a
#' \eqn{1/\xi}, so the shape direction needs no series at the exponential
#' limit.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, GPDDistrib) <- function(distrib, y, theta,
                                                    scale = c("parameter", "link"),
                                                    ...) {
  sg <- theta[[1]]
  xi <- theta[[2]]
  z <- y / sg
  t <- 1 + xi * z
  list(sigma = (xi + 1) * (t - xi * z) / (sg^2 * t^2) + 0 * y,
       xi = -1 / (sg * t) + (xi + 1) * y / (sg^2 * t^2))
}


# --- the response derivatives two families were missing --------------------

#' @title Generalized Pareto Response Derivatives
#' @name distrib_grad_y.GPDDistrib
#' @description
#' Closed form: \eqn{\ell^{(y)} = -(\xi+1)/(\sigma t)} and
#' \eqn{\ell^{(yy)} = \xi(\xi+1)/(\sigma t)^2}, with
#' \eqn{t = 1 + \xi y/\sigma}.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_grad_y, GPDDistrib) <- function(distrib, y, theta, ...) {
  sg <- theta[[1]]
  xi <- theta[[2]]
  -(xi + 1) / (sg * (1 + xi * y / sg)) + 0 * y
}

#' @rdname distrib_grad_y.GPDDistrib
#' @name distrib_hess_y.GPDDistrib
#' @keywords internal
S7::method(distrib_hess_y, GPDDistrib) <- function(distrib, y, theta, ...) {
  sg <- theta[[1]]
  xi <- theta[[2]]
  t <- 1 + xi * y / sg
  xi * (xi + 1) / (sg * t)^2 + 0 * y
}

#' @title Generalized Gamma Response Derivatives
#' @name distrib_grad_y.GenGamma1Distrib
#' @description
#' Closed form. With \eqn{w = (y/a)^p},
#' \eqn{\ell^{(y)} = ((d-1) - pw)/y} and
#' \eqn{\ell^{(yy)} = (pw(1-p) - (d-1))/y^2}.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_grad_y, GenGamma1Distrib) <- function(distrib, y, theta, ...) {
  w <- (y / theta[[1]])^theta[[3]]
  ((theta[[2]] - 1) - theta[[3]] * w) / y
}

#' @rdname distrib_grad_y.GenGamma1Distrib
#' @name distrib_hess_y.GenGamma1Distrib
#' @keywords internal
S7::method(distrib_hess_y, GenGamma1Distrib) <- function(distrib, y, theta, ...) {
  p <- theta[[3]]
  w <- (y / theta[[1]])^p
  (p * w * (1 - p) - (theta[[2]] - 1)) / y^2
}


# --- families written as a map of another ----------------------------------
#
# A derivative in the response does not see the parameters, so grad_y and
# hess_y are the parent's at the mapped parameter and nothing is chained. The
# mixed block is one first-order chain, the response derivative of the score:
# d^2 l / dy dpsi_i = sum_k (d^2 l / dy dtheta_k) h^k_i.

#' Mixed Derivatives Through a Map
#'
#' @description
#' The parent's mixed block carried onto a new parametrization by the
#' first-order chain rule, which is all that is needed: the response
#' derivative does not interact with a reparametrization of the parameters.
#'
#' @param distrib The distribution in the new parametrization.
#' @param parent The parent distribution.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param y A numeric vector of observations.
#'
#' @return A named list with one numeric vector per new parameter.
#'
#' @seealso \code{\link{distrib_cross_y}}, \code{\link{reparam_tables}}
#' @keywords internal
mapped_cross_y <- function(distrib, parent, th_par, maps, y) {
  cy <- distrib_cross_y(parent, y, th_par)
  new_params <- distrib@params
  zero <- 0 * cy[[1L]]
  out <- lapply(seq_along(new_params), function(i) {
    s <- zero
    for (k in seq_along(cy)) {
      v <- maps[[k]][[as.character(i)]]
      if (!is.null(v)) s <- s + cy[[k]] * v
    }
    s
  })
  stats::setNames(out, new_params)
}

#' @title Mixed Derivatives of a Reparametrized Distribution
#' @name distrib_cross_y.ReparamContinuousDistrib
#' @description
#' The parent's mixed block carried by the first-order chain rule on the map.
#' @param distrib A reparametrized distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the new parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, ReparamContinuousDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    mapped_cross_y(distrib, distrib@parent_distrib,
                   reparam_theta(distrib, theta),
                   reparam_tables(distrib, theta), y)
  }

#' @title Inverse Gaussian Response and Mixed Derivatives in Mean and Shape
#' @name distrib_grad_y.InvGauss2Distrib
#' @description
#' The dispersion parametrization's, at \eqn{\phi = 1/\lambda}: a derivative
#' in the response does not see the parameters, so the response derivatives
#' are the parent's unchanged, and the mixed block is one first-order chain.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{lambda}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_grad_y, InvGauss2Distrib) <- function(distrib, y, theta, ...) {
  distrib_grad_y(invgauss1_distrib(), y,
                 list(mu = theta[[1]], phi = 1 / theta[[2]]))
}

#' @rdname distrib_grad_y.InvGauss2Distrib
#' @name distrib_hess_y.InvGauss2Distrib
#' @keywords internal
S7::method(distrib_hess_y, InvGauss2Distrib) <- function(distrib, y, theta, ...) {
  distrib_hess_y(invgauss1_distrib(), y,
                 list(mu = theta[[1]], phi = 1 / theta[[2]]))
}

#' @rdname distrib_grad_y.InvGauss2Distrib
#' @name distrib_cross_y.InvGauss2Distrib
#' @param scale Handled by the generic before dispatch.
#' @keywords internal
S7::method(distrib_cross_y, InvGauss2Distrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    mapped_cross_y(distrib, invgauss1_distrib(),
                   list(mu = theta[[1]], phi = 1 / theta[[2]]),
                   md_invgauss2(theta), y)
  }

#' @title von Mises Response and Mixed Derivatives in the Resultant Length
#' @name distrib_grad_y.VonMises2Distrib
#' @description
#' The concentration parametrization's, at \eqn{\kappa = A^{-1}(\rho)}. The
#' response derivatives are the parent's unchanged; the mixed block is the
#' parent's carried by \eqn{\kappa'(\rho)}, the location passing through.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{rho}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_grad_y, VonMises2Distrib) <- function(distrib, y, theta, ...) {
  p <- vm2_parts(theta)
  distrib_grad_y(vonmises1_distrib(), y,
                 list(mu = theta[[1]], kappa = p$kappa))
}

#' @rdname distrib_grad_y.VonMises2Distrib
#' @name distrib_hess_y.VonMises2Distrib
#' @keywords internal
S7::method(distrib_hess_y, VonMises2Distrib) <- function(distrib, y, theta, ...) {
  p <- vm2_parts(theta)
  distrib_hess_y(vonmises1_distrib(), y,
                 list(mu = theta[[1]], kappa = p$kappa))
}

#' @rdname distrib_grad_y.VonMises2Distrib
#' @name distrib_cross_y.VonMises2Distrib
#' @param scale Handled by the generic before dispatch.
#' @keywords internal
S7::method(distrib_cross_y, VonMises2Distrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    p <- vm2_parts(theta)
    one <- rep_len(1, max(lengths(theta[1:2])))
    mapped_cross_y(distrib, vonmises1_distrib(),
                   list(mu = theta[[1]], kappa = p$kappa),
                   list(list("1" = one), list("2" = p$kd$d1)), y)
  }
