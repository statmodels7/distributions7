#' @include cross_derivatives_families.R reparam_maps.R
#' @include vonmises1_distrib.R gengamma1_distrib.R skewnormal2_distrib.R
#' @include gaussian2_distrib.R gaussian3_distrib.R lognormal1_distrib.R
#' @include invgauss1_distrib.R beta1_distrib.R beta2_distrib.R
#' @include chisq_distrib.R gamma1_distrib.R gamma2_distrib.R
NULL

# The families whose response derivative is a short expression in the
# parameters. For these the mixed block is read straight off l_y: each
# parameter enters it through one or two coefficients, so differentiating is
# a line apiece and no identity is needed.

#' @title von Mises Mixed Derivatives
#' @name distrib_cross_y.VonMises1Distrib
#' @description
#' Closed form. With \eqn{\ell^{(y)} = -\kappa\sin(y-\mu)}, the location
#' component is \eqn{\kappa\cos(y-\mu)} and the concentration one
#' \eqn{-\sin(y-\mu)}.
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, VonMises1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
  d <- y - theta[[1]]
  list(mu = theta[[2]] * cos(d), kappa = -sin(d))
}

#' @title Generalized Gamma Mixed Derivatives
#' @name distrib_cross_y.GenGamma1Distrib
#' @description
#' Closed form. With \eqn{w = (y/a)^p} and \eqn{L = \log(y/a)},
#' \eqn{\ell^{(y)} = ((d-1) - pw)/y}, so the components are \eqn{p^2w/(ay)},
#' \eqn{1/y} and \eqn{-w(1 + pL)/y}.
#' @param distrib A `GenGamma1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `a`, `d` and `p`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
  a <- theta[[1]]
  pp <- theta[[3]]
  r <- y / a
  w <- r^pp
  list(a = pp^2 * w / (a * y),
       d = 1 / y + 0 * w,
       p = -w * (1 + pp * base::log(r)) / y)
}

#' @title Gaussian Mixed Derivatives in Mean and Variance
#' @name distrib_cross_y.Gaussian2Distrib
#' @description
#' Closed form: \eqn{\ell^{(y)} = -(y-\mu)/\sigma^2} gives \eqn{1/\sigma^2}
#' and \eqn{(y-\mu)/\sigma^4}.
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `sigma2`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Gaussian2Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
  v <- theta[[2]]
  list(mu = rep_len(1 / v, length(y)), sigma2 = (y - theta[[1]]) / v^2)
}

#' @title Gaussian Mixed Derivatives in Mean and Precision
#' @name distrib_cross_y.Gaussian3Distrib
#' @description
#' Closed form: \eqn{\ell^{(y)} = -\tau(y-\mu)} gives \eqn{\tau} and
#' \eqn{-(y-\mu)}.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `tau`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Gaussian3Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
  list(mu = rep_len(theta[[2]], length(y)), tau = -(y - theta[[1]]))
}

#' @title Lognormal Mixed Derivatives
#' @name distrib_cross_y.Lognormal1Distrib
#' @description
#' Closed form. On the log scale the family is location-scale, and
#' \eqn{\ell^{(y)} = -1/y - (\log y - \mu)/(\sigma^2 y)} gives
#' \eqn{1/(\sigma^2 y)} and \eqn{(\log y - \mu)/(\sigma^4 y)}.
#' @param distrib A `Lognormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `sigma2`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Lognormal1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"),
                                                           ...) {
  v <- theta[[2]]
  r <- base::log(y) - theta[[1]]
  list(mu = 1 / (v * y), sigma2 = r / (v^2 * y))
}

#' @title Inverse Gaussian Mixed Derivatives
#' @name distrib_cross_y.InvGauss1Distrib
#' @description
#' Closed form. With
#' \eqn{\ell^{(y)} = -3/(2y) - 1/(2\phi\mu^2) + 1/(2\phi y^2)}, the components
#' are \eqn{1/(\phi\mu^3)} and \eqn{1/(2\phi^2\mu^2) - 1/(2\phi^2 y^2)}.
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `phi`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, InvGauss1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
  mu <- theta[[1]]
  ph <- theta[[2]]
  list(mu = rep_len(1 / (ph * mu^3), length(y)),
       phi = 1 / (2 * ph^2 * mu^2) - 1 / (2 * ph^2 * y^2))
}

#' @title Beta Mixed Derivatives in Mean and Precision
#' @name distrib_cross_y.Beta1Distrib
#' @description
#' Closed form. The shapes are \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi},
#' and \eqn{\ell^{(y)} = (a-1)/y - (b-1)/(1-y)}, so the components are
#' \eqn{\phi/y + \phi/(1-y)} and \eqn{\mu/y - (1-\mu)/(1-y)}.
#' @param distrib A `Beta1Distrib` object.
#' @param y A numeric vector in the unit interval.
#' @param theta A list containing `mu` and `phi`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Beta1Distrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
  mu <- theta[[1]]
  ph <- theta[[2]]
  list(mu = ph / y + ph / (1 - y),
       phi = mu / y - (1 - mu) / (1 - y))
}

#' @title Beta Mixed Derivatives in the Shapes
#' @name distrib_cross_y.Beta2Distrib
#' @description
#' Closed form: \eqn{\ell^{(y)} = (\alpha-1)/y - (\beta-1)/(1-y)} gives
#' \eqn{1/y} and \eqn{-1/(1-y)}.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector in the unit interval.
#' @param theta A list containing `alpha` and `beta`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Beta2Distrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
  list(alpha = 1 / y, beta = -1 / (1 - y))
}

#' @title Chi-Squared Mixed Derivative
#' @name distrib_cross_y.ChisqDistrib
#' @description
#' Closed form: the degrees of freedom are the mean, and
#' \eqn{\ell^{(y)} = (\mu/2 - 1)/y - 1/2} gives \eqn{1/(2y)}.
#' @param distrib A `ChisqDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector.
#' @keywords internal
S7::method(distrib_cross_y, ChisqDistrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
  list(mu = 1 / (2 * y))
}

#' @title Gamma Mixed Derivatives in Mean and Dispersion
#' @name distrib_cross_y.Gamma1Distrib
#' @description
#' Closed form. The shape is \eqn{1/\phi} and the rate \eqn{1/(\phi\mu)}, so
#' the mean component is \eqn{1/(\phi\mu^2)} and the dispersion one
#' \eqn{1/(\phi^2\mu) - 1/(\phi^2 y)}.
#' @param distrib A `Gamma1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `phi`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Gamma1Distrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"),
                                                       ...) {
  mu <- theta[[1]]
  ph <- theta[[2]]
  list(mu = rep_len(1 / (ph * mu^2), length(y)),
       phi = 1 / (ph^2 * mu) - 1 / (ph^2 * y))
}

#' @title Gamma Mixed Derivatives in Mean and Variance
#' @name distrib_cross_y.Gamma2Distrib
#' @description
#' Closed form: \eqn{\ell^{(y)} = (\mu^2/\sigma^2 - 1)/y - \mu/\sigma^2} gives
#' \eqn{2\mu/(\sigma^2 y) - 1/\sigma^2} and
#' \eqn{-\mu^2/(\sigma^4 y) + \mu/\sigma^4}.
#' @param distrib A `Gamma2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `sigma2`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, Gamma2Distrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"),
                                                       ...) {
  mu <- theta[[1]]
  v <- theta[[2]]
  list(mu = 2 * mu / (v * y) - 1 / v,
       sigma2 = -mu^2 / (v^2 * y) + mu / v^2)
}

#' @title Skew Normal Mixed Derivatives in the Centered Parametrization
#' @name distrib_cross_y.SkewNormal2Distrib
#' @description
#' The direct parametrization's mixed block carried by the first-order chain
#' rule on the centered-to-direct map.
#' @param distrib A `SkewNormal2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `gamma1`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  mapped_cross_y(distrib, skewnormal1_distrib(), sn2_theta(theta),
                 md_skewnormal2(theta[1:3]), y)
}
