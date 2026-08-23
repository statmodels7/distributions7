#' @include distrib.R generics.R utility_functions.R y_derivatives.R cross_derivatives.R gaussian1_distrib.R student_t1_distrib.R
NULL

# Two derivatives with respect to the response, one with respect to each
# parameter. This is what a marginal criterion needs of a penalty: a penalty is
# a negative log-density evaluated at the coefficients, so its second
# derivative in the coefficients carries the density's second derivative in the
# response, and estimating the hyperparameters by that criterion asks how THAT
# moves when the hyperparameters do.
#
# Continuous distributions only, like the other response derivatives.

#' Mixed Second-Response Parameter Derivatives of the Log-Density
#'
#' @description
#' Computes \eqn{\partial^3 \ell / \partial y^2\, \partial \theta_i}, one
#' component per parameter, each a vector along `y`.
#'
#' @details
#' [distrib_hess_y()] says how curved the log-density is in the
#' response; this says how that curvature moves with each parameter. It
#' completes the surface a marginal likelihood needs: where a penalty is
#' \eqn{-\log f(D\beta;\theta)}, its Hessian in the coefficients is
#' \eqn{-D'\mathrm{diag}(\ell^{(yy)})D} and the derivative of that in
#' \eqn{\theta} is this component placed the same way.
#'
#' On the link scale the component for \eqn{\eta_i} is the parameter-scale
#' component multiplied by \eqn{h_i'(\eta_i)}: the response derivatives are
#' untouched by a reparametrization of \eqn{\theta}, so only the first-order
#' diagonal chain rule enters, exactly as for [distrib_cross_y()].
#'
#' Distributions with a closed form provide it directly; the others fall back
#' to one central difference of [distrib_hess_y()] in each parameter
#' (see [numerical_cross2_y()]). The reference is the analytic
#' response Hessian, so a family with one pays for exactly one difference.
#'
#' @param distrib A distribution object inheriting from the `distrib`
#'   class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution
#'   parameters. Each parameter must have length 1 or `length(y)`.
#' @inheritParams distrib_gradient
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' distrib_cross2_y(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#'
#' @seealso [distrib_hess_y()], [distrib_cross_y()]
#' @export
distrib_cross2_y <- S7::new_generic("distrib_cross2_y", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") to_link_scale(distrib, theta, list(res), 1L) else res
})


#' Numerical Mixed Second-Response Parameter Derivatives
#'
#' @description
#' Computes \eqn{\partial^3 \ell / \partial y^2\, \partial \theta_i} by one
#' central difference of [distrib_hess_y()] in each parameter.
#' Powers the default [distrib_cross2_y()] method for continuous
#' distributions without a closed form.
#'
#' @details
#' The reference is the response Hessian, not the log-density, so a
#' distribution with an analytical `distrib_hess_y` pays for exactly one
#' finite-difference layer. The two differences act on different variables and
#' therefore compose into a single mixed stencil rather than into the nested
#' differencing of one variable that the package forbids.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param h_rel The relative step.
#' @param which Optional subset of parameters.
#'
#' @return A named list with one numeric vector per parameter.
#'
#' @examples
#' numerical_cross2_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_cross2_y <- function(distrib, y, theta,
                               h_rel = .Machine$double.eps^(1 / 3),
                               which = NULL) {
  params <- distrib@params
  keep <- if (is.null(which)) params else which
  out <- vector("list", length(keep))
  names(out) <- keep
  for (j in match(keep, params)) {
    p <- params[j]
    h <- fd_steps(theta[[j]], distrib@params_bounds[[p]], h_rel)
    th_up <- theta
    th_dn <- theta
    th_up[[j]] <- theta[[j]] + h
    th_dn[[j]] <- theta[[j]] - h
    out[[p]] <- (distrib_hess_y(distrib, y, th_up) -
      distrib_hess_y(distrib, y, th_dn)) / (2 * h)
  }
  out
}


#' @title Default Mixed Second-Response Derivatives for Continuous Distributions
#' @name distrib_cross2_y.continuous_distrib
#' @description Fallback: one central difference of
#'   [distrib_hess_y()] in each parameter (see
#'   [numerical_cross2_y()]).
#' @param distrib A `continuous_distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross2_y, continuous_distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"),
                                                             ...) {
  numerical_cross2_y(distrib, y, theta)
}


#' @title Gaussian Mixed Second-Response Derivatives
#' @name distrib_cross2_y.Gaussian1Distrib
#' @description Closed form: \eqn{\ell^{(yy)} = -1/\sigma^2} does not depend on
#'   the location, so \eqn{\partial^3\ell/\partial y^2\partial\mu = 0} and
#'   \eqn{\partial^3\ell/\partial y^2\partial\sigma = 2/\sigma^3}.
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `sigma`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross2_y, Gaussian1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"),
                                                           ...) {
  sigma <- theta[[2]]
  n <- length(y)
  list(mu = rep(0, length.out = n),
       sigma = rep(2 / sigma^3, length.out = n))
}


#' @title Student's t Mixed Second-Response Derivatives
#' @name distrib_cross2_y.StudentT1Distrib
#' @description
#' Closed form. With \eqn{z = (y-\mu)/\sigma} and
#' \eqn{Q(z) = (\nu - z^2)/(\nu + z^2)^2}, the response Hessian is
#' \eqn{\ell^{(yy)} = -(\nu+1)Q/\sigma^2}, and with
#' \eqn{Q' = 2z(z^2 - 3\nu)/(\nu+z^2)^3} and
#' \eqn{\partial Q/\partial\nu = (3z^2-\nu)/(\nu+z^2)^3},
#' \deqn{\partial_\mu \ell^{(yy)} = (\nu+1)Q'/\sigma^3, \quad
#'   \partial_\sigma \ell^{(yy)} = (\nu+1)(2Q + zQ')/\sigma^3, \quad
#'   \partial_\nu \ell^{(yy)} = -(Q + (\nu+1)\partial_\nu Q)/\sigma^2.}
#' @param distrib A `StudentT1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `nu`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross2_y, StudentT1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"),
                                                           ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  nu <- theta[[3]]
  z <- (y - mu) / sigma
  d <- nu + z^2
  Q <- (nu - z^2) / d^2
  Qz <- 2 * z * (z^2 - 3 * nu) / d^3
  Qnu <- (3 * z^2 - nu) / d^3
  list(
    mu = (nu + 1) * Qz / sigma^3,
    sigma = (nu + 1) * (2 * Q + z * Qz) / sigma^3,
    nu = -(Q + (nu + 1) * Qnu) / sigma^2
  )
}
