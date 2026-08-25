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

#' @title Mixed Second-Response Parameter Derivatives of the Log-Density
#'
#' @description
#' Computes \eqn{\partial^3 \ell / \partial y^2\, \partial \theta_i}, one
#' component per parameter, each a vector along `y`. Where [distrib_hess_y()]
#' says how curved the log-density is in the RESPONSE, this says how that
#' curvature moves with each parameter.
#'
#' @details
#' # What consumes it
#'
#' It completes the surface a marginal likelihood needs. A penalty is
#' \eqn{-\log f(D\beta;\theta)}, so its Hessian in the coefficients is
#' \eqn{-D'\mathrm{diag}(\ell^{(yy)})D}, and the derivative of that in
#' \eqn{\theta} is this component placed the same way. Differentiating the
#' Laplace approximation once needs it; differentiating twice needs
#' [distrib_hess_y_hess()], which is this quantity one order further in the
#' parameters.
#'
#' # The link scale
#'
#' The component for \eqn{\eta_i} is the parameter-scale component multiplied
#' by \eqn{h_i'(\eta_i)}. The response derivatives are untouched by a
#' reparametrization of \eqn{\theta}, so only the first-order diagonal chain
#' rule enters, exactly as for [distrib_cross_y()].
#'
#' # Where the numbers come from
#'
#' A distribution with a closed form provides it. The rest take one central
#' difference of [distrib_hess_y()] in each parameter, through
#' [numerical_cross2_y()]. The reference is the analytic response Hessian, so a
#' family carrying one pays for exactly one difference.
#'
#' Continuous distributions only, as with every response derivative: a discrete
#' family has no method and the call raises with its class named.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list, or named numeric vector, of distribution
#'   parameters. Each must have length 1 or `length(y)`.
#' @inheritParams distrib_gradient
#' @param ... Passed to the method.
#'
#' @return A named list with one numeric vector per parameter, each as long as
#'   `y`, keyed by `distrib@params`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one
#' to the other. \eqn{\ell^{(yy)}} is \eqn{\partial^2\ell/\partial y^2}.
#'
#' @seealso [distrib_hess_y()] for the quantity being differentiated,
#'   [distrib_cross_y()] for the first-response counterpart,
#'   [distrib_hess_y_hess()] for the next order in \eqn{\theta}, and
#'   [numerical_cross2_y()] for the fallback.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' distrib_cross2_y(d, y, theta)
#'
#' # The gaussian's response curvature is -1 / sigma^2, which carries no
#' # location, so the mu component is exactly zero.
#' c(sigma = 2 / 1.3^3)
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::grad(f, c(0.4, 1.3))
#'
#' # On the link scale, sigma riding a log link, the component is multiplied
#' # by h' = sigma and nothing else.
#' c(link = distrib_cross2_y(d, y, theta, scale = "link")$sigma[1],
#'   hand = distrib_cross2_y(d, y, theta)$sigma[1] * 1.3)
#'
#' # A discrete family has no method at all.
#' try(distrib_cross2_y(poisson_distrib(), 1:3, list(mu = 2)))
#'
#' @export
distrib_cross2_y <- S7::new_generic("distrib_cross2_y", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") to_link_scale(distrib, theta, list(res), 1L) else res
})


#' @title Numerical Mixed Second-Response Parameter Derivatives
#'
#' @description
#' Computes \eqn{\partial^3 \ell / \partial y^2\, \partial \theta_i} by one
#' central difference of [distrib_hess_y()] in each parameter. It is what the
#' default [distrib_cross2_y()] method runs for a continuous family with no
#' closed form of its own.
#'
#' @details
#' The reference is the response HESSIAN, one order up from the log-density, so
#' a distribution carrying an analytic `distrib_hess_y` pays for exactly one
#' finite-difference layer. Where that Hessian is itself a fallback the two
#' differences act on different variables and compose into a single mixed
#' stencil, rather than into the nested differencing of one variable the
#' package forbids.
#'
#' The step is [fd_steps()]'s, so a parameter near a bound is differenced
#' inward instead of across it.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, the point to differentiate at.
#' @param h_rel The relative step, defaulting to `.Machine$double.eps^(1/3)`,
#'   the optimal exponent for a central first difference.
#' @param which An optional character vector naming a subset of
#'   `distrib@params`. The default, `NULL`, computes every component; a subset
#'   returns only those, keyed the same way, and costs one difference each.
#'
#' @return A named list with one numeric vector per requested parameter, each
#'   as long as `y`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one
#' to the other. \eqn{\ell^{(yy)}} is \eqn{\partial^2\ell/\partial y^2}.
#'
#' @seealso [distrib_cross2_y()], the generic it serves, [distrib_hess_y()] for
#'   the quantity it differences, and [fd_steps()] for the step rule.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#'
#' num <- numerical_cross2_y(d, y, theta)
#' num$sigma
#'
#' # Against the gaussian's own closed form, which it does not need this for.
#' max(abs(unlist(num) - unlist(distrib_cross2_y(d, y, theta))))
#'
#' # One parameter at a time, when only one is wanted.
#' numerical_cross2_y(d, y, theta, which = "sigma")
#'
#' # A family that does take this route, checked against numDeriv.
#' dg <- gamma2_distrib()
#' yg <- c(0.5, 1, 2)
#' vapply(distrib_cross2_y(dg, yg, list(mu = 2, sigma2 = 1)),
#'        function(z) z[1], numeric(1))
#' f <- function(v) distrib_hess_y(dg, yg[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::grad(f, c(2, 1))
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
#'
#' @description
#' Falls back to one central difference of [distrib_hess_y()] in each
#' parameter, through [numerical_cross2_y()]. Registering it on
#' `continuous_distrib` gives the quantity to every continuous family, whether
#' or not it writes one out.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`.
#'
#' @seealso [distrib_cross2_y()] for the generic, [numerical_cross2_y()], which
#'   does the work, and [distrib_cross2_y.Gaussian1Distrib()] for a family that
#'   overrides it.
#'
#' @keywords internal
#'
#' @examples
#' # The gamma writes no closed form, so this method answers for it.
#' d <- gamma2_distrib()
#' y <- c(0.5, 1, 2)
#' theta <- list(mu = 2, sigma2 = 1)
#' attr(S7::method(distrib_cross2_y, S7::S7_class(d)), "signature")[[1]]
#'
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::grad(f, c(2, 1))
S7::method(distrib_cross2_y, continuous_distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"),
                                                             ...) {
  numerical_cross2_y(distrib, y, theta)
}


#' @title Gaussian Mixed Second-Response Derivatives
#' @name distrib_cross2_y.Gaussian1Distrib
#'
#' @description
#' Closed form, and half of it zero. The gaussian's response curvature is
#' \eqn{\ell^{(yy)} = -1/\sigma^2}, which carries no location, so
#' \deqn{\frac{\partial^3\ell}{\partial y^2\,\partial\mu} = 0, \qquad
#'   \frac{\partial^3\ell}{\partial y^2\,\partial\sigma} = \frac{2}{\sigma^3}.}
#' Neither component varies with the data.
#'
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with `mu` and `sigma`.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors of length `length(y)`, keyed
#'   `mu` and `sigma`, the first exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one
#' to the other. \eqn{\ell^{(yy)}} is \eqn{\partial^2\ell/\partial y^2}.
#'
#' @seealso [distrib_cross2_y()] for the generic,
#'   [distrib_hess_y_hess.Gaussian1Distrib()] for the next order in
#'   \eqn{\theta}, and [distrib_cross2_y.StudentT1Distrib()], where nothing
#'   vanishes.
#'
#' @keywords internal
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' distrib_cross2_y(d, y, theta)
#'
#' # The one surviving component, written out.
#' 2 / 1.3^3
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::grad(f, c(0.4, 1.3))
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
#'
#' @description
#' Closed form in all three parameters. With \eqn{z = (y-\mu)/\sigma} and
#' \eqn{Q(z) = (\nu - z^2)/(\nu + z^2)^2} the response curvature is
#' \eqn{\ell^{(yy)} = -(\nu+1)Q/\sigma^2}, so writing
#' \eqn{Q' = 2z(z^2 - 3\nu)/(\nu+z^2)^3} and
#' \eqn{\partial Q/\partial\nu = (3z^2-\nu)/(\nu+z^2)^3},
#' \deqn{\frac{\partial\ell^{(yy)}}{\partial\mu} = \frac{(\nu+1)Q'}{\sigma^3},
#'   \qquad \frac{\partial\ell^{(yy)}}{\partial\sigma}
#'     = \frac{(\nu+1)(2Q + zQ')}{\sigma^3}, \qquad
#'   \frac{\partial\ell^{(yy)}}{\partial\nu}
#'     = -\frac{Q + (\nu+1)\partial_\nu Q}{\sigma^2}.}
#'
#' @details
#' Nothing vanishes here, unlike the gaussian's, and every component varies
#' with the data: the t's response curvature is a function of \eqn{z} rather
#' than a constant, and it is that dependence which makes the family
#' redescending. The three formulas share \eqn{Q}, \eqn{Q'} and
#' \eqn{\partial_\nu Q}, so each is evaluated once per observation.
#'
#' @param distrib A `StudentT1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with `mu`, `sigma` and `nu`.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors of length `length(y)`, keyed
#'   `mu`, `sigma` and `nu`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one
#' to the other. \eqn{\ell^{(yy)}} is \eqn{\partial^2\ell/\partial y^2}.
#'
#' @seealso [distrib_cross2_y()] for the generic, [distrib_hess_y()] for the
#'   quantity being differentiated, and [distrib_cross2_y.Gaussian1Distrib()],
#'   the limit as \eqn{\nu} grows.
#'
#' @keywords internal
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3, nu = 6)
#' distrib_cross2_y(d, y, theta)
#'
#' # Against a numerical derivative of the response Hessian.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2],
#'                                               nu = v[3]))
#' numDeriv::grad(f, c(0.4, 1.3, 6))
#'
#' # Every component varies with y, where the gaussian's do not.
#' distrib_cross2_y(d, y, theta)$mu
#' distrib_cross2_y(gaussian1_distrib(), y, theta[1:2])$mu
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
