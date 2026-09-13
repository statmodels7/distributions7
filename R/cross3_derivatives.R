#' @include distrib.R generics.R utility_functions.R y_derivatives.R y_higher.R cross2_derivatives.R
NULL

# Three derivatives with respect to the response, one with respect to each
# parameter. A penalty is a negative log-density evaluated at the coefficients,
# so its Hessian in the coefficients carries the density's second response
# derivative and the movement of that Hessian with the coefficients carries the
# third. The second derivative of a marginal criterion in the hyperparameters
# asks how THAT movement changes when a hyperparameter does, which is this
# quantity.
#
# Continuous distributions only, like the other response derivatives.

#' @title Mixed Third-Response Parameter Derivatives of the Log-Density
#'
#' @description
#' Computes \eqn{\partial^4 \ell / \partial y^3\, \partial \theta_i}, one
#' component per parameter, each a vector along `y`. It is
#' [distrib_deriv3_y()] differentiated once in each parameter.
#'
#' @details
#' # What consumes it
#'
#' A penalty \eqn{-\log f(D\beta;\theta)} has the Hessian
#' \eqn{S = -D'\mathrm{diag}(\ell^{(yy)})D} in the coefficients, whose
#' derivative along a direction \eqn{v} is
#' \eqn{-D'\mathrm{diag}(\ell^{(yyy)} \odot Dv)D}. The second derivative of a
#' marginal criterion in the hyperparameters reads how that matrix moves with
#' each hyperparameter, \eqn{-D'\mathrm{diag}(\partial_{\theta_i}\ell^{(yyy)}
#' \odot Dv)D}, and this is the vector placed on its diagonal.
#'
#' # Where the numbers come from
#'
#' For a family whose response enters only as \eqn{y - \mu} the response
#' derivative of order three is \eqn{-\partial^3\ell/\partial\mu^3}, so
#' \deqn{\frac{\partial^4\ell}{\partial y^3\,\partial\theta_i} =
#'   -\frac{\partial^4\ell}{\partial\mu^3\,\partial\theta_i},}
#' a component of [distrib_deriv4()] with a sign. The fourteen location
#' families [distrib_deriv3_y()] serves by the same identity are served this
#' way. Every other continuous family takes one central difference of
#' [distrib_deriv3_y()] in each parameter, through [numerical_cross3_y()].
#'
#' # The link scale
#'
#' The component for \eqn{\eta_i} is the parameter-scale component multiplied
#' by \eqn{h_i'(\eta_i)}: the response derivatives are untouched by a
#' reparametrization of \eqn{\theta}, so only the first-order diagonal chain
#' rule enters, as for [distrib_cross_y()] and [distrib_cross2_y()].
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
#' to the other. \eqn{\ell^{(yyy)}} is \eqn{\partial^3\ell/\partial y^3}.
#'
#' @seealso [distrib_deriv3_y()] for the quantity being differentiated,
#'   [distrib_cross2_y()] for the order below, and [numerical_cross3_y()] for
#'   the fallback.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3, nu = 6)
#' distrib_cross3_y(d, y, theta)
#'
#' # Against a numerical derivative of the analytic third response derivative.
#' f <- function(v) distrib_deriv3_y(d, y[1], list(mu = 0.4, sigma = v, nu = 6))
#' numDeriv::grad(f, 1.3)
#'
#' # A family whose response is not a location takes the difference.
#' distrib_cross3_y(gamma1_distrib(), c(0.5, 1, 2), list(mu = 1.5, phi = 0.8))
#'
#' @export
distrib_cross3_y <- S7::new_generic("distrib_cross3_y", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") to_link_scale(distrib, theta, list(res), 1L) else res
})


#' @title Numerical Mixed Third-Response Parameter Derivatives
#'
#' @description
#' Computes \eqn{\partial^4 \ell / \partial y^3\, \partial \theta_i} by one
#' central difference of [distrib_deriv3_y()] in each parameter. It is what the
#' default [distrib_cross3_y()] method runs for a continuous family that is not
#' a location family.
#'
#' @details
#' A family carrying an analytic `distrib_deriv3_y` pays for exactly one
#' difference. Where that derivative is itself the base class's stencil on the
#' log-density, the two differences act on different variables and compose
#' into one mixed stencil. The step is chosen by [fd_stable_step()], so a
#' parameter near a bound is differenced on the side the bound allows.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, the point to differentiate at.
#' @param h_rel The relative step, defaulting to `.Machine$double.eps^(1/3)`,
#'   the optimal exponent for a central first difference.
#' @param which An optional character vector naming a subset of
#'   `distrib@params`. The default, `NULL`, computes every component.
#'
#' @return A named list with one numeric vector per requested parameter, each
#'   as long as `y`.
#'
#' @seealso [distrib_cross3_y()], the generic it serves, and
#'   [distrib_deriv3_y()] for the quantity it differences.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3, nu = 6)
#'
#' # Against the identity the location family uses.
#' max(abs(unlist(numerical_cross3_y(d, y, theta)) -
#'         unlist(distrib_cross3_y(d, y, theta))))
#'
#' numerical_cross3_y(d, y, theta, which = "nu")
#'
#' @export
numerical_cross3_y <- function(distrib, y, theta,
                               h_rel = .Machine$double.eps^(1 / 3),
                               which = NULL) {
  params <- distrib@params
  keep <- if (is.null(which)) params else which
  out <- vector("list", length(keep))
  names(out) <- keep
  for (j in match(keep, params)) {
    p <- params[j]
    quotient <- function(h) {
      th_up <- th_dn <- theta
      th_up[[j]] <- theta[[j]] + h
      th_dn[[j]] <- theta[[j]] - h
      (distrib_deriv3_y(distrib, y, th_up) -
        distrib_deriv3_y(distrib, y, th_dn)) / (2 * h)
    }
    out[[p]] <- fd_stable_step(quotient, theta[[j]],
                               distrib@params_bounds[[p]], h_rel)$value
  }
  out
}


#' @title Default Mixed Third-Response Derivatives for Continuous Distributions
#' @name distrib_cross3_y.continuous_distrib
#'
#' @description
#' Falls back to one central difference of [distrib_deriv3_y()] in each
#' parameter, through [numerical_cross3_y()].
#'
#' @param distrib A `continuous_distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic after
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`.
#'
#' @seealso [distrib_cross3_y()] for the generic, [numerical_cross3_y()], which
#'   does the work, and [loc_cross3_y()] for the location families.
#'
#' @keywords internal
S7::method(distrib_cross3_y, continuous_distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"),
                                                             ...) {
  numerical_cross3_y(distrib, y, theta)
}


#' The Mixed Third-Response Derivative of a Location Family
#'
#' @description
#' Reads \eqn{\partial^4\ell/\partial y^3\,\partial\theta_i} as
#' \eqn{-\partial^4\ell/\partial\mu^3\,\partial\theta_i} from
#' [distrib_deriv4()], for a family whose response enters only as
#' \eqn{y - \mu}.
#'
#' @details
#' The identity is exact: each derivative in the response is minus the
#' derivative in the location, and a derivative in any parameter commutes with
#' both. The component for \eqn{\theta_i = \mu} is therefore
#' \eqn{-\partial^4\ell/\partial\mu^4}. The keys are generated from the
#' parameter indices, as [deriv_names()] generates them, rather than written.
#'
#' @param distrib A location-family distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic after
#'   dispatch.
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per parameter, each as long as
#'   `y`.
#'
#' @section Registered on:
#' The fourteen families [distrib_deriv3_y()] serves through
#' [loc_deriv_y_k()]: `Gaussian1Distrib`, `Gaussian2Distrib`,
#' `Gaussian3Distrib`, `CauchyDistrib`, `LogisticDistrib`, `LaplaceDistrib`,
#' `Laplace2Distrib`, `EnetDistrib`, `PseudoHuberDistrib`, `StudentT1Distrib`,
#' `SkewNormal1Distrib`, `SkewNormal2Distrib`, `SkewTDistrib`,
#' `GumbelDistrib`.
#'
#' @aliases distrib_cross3_y.Gaussian1Distrib distrib_cross3_y.Gaussian2Distrib
#' @aliases distrib_cross3_y.Gaussian3Distrib distrib_cross3_y.CauchyDistrib
#' @aliases distrib_cross3_y.LogisticDistrib distrib_cross3_y.LaplaceDistrib
#' @aliases distrib_cross3_y.Laplace2Distrib distrib_cross3_y.EnetDistrib
#' @aliases distrib_cross3_y.PseudoHuberDistrib distrib_cross3_y.StudentT1Distrib
#' @aliases distrib_cross3_y.SkewNormal1Distrib distrib_cross3_y.SkewNormal2Distrib
#' @aliases distrib_cross3_y.SkewTDistrib distrib_cross3_y.GumbelDistrib
#'
#' @seealso [distrib_cross3_y()], [loc_deriv_y_k()].
#'
#' @keywords internal
loc_cross3_y <- function(distrib, y, theta, scale = c("parameter", "link"),
                         ...) {
  params <- distrib@params
  im <- match("mu", params)
  d4 <- distrib_deriv4(distrib, y, theta)
  n <- length(y)
  out <- lapply(seq_along(params), function(j) {
    key <- paste(params[sort(c(im, im, im, j))], collapse = "_")
    rep_len(-d4[[key]], n)
  })
  names(out) <- params
  out
}

for (.cls in list(Gaussian1Distrib, Gaussian2Distrib, Gaussian3Distrib,
                  CauchyDistrib, LogisticDistrib, LaplaceDistrib,
                  Laplace2Distrib, EnetDistrib, PseudoHuberDistrib,
                  StudentT1Distrib, SkewNormal1Distrib,
                  SkewNormal2Distrib, SkewTDistrib, GumbelDistrib)) {
  S7::method(distrib_cross3_y, .cls) <- loc_cross3_y
}
rm(.cls)
