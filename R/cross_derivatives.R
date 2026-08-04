#' @include distrib.R generics.R utility_functions.R y_derivatives.R gaussian_distrib.R student_t_distrib.R truncated.R
NULL

# Mixed second derivatives of the log-density: one derivative with respect to
# the response, one with respect to each parameter. Together with
# distrib_hess_y and distrib_hessian these complete the joint Hessian of a
# quantity that depends on both, which is what estimating the coefficients and
# the hyperparameters of a penalty together requires: the penalty is a negative
# log-density evaluated at the coefficients, and this block is its off-diagonal.
#
# Defined for continuous distributions, like the other response derivatives;
# a discrete distribution has no derivative in y to cross.

#' Mixed Response-Parameter Derivatives of the Log-Density
#'
#' @description
#' Computes the mixed second derivatives
#' \eqn{\partial^2 \ell / \partial y\, \partial \theta_i}, one component per
#' parameter, each a vector along \code{y}.
#'
#' @details
#' This is the off-diagonal block of the joint Hessian in \eqn{(y, \theta)},
#' whose diagonal blocks are \code{\link{distrib_hess_y}} and
#' \code{\link{distrib_hessian}}. It is what joint estimation over both
#' arguments needs, and what the gradient of a profiled objective needs
#' through the implicit function theorem.
#'
#' On the link scale the component for \eqn{\eta_i} is the parameter-scale
#' component multiplied by \eqn{h_i'(\eta_i)}: the response derivative is
#' untouched by a reparametrisation of \eqn{\theta}, so only the first-order
#' diagonal chain rule enters, exactly as for \code{\link{distrib_gradient}}.
#'
#' Distributions with a closed form provide it directly; the others fall back
#' to one central difference of \code{\link{distrib_grad_y}} in each
#' parameter (see \code{\link{numerical_cross_y}}).
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @inheritParams distrib_gradient
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   \code{distrib@params}.
#'
#' @examples
#' d <- gaussian_distrib()
#' distrib_cross_y(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#'
#' @export
distrib_cross_y <- S7::new_generic("distrib_cross_y", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") to_link_scale(distrib, theta, list(res), 1L) else res
})

#' Numerical Mixed Response-Parameter Derivatives
#'
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y\, \partial \theta_i} by one
#' central difference of \code{\link{distrib_grad_y}} in each parameter.
#' Powers the default \code{\link{distrib_cross_y}} method for continuous
#' distributions without a closed form.
#'
#' @details
#' The reference is the response gradient, not the log-density, so that a
#' distribution with an analytical \code{distrib_grad_y} pays for exactly one
#' finite-difference layer. When the response gradient is itself the
#' finite-difference fallback, the composition is the four-point mixed stencil
#' on the log-density -- the two differences act on different variables, so
#' they commute into a single stencil rather than compounding the way nested
#' differences in the same variable do.
#'
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param h_rel Numeric. Relative finite-difference step. Defaults to
#'   \code{.Machine$double.eps^(1/3)}.
#'
#' @return A named list with one numeric vector per parameter.
#'
#' @seealso \code{\link{numerical_grad_y}}, \code{\link{distrib_cross_y}}
#' @examples
#' numerical_cross_y(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_cross_y <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3)) {
  params <- distrib@params
  out <- vector("list", length(params))
  names(out) <- params
  for (j in seq_along(params)) {
    p <- params[j]
    h <- fd_steps(theta[[j]], distrib@params_bounds[[p]], h_rel)
    th_up <- theta
    th_dn <- theta
    th_up[[j]] <- theta[[j]] + h
    th_dn[[j]] <- theta[[j]] - h
    out[[p]] <- (distrib_grad_y(distrib, y, th_up) -
      distrib_grad_y(distrib, y, th_dn)) / (2 * h)
  }
  out
}

#' @title Default Mixed Derivatives for Continuous Distributions
#' @name distrib_cross_y.continuous_distrib
#' @description Fallback: one central difference of \code{\link{distrib_grad_y}}
#'   in each parameter (see \code{\link{numerical_cross_y}}).
#' @param distrib A \code{continuous_distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, continuous_distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  numerical_cross_y(distrib, y, theta)
}

#' @title Gaussian Mixed Derivatives
#' @name distrib_cross_y.GaussianDistrib
#' @description Closed form: with \eqn{r = y - \mu},
#'   \eqn{\partial^2 \ell / \partial y\, \partial \mu = 1/\sigma^2} and
#'   \eqn{\partial^2 \ell / \partial y\, \partial \sigma = 2r/\sigma^3}.
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with components \code{mu} and \code{sigma}.
#' @keywords internal
S7::method(distrib_cross_y, GaussianDistrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  r <- y - mu
  list(
    mu = rep(1 / sigma^2, length.out = length(y)),
    sigma = 2 * r / sigma^3
  )
}

#' @title Student's t Mixed Derivatives
#' @name distrib_cross_y.StudentTDistrib
#' @description Closed form: with \eqn{r = y - \mu} and
#'   \eqn{D = \nu\sigma^2 + r^2},
#'   \eqn{\partial^2 \ell / \partial y\, \partial \mu = (\nu+1)(\nu\sigma^2 - r^2)/D^2},
#'   \eqn{\partial^2 \ell / \partial y\, \partial \sigma = 2\nu\sigma(\nu+1)\, r/D^2},
#'   \eqn{\partial^2 \ell / \partial y\, \partial \nu = -r\,(r^2 - \sigma^2)/D^2}.
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with components \code{mu}, \code{sigma} and \code{nu}.
#' @keywords internal
S7::method(distrib_cross_y, StudentTDistrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  nu <- theta[[3]]
  r <- y - mu
  d <- nu * sigma^2 + r^2
  list(
    mu = (nu + 1) * (nu * sigma^2 - r^2) / d^2,
    sigma = 2 * nu * sigma * (nu + 1) * r / d^2,
    nu = -r * (r^2 - sigma^2) / d^2
  )
}

#' @title Mixed Derivatives of a Truncated Distribution
#' @name distrib_cross_y.TruncatedContinuousDistrib
#' @description The parent's mixed derivatives, unchanged: the truncated
#'   log-density is the parent's minus \eqn{\log Z(\theta)}, and the
#'   normalising constant does not depend on \eqn{y}, so it vanishes from any
#'   derivative that involves the response.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross_y, TruncatedContinuousDistrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    ...) {
  distrib_cross_y(distrib@parent_distrib, y, theta, scale = "parameter", ...)
}
