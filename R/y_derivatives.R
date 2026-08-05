#' @include distrib.R generics.R
NULL

# Derivatives of the log-density with respect to the response y (not the
# parameters). Registered on continuous_distrib; distributions with a closed
# form override these. For discrete distributions the derivative w.r.t. y is not
# defined, so no default method is provided.

#' Finite-Difference Steps That Respect the Support
#'
#' @description
#' The step for a central difference in the response: scaled by \eqn{|y|} and
#' shrunk so that \eqn{y \pm h} stays strictly inside the distribution's support.
#'
#' @details
#' The response counterpart of \code{\link{fd_steps}}, and clamped for the same
#' reason: a Gamma observation close to zero would otherwise be differentiated
#' using a point outside the support, where the density is not defined.
#'
#' @param y A numeric vector of observations.
#' @param bounds A length-2 numeric vector, the distribution's support.
#' @param h_rel The relative step size.
#'
#' @return A numeric vector of steps, the same length as \code{y}.
#'
#' @seealso \code{\link{fd_steps}}, \code{\link{numerical_grad_y}}
#' @keywords internal
fd_steps_y <- function(y, bounds, h_rel) {
  h <- h_rel * pmax(1, abs(y))
  if (is.finite(bounds[1])) h <- pmin(h, 0.49 * (y - bounds[1]))
  if (is.finite(bounds[2])) h <- pmin(h, 0.49 * (bounds[2] - y))
  h
}

#' Numerical Gradient of the Log-Density with Respect to the Response
#'
#' @description
#' Computes \eqn{\partial \ell / \partial y} by central finite differences of
#' \code{distrib_pdf(..., log = TRUE)}. Powers the default \code{\link{distrib_grad_y}}
#' method for continuous distributions without a closed form.
#'
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param h_rel Numeric. Relative finite-difference step. Defaults to
#'   \code{.Machine$double.eps^(1/3)}.
#' @return A numeric vector of the same length as \code{y}.
#' @seealso \code{\link{numerical_hess_y}}, \code{\link{distrib_grad_y}}
#' @examples
#' numerical_grad_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_grad_y <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3)) {
  h <- fd_steps_y(y, distrib@bounds, h_rel)
  (distrib_pdf(distrib, y + h, theta, log = TRUE) -
    distrib_pdf(distrib, y - h, theta, log = TRUE)) / (2 * h)
}

#' Numerical Second Derivative of the Log-Density with Respect to the Response
#'
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2} by a central three-point stencil of
#' \code{distrib_pdf(..., log = TRUE)}. Powers the default \code{\link{distrib_hess_y}}
#' method for continuous distributions without a closed form.
#'
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param h_rel Numeric. Relative finite-difference step. Defaults to
#'   \code{.Machine$double.eps^(1/4)}.
#' @return A numeric vector of the same length as \code{y}.
#' @seealso \code{\link{numerical_grad_y}}, \code{\link{distrib_hess_y}}
#' @examples
#' numerical_hess_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_hess_y <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 4)) {
  h <- fd_steps_y(y, distrib@bounds, h_rel)
  (distrib_pdf(distrib, y + h, theta, log = TRUE) -
    2 * distrib_pdf(distrib, y, theta, log = TRUE) +
    distrib_pdf(distrib, y - h, theta, log = TRUE)) / (h^2)
}

#' @title Default Response Gradient for Continuous Distributions
#' @name distrib_grad_y.continuous_distrib
#' @description Fallback: \eqn{\partial \ell / \partial y} via finite differences (see \code{\link{numerical_grad_y}}).
#' @param distrib A \code{continuous_distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_grad_y, continuous_distrib) <- function(distrib, y, theta) {
  numerical_grad_y(distrib, y, theta)
}

#' @title Default Response Hessian for Continuous Distributions
#' @name distrib_hess_y.continuous_distrib
#' @description Fallback: \eqn{\partial^2 \ell / \partial y^2} via finite differences (see \code{\link{numerical_hess_y}}).
#' @param distrib A \code{continuous_distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_hess_y, continuous_distrib) <- function(distrib, y, theta) {
  numerical_hess_y(distrib, y, theta)
}
