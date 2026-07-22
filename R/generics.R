#' @title Distribution Generics
#' @description A collection of S7 generic functions for mathematical and statistical
#' operations on probability distributions.
#' @import S7
#' @name distrib_generics
NULL

# All generics normalize `theta` before dispatch (see align_theta() in
# utility_functions.R): a named theta is reordered to match distrib@params, so
# methods can safely access parameters by position. The derivative generics
# additionally validate parameter lengths against `y` and recycle a scalar `y`,
# since the Rcpp kernels iterate over length(y) and would otherwise silently
# truncate (or read out of bounds) with mismatched lengths.

#' Probability Density Function
#'
#' @description Evaluates the probability density function (PDF) or probability mass function (PMF).
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method (e.g., \code{log}).
#' @export
distrib_pdf <- S7::new_generic("distrib_pdf", "distrib", function(distrib, y, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Cumulative Distribution Function
#'
#' @description Evaluates the cumulative distribution function (CDF) for a given distribution.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param q A numeric vector of quantiles.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method (e.g., \code{lower.tail}, \code{log.p}).
#' @export
distrib_cdf <- S7::new_generic("distrib_cdf", "distrib", function(distrib, q, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Quantile Function
#'
#' @description Evaluates the quantile function for a given distribution.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param p A numeric vector of probabilities.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method (e.g., \code{lower.tail}, \code{log.p}).
#' @export
distrib_quantile <- S7::new_generic("distrib_quantile", "distrib", function(distrib, p, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Random Number Generator
#'
#' @description Generates random variates from the given distribution.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param n Number of observations to generate.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method.
#' @export
distrib_rng <- S7::new_generic("distrib_rng", "distrib", function(distrib, n, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

# Shared validation for the derivative generics: aligns theta, checks that all
# parameter lengths are 1 or n, and recycles a scalar y up to n when theta is
# vectorized. Returns list(y = ..., theta = ...).
check_derivative_args <- function(distrib, y, theta) {
  theta <- align_theta(distrib, theta)
  n <- max(length(y), lengths(theta[seq_len(distrib@n_params)]))
  check_params_dim(theta[seq_len(distrib@n_params)], n = n)
  if (length(y) != n) {
    if (length(y) != 1) {
      stop(sprintf("'y' must have length 1 or %d, not %d.", n, length(y)), call. = FALSE)
    }
    y <- rep(y, n)
  }
  list(y = y, theta = theta)
}

#' Analytical Gradient
#'
#' @description Computes the analytical first derivatives of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @param scale Either \code{"parameter"} (default) for derivatives with respect to
#'   the parameters \eqn{\theta} on their natural (constrained) scale, or
#'   \code{"link"} for derivatives with respect to the unconstrained linear
#'   predictors \eqn{\eta = g(\theta)} defined by \code{distrib@link_params}.
#'   See \code{\link{link_scale_derivatives}}.
#' @param ... Additional arguments passed to the specific method.
#' @export
distrib_gradient <- S7::new_generic("distrib_gradient", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") to_link_scale(distrib, theta, list(res), 1L) else res
})

#' Analytical Hessian
#'
#' @description Computes the analytical observed second derivatives (Hessian matrix) of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @inheritParams distrib_gradient
#' @param ... Additional arguments passed to the specific method.
#' @export
distrib_hessian <- S7::new_generic("distrib_hessian", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    nat1 <- distrib_gradient(distrib, y, theta)
    to_link_scale(distrib, theta, list(nat1, res), 2L)
  } else {
    res
  }
})

#' Analytical Expected Hessian
#'
#' @description Computes the analytical expected second derivatives of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @inheritParams distrib_gradient
#' @param approx How the expectation is approximated \strong{when the distribution
#'   has no closed-form expected Hessian}; ignored otherwise. One of
#'   \code{"bartlett"} (default, equivalently \code{"opg"}: the outer product of
#'   the score), \code{"integrate"} (quadrature/summation of the observed Hessian)
#'   or \code{"mc"} (Monte Carlo). See \code{\link{expected_derivative_methods}}
#'   for the accuracy/speed trade-offs.
#' @param nsim Monte Carlo sample size used when \code{approx = "mc"}. Defaults to 10000.
#' @param ... Additional arguments passed to the specific method.
#' @details
#' On the link scale the first-order term of the chain rule drops out because the
#' score has zero expectation, so the expected Hessian transforms as the simple
#' congruence \eqn{\mathrm{diag}(h')\, \mathbb{E}[H]\, \mathrm{diag}(h')} with
#' \eqn{h' = dg^{-1}/d\eta}.
#' @export
distrib_expected_hessian <- S7::new_generic("distrib_expected_hessian", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    zero1 <- stats::setNames(
      lapply(distrib@params, function(p) rep(0, length(y))), distrib@params
    )
    to_link_scale(distrib, theta, list(zero1, res), 2L)
  } else {
    res
  }
})

#' Analytical Third-Order Derivatives
#'
#' @description
#' Computes the unique third-order partial derivatives of the log-likelihood with
#' respect to the distribution's parameters. Distributions with a closed-form
#' implementation provide it directly (in C++); the others fall back to finite
#' differences of the Hessian (see \code{\link{numerical_deriv3}}).
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives
#'   (\eqn{\mathbb{E}[\partial^3 \ell]}) instead of the observed ones. Defaults to \code{FALSE}.
#' @inheritParams distrib_gradient
#' @param approx How the expectation is approximated when \code{expected = TRUE}
#'   and the distribution has no closed-form expected third derivatives; ignored
#'   otherwise. One of \code{"integrate"} (default), \code{"bartlett"} (equivalently
#'   \code{"opg"}) or \code{"mc"}. See \code{\link{expected_derivative_methods}}.
#' @param nsim Monte Carlo sample size used when \code{approx = "mc"}. Defaults to 10000.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of derivative-component vectors, keyed as in
#'   \code{\link{deriv_names}(distrib@params, 3)} (e.g. \code{"mu_mu_sigma"}).
#' @export
distrib_deriv3 <- S7::new_generic("distrib_deriv3", "distrib", function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    nat <- link_scale_lower_orders(distrib, y, theta, expected, 3L)
    to_link_scale(distrib, theta, c(nat, list(res)), 3L)
  } else {
    res
  }
})

#' Analytical Fourth-Order Derivatives
#'
#' @description
#' Computes the unique fourth-order partial derivatives of the log-likelihood with
#' respect to the distribution's parameters. Distributions with a closed-form
#' implementation provide it directly (in C++); the others fall back to finite
#' differences of the Hessian (see \code{\link{numerical_deriv4}}).
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives
#'   (\eqn{\mathbb{E}[\partial^4 \ell]}) instead of the observed ones. Defaults to \code{FALSE}.
#' @inheritParams distrib_gradient
#' @param approx How the expectation is approximated when \code{expected = TRUE}
#'   and the distribution has no closed-form expected fourth derivatives; ignored
#'   otherwise. One of \code{"integrate"} (default), \code{"bartlett"} (equivalently
#'   \code{"opg"}) or \code{"mc"}. See \code{\link{expected_derivative_methods}}.
#' @param nsim Monte Carlo sample size used when \code{approx = "mc"}. Defaults to 10000.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of derivative-component vectors, keyed as in
#'   \code{\link{deriv_names}(distrib@params, 4)} (e.g. \code{"mu_mu_sigma_sigma"}).
#' @export
distrib_deriv4 <- S7::new_generic("distrib_deriv4", "distrib", function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    nat <- link_scale_lower_orders(distrib, y, theta, expected, 4L)
    to_link_scale(distrib, theta, c(nat, list(res)), 4L)
  } else {
    res
  }
})

#' Gradient of the Log-Density with Respect to the Response
#'
#' @description
#' Computes the first derivative of the log-density with respect to the random
#' variable \eqn{y} (as opposed to the parameters), \eqn{\partial \ell / \partial y}.
#' This is defined for continuous distributions; distributions with a closed form
#' provide it directly, the others fall back to finite differences of the log-density.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @param ... Additional arguments passed to the specific method.
#' @return A numeric vector of the same length as \code{y}.
#' @export
distrib_grad_y <- S7::new_generic("distrib_grad_y", "distrib", function(distrib, y, theta, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  S7::S7_dispatch()
})

#' Second Derivative of the Log-Density with Respect to the Response
#'
#' @description
#' Computes the second derivative of the log-density with respect to the random
#' variable \eqn{y}, \eqn{\partial^2 \ell / \partial y^2}. Defined for continuous
#' distributions; distributions with a closed form provide it directly, the others
#' fall back to finite differences of the log-density.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @param ... Additional arguments passed to the specific method.
#' @return A numeric vector of the same length as \code{y}.
#' @export
distrib_hess_y <- S7::new_generic("distrib_hess_y", "distrib", function(distrib, y, theta, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  S7::S7_dispatch()
})

#' Generate Random Parameters
#'
#' @description Generates sensible random parameters for a distribution based on its mathematical domain.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method.
#' @export
generate_random_theta <- S7::new_generic("generate_random_theta", "distrib")
