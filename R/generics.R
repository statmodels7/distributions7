#' @title Distribution Generics
#' @description A collection of S7 generic functions for mathematical and statistical
#' operations on probability distributions.
#' @import S7
#' @name distrib_generics
#' @return Nothing. This page is the index of the generics documented on their
#'   own pages; the value returned is theirs.
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
#'
#' @details
#' For a continuous family \eqn{f(y; \theta)} is the density with respect to
#' the Lebesgue measure and for a discrete one the mass
#' \eqn{f(y; \theta) = P(Y = y)}; \code{log = TRUE} returns
#' \eqn{\log f(y; \theta)}, which is the quantity every derivative generic
#' differentiates. This is the only method a distribution has to supply:
#' every other quantity has a numerical fallback derived from it.
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method (e.g., \code{log}).
#' @return A numeric vector of density values, one per observation.
#' @examples
#' distrib_pdf(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' distrib_pdf(poisson_distrib(), 0:3, list(mu = 2), log = TRUE)
#' @export
distrib_pdf <- S7::new_generic("distrib_pdf", "distrib", function(distrib, y, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Cumulative Distribution Function
#'
#' @description Evaluates the cumulative distribution function (CDF) for a given distribution.
#'
#' @details
#' \deqn{F(q; \theta) = P(Y \le q),}
#'
#' the integral of the density up to \eqn{q} for a continuous family and the
#' sum of the mass over the support points at or below \eqn{q} for a discrete
#' one. \code{lower.tail = FALSE} returns \eqn{1 - F(q; \theta)} and
#' \code{log.p = TRUE} its logarithm, both computed on the log scale where a
#' family provides one. Without a method the value comes from quadrature of
#' the density.
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param q A numeric vector of quantiles.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method (e.g., \code{lower.tail}, \code{log.p}).
#' @return A numeric vector of cumulative probabilities.
#' @examples
#' distrib_cdf(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' distrib_cdf(poisson_distrib(), 0:3, list(mu = 2), lower.tail = FALSE)
#' @export
distrib_cdf <- S7::new_generic("distrib_cdf", "distrib", function(distrib, q, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Quantile Function
#'
#' @description Evaluates the quantile function for a given distribution.
#'
#' @details
#' The generalized inverse of the distribution function,
#'
#' \deqn{Q(p; \theta) = \inf\{y : F(y; \theta) \ge p\},}
#'
#' which for a continuous strictly increasing \eqn{F} is the ordinary inverse
#' and for a discrete family the smallest support point whose cumulative mass
#' reaches \eqn{p}. Without a method the value comes from root-finding on
#' \code{\link{distrib_cdf}} in the continuous case and from an exact table
#' lookup in the discrete one.
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param p A numeric vector of probabilities.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of \code{distrib@params}.
#' @param ... Additional arguments passed to the specific method (e.g., \code{lower.tail}, \code{log.p}).
#' @return A numeric vector of quantiles.
#' @examples
#' distrib_quantile(gaussian1_distrib(), c(0.025, 0.5, 0.975), list(mu = 0, sigma = 1))
#' distrib_quantile(poisson_distrib(), c(0.1, 0.9), list(mu = 2))
#' @export
distrib_quantile <- S7::new_generic("distrib_quantile", "distrib", function(distrib, p, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Atoms of a Distribution
#'
#' @description
#' Returns the locations and probabilities of the point masses a distribution
#' places on individual values --- the discrete part of a \emph{mixed}
#' distribution, one that is neither purely continuous nor purely discrete.
#' \code{\link{zero_adjusted}} applied to a continuous distribution builds exactly
#' such an object: a point mass at zero next to a density.
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param theta A named list (or named numeric vector) of distribution parameters,
#'   with scalar entries.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A list with components \code{y} (the locations) and \code{p} (their
#'   probabilities), both numeric vectors of the same length, possibly of length
#'   zero.
#'
#' @details
#' The default returns no atoms, which is the right answer for every ordinary
#' distribution: a continuous one has none, and a discrete one is made of nothing
#' else, so listing them would be pointless. The generic exists for the case in
#' between, where a routine written for densities has to be told that part of the
#' mass is not in the integral --- \code{\link{check_distrib}} uses it to know
#' that the density is expected to integrate to \eqn{1 - \sum p} rather than 1,
#' and to keep its finite differences away from the jumps.
#'
#' @examples
#' distrib_atoms(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#' distrib_atoms(zero_adjusted(gamma2_distrib()), list(mu = 2, sigma2 = 1, za = 0.3))
#'
#' @export
distrib_atoms <- S7::new_generic("distrib_atoms", "distrib", function(distrib, theta, ...) {
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
#' @return A numeric vector of \code{n} draws for a univariate distribution, and
#'   an \eqn{n \times p} matrix for a multivariate one.
#' @examples
#' set.seed(1)
#' distrib_rng(gaussian1_distrib(), 5, list(mu = 0, sigma = 1))
#' distrib_rng(mvgaussian_distrib(2), 3, list(mu1 = 0, mu2 = 0,
#'   sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5))
#' @export
distrib_rng <- S7::new_generic("distrib_rng", "distrib", function(distrib, n, theta, ...) {
  theta <- align_theta(distrib, theta)
  # Recycle per-observation parameters up to n, as the base generators do
  # (rnorm(4, c(0, 10)) draws two from each mean). Methods are free to index
  # theta by draw; without this a shorter parameter vector either silently
  # produces NAs, when the method subsets it by a longer logical, or is rejected
  # outright, when the method delegates to distrib_quantile. Scalars are left
  # alone so that the C++ kernels keep their fast path.
  idx <- seq_len(distrib@n_params)
  lens <- lengths(theta[idx])
  if (any(lens > 1L) && !all(lens %in% c(1L, n))) {
    theta[idx] <- lapply(theta[idx], function(x) if (length(x) == 1L) x else rep_len(x, n))
  }
  S7::S7_dispatch()
})

#' Shared Argument Handling for the Derivative Generics
#'
#' @description
#' Aligns \code{theta}, checks that every parameter has length 1 or \eqn{n}, and
#' recycles a scalar \code{y} up to \eqn{n} when \code{theta} is vectorized.
#'
#' @details
#' An empty \code{y} is allowed through untouched, giving empty derivatives, the
#' way \code{dnorm(numeric(0))} gives \code{numeric(0)} and the way
#' \code{\link{distrib_pdf}} already behaves. Without the special case the
#' recycling check below rejects it with the nonsensical message
#' \code{"'y' must have length 1 or 1, not 0"}.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#'
#' @return A list with elements \code{y} and \code{theta}, both conformable.
#'
#' @seealso \code{\link{align_theta}}, \code{\link{check_params_dim}}
#' @keywords internal
check_derivative_args <- function(distrib, y, theta) {
  theta <- align_theta(distrib, theta)
  pars <- theta[seq_len(distrib@n_params)]

  # An empty y yields empty derivatives, the way dnorm(numeric(0)) yields
  # numeric(0) and the way distrib_pdf already behaves. Without this the
  # recycling check below rejects it with the nonsensical message
  # "'y' must have length 1 or 1, not 0".
  if (length(y) == 0L && all(lengths(pars) <= 1L)) {
    return(list(y = y, theta = theta))
  }

  # A multivariate response is a matrix, so the number of observations is its
  # row count rather than its length, and there is nothing to recycle: the
  # parameters of a multivariate distribution are scalars for the whole sample
  # (see multivariate_distrib), so the conformability question does not arise.
  if (S7::S7_inherits(distrib, multivariate_distrib)) {
    return(list(y = y, theta = theta))
  }

  n <- max(length(y), lengths(pars))
  check_params_dim(pars, n = n)
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
#' @return A named list with one numeric vector per parameter, keyed by
#'   \code{distrib@params}.
#' @examples
#' d <- gaussian1_distrib()
#' distrib_gradient(d, c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' # the same score with respect to the unconstrained parameters
#' distrib_gradient(d, c(-1, 0, 1), list(mu = 0, sigma = 1), scale = "link")
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
#'
#' @details
#' With \eqn{l(\theta; y) = \log f(y; \theta)}, one entry per unordered pair
#' of parameters,
#'
#' \deqn{l^{(ij)} = \frac{\partial^{2} l}{\partial \theta_i \partial \theta_j},}
#'
#' evaluated at the observed \eqn{y} (see
#' \code{\link{distrib_expected_hessian}} for its expectation). On the link
#' scale the reparametrization \eqn{\theta_i = h_i(\eta_i)} is diagonal, so
#' the second-order chain rule carries a first-order term on the diagonal
#' alone:
#'
#' \deqn{\frac{\partial^{2} l}{\partial \eta_i \partial \eta_j}
#'   = l^{(ij)} h_i'(\eta_i) h_j'(\eta_j)
#'   + \delta_{ij}\, l^{(i)} h_i''(\eta_i).}
#'
#' The transformation is applied in the generic, so a method always returns
#' the parameter scale.
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(y)}.
#' @inheritParams distrib_gradient
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of numeric vectors, keyed as
#'   \code{\link{hess_names}(distrib@params)} (e.g. \code{"mu_sigma"}).
#' @examples
#' d <- gaussian1_distrib()
#' distrib_hessian(d, c(-1, 0, 1), list(mu = 0, sigma = 1))
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
#' @return A named list of numeric vectors, keyed as
#'   \code{\link{hess_names}(distrib@params)}, holding the expected second
#'   derivatives, that is minus the Fisher information.
#' @examples
#' d <- gaussian1_distrib()
#' distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1))
#'
#' # a family with no closed form uses the strategy named by 'approx'
#' distrib_expected_hessian(
#'   pseudohuber_distrib(), 0, list(mu = 0, sigma = 1, nu = 1),
#'   approx = "integrate"
#' )
#' @export
distrib_expected_hessian <- S7::new_generic("distrib_expected_hessian", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    zero1 <- stats::setNames(
      # n_obs(), not length(y): for a matrix response length(y) counts entries,
      # so this vector came out n*p long. Recycled against the p-long
      # components of the expected Hessian it inflated every DIAGONAL entry of
      # the information by a factor of p -- the first-order term of the
      # order-2 chain rule appears only there -- and every standard error of a
      # multivariate fit came out a factor of sqrt(p) too small.
      lapply(distrib@params, function(p) rep(0, n_obs(distrib, y))), distrib@params
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
#' @examples
#' distrib_deriv3(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
#' @examples
#' distrib_deriv4(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
#' @examples
#' distrib_grad_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
#' @examples
#' distrib_hess_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' @export
distrib_hess_y <- S7::new_generic("distrib_hess_y", "distrib", function(distrib, y, theta, ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  S7::S7_dispatch()
})

#' Gradient of the Log Distribution Function
#'
#' @description
#' Computes the first derivatives, with respect to the parameters, of
#' \eqn{\log F(q;\theta)} --- or of \eqn{\log(1 - F(q;\theta))} when
#' \code{lower.tail = FALSE}.
#'
#' These are what a \strong{censored} observation contributes to the score. An
#' observation known only to be at most \eqn{q} contributes \eqn{\log F(q)}, one
#' known only to exceed \eqn{q} contributes \eqn{\log(1-F(q))}, and an interval
#' censored one contributes \eqn{\log(F(b) - F(a))}, which is assembled from the
#' unlogged derivatives (\code{log = FALSE}) at the two endpoints. They are also
#' what the delta method needs for the standard error of a quantile residual.
#'
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param q A numeric vector of quantiles.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or \code{length(q)}.
#' @param lower.tail Logical; if \code{TRUE} (default), derivatives of
#'   \eqn{\log F(q)}, otherwise of \eqn{\log(1 - F(q))}.
#' @param log Logical; if \code{TRUE} (default), derivatives of the \emph{log}
#'   tail probability. With \code{FALSE} the derivatives of the probability
#'   itself are returned, which is what interval censoring and the truncation
#'   constant are built from.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list with one numeric vector per parameter.
#'
#' @details
#' The mathematics is one exchange of derivative and integral. Since the region
#' of integration does not depend on \eqn{\theta},
#' \deqn{\frac{\partial F(q;\theta)}{\partial\theta_i}
#'   = \int_{-\infty}^{q}\frac{\partial f}{\partial\theta_i}
#'   = \int_{-\infty}^{q} f\,\ell^{(i)}
#'   = F(q)\;\mathbb{E}\!\left[\ell^{(i)} \mid Y \leq q\right],}
#' so the gradient of the log distribution function is a \emph{partial mean of
#' the score}. For a discrete distribution the integral is a finite sum and the
#' identity is exact, which is how the default method computes it there; for a
#' continuous one the default differentiates \code{\link{distrib_cdf}}
#' numerically, and distributions with a closed form register it directly.
#'
#' @examples
#' d <- gaussian1_distrib()
#' theta <- list(mu = 0, sigma = 1)
#'
#' # what a right-censored observation at q = 1 contributes to the score
#' distrib_grad_cdf(d, 1, theta, lower.tail = FALSE)
#'
#' @seealso \code{\link{distrib_hess_cdf}}, \code{\link{distrib_gradient}}
#' @export
distrib_grad_cdf <- S7::new_generic("distrib_grad_cdf", "distrib", function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  args <- check_derivative_args(distrib, q, theta)
  q <- args$y
  theta <- args$theta
  S7::S7_dispatch()
})

#' Second Derivatives of the Log Distribution Function
#'
#' @description
#' Computes the second derivatives, with respect to the parameters, of
#' \eqn{\log F(q;\theta)} --- or of \eqn{\log(1 - F(q;\theta))} when
#' \code{lower.tail = FALSE}. Together with \code{\link{distrib_grad_cdf}} these
#' give the observed information of a censored observation.
#'
#' @inheritParams distrib_grad_cdf
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list keyed as \code{\link{hess_names}(distrib@params)}.
#'
#' @details
#' By the same exchange as in \code{\link{distrib_grad_cdf}}, and using
#' \eqn{\partial_{ij} f / f = \ell^{(ij)} + \ell^{(i)}\ell^{(j)}},
#' \deqn{\frac{\partial^{2} F(q)}{\partial\theta_i\partial\theta_j}
#'   = F(q)\;\mathbb{E}\!\left[\ell^{(ij)} + \ell^{(i)}\ell^{(j)} \mid Y \leq q\right],}
#' and the log scale follows from
#' \eqn{\partial_{ij}\log P = \partial_{ij}P/P - (\partial_i P/P)(\partial_j P/P)}.
#'
#' @examples
#' distrib_hess_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1))
#'
#' @seealso \code{\link{distrib_grad_cdf}}, \code{\link{distrib_hessian}}
#' @export
distrib_hess_cdf <- S7::new_generic("distrib_hess_cdf", "distrib", function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  args <- check_derivative_args(distrib, q, theta)
  q <- args$y
  theta <- args$theta
  S7::S7_dispatch()
})

#' Generate Random Parameters
#'
#' @description Generates sensible random parameters for a distribution based on its mathematical domain.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of parameter values, one per element of
#'   \code{distrib@params}, each inside that parameter's bounds.
#' @examples
#' set.seed(1)
#' generate_random_theta(gamma2_distrib())
#' @export
generate_random_theta <- S7::new_generic("generate_random_theta", "distrib")
