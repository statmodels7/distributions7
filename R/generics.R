#' @title Distribution Generics
#' @name distrib_generics
#'
#' @description
#' The S7 generics that apply to every `distrib` object, grouped by what they
#' compute. This page is the index; each generic is documented on its own page,
#' and the method of each family for it on a third.
#'
#' The only compulsory method is [distrib_pdf()]. Every other quantity has a
#' fallback: the distribution function by quadrature (continuous) or by
#' summation (discrete), the quantile by root finding or by a lookup in the
#' cumulative table, random draws by inverse transform or by the
#' ratio-of-uniforms method, the derivatives of the log-density by finite
#' differences, and the expected derivatives by the strategies of
#' [expected_derivative_methods()]. A family defined by its density alone
#' therefore has every quantity, and a family that writes a quantity out
#' overrides the fallback by ordinary S7 dispatch.
#'
#' @details
#' # The surface
#' \describe{
#'   \item{probability}{[distrib_pdf()], [distrib_cdf()],
#'     [distrib_quantile()], [distrib_rng()], [distrib_atoms()].}
#'   \item{derivatives in the parameters}{[distrib_gradient()],
#'     [distrib_hessian()], [distrib_expected_hessian()], [distrib_deriv3()],
#'     [distrib_deriv4()], [distrib_deriv5()], and the derivatives of the
#'     expected information [distrib_dexpected_hessian()],
#'     [distrib_d2expected_hessian()], [distrib_d3expected_hessian()] and
#'     [distrib_d4expected_hessian()].}
#'   \item{derivatives in the response}{[distrib_grad_y()],
#'     [distrib_hess_y()], [distrib_deriv3_y()], [distrib_deriv4_y()].
#'     A discrete family has no method for them, its support having no
#'     continuous direction to differentiate along, so a call signals an
#'     error.}
#'   \item{mixed}{[distrib_cross_y()], [distrib_cross2_y()],
#'     [distrib_cross3_y()], [distrib_grad_y_hess()],
#'     [distrib_hess_y_hess()].}
#'   \item{of the distribution function}{[distrib_grad_cdf()],
#'     [distrib_hess_cdf()], [distrib_deriv3_cdf()], [distrib_deriv4_cdf()],
#'     which a censored likelihood needs.}
#'   \item{moments}{[mean()], [variance()], [std_dev()], [skewness()],
#'     [kurtosis()], [moment()], [expectation()].}
#'   \item{fitting}{[distrib_start()], and [fit_distrib()], which is a function
#'     rather than a generic.}
#'   \item{multivariate}{[mv_location()], [mv_sigma()], [mv_marginal()],
#'     [mv_support()], [mv_reference_draw()], [mv_derived()].}
#' }
#'
#' # The arguments theta and scale
#' `theta` is normalized before dispatch by `align_theta()`, which reorders it
#' by name, removes the names carried by the values and validates it against
#' `params_bounds`, treated as **open** intervals. The derivative generics in
#' the parameters of orders one to four also take `scale`, which the generic
#' applies **after** dispatch, so a method returns the parameter scale and
#' never reads it; see [link_scale_derivatives()]. [distrib_deriv5()] is the
#' exception: its methods read `scale` and apply it themselves.
#'
#' @return Nothing. This page is the index of the generics documented on their
#'   own pages; the value returned is theirs.
#'
#' @examples
#' # The Gaussian writes these quantities out; a family defined by its
#' # density alone would reach the fallbacks for the same calls.
#' d <- gaussian1_distrib()
#' th <- list(mu = 0, sigma = 1)
#'
#' distrib_pdf(d, 0, th)
#' distrib_cdf(d, 0, th)
#' distrib_quantile(d, 0.975, th)
#' distrib_gradient(d, c(-1, 1), th)
#'
#' # theta is reordered by name before dispatch, so the order it is written in
#' # does not matter.
#' identical(distrib_pdf(d, 1, list(mu = 0, sigma = 2)),
#'           distrib_pdf(d, 1, list(sigma = 2, mu = 0)))
#'
#' # The bounds are OPEN, so a scale of exactly zero is rejected rather than
#' # returning an infinite density.
#' try(distrib_pdf(d, 1, list(mu = 0, sigma = 0)))
#'
#' # A discrete family has no response derivatives, so this signals an error.
#' try(distrib_grad_y(poisson_distrib(), 2, list(mu = 3)))
#'
#' @seealso [distrib()] for the base class and its properties;
#'   [link_scale_derivatives()] for the `scale` argument;
#'   [expected_derivative_methods()] for the `approx` argument;
#'   [check_distrib()], which exercises the whole surface;
#'   the `defining-a-distribution` vignette.
#' @import S7
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
#' \eqn{f(y; \theta) = P(Y = y)}; `log = TRUE` returns
#' \eqn{\log f(y; \theta)}, the quantity that the derivative generics in the
#' parameters and in the response differentiate. This is the only method that
#' a distribution has to supply: every other quantity has a numerical fallback
#' derived from it.
#'
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of `distrib@params`.
#' @param ... Additional arguments passed to the specific method (e.g., `log`).
#' @return A numeric vector of density values, one per observation.
#' @examples
#' distrib_pdf(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' distrib_pdf(poisson_distrib(), 0:3, list(mu = 2), log = TRUE)
#' @seealso [distrib_cdf()], [distrib_quantile()], [distrib_rng()], [distrib_atoms()]
#' @export
distrib_pdf <- S7::new_generic("distrib_pdf", "distrib", function(distrib, y, theta, ...) {
  theta <- align_theta(distrib, theta)
  y <- recycle_point(distrib, y, theta)
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
#' one. `lower.tail = FALSE` returns \eqn{1 - F(q; \theta)} and
#' `log.p = TRUE` its logarithm, both computed on the log scale where a
#' family provides one. Without a method the value comes from quadrature of
#' the density for a continuous family and from summation of the mass for a
#' discrete one.
#'
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param q A numeric vector of quantiles.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of `distrib@params`.
#' @param ... Additional arguments passed to the specific method (e.g., `lower.tail`, `log.p`).
#' @return A numeric vector of cumulative probabilities.
#' @examples
#' distrib_cdf(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' distrib_cdf(poisson_distrib(), 0:3, list(mu = 2), lower.tail = FALSE)
#' @seealso [distrib_pdf()], [distrib_quantile()], [distrib_rng()], [distrib_atoms()]
#' @export
distrib_cdf <- S7::new_generic("distrib_cdf", "distrib", function(distrib, q, theta, ...) {
  theta <- align_theta(distrib, theta)
  q <- recycle_point(distrib, q, theta)
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
#' [distrib_cdf()] in the continuous case and from an exact table
#' lookup in the discrete one.
#'
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param p A numeric vector of probabilities.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of `distrib@params`.
#' @param ... Additional arguments passed to the specific method (e.g., `lower.tail`, `log.p`).
#' @return A numeric vector of quantiles.
#' @examples
#' distrib_quantile(gaussian1_distrib(), c(0.025, 0.5, 0.975), list(mu = 0, sigma = 1))
#' distrib_quantile(poisson_distrib(), c(0.1, 0.9), list(mu = 2))
#' @seealso [distrib_pdf()], [distrib_cdf()], [distrib_rng()], [distrib_atoms()]
#' @export
distrib_quantile <- S7::new_generic("distrib_quantile", "distrib", function(distrib, p, theta, ...) {
  theta <- align_theta(distrib, theta)
  p <- recycle_point(distrib, p, theta)
  S7::S7_dispatch()
})

#' Atoms of a Distribution
#'
#' @description
#' Returns the locations and probabilities of the point masses that a
#' distribution places on individual values. These form the discrete part of a
#' *mixed* distribution, one that is neither purely continuous nor purely
#' discrete.
#' [zero_adjusted()] applied to a continuous distribution builds exactly
#' such an object: a point mass at zero next to a density.
#'
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param theta A named list (or named numeric vector) of distribution parameters,
#'   with scalar entries.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A list with components `y` (the locations) and `p` (their
#'   probabilities), both numeric vectors of the same length, possibly of length
#'   zero.
#'
#' @details
#' The default returns no atoms, which is the right answer for every ordinary
#' distribution: a continuous one has none, and a discrete one is made of nothing
#' else, so listing them would be pointless. The generic exists for the case in
#' between, where a routine written for densities has to be told that part of the
#' mass is not in the integral --- [check_distrib()] uses it to know
#' that the density is expected to integrate to \eqn{1 - \sum p} rather than 1,
#' and to keep its finite differences away from the jumps.
#'
#' @examples
#' distrib_atoms(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#' distrib_atoms(zero_adjusted(gamma2_distrib()), list(mu = 2, sigma2 = 1, za = 0.3))
#'
#' @seealso [distrib_pdf()], [distrib_cdf()], [distrib_quantile()], [distrib_rng()]
#' @export
distrib_atoms <- S7::new_generic("distrib_atoms", "distrib", function(distrib, theta, ...) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' Random Number Generator
#'
#' @description Generates random variates from the given distribution.
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param n Number of observations to generate.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   If unnamed, parameters are taken in the order of `distrib@params`.
#' @param ... Additional arguments passed to the specific method.
#' @return A numeric vector of `n` draws for a univariate distribution, and
#'   an \eqn{n \times p} matrix for a multivariate one.
#' @examples
#' set.seed(1)
#' distrib_rng(gaussian1_distrib(), 5, list(mu = 0, sigma = 1))
#' distrib_rng(mvgaussian1_distrib(2), 3, list(mu1 = 0, mu2 = 0,
#'   sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5))
#' @seealso [distrib_pdf()], [distrib_cdf()], [distrib_quantile()], [distrib_atoms()]
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
#' Aligns `theta`, checks that every parameter has length 1 or \eqn{n}, and
#' recycles a scalar `y` up to \eqn{n} when `theta` is vectorized.
#'
#' @details
#' An empty `y` is allowed through untouched, giving empty derivatives, the
#' way `dnorm(numeric(0))` gives `numeric(0)` and the way
#' [distrib_pdf()] already behaves. Without the special case the
#' recycling check below rejects it with the nonsensical message
#' `"'y' must have length 1 or 1, not 0"`.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#'
#' @return A list with elements `y` and `theta`, both conformable.
#'
#' @seealso [align_theta()], [check_params_dim()]
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

#' Recycle a Single Point Against Parameters That Vary by Observation
#'
#' @description
#' Returns `x` repeated to the common length of the parameters when `x` has
#' length one and some parameter is longer, and `x` unchanged otherwise.
#'
#' @details
#' [distrib_pdf()], [distrib_cdf()] and [distrib_quantile()] call it before
#' dispatch, in the way that [check_derivative_args()] recycles the response
#' for the derivative generics. A compiled kernel sizes its output by the
#' length of the response, so a single point evaluated at parameters of length
#' \eqn{n} has to be repeated \eqn{n} times to give one value per observation.
#' This matters wherever a wrapper evaluates its parent at a fixed point, as
#' [zero_inflated()] does for \eqn{f(0)}. A multivariate response is a matrix
#' and is left alone.
#'
#' @param distrib A `distrib` object.
#' @param x The response, quantile or probability, a numeric vector.
#' @param theta The aligned parameter list.
#'
#' @return A numeric vector.
#'
#' @seealso [check_derivative_args()], the derivative generics' counterpart.
#' @keywords internal
recycle_point <- function(distrib, x, theta) {
  if (length(x) != 1L || S7::S7_inherits(distrib, multivariate_distrib)) return(x)
  n <- max(lengths(theta[seq_len(distrib@n_params)]), 1L)
  if (n > 1L) rep(x, n) else x
}

#' Analytical Gradient
#'
#' @description Computes the analytical first derivatives of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @param scale Either `"parameter"` (default) for derivatives with respect to
#'   the parameters \eqn{\theta} on their natural (constrained) scale, or
#'   `"link"` for derivatives with respect to the unconstrained linear
#'   predictors \eqn{\eta = g(\theta)} defined by `distrib@link_params`.
#'   See [link_scale_derivatives()].
#' @param ... Additional arguments passed to the specific method.
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`.
#' @examples
#' d <- gaussian1_distrib()
#' distrib_gradient(d, c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' # the same score with respect to the unconstrained parameters
#' distrib_gradient(d, c(-1, 0, 1), list(mu = 0, sigma = 1), scale = "link")
#' @seealso [distrib_hessian()], [distrib_expected_hessian()], [distrib_deriv3()], [distrib_deriv4()]
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
#' [distrib_expected_hessian()] for its expectation). On the link
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
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @inheritParams distrib_gradient
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of numeric vectors, keyed as
#'   [`hess_names(distrib@params)`][hess_names] (e.g. `"mu_sigma"`).
#' @examples
#' d <- gaussian1_distrib()
#' distrib_hessian(d, c(-1, 0, 1), list(mu = 0, sigma = 1))
#' @seealso [distrib_gradient()], [distrib_expected_hessian()], [distrib_deriv3()], [distrib_deriv4()]
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
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @inheritParams distrib_gradient
#' @param approx How the expectation is approximated **when the distribution
#'   has no closed-form expected Hessian**; ignored otherwise. One of
#'   `"opg"` (the default: the outer product of the score, read at each
#'   observation without taking an expectation), `"bartlett"` (the second
#'   Bartlett identity, with the expectation taken by summation or
#'   quadrature), `"integrate"` (quadrature or summation of the observed
#'   Hessian) or `"mc"` (Monte Carlo). See [expected_derivative_methods()] for
#'   the accuracy and speed of each.
#' @param nsim Monte Carlo sample size used when `approx = "mc"`. Defaults to 10000.
#' @param ... Additional arguments passed to the specific method.
#' @details
#' On the link scale the first-order term of the chain rule drops out because the
#' score has zero expectation, so the expected Hessian transforms as the simple
#' congruence \eqn{\mathrm{diag}(h')\, \mathbb{E}[H]\, \mathrm{diag}(h')} with
#' \eqn{h' = dg^{-1}/d\eta}.
#' @return A named list of numeric vectors, keyed as
#'   [`hess_names(distrib@params)`][hess_names], holding the expected second
#'   derivatives, that is minus the Fisher information.
#' @examples
#' d <- gaussian1_distrib()
#' distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1))
#'
#' # the pseudo-Huber computes its expected information itself, so 'approx'
#' # has no effect here; it applies to a family without a closed form, such as
#' # one defined by its density alone
#' distrib_expected_hessian(
#'   pseudohuber_distrib(), 0, list(mu = 0, sigma = 1, nu = 1),
#'   approx = "integrate"
#' )
#' @seealso [distrib_gradient()], [distrib_hessian()], [distrib_deriv3()], [distrib_deriv4()]
#' @export
distrib_expected_hessian <- S7::new_generic("distrib_expected_hessian", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
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
#' differences of the highest order the family implements itself (see
#' [numerical_deriv3()]).
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @param expected Logical; if `TRUE`, returns the expected third derivatives
#'   (\eqn{\mathbb{E}[\partial^3 \ell]}) instead of the observed ones. Defaults to `FALSE`.
#' @inheritParams distrib_gradient
#' @param approx How the expectation is approximated when `expected = TRUE`
#'   and the distribution has no closed-form expected third derivatives; ignored
#'   otherwise. One of `"integrate"` (default), `"bartlett"` (equivalently
#'   `"opg"`) or `"mc"`. See [expected_derivative_methods()].
#' @param nsim Monte Carlo sample size used when `approx = "mc"`. Defaults to 10000.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of derivative-component vectors, keyed as in
#'   [`deriv_names(distrib@params, 3)`][deriv_names] (e.g. `"mu_mu_sigma"`).
#' @examples
#' distrib_deriv3(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' @seealso [distrib_gradient()], [distrib_hessian()], [distrib_expected_hessian()], [distrib_deriv4()]
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
#' differences of the highest order the family implements itself, up to the
#' Hessian (see [numerical_deriv4()]).
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @param expected Logical; if `TRUE`, returns the expected fourth derivatives
#'   (\eqn{\mathbb{E}[\partial^4 \ell]}) instead of the observed ones. Defaults to `FALSE`.
#' @inheritParams distrib_gradient
#' @param approx How the expectation is approximated when `expected = TRUE`
#'   and the distribution has no closed-form expected fourth derivatives; ignored
#'   otherwise. One of `"integrate"` (default), `"bartlett"` (equivalently
#'   `"opg"`) or `"mc"`. See [expected_derivative_methods()].
#' @param nsim Monte Carlo sample size used when `approx = "mc"`. Defaults to 10000.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of derivative-component vectors, keyed as in
#'   [`deriv_names(distrib@params, 4)`][deriv_names] (e.g. `"mu_mu_sigma_sigma"`).
#' @examples
#' distrib_deriv4(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' @seealso [distrib_gradient()], [distrib_hessian()], [distrib_expected_hessian()], [distrib_deriv3()]
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

#' Fifth-Order Derivatives
#'
#' @description
#' The unique fifth-order partial derivatives of the log-likelihood with
#' respect to the distribution's parameters, on the parameter scale or on the
#' link scale.
#'
#' @details
#' The order is needed because each order of differentiating the predictor of
#' a score-driven filter through its own recursion requires one more order of
#' the family: the curvature reaches the third order, the directional third
#' derivative the fourth, and the outer Hessian of a model with such a term the
#' fifth. The fifth order is therefore provided by the family, like the lower
#' orders.
#'
#' The families that register their own method ([gengamma1_distrib()],
#' [gengamma2_distrib()], [gpd_distrib()], [lognormal2_distrib()],
#' [pseudohuber2_distrib()], [skewnormal2_distrib()], [skewt_distrib()],
#' [student_t1_distrib()], [student_t2_distrib()] and [weibull3_distrib()])
#' compute the fifth order in closed form with compiled kernels. Every other
#' family reaches [numerical_deriv5()], which applies **one** central
#' difference to the fourth order. Where the fourth order is analytic this is
#' an ordinary numerical derivative; where the fourth order is itself a
#' fallback the result is a difference of a difference, with the accuracy
#' described on [has_exact_deriv4()].
#'
#' Unlike the generics of lower order this one carries no `expected` argument.
#' The criterion of a model with a filter reads the observed information (a
#' filter has no expected one), and an expectation of the fifth order would
#' sit on top of a numerical quantity. [expected_derivative()] is generic in
#' the order, so the expected fifth order can be added when a caller needs
#' it.
#'
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution
#'   parameters. Each parameter must have length 1 or `length(y)`.
#' @param scale `"parameter"` (the default) for derivatives with respect to
#'   the parameters, `"link"` for derivatives with respect to the
#'   unconstrained predictors. In the numerical default the link scale is
#'   obtained by differentiating the order-4 component **already on the link
#'   scale**, which needs neither the fifth derivative of a link nor an
#'   order-5 entry in [bell_partial()]. A family with a closed-form fifth
#'   order composes it with the lower orders through [to_link_scale()], which
#'   uses both.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list of derivative-component vectors, each of length
#'   `length(y)`, keyed as
#'   [`deriv_names(distrib@params, 5)`][deriv_names] gives them (e.g.
#'   `"mu_mu_sigma_sigma_sigma"`).
#'
#' @examples
#' distrib_deriv5(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @seealso [distrib_deriv4()] for the order below, [numerical_deriv5()] for
#'   the differencing, and [has_exact_deriv4()] for the predicate that says
#'   whether the quantity being differenced is itself analytic.
#' @export
distrib_deriv5 <- S7::new_generic("distrib_deriv5", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  S7::S7_dispatch()
})

#' Whether a Family Supplies Its Own Fourth Derivative
#'
#' @description
#' Whether the family supplies its own fourth-order derivatives, rather than
#' reaching one of the base classes' numerical fallbacks.
#'
#' @details
#' The numerical fifth order is one central difference of the fourth, so its
#' accuracy depends on the quantity that it differences. Where the fourth order
#' is the family's own, the result is an ordinary numerical derivative of a
#' quantity that the family computes directly. Where the fourth order is the
#' numerical fallback of the package, the fifth order is a difference of a
#' difference, and no choice of tolerance makes it reliable.
#'
#' A family in the second case still receives a value, since a distribution
#' needs only [distrib_pdf()], and the value comes back without a warning. Its
#' error can reach the first significant digit, so it is not suitable for
#' building on. [check_distrib()] emits no order-5 row for such a family,
#' because the check does not apply to it.
#'
#' The default method reads the class that owns the registered
#' [distrib_deriv4()] method, the documented S7 route. A wrapper overrides it
#' and returns the answer of its parent instead: the fourth derivative of a
#' wrapper is a partition sum over the first four derivatives of the parent,
#' so it is the wrapper's own exactly when those of the parent are.
#'
#' The predicate is one logical for the whole family. A family may own its
#' fourth-order method and still build part of it from single stencils on
#' analytic quantities; it returns `TRUE`, and the order-5 row of
#' [check_distrib()] is then emitted for it. That row appears only when `5` is
#' among the requested `orders`, which default to `1:4`.
#'
#' @param x A distribution object.
#' @param ... Unused.
#'
#' @return A single logical.
#'
#' @examples
#' has_exact_deriv4(gaussian1_distrib())
#' has_exact_deriv4(truncated(gaussian1_distrib(), lower = 0))
#'
#' @seealso [distrib_deriv5()], the consumer; [expected_hessian_exact()],
#'   the predicate this is modeled on.
#' @export
has_exact_deriv4 <- S7::new_generic("has_exact_deriv4", "x")


#' Gradient of the Log-Density with Respect to the Response
#'
#' @description
#' Computes the first derivative of the log-density with respect to the random
#' variable \eqn{y} (as opposed to the parameters), \eqn{\partial \ell / \partial y}.
#' This is defined for continuous distributions; distributions with a closed form
#' provide it directly, the others fall back to finite differences of the log-density.
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @param ... Additional arguments passed to the specific method.
#' @return A numeric vector of the same length as `y`.
#' @examples
#' distrib_grad_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' @seealso [distrib_hess_y()], [distrib_cross_y()]
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
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(y)`.
#' @param ... Additional arguments passed to the specific method.
#' @return A numeric vector of the same length as `y`.
#' @examples
#' distrib_hess_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#' @seealso [distrib_grad_y()], [distrib_cross_y()]
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
#' `lower.tail = FALSE`.
#'
#' These are what a **censored** observation contributes to the score. An
#' observation known only to be at most \eqn{q} contributes \eqn{\log F(q)}, one
#' known only to exceed \eqn{q} contributes \eqn{\log(1-F(q))}, and an interval
#' censored one contributes \eqn{\log(F(b) - F(a))}, which is assembled from the
#' unlogged derivatives (`log = FALSE`) at the two endpoints. They are also
#' what the delta method needs for the standard error of a quantile residual.
#'
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param q A numeric vector of quantiles.
#' @param theta A named list (or named numeric vector) of distribution parameters.
#'   Each parameter must have length 1 or `length(q)`.
#' @param lower.tail Logical; if `TRUE` (default), derivatives of
#'   \eqn{\log F(q)}, otherwise of \eqn{\log(1 - F(q))}.
#' @param log Logical; if `TRUE` (default), derivatives of the *log*
#'   tail probability. With `FALSE` the derivatives of the probability
#'   itself are returned, which interval censoring and the truncation
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
#' so the gradient of the log distribution function is a *partial mean of
#' the score*. For a discrete distribution the integral is a finite sum and the
#' identity is exact, which is how the default method computes it there; for a
#' continuous one the default differentiates [distrib_cdf()]
#' numerically, and distributions with a closed form register it directly.
#'
#' @examples
#' d <- gaussian1_distrib()
#' theta <- list(mu = 0, sigma = 1)
#'
#' # what a right-censored observation at q = 1 contributes to the score
#' distrib_grad_cdf(d, 1, theta, lower.tail = FALSE)
#'
#' @seealso [distrib_hess_cdf()], [distrib_gradient()]
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
#' `lower.tail = FALSE`. Together with [distrib_grad_cdf()] these
#' give the observed information of a censored observation.
#'
#' @inheritParams distrib_grad_cdf
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @details
#' By the same exchange as in [distrib_grad_cdf()], and using
#' \eqn{\partial_{ij} f / f = \ell^{(ij)} + \ell^{(i)}\ell^{(j)}},
#' \deqn{\frac{\partial^{2} F(q)}{\partial\theta_i\partial\theta_j}
#'   = F(q)\;\mathbb{E}\!\left[\ell^{(ij)} + \ell^{(i)}\ell^{(j)} \mid Y \leq q\right],}
#' and the log scale follows from
#' \eqn{\partial_{ij}\log P = \partial_{ij}P/P - (\partial_i P/P)(\partial_j P/P)}.
#'
#' @examples
#' distrib_hess_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1))
#'
#' @seealso [distrib_grad_cdf()], [distrib_hessian()]
#' @export
distrib_hess_cdf <- S7::new_generic("distrib_hess_cdf", "distrib", function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  args <- check_derivative_args(distrib, q, theta)
  q <- args$y
  theta <- args$theta
  S7::S7_dispatch()
})

#' Whether a Family's Expected Information Is Written Out
#'
#' @description
#' `TRUE` when the family computes its expected information exactly, in
#' closed form or by a quadrature of its own to working precision, and `FALSE`
#' when a call to [distrib_expected_hessian()] reaches a fallback and the
#' answer is therefore an approximation. It is a generic so that a family
#' whose registered method is not what its owning class suggests can report
#' this; [has_exact_expected_hessian()] calls it and
#' [expected_hessian_exact.distrib()] is the default.
#'
#' @details
#' The predicate is what a consumer reads to decide which information to
#' report. Where it returns `TRUE` the expected information costs one
#' evaluation and has the smaller variance, so it is the better matrix to
#' invert; where it returns `FALSE` the choice is between an approximation and
#' the observed Hessian, which every family has and which is exact. That is
#' the rule that [fit_distrib()] follows for its standard errors and that
#' \pkg{statmodels7} follows in `vcov()`.
#'
#' Every shipped univariate family returns `TRUE`. The Poisson-inverse
#' Gaussian families sum their expected information exactly over the support,
#' and the skew normals, the skew t and the pseudo-Huber evaluate one
#' quadrature over the standardized response per distinct shape; see
#' [loc_scale_expected()]. Among the wrappers, [truncated()] returns `FALSE`.
#'
#' Reading the class that owns the registered method is not sufficient on its
#' own, which is why the generic exists: a method registered on the class of a
#' family may chain onto a parent that approximates, and a family that does so
#' returns the answer of its parent.
#'
#' @param x An object inheriting from class `"distrib"`.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' # A Gaussian writes its information out; a truncated Gaussian does not.
#' expected_hessian_exact(gaussian1_distrib())
#' expected_hessian_exact(truncated(gaussian1_distrib(), lower = 0))
#'
#' @seealso [has_exact_expected_hessian()], the wrapper that calls it;
#'   [distrib_expected_hessian()] for the quantity itself, and
#'   [expected_by_opg()] for the fallback a `FALSE` answer reaches.
#'
#' @export
expected_hessian_exact <- S7::new_generic("expected_hessian_exact", "x")

#' The Kink of a Family
#'
#' @description
#' The one non-smooth piece of a log-density, as the composition
#' \eqn{c(\theta)\,\phi(v(y,\theta))} of [kink_spec()], or `NULL` for a family
#' that is smooth in every parameter. [params_order()] reads the order of
#' differentiability off the answer.
#'
#' @details
#' Three families declare one ([laplace_distrib()], [laplace2_distrib()] and
#' [enet_distrib()]), all for the same reason: an absolute value of the
#' residual in the log-density. Every other family inherits the base method,
#' which returns `NULL`.
#'
#' A wrapper does not propagate the declaration today, so a wrapper of a kinked
#' family returns `NULL` here while `params_smooth` still records the kink;
#' [params_order()] reports `NA` there rather than `Inf`, which says the order
#' has not been established rather than that the parameter is smooth.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#'
#' @return A [kink_spec()], or `NULL` when the family is smooth.
#'
#' @examples
#' kink_decomposition(gaussian1_distrib())
#' kink_decomposition(laplace_distrib())
#'
#' @seealso [kink_spec()] for the object, [params_order()] for the order
#'   deduced from it, [check_kink()] for the validator, and
#'   [param_smoothness()] for the logical that this generalizes.
#' @export
kink_decomposition <- S7::new_generic("kink_decomposition", "distrib",
                                      function(distrib) S7::S7_dispatch())

#' Generate Random Parameters
#'
#' @description Generates sensible random parameters for a distribution based on its mathematical domain.
#' @param distrib A distribution object inheriting from the `distrib` class.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of parameter values, one per element of
#'   `distrib@params`, each inside that parameter's bounds.
#' @examples
#' set.seed(1)
#' generate_random_theta(gamma2_distrib())
#' @seealso [deriv_names()], [hess_names()], [expand_params()], [transpose_params()], [check_params_dim()], [check_theta_bounds()], [param_smoothness()]
#' @export
generate_random_theta <- S7::new_generic("generate_random_theta", "distrib")
