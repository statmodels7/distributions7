#' @include distrib.R generics.R
NULL

#' @title Weibull Distribution Class, Scale and Shape
#' @name Weibull1Distrib
#'
#' @description
#' The S7 class of the Weibull family parametrized by a scale \eqn{\mu > 0} and
#' a shape \eqn{\sigma > 0}, with density
#' \eqn{f(y) = (\sigma/\mu)(y/\mu)^{\sigma-1}\exp\{-(y/\mu)^{\sigma}\}} on the
#' positive half line. It inherits from `continuous_distrib`, so it answers
#' every generic of the `distrib` contract; the methods listed below are
#' registered on it directly and everything else comes from the parent.
#'
#' Build one with [weibull1_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Weibull1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [weibull1_distrib()] they hold `"weibull1"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "sigma")`, the interpretations
#'   `c(mu = "scale", sigma = "shape")`, `2`, the domain \eqn{(0, \infty)} for
#'   both parameters, and the two links.
#'
#' @seealso [weibull1_distrib()] to build one;
#'   [weibull3_distrib()] for the same law parametrized by a quantile;
#'   [gumbel_distrib()], which is this family on a negative log scale;
#'   [distrib_pdf.Weibull1Distrib()] and
#'   [distrib_gradient.Weibull1Distrib()] for the closed forms it supplies.
#'
#' @section Methods:
#' Registered on this class in this file:
#'   [`distrib_cdf()`][distrib_cdf.Weibull1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Weibull1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Weibull1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Weibull1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Weibull1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Weibull1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Weibull1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Weibull1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Weibull1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Weibull1Distrib],
#'   [`distrib_rng()`][distrib_rng.Weibull1Distrib].
#'
#' Four more are registered on the class from other files:
#'   [`distrib_cross_y()`][distrib_cross_y.Weibull1Distrib] in
#'   `cross_derivatives_families.R`,
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.Weibull1Distrib] in
#'   `cdf_survival_higher.R`, and the four moments
#'   [`mean()`][mean.Weibull1Distrib], [`variance()`][variance.Weibull1Distrib],
#'   [`skewness()`][skewness.Weibull1Distrib] and
#'   [`kurtosis()`][kurtosis.Weibull1Distrib] in `moments.R`.
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- weibull1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # Both parameters are positive, so both ride a log by default.
#' vapply(d@link_params, function(l) l@link_name, character(1))
Weibull1Distrib <- S7::new_class("Weibull1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' The Pieces a Weibull Evaluates From
#'
#' @description
#' Returns the standardized variable \eqn{z = y/\mu}, its logarithm and the
#' quantity \eqn{u = z^{\sigma}}. Every derivative of the Weibull log-density,
#' in the parameters and in the response, is a polynomial in \eqn{u} and
#' \eqn{u \log z}, so the three are computed once and shared.
#'
#' @details
#' \eqn{u} is the substitution that makes the family tractable. Under the model
#' \eqn{u \sim \mathrm{Exp}(1)} whatever the parameters are, so an expectation
#' of any polynomial in \eqn{u} and \eqn{\log u} is a derivative of the gamma
#' function at 2, and [distrib_expected_hessian.Weibull1Distrib()] rests on
#' that. The logarithm is formed before the power so that \eqn{u} is computed
#' as \eqn{\exp(\sigma \log z)}, which stays representable for a large shape at
#' a moderate \eqn{z}.
#'
#' @param y A numeric vector of observations, each positive. A non-positive
#'   value gives `NaN` for `lz` and propagates.
#' @param mu The scale parameter, a numeric vector of length 1 or of the length
#'   of `y`, strictly positive.
#' @param sigma The shape parameter, a numeric vector of length 1 or of the
#'   length of `y`, strictly positive.
#'
#' @return A named list of three numeric vectors: `z`, the ratio \eqn{y/\mu};
#'   `lz`, its logarithm; and `u`, the power \eqn{z^{\sigma}}. Each has the
#'   recycled length of the inputs.
#'
#' @keywords internal
weibull_pieces <- function(y, mu, sigma) {
  z <- y / mu
  lz <- log(z)
  list(z = z, lz = lz, u = exp(sigma * lz))
}

#' @title Weibull Probability Density Function
#' @name distrib_pdf.Weibull1Distrib
#' @description
#' Computes the Weibull density
#' \deqn{f(y; \mu, \sigma) = \dfrac{\sigma}{\mu}
#'   \left(\dfrac{y}{\mu}\right)^{\sigma - 1}
#'   \exp\left\{-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\}, \qquad y > 0,}
#' by calling [stats::dweibull()] at `shape = sigma` and `scale = mu`, so the
#' accuracy and the underflow behavior are R's own. Outside the support the
#' density is 0. With `log = TRUE` the logarithm is formed inside `dweibull()`
#' and stays finite far into the upper tail, where the density itself
#' underflows.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of observations. A value at or below zero is
#'   outside the support and gives a density of 0.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive; a zero or negative value gives
#'   `NaN` with a warning from [stats::dweibull()].
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.Weibull1Distrib()] for the distribution function,
#'   [distrib_gradient.Weibull1Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf()] for the generic and
#'   [weibull1_distrib()] for the family.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#'
#' # The method is stats::dweibull with the shape and the scale swapped into
#' # this parametrization's order.
#' all.equal(distrib_pdf(d, y, list(mu = 2, sigma = 1.5)),
#'           dweibull(y, shape = 1.5, scale = 2))
#'
#' # Shape 1 is the exponential with mean mu, and shape 2 the Rayleigh.
#' all.equal(distrib_pdf(d, y, list(mu = 2, sigma = 1)), dexp(y, rate = 1 / 2))
#'
#' # Below the support the density is zero; in the far upper tail it
#' # underflows and its logarithm does not.
#' distrib_pdf(d, c(-1, 0), list(mu = 2, sigma = 1.5))
#' distrib_pdf(d, 60, list(mu = 2, sigma = 1.5))
#' distrib_pdf(d, 60, list(mu = 2, sigma = 1.5), log = TRUE)
S7::method(distrib_pdf, Weibull1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dweibull(x = y, shape = theta[[2]], scale = theta[[1]], log = log)
}

#' @title Weibull Cumulative Distribution Function
#' @name distrib_cdf.Weibull1Distrib
#' @description
#' Computes the Weibull distribution function
#' \deqn{F(q; \mu, \sigma) = 1 - \exp\left\{-(q/\mu)^{\sigma}\right\}}
#' by calling [stats::pweibull()]. The survival function is elementary,
#' \eqn{1 - F(q) = \exp\{-(q/\mu)^{\sigma}\}}, so `lower.tail = FALSE` with
#' `log.p = TRUE` returns \eqn{-(q/\mu)^{\sigma}} exactly and never forms the
#' difference. That combination is what a right-censored observation far out in
#' the tail contributes to a log-likelihood.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param q A numeric vector of quantiles. Below the support the probability is
#'   0 in the lower tail and 1 in the upper.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Weibull1Distrib()] for the inverse,
#'   [distrib_pdf.Weibull1Distrib()] for the density,
#'   [distrib_grad_cdf.Weibull1Distrib()] for the derivatives of this function
#'   in the parameters, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' th <- list(mu = 2, sigma = 1.5)
#'
#' # The method is stats::pweibull at this parametrization.
#' all.equal(distrib_cdf(d, c(0.5, 1.2, 3.0), th),
#'           pweibull(c(0.5, 1.2, 3.0), shape = 1.5, scale = 2))
#'
#' # The log survival function is the exponent itself, with no cancellation.
#' all.equal(distrib_cdf(d, 40, th, lower.tail = FALSE, log.p = TRUE),
#'           -(40 / 2)^1.5)
#'
#' # There the survival probability itself has underflowed to zero.
#' distrib_cdf(d, 40, th, lower.tail = FALSE)
S7::method(distrib_cdf, Weibull1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  stats::pweibull(
    q = q, shape = theta[[2]], scale = theta[[1]],
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @title Weibull Quantile Function
#' @name distrib_quantile.Weibull1Distrib
#' @description
#' Computes the Weibull quantile function
#' \deqn{Q(p; \mu, \sigma) = \mu \left\{-\log(1 - p)\right\}^{1/\sigma}}
#' by calling [stats::qweibull()]. The distribution function is strictly
#' increasing on \eqn{(0, \infty)}, so the inverse is exact and unique; the
#' root-finding fallback the base class supplies is bypassed. `Q(0)` is 0 and
#' `Q(1)` is `Inf`.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN` with
#'   a warning from [stats::qweibull()].
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is read as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty]}, of length
#'   `max(length(p), length(mu), length(sigma))`.
#'
#' @seealso [distrib_cdf.Weibull1Distrib()] for the function inverted here,
#'   [distrib_rng.Weibull1Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' th <- list(mu = 2, sigma = 1.5)
#'
#' # The quartiles, and the round trip back through the distribution function.
#' q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#' q
#' all.equal(distrib_cdf(d, q, th), c(0.25, 0.5, 0.75))
#'
#' # The median is mu times (log 2)^(1/sigma), whatever the shape.
#' all.equal(distrib_quantile(d, 0.5, th), 2 * log(2)^(1 / 1.5))
S7::method(distrib_quantile, Weibull1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  stats::qweibull(
    p = p, shape = theta[[2]], scale = theta[[1]],
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @title Weibull Random Number Generator
#' @name distrib_rng.Weibull1Distrib
#' @description
#' Draws `n` independent Weibull variates by calling [stats::rweibull()], which
#' inverts the distribution function at a uniform draw. The inversion is exact
#' because the quantile function is elementary, so the generalized
#' ratio-of-uniforms fallback the base class supplies is bypassed. The draws
#' depend on `.Random.seed` in the usual way.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. Both
#'   must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` positive draws.
#'
#' @seealso [distrib_quantile.Weibull1Distrib()] for the function inverted,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#'
#' # Same generator as stats::rweibull, so the same seed gives the same draws.
#' set.seed(3)
#' a <- distrib_rng(d, 3, list(mu = 2, sigma = 1.5))
#' set.seed(3)
#' identical(a, rweibull(3, shape = 1.5, scale = 2))
#'
#' # The sample mean recovers mu * gamma(1 + 1/sigma), not mu itself.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, list(mu = 2, sigma = 1.5))
#' c(sample = mean(z), theoretical = 2 * gamma(1 + 1 / 1.5), scale = 2)
S7::method(distrib_rng, Weibull1Distrib) <- function(distrib, n, theta, ...) {
  stats::rweibull(n = n, shape = theta[[2]], scale = theta[[1]])
}

#' @title Weibull Score
#' @name distrib_gradient.Weibull1Distrib
#' @description
#' Computes the first derivatives of the Weibull log-density with respect to
#' the scale \eqn{\mu} and the shape \eqn{\sigma}, one value per observation,
#' in closed form. Writing \eqn{z = y/\mu} and \eqn{u = z^{\sigma}},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1),
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}
#'         + (1 - u)\log z.}
#' Both components are polynomials in \eqn{u} and \eqn{u \log z}, so every
#' higher order is elementary as well.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning, giving \eqn{\partial \ell / \partial \eta_j
#' = h_j'(\eta_j)\, \partial \ell / \partial \theta_j}. This method always
#' returns the parameter scale; the transformation happens in the generic.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of positive observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each of
#'   length `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the scale
#' and \eqn{\sigma > 0} the shape. \eqn{\eta_j} is the coordinate of parameter
#' \eqn{j} on the unconstrained scale of its link, and \eqn{h_j' = \partial
#' \theta_j / \partial \eta_j} the chain-rule factor onto it.
#'
#' @seealso [distrib_hessian.Weibull1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Weibull1Distrib()] for their expectation,
#'   [distrib_grad_y.Weibull1Distrib()] for the derivative in the response,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' u <- (y / 2)^1.5
#' all.equal(g$mu, 1.5 * (u - 1) / 2)
#' all.equal(g$sigma, 1 / 1.5 + (1 - u) * log(y / 2))
#'
#' # The summed score vanishes at the maximum likelihood estimate.
#' set.seed(5)
#' z <- distrib_rng(d, 500, th)
#' mle <- coef(fit_distrib(d, z))
#' vapply(distrib_gradient(d, z, as.list(mle)), sum, numeric(1))
#'
#' # On the link scale both components are multiplied by h' = theta, both
#' # parameters riding a log by default.
#' distrib_gradient(d, y, th, scale = "link")$mu / g$mu
S7::method(distrib_gradient, Weibull1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  p <- weibull_pieces(y, mu, sigma)
  list(
    mu = sigma * (p$u - 1) / mu,
    sigma = 1 / sigma + (1 - p$u) * p$lz
  )
}

#' @title Weibull Observed Hessian
#' @name distrib_hessian.Weibull1Distrib
#' @description
#' Computes the three distinct second derivatives of the Weibull log-density
#' with respect to the scale \eqn{\mu} and the shape \eqn{\sigma}, one value
#' per observation, in closed form. With \eqn{z = y/\mu} and
#' \eqn{u = z^{\sigma}},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2}
#'         = \dfrac{\sigma}{\mu^2}\left\{1 - (1 + \sigma) u\right\},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = -\dfrac{1}{\sigma^2} - u (\log z)^2,}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{1}{\mu}\left(u - 1 + \sigma u \log z\right).}
#' All three carry the data through \eqn{u}, so none is free of \eqn{y} and the
#' observed matrix differs from its expectation at every observation.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of positive observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
#'   The three name the distinct entries of a symmetric \eqn{2 \times 2} matrix
#'   per observation.
#'
#' @seealso [distrib_gradient.Weibull1Distrib()] for the score,
#'   [distrib_expected_hessian.Weibull1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Weibull1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#' h <- distrib_hessian(d, y, th)
#'
#' # The three closed forms, written out.
#' u <- (y / 2)^1.5; lz <- log(y / 2)
#' all.equal(h$mu_mu, 1.5 * (1 - 2.5 * u) / 4)
#' all.equal(h$sigma_sigma, -1 / 1.5^2 - u * lz^2)
#' all.equal(h$mu_sigma, (u - 1 + 1.5 * u * lz) / 2)
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-5
#' up <- distrib_gradient(d, y, list(mu = 2 + eps, sigma = 1.5))$mu
#' dn <- distrib_gradient(d, y, list(mu = 2 - eps, sigma = 1.5))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#'
#' # The curvature in mu is positive wherever u < 1/(1 + sigma), that is
#' # below the 33rd percentile at this shape, so the observed information is
#' # not positive definite at every observation while its expectation is.
#' h$mu_mu
S7::method(distrib_hessian, Weibull1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  p <- weibull_pieces(y, mu, sigma)
  list(
    mu_mu = sigma * (1 - (1 + sigma) * p$u) / mu^2,
    sigma_sigma = -1 / sigma^2 - p$u * p$lz^2,
    mu_sigma = (p$u - 1 + sigma * p$u * p$lz) / mu
  )
}

#' @title Weibull Expected Hessian
#' @name distrib_expected_hessian.Weibull1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'         = -\dfrac{\sigma^2}{\mu^2},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right]
#'         = -\dfrac{1}{\sigma^2}\left\{(1-\gamma)^2 + \dfrac{\pi^2}{6}\right\},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}\right]
#'         = \dfrac{1 - \gamma}{\mu},}
#' with \eqn{\gamma \approx 0.5772} the Euler-Mascheroni constant. The mixed
#' entry does not vanish, so the scale and the shape are not orthogonal and
#' their estimates are asymptotically correlated.
#'
#' Because a closed form exists, `approx` and `nsim` are ignored: every
#' strategy returns the same three numbers.
#'
#' @details
#' # Why the expectations are gamma derivatives
#'
#' Under the model \eqn{u = (Y/\mu)^{\sigma}} is standard exponential whatever
#' the parameters are, so every expectation the Hessian needs is a moment of
#' \eqn{u} against a power of \eqn{\log u}, and each of those is a derivative
#' of \eqn{\Gamma} at 2:
#' \deqn{\mathbb{E}[u] = 1, \qquad
#'       \mathbb{E}[u \log u] = 1 - \gamma, \qquad
#'       \mathbb{E}[u (\log u)^2] = (1-\gamma)^2 + \dfrac{\pi^2}{6} - 1.}
#' Substituting \eqn{\log z = (\log u)/\sigma} into the three observed
#' components gives the display above. The Gumbel family shares the
#' substitution, since \eqn{\exp(-\text{Gumbel})} is Weibull, so
#' [distrib_expected_hessian.GumbelDistrib()] rests on the same three moments.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of observations. Only its length is read, the
#'   expectation not depending on the data; the values are ignored.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being exact. Accepted so that
#'   the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length `length(y)` and each constant along it.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the scale,
#' \eqn{\sigma > 0} the shape and \eqn{\gamma} the Euler-Mascheroni constant,
#' `-digamma(1)`.
#'
#' @seealso [distrib_hessian.Weibull1Distrib()] for the quantity this is the
#'   expectation of, [distrib_gradient.Weibull1Distrib()] for the score, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#' eh <- distrib_expected_hessian(d, y, th)
#' vapply(eh, function(z) z[1], numeric(1))
#'
#' # The closed forms, written out.
#' eg <- -digamma(1)
#' c(mu_mu = -1.5^2 / 4,
#'   sigma_sigma = -((1 - eg)^2 + pi^2 / 6) / 1.5^2,
#'   mu_sigma = (1 - eg) / 2)
#'
#' # Averaging the observed Hessian over draws reaches the same three numbers.
#' set.seed(4)
#' z <- distrib_rng(d, 2e5, th)
#' vapply(distrib_hessian(d, z, th), mean, numeric(1))
#'
#' # The strategy argument is inert, the expectation being exact.
#' identical(eh, distrib_expected_hessian(d, y, th, approx = "mc", nsim = 50))
S7::method(distrib_expected_hessian, Weibull1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  n <- length(y)
  eg <- -digamma(1)   # the Euler-Mascheroni constant
  list(
    mu_mu = rep(-sigma^2 / mu^2, length.out = n),
    sigma_sigma = rep(-((1 - eg)^2 + pi^2 / 6) / sigma^2, length.out = n),
    mu_sigma = rep((1 - eg) / mu, length.out = n)
  )
}

#' @title Weibull Third-Order Derivatives
#' @name distrib_deriv3.Weibull1Distrib
#' @description
#' Computes the four distinct third derivatives of the Weibull log-density in
#' \eqn{\mu} and \eqn{\sigma}, in closed form, observed at the data or expected
#' under the model. With \eqn{z = y/\mu}, \eqn{u = z^{\sigma}} and
#' \eqn{L = \log z}, each component is a polynomial in \eqn{u}, \eqn{uL},
#' \eqn{uL^2} and \eqn{uL^3} with coefficients rational in \eqn{\mu} and
#' \eqn{\sigma}; the pure-\eqn{\mu} component, for instance, is
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3}
#'   = \dfrac{\sigma}{\mu^3}\left\{-2 + (1+\sigma)(2+\sigma) u\right\}.}
#' The expected values replace each \eqn{\mathbb{E}[u L^k]} by
#' \eqn{\Gamma^{(k)}(2)/\sigma^k}, so `expected = TRUE` is exact here and needs
#' no quadrature. The arithmetic runs in a compiled kernel decomposed over the
#' elements of the output, so the result does not depend on the thread count.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of positive observations. With `expected = TRUE`
#'   only its length is read.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, both branches being exact. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`. The four name the distinct
#'   entries of a symmetric third-order array over two parameters.
#'
#' @seealso [distrib_hessian.Weibull1Distrib()] for the order below,
#'   [distrib_deriv4.Weibull1Distrib()] for the order above, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # The pure-mu component, written out.
#' u <- (y / 2)^1.5
#' all.equal(d3$mu_mu_mu, 1.5 * (-2 + 2.5 * 3.5 * u) / 8)
#'
#' # A central difference of the Hessian reproduces the mixed component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 2 + eps, sigma = 1.5))$mu_sigma
#' dn <- distrib_hessian(d, y, list(mu = 2 - eps, sigma = 1.5))$mu_sigma
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_sigma, tolerance = 1e-6)
#'
#' # The expected branch is exact, and averaging the observed one over draws
#' # reaches it.
#' set.seed(6)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(expected = vapply(distrib_deriv3(d, y, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       averaged = vapply(distrib_deriv3(d, z, th), mean, numeric(1)))
S7::method(distrib_deriv3, Weibull1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) weibull_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else weibull_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Weibull Fourth-Order Derivatives
#' @name distrib_deriv4.Weibull1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the Weibull log-density in
#' \eqn{\mu} and \eqn{\sigma}, in closed form, observed at the data or expected
#' under the model. The construction is the one
#' [distrib_deriv3.Weibull1Distrib()] describes carried one order further: with
#' \eqn{z = y/\mu}, \eqn{u = z^{\sigma}} and \eqn{L = \log z}, each component
#' is a polynomial in \eqn{u L^k} for \eqn{k \le 4}, and the pure-\eqn{\mu}
#' component is
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4}
#'   = \dfrac{\sigma}{\mu^4}\left\{6 - (1+\sigma)(2+\sigma)(3+\sigma) u\right\}.}
#' The expected values replace each \eqn{\mathbb{E}[u L^k]} by
#' \eqn{\Gamma^{(k)}(2)/\sigma^k}, so `expected = TRUE` is exact and needs no
#' quadrature. The arithmetic runs in a compiled kernel decomposed over the
#' elements of the output, so the result does not depend on the thread count.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of positive observations. With `expected = TRUE`
#'   only its length is read.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, both branches being exact. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_sigma`, `mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
#'   `sigma_sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @seealso [distrib_deriv3.Weibull1Distrib()] for the order below,
#'   [distrib_hessian.Weibull1Distrib()] for the second order, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The pure-mu component, written out.
#' u <- (y / 2)^1.5
#' all.equal(d4$mu_mu_mu_mu, 1.5 * (6 - 2.5 * 3.5 * 4.5 * u) / 16)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-4
#' up <- distrib_deriv3(d, y, list(mu = 2 + eps, sigma = 1.5))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 2 - eps, sigma = 1.5))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#'
#' # The expected branch is exact, and averaging the observed one over draws
#' # reaches it.
#' set.seed(8)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(expected = vapply(distrib_deriv4(d, y, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       averaged = vapply(distrib_deriv4(d, z, th), mean, numeric(1)))
S7::method(distrib_deriv4, Weibull1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) weibull_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else weibull_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Weibull First Derivative in the Response
#' @name distrib_grad_y.Weibull1Distrib
#' @description
#' Computes \eqn{\partial \ell / \partial y}, the derivative of the Weibull
#' log-density with respect to the response, in closed form. With
#' \eqn{u = (y/\mu)^{\sigma}},
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\sigma - 1 - \sigma u}{y}.}
#' At \eqn{\sigma > 1} it vanishes at the mode
#' \eqn{y = \mu\{(\sigma-1)/\sigma\}^{1/\sigma}}; at \eqn{\sigma \le 1} the
#' density is decreasing on the whole support and the derivative is negative
#' everywhere. This quantity is what a quantile residual's delta-method
#' standard error and a change of variable in the response both need.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of positive observations. At or below zero the
#'   result is not defined and propagates as `NaN` or `Inf`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Weibull1Distrib()] for the second derivative in the
#'   response, [distrib_cross_y.Weibull1Distrib()] for the mixed derivative in
#'   the response and the parameters,
#'   [distrib_gradient.Weibull1Distrib()] for the derivatives in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#'
#' # The closed form, written out.
#' u <- (y / 2)^1.5
#' all.equal(distrib_grad_y(d, y, th), (1.5 - 1 - 1.5 * u) / y)
#'
#' # It is the derivative of the log-density, so a central difference of the
#' # log-density in y reproduces it.
#' eps <- 1e-6
#' all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
#'            distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
#'           distrib_grad_y(d, y, th), tolerance = 1e-6)
#'
#' # It vanishes at the mode, which exists because the shape exceeds one.
#' mode <- 2 * ((1.5 - 1) / 1.5)^(1 / 1.5)
#' c(mode = mode, deriv = distrib_grad_y(d, mode, th))
S7::method(distrib_grad_y, Weibull1Distrib) <- function(distrib, y, theta, ...) {
  sigma <- theta[[2]]
  p <- weibull_pieces(y, theta[[1]], sigma)
  (sigma - 1 - sigma * p$u) / y
}

#' @title Weibull Second Derivative in the Response
#' @name distrib_hess_y.Weibull1Distrib
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2}, the second derivative of the
#' Weibull log-density with respect to the response, in closed form. With
#' \eqn{u = (y/\mu)^{\sigma}},
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2}
#'   = -\dfrac{(\sigma - 1)(1 + \sigma u)}{y^2}.}
#' The sign is the shape's: at \eqn{\sigma > 1} the log-density is concave on
#' the whole support, at \eqn{\sigma < 1} convex, and at \eqn{\sigma = 1}, the
#' exponential case, exactly zero, the log-density being linear in \eqn{y}
#' there.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param y A numeric vector of positive observations. At or below zero the
#'   result is not defined and propagates as `NaN` or `Inf`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_grad_y.Weibull1Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Weibull1Distrib()] for the second derivatives
#'   in the parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- weibull1_distrib()
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#'
#' # The closed form, written out.
#' u <- (y / 2)^1.5
#' all.equal(distrib_hess_y(d, y, th), -(1.5 - 1) * (1 + 1.5 * u) / y^2)
#'
#' # A central difference of the first derivative reproduces it.
#' eps <- 1e-5
#' all.equal((distrib_grad_y(d, y + eps, th) -
#'            distrib_grad_y(d, y - eps, th)) / (2 * eps),
#'           distrib_hess_y(d, y, th), tolerance = 1e-6)
#'
#' # At shape one the log-density is linear in y, so the curvature is exactly
#' # zero; at a smaller shape it turns positive.
#' rbind(shape_1.0 = distrib_hess_y(d, y, list(mu = 2, sigma = 1)),
#'       shape_0.6 = distrib_hess_y(d, y, list(mu = 2, sigma = 0.6)))
S7::method(distrib_hess_y, Weibull1Distrib) <- function(distrib, y, theta, ...) {
  sigma <- theta[[2]]
  p <- weibull_pieces(y, theta[[1]], sigma)
  -(sigma - 1) * (1 + sigma * p$u) / y^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Weibull Distribution, Scale and Shape
#'
#' @description
#' Builds the distribution object for the Weibull family parametrized by a
#' scale \eqn{\mu > 0} and a shape \eqn{\sigma > 0}, with density
#' \eqn{f(y) = (\sigma/\mu)(y/\mu)^{\sigma-1}\exp\{-(y/\mu)^{\sigma}\}} on
#' \eqn{y > 0}. The returned object carries closed-form derivatives of the
#' log-density to fourth order, observed and expected, in the parameters and in
#' the response, and closed-form moments, so every generic of the toolkit
#' answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. Both default to the logarithm,
#' both parameters being positive.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the scale
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#' @param link_sigma A `link` object from `linkfunctions7` for the shape
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, \infty)} is
#' \deqn{f(y; \mu, \sigma) = \dfrac{\sigma}{\mu}
#'   \left(\dfrac{y}{\mu}\right)^{\sigma - 1}
#'   \exp\left\{-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\},}
#' the distribution function \eqn{F(q) = 1 - \exp\{-(q/\mu)^{\sigma}\}} and the
#' quantile function \eqn{Q(p) = \mu\{-\log(1-p)\}^{1/\sigma}}. This is `WEI` in
#' \pkg{gamlss}.
#'
#' **\eqn{\mu} is the scale and not the mean.** The mean is
#' \eqn{\mu\,\Gamma(1 + 1/\sigma)}, which involves the shape, so a mean
#' parametrization would make every derivative a derivative of the gamma
#' function and of its inverse. The scale-shape form keeps the whole family
#' elementary; [mean.Weibull1Distrib()] reports the mean, and
#' [weibull3_distrib()] is the parametrization by a quantile for a reader who
#' wants a location that is one.
#'
#' # Derivatives
#'
#' Write \eqn{z = y/\mu} and \eqn{u = z^{\sigma}}. The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1),
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}
#'         + (1 - u)\log z,}
#' the observed Hessian
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2}
#'         = \dfrac{\sigma}{\mu^2}\left\{1 - (1 + \sigma) u\right\}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = -\dfrac{1}{\sigma^2} - u (\log z)^2, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{u - 1 + \sigma u \log z}{\mu},}
#' and its expectation
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'         = -\dfrac{\sigma^2}{\mu^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right]
#'         = -\dfrac{(1-\gamma)^2 + \pi^2/6}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right]
#'         = \dfrac{1 - \gamma}{\mu},}
#' with \eqn{\gamma} the Euler-Mascheroni constant. The mixed entry is not
#' zero, so the scale and the shape are asymptotically correlated, unlike the
#' mean and the standard deviation of a Gaussian.
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.Weibull1Distrib()] and [distrib_deriv4.Weibull1Distrib()],
#' as are the derivatives in the response,
#' [distrib_grad_y.Weibull1Distrib()] and [distrib_hess_y.Weibull1Distrib()],
#' and the mixed derivative [distrib_cross_y.Weibull1Distrib()].
#'
#' # Why every expectation is elementary
#'
#' Under the model \eqn{u = (Y/\mu)^{\sigma}} is standard exponential whatever
#' the parameters are. Every expectation any order needs is therefore a moment
#' of \eqn{u} against a power of \eqn{\log u}, and each of those is a
#' derivative of \eqn{\Gamma} at 2: \eqn{\mathbb{E}[u] = 1},
#' \eqn{\mathbb{E}[u\log u] = 1 - \gamma},
#' \eqn{\mathbb{E}[u(\log u)^2] = (1-\gamma)^2 + \pi^2/6 - 1}. The Gumbel
#' family shares the substitution, \eqn{\exp(-\text{Gumbel})} being Weibull, so
#' [gumbel_distrib()] rests on the same three moments.
#'
#' # Moments and shape
#'
#' With \eqn{g_k = \Gamma(1 + k/\sigma)}, the mean is \eqn{\mu g_1} and the
#' variance \eqn{\mu^2(g_2 - g_1^2)}; the skewness and the excess kurtosis
#' follow from \eqn{g_3} and \eqn{g_4} and do not depend on \eqn{\mu}, the
#' scale being a multiplier. The hazard \eqn{f/(1-F)} is
#' \eqn{(\sigma/\mu)(y/\mu)^{\sigma-1}}, increasing for \eqn{\sigma > 1} and
#' decreasing for \eqn{\sigma < 1}. That monotone hazard is why the family is
#' used in survival work.
#'
#' Two shapes are named families: \eqn{\sigma = 1} is the exponential with mean
#' \eqn{\mu}, and \eqn{\sigma = 2} the Rayleigh. Both are exact, not limiting,
#' and the example checks the first.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' estimate is closed form: the shape solves
#' \eqn{1/\hat\sigma + \overline{\log y} = \sum y_i^{\hat\sigma}\log y_i /
#' \sum y_i^{\hat\sigma}} by iteration, and the scale follows from it as
#' \eqn{\hat\mu = (n^{-1}\sum y_i^{\hat\sigma})^{1/\hat\sigma}}. The example
#' checks that the fit satisfies the second of these.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the scale,
#' \eqn{\sigma > 0} the shape and \eqn{\gamma} the Euler-Mascheroni constant.
#' \eqn{\eta} is a parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Weibull1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"weibull1"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "sigma")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 21.
#' Wiley, New York.
#'
#' Rigby, R. A. and Stasinopoulos, D. M. (2005).
#' Generalized additive models for location, scale and shape.
#' *Journal of the Royal Statistical Society, Series C*, **54**(3), 507-554.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dweibull pweibull qweibull rweibull
#'
#' @examples
#' d <- weibull1_distrib()
#' d
#'
#' # The density is stats::dweibull with the arguments in this order.
#' y <- c(0.5, 1.2, 3.0)
#' th <- list(mu = 2, sigma = 1.5)
#' all.equal(distrib_pdf(d, y, th), dweibull(y, shape = 1.5, scale = 2))
#'
#' # The scale is not the mean: the mean carries a gamma function of the shape.
#' c(scale = th$mu, mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # Shape one is exactly the exponential with mean mu.
#' all.equal(distrib_pdf(d, y, list(mu = 2, sigma = 1)), dexp(y, rate = 1 / 2))
#'
#' # The hazard increases with y above shape one and decreases below it.
#' haz <- function(s) {
#'   th <- list(mu = 2, sigma = s)
#'   distrib_pdf(d, c(0.5, 1, 2, 4), th) /
#'     distrib_cdf(d, c(0.5, 1, 2, 4), th, lower.tail = FALSE)
#' }
#' rbind(shape_0.5 = haz(0.5), shape_2.5 = haz(2.5))
#'
#' # Fitting recovers the parameters, and the scale satisfies the profile
#' # identity given the fitted shape.
#' set.seed(9)
#' z <- distrib_rng(d, 2000, list(mu = 2, sigma = 1.5))
#' fit <- fit_distrib(d, z)
#' cf <- coef(fit)
#' c(cf, profile_mu = mean(z^cf[["sigma"]])^(1 / cf[["sigma"]]))
#'
#' @seealso
#' [weibull3_distrib()] for the same law parametrized by a quantile;
#' [gumbel_distrib()], which is this family after a negative logarithm;
#' [exponential_distrib()] and [gengamma1_distrib()] for the special case and
#' the generalization; [gpd_distrib()] for another extreme-value family;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [Weibull1Distrib] for the class.
#' @export
weibull1_distrib <- function(link_mu = log_link(), link_sigma = log_link()) {
  Weibull1Distrib(
    distrib_name = "weibull1",
    dimension = "univariate",
    bounds = c(0, Inf),

    params = c("mu", "sigma"),
    params_interpretation = c(mu = "scale", sigma = "shape"),
    n_params = 2,

    params_bounds = list(
      mu = c(0, Inf),
      sigma = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma
    )
  )
}
