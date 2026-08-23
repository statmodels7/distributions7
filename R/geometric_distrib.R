#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Geometric Distribution
#' @name GeometricDistrib
#'
#' @description A subclass of `discrete_distrib` representing the
#'   geometric distribution on \eqn{\{0, 1, 2, \dots\}} in its mean
#'   parametrization.
#' @inheritParams distrib
#' @return An object of class `GeometricDistrib`.
#' @seealso [geometric_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.GeometricDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.GeometricDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.GeometricDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GeometricDistrib],
#'   [`distrib_gradient()`][distrib_gradient.GeometricDistrib],
#'   [`distrib_hessian()`][distrib_hessian.GeometricDistrib],
#'   [`distrib_pdf()`][distrib_pdf.GeometricDistrib],
#'   [`distrib_quantile()`][distrib_quantile.GeometricDistrib],
#'   [`distrib_rng()`][distrib_rng.GeometricDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
GeometricDistrib <- S7::new_class("GeometricDistrib", parent = discrete_distrib)

#' The Success Probability Behind a Geometric Mean
#'
#' @description
#' Converts the mean into the success probability the base R functions take.
#'
#' @details
#' The mean number of failures before the first success is
#' \eqn{(1-p)/p}, so \eqn{p = 1/(1+\mu)}. Writing it once keeps the four
#' functions that call \pkg{stats} from each repeating the algebra.
#'
#' @param mu The mean, a positive numeric vector.
#'
#' @return A numeric vector of probabilities in \eqn{(0, 1)}.
#'
#' @seealso [geometric_distrib()]
#'
#' @keywords internal
geom_prob <- function(mu) 1 / (1 + mu)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Geometric Probability Mass Function
#' @name distrib_pdf.GeometricDistrib
#' @description
#' \deqn{P(Y = y; \mu) = \dfrac{1}{1+\mu}
#'       \left(\dfrac{\mu}{1+\mu}\right)^{y}, \qquad y = 0, 1, 2, \dots}
#' @param distrib A `GeometricDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param log Logical; if `TRUE`, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso [geometric_distrib()]
S7::method(distrib_pdf, GeometricDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dgeom(y, prob = geom_prob(theta[[1]]), log = log)
}

#' @title Geometric Cumulative Distribution Function
#' @name distrib_cdf.GeometricDistrib
#' @description
#' \deqn{F(q; \mu) = 1 - \left(\dfrac{\mu}{1+\mu}\right)^{\lfloor q \rfloor + 1}}
#' @param distrib A `GeometricDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter `mu`.
#' @param lower.tail Logical; if `TRUE` (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [geometric_distrib()]
S7::method(distrib_cdf, GeometricDistrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE) {
  stats::pgeom(q, prob = geom_prob(theta[[1]]), lower.tail = lower.tail,
               log.p = log.p)
}

#' @title Geometric Quantile Function
#' @name distrib_quantile.GeometricDistrib
#' @description The generalized inverse
#'   \eqn{Q(p) = \min\{y \in \mathbb{N}_0 : F(y) \ge p\}}.
#' @param distrib A `GeometricDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter `mu`.
#' @param lower.tail Logical; if `TRUE` (default), `p` is \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, `p` is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso [geometric_distrib()]
S7::method(distrib_quantile, GeometricDistrib) <- function(distrib, p, theta,
                                                           lower.tail = TRUE,
                                                           log.p = FALSE) {
  stats::qgeom(p, prob = geom_prob(theta[[1]]), lower.tail = lower.tail,
               log.p = log.p)
}

#' @title Geometric Random Generation
#' @name distrib_rng.GeometricDistrib
#' @description Draws through [stats::rgeom()] at
#'   \eqn{p = 1/(1+\mu)}.
#' @param distrib A `GeometricDistrib` object.
#' @param n The number of draws.
#' @param theta A list containing the parameter `mu`.
#' @return A numeric vector of length `n`.
#' @seealso [geometric_distrib()]
S7::method(distrib_rng, GeometricDistrib) <- function(distrib, n, theta) {
  stats::rgeom(n, prob = geom_prob(theta[[1]]))
}

#' @title Geometric Analytical Gradient
#' @name distrib_gradient.GeometricDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1+\mu)}}
#' the deviation from the mean over the variance.
#' @param distrib A `GeometricDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list with the `mu` component.
#' @seealso [geometric_distrib()]
S7::method(distrib_gradient, GeometricDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  geometric_gradient_cpp(y, theta[[1]], threads)
}

#' @title Geometric Analytical Observed Hessian
#' @name distrib_hessian.GeometricDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} =
#'       -\dfrac{y}{\mu^2} + \dfrac{y+1}{(1+\mu)^2}}
#' @param distrib A `GeometricDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list with the `mu_mu` component.
#' @seealso [geometric_distrib()]
S7::method(distrib_hessian, GeometricDistrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ..., threads = 1L) {
  geometric_hessian_cpp(y, theta[[1]], threads)
}

#' @title Geometric Analytical Expected Hessian
#' @name distrib_expected_hessian.GeometricDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'       = -\dfrac{1}{\mu(1+\mu)}}
#' the reciprocal of the variance, as it must be for a family written in its
#' mean.
#' @param distrib A `GeometricDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the `mu_mu` component.
#' @seealso [geometric_distrib()]
S7::method(distrib_expected_hessian, GeometricDistrib) <- function(distrib, y, theta,
                                                                   scale = c("parameter", "link"),
                                                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                                                   nsim = 10000, ..., threads = 1L) {
  geometric_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Geometric Analytical Third-Order Derivative
#' @name distrib_deriv3.GeometricDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu)} = 2\left(\dfrac{y}{\mu^3}
#'       - \dfrac{y+1}{(1+\mu)^3}\right), \qquad
#'       \mathbb{E}[\ell^{(\mu\mu\mu)}] =
#'       2\left(\dfrac{1}{\mu^2} - \dfrac{1}{(1+\mu)^2}\right)}
#' @param distrib A `GeometricDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param expected Logical; if `TRUE`, returns the expected derivative.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the `mu_mu_mu` component.
#' @seealso [geometric_distrib()]
S7::method(distrib_deriv3, GeometricDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ..., threads = 1L) {
  if (expected) geometric_deriv3_expected_cpp(y, theta[[1]], threads)
  else geometric_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Geometric Analytical Fourth-Order Derivative
#' @name distrib_deriv4.GeometricDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu\mu)} = -6\left(\dfrac{y}{\mu^4}
#'       - \dfrac{y+1}{(1+\mu)^4}\right), \qquad
#'       \mathbb{E}[\ell^{(\mu\mu\mu\mu)}] =
#'       -6\left(\dfrac{1}{\mu^3} - \dfrac{1}{(1+\mu)^3}\right)}
#' @param distrib A `GeometricDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param expected Logical; if `TRUE`, returns the expected derivative.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the `mu_mu_mu_mu` component.
#' @seealso [geometric_distrib()]
S7::method(distrib_deriv4, GeometricDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ..., threads = 1L) {
  if (expected) geometric_deriv4_expected_cpp(y, theta[[1]], threads)
  else geometric_deriv4_cpp(y, theta[[1]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Geometric Distribution Object
#'
#' @description
#' Creates a distribution object for the geometric distribution on
#' \eqn{\{0, 1, 2, \dots\}}, parametrized by its mean \eqn{\mu}.
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The distribution counts the failures before the first success, so the
#' support includes zero and the success probability is \eqn{p = 1/(1+\mu)}.
#'
#' **Probability mass function:**
#' \deqn{P(Y = y; \mu) = \dfrac{1}{1+\mu}\left(\dfrac{\mu}{1+\mu}\right)^{y}}
#'
#' **Score, observed and expected Hessian:**
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\mu(1+\mu)},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}
#'       + \dfrac{y+1}{(1+\mu)^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'       = -\dfrac{1}{\mu(1+\mu)}}
#'
#' Every order is the same expression evaluated at \eqn{\mu} and at
#' \eqn{1+\mu},
#' \deqn{\ell^{(j)} = (-1)^{j-1}(j-1)!\left(\dfrac{y}{\mu^{j}}
#'       - \dfrac{y+1}{(1+\mu)^{j}}\right), \qquad
#'       \mathbb{E}[\ell^{(j)}] = (-1)^{j-1}(j-1)!\left(\mu^{1-j}
#'       - (1+\mu)^{1-j}\right)}
#' so the expected orders are closed form and vanish at \eqn{j = 1}.
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{\mu(1+\mu)}, so the family
#' is overdispersed relative to the Poisson at every mean.
#'
#' **Parameter domains:**
#'
#' - \eqn{\mu \in (0, +\infty)}
#'
#' The family is the negative binomial at \eqn{\theta = 1}, so
#' `fixed(negbin2_distrib(), theta = 1)` describes the same law and is used
#' in the tests as an independent implementation.
#'
#' @return An S7 object of class `GeometricDistrib`.
#'
#' @seealso [negbin2_distrib()], [poisson_distrib()],
#'   [exponential_distrib()]
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dgeom pgeom qgeom rgeom
#' @examples
#' d <- geometric_distrib()
#' d@params
#'
#' theta <- list(mu = 2)
#' distrib_pdf(d, 0:4, theta)
#' distrib_gradient(d, 0:4, theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' @export
geometric_distrib <- function(link_mu = log_link()) {
  GeometricDistrib(
    distrib_name = "geometric", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = link_mu)
  )
}
