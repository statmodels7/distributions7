#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Exponential Distribution
#' @name ExponentialDistrib
#'
#' @description A subclass of `continuous_distrib` representing the
#'   exponential distribution in its mean parametrization.
#' @inheritParams distrib
#' @return An object of class `ExponentialDistrib`.
#' @seealso [exponential_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.ExponentialDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.ExponentialDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.ExponentialDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ExponentialDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ExponentialDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.ExponentialDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ExponentialDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.ExponentialDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ExponentialDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ExponentialDistrib],
#'   [`distrib_rng()`][distrib_rng.ExponentialDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
ExponentialDistrib <- S7::new_class("ExponentialDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Exponential Density
#' @name distrib_pdf.ExponentialDistrib
#' @description
#' \deqn{f(y; \mu) = \dfrac{1}{\mu} e^{-y/\mu}, \qquad y > 0}
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [exponential_distrib()]
S7::method(distrib_pdf, ExponentialDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dexp(y, rate = 1 / theta[[1]], log = log)
}

#' @title Exponential Distribution Function
#' @name distrib_cdf.ExponentialDistrib
#' @description
#' \deqn{F(q; \mu) = 1 - e^{-q/\mu}}
#' @param distrib An `ExponentialDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter `mu`.
#' @param lower.tail Logical; if `TRUE` (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [exponential_distrib()]
S7::method(distrib_cdf, ExponentialDistrib) <- function(distrib, q, theta,
                                                        lower.tail = TRUE,
                                                        log.p = FALSE) {
  stats::pexp(q, rate = 1 / theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Exponential Quantile Function
#' @name distrib_quantile.ExponentialDistrib
#' @description
#' \deqn{Q(p; \mu) = -\mu \log(1 - p)}
#' @param distrib An `ExponentialDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter `mu`.
#' @param lower.tail Logical; if `TRUE` (default), `p` is \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, `p` is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso [exponential_distrib()]
S7::method(distrib_quantile, ExponentialDistrib) <- function(distrib, p, theta,
                                                             lower.tail = TRUE,
                                                             log.p = FALSE) {
  stats::qexp(p, rate = 1 / theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Exponential Random Generation
#' @name distrib_rng.ExponentialDistrib
#' @description Draws from the exponential distribution through
#'   [stats::rexp()].
#' @param distrib An `ExponentialDistrib` object.
#' @param n The number of draws.
#' @param theta A list containing the parameter `mu`.
#' @return A numeric vector of length `n`.
#' @seealso [exponential_distrib()]
S7::method(distrib_rng, ExponentialDistrib) <- function(distrib, n, theta) {
  stats::rexp(n, rate = 1 / theta[[1]])
}

#' @title Exponential Analytical Gradient
#' @name distrib_gradient.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}}
#' the score of a one-parameter family written as the deviation from the mean
#' over the variance.
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list with the `mu` component.
#' @seealso [exponential_distrib()]
S7::method(distrib_gradient, ExponentialDistrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ..., threads = 1L) {
  exponential_gradient_cpp(y, theta[[1]], threads)
}

#' @title Exponential Analytical Observed Hessian
#' @name distrib_hessian.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} - \dfrac{2y}{\mu^3}}
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list with the `mu_mu` component.
#' @seealso [exponential_distrib()]
S7::method(distrib_hessian, ExponentialDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  exponential_hessian_cpp(y, theta[[1]], threads)
}

#' @title Exponential Analytical Expected Hessian
#' @name distrib_expected_hessian.ExponentialDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu^2}}
#' obtained from the observed form by \eqn{\mathbb{E}[y] = \mu}.
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the `mu_mu` component.
#' @seealso [exponential_distrib()]
S7::method(distrib_expected_hessian, ExponentialDistrib) <- function(distrib, y, theta,
                                                                     scale = c("parameter", "link"),
                                                                     approx = c("bartlett", "integrate", "mc", "opg"),
                                                                     nsim = 10000, ..., threads = 1L) {
  exponential_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Exponential Analytical Third-Order Derivative
#' @name distrib_deriv3.ExponentialDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu)} = -\dfrac{2}{\mu^3} + \dfrac{6y}{\mu^4},
#'       \qquad \mathbb{E}[\ell^{(\mu\mu\mu)}] = \dfrac{4}{\mu^3}}
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param expected Logical; if `TRUE`, returns the expected derivative.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the `mu_mu_mu` component.
#' @seealso [exponential_distrib()]
S7::method(distrib_deriv3, ExponentialDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ..., threads = 1L) {
  if (expected) exponential_deriv3_expected_cpp(y, theta[[1]], threads)
  else exponential_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Exponential Analytical Fourth-Order Derivative
#' @name distrib_deriv4.ExponentialDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu\mu)} = \dfrac{6}{\mu^4} - \dfrac{24y}{\mu^5},
#'       \qquad \mathbb{E}[\ell^{(\mu\mu\mu\mu)}] = -\dfrac{18}{\mu^4}}
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param expected Logical; if `TRUE`, returns the expected derivative.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the `mu_mu_mu_mu` component.
#' @seealso [exponential_distrib()]
S7::method(distrib_deriv4, ExponentialDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ..., threads = 1L) {
  if (expected) exponential_deriv4_expected_cpp(y, theta[[1]], threads)
  else exponential_deriv4_cpp(y, theta[[1]], threads)
}

#' @title Exponential Response Gradient
#' @name distrib_grad_y.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{1}{\mu}}
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [exponential_distrib()]
S7::method(distrib_grad_y, ExponentialDistrib) <- function(distrib, y, theta, ...) {
  rep(-1 / theta[[1]], length.out = length(y))
}

#' @title Exponential Response Hessian
#' @name distrib_hess_y.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = 0}
#' the log-density being linear in the response.
#' @param distrib An `ExponentialDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param ... Unused.
#' @return A numeric vector of zeros.
#' @seealso [exponential_distrib()]
S7::method(distrib_hess_y, ExponentialDistrib) <- function(distrib, y, theta, ...) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Exponential Distribution Object
#'
#' @description
#' Creates a distribution object for the exponential distribution parametrized
#' by its mean \eqn{\mu}.
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The exponential distribution models waiting times on \eqn{y > 0}, with mean
#' \eqn{\mu} and variance \eqn{\mu^2}: the coefficient of variation is one, and
#' fixing it is what distinguishes the family from the Gamma.
#'
#' **Density:** \deqn{f(y; \mu) = \dfrac{1}{\mu} e^{-y/\mu}}
#'
#' **Distribution function:** \deqn{F(q; \mu) = 1 - e^{-q/\mu}}
#'
#' **Score, observed and expected Hessian:**
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} - \dfrac{2y}{\mu^3},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu^2}}
#'
#' Every order follows the same pattern, the log-density being a logarithm plus
#' a reciprocal:
#' \deqn{\ell^{(k)} = \dfrac{(-1)^k (k-1)!}{\mu^k}
#'       + \dfrac{(-1)^{k+1} k!\, y}{\mu^{k+1}}, \qquad
#'       \mathbb{E}[\ell^{(k)}] = \dfrac{(-1)^k (k-1)! (1-k)}{\mu^k}}
#' so the expected orders are closed form as well, and vanish at \eqn{k = 1}.
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{\mu^2}, skewness 2, excess
#' kurtosis 6.
#'
#' **Parameter domains:**
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#' }
#'
#' The family is the Weibull with unit shape, so
#' `fixed(weibull1_distrib(), sigma = 1)` describes the same law and is
#' used in the tests as an independent implementation. It is **not** a
#' Gamma with a fixed parameter: this package writes the Gamma in
#' \eqn{(\mu, \sigma^2)}, whose shape is \eqn{\mu^2/\sigma^2}, so unit shape is
#' the relation \eqn{\sigma^2 = \mu^2} between two parameters rather than a
#' value one of them can be held at.
#'
#' @return An S7 object of class `ExponentialDistrib`.
#'
#' @seealso [gamma2_distrib()], [weibull1_distrib()],
#'   [geometric_distrib()]
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dexp pexp qexp rexp
#' @examples
#' d <- exponential_distrib()
#' d@params
#'
#' theta <- list(mu = 2)
#' distrib_pdf(d, c(0.5, 1, 3), theta)
#' distrib_gradient(d, c(0.5, 1, 3), theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' @export
exponential_distrib <- function(link_mu = log_link()) {
  ExponentialDistrib(
    distrib_name = "exponential", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = link_mu)
  )
}
