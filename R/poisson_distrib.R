#' @include distrib.R generics.R
NULL

#' @title S7 Class for Poisson Distribution
#' @name PoissonDistrib
#'
#' @description A subclass of `discrete_distrib` representing the Poisson distribution.
#' @inheritParams distrib
#' @return An object of class `PoissonDistrib`.
#' @seealso [poisson_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.PoissonDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.PoissonDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.PoissonDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.PoissonDistrib],
#'   [`distrib_gradient()`][distrib_gradient.PoissonDistrib],
#'   [`distrib_hessian()`][distrib_hessian.PoissonDistrib],
#'   [`distrib_pdf()`][distrib_pdf.PoissonDistrib],
#'   [`distrib_quantile()`][distrib_quantile.PoissonDistrib],
#'   [`distrib_rng()`][distrib_rng.PoissonDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
PoissonDistrib <- S7::new_class("PoissonDistrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Poisson Probability Mass Function
#' @name distrib_pdf.PoissonDistrib
#' @description
#' Computes the probability mass function for the Poisson distribution:
#' \deqn{P(Y=y; \mu) = \dfrac{\mu^y e^{-\mu}}{y!}}
#'
#' @param distrib A `PoissonDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param log Logical; if `TRUE`, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso [poisson_distrib()]
S7::method(distrib_pdf, PoissonDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dpois(
    x = y,
    lambda = theta[[1]],
    log = log
  )
}

#' @title Poisson Cumulative Distribution Function
#' @name distrib_cdf.PoissonDistrib
#' @description
#' Computes the cumulative distribution function for the Poisson distribution:
#' \deqn{F(q; \mu) = \sum_{k=0}^{\lfloor q \rfloor} \dfrac{\mu^k e^{-\mu}}{k!}}
#'
#' @param distrib A `PoissonDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter `mu`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [poisson_distrib()]
S7::method(distrib_cdf, PoissonDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::ppois(
    q = q,
    lambda = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Poisson Quantile Function
#' @name distrib_quantile.PoissonDistrib
#' @description
#' Computes the quantile function for the Poisson distribution, defined as the
#' generalized inverse of the CDF:
#' \deqn{Q(p; \mu) = \min\left\{y \in \mathbb{N}_0 : F(y; \mu) \ge p\right\}}
#'
#' @param distrib A `PoissonDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter `mu`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso [poisson_distrib()]
S7::method(distrib_quantile, PoissonDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qpois(
    p = p,
    lambda = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Poisson Random Number Generator
#' @name distrib_rng.PoissonDistrib
#' @description
#' Generates random numbers from the Poisson distribution.
#'
#' @param distrib A `PoissonDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameter `mu`.
#' @return A numeric vector of random draws.
#' @seealso [poisson_distrib()]
S7::method(distrib_rng, PoissonDistrib) <- function(distrib, n, theta) {
  stats::rpois(
    n = n,
    lambda = theta[[1]]
  )
}

#' @title Poisson Analytical Gradient
#' @name distrib_gradient.PoissonDistrib
#' @description
#' Computes the analytical gradient (first derivative) of the Poisson log-probability 
#' with respect to the parameter \eqn{\mu}.
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu}}
#'
#' @param distrib A `PoissonDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vector of first derivatives.
#' @seealso [poisson_distrib()]
S7::method(distrib_gradient, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  poisson_gradient_cpp(y, theta[[1]], threads)
}

#' @title Poisson Analytical Observed Hessian
#' @name distrib_hessian.PoissonDistrib
#' @description
#' Computes the analytical observed Hessian (second derivative) of the Poisson log-probability 
#' with respect to the parameter \eqn{\mu}.
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}}
#'
#' @param distrib A `PoissonDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vector of second derivatives.
#' @seealso [poisson_distrib()]
S7::method(distrib_hessian, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  poisson_hessian_cpp(y, theta[[1]], threads)
}

#' @title Poisson Analytical Expected Hessian
#' @name distrib_expected_hessian.PoissonDistrib
#' @description
#' Computes the analytical expected Hessian of the Poisson log-probability 
#' with respect to the parameter \eqn{\mu}.
#'
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu}}
#'
#' @param distrib A `PoissonDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vector of expected second derivatives.
#' @seealso [poisson_distrib()]
S7::method(distrib_expected_hessian, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  poisson_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Poisson Analytical Third-Order Derivatives
#' @name distrib_deriv3.PoissonDistrib
#' @description Closed-form third-order derivative of the Poisson log-mass (observed, or expected when `expected = TRUE`).
#' @param distrib A `PoissonDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param expected Logical; if `TRUE`, returns the expected third derivative.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list with the `mu_mu_mu` component.
#' @seealso [poisson_distrib()]
S7::method(distrib_deriv3, PoissonDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) poisson_deriv3_expected_cpp(y, theta[[1]], threads)
  else poisson_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Poisson Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.PoissonDistrib
#' @description Closed-form fourth-order derivative of the Poisson log-mass (observed, or expected when `expected = TRUE`).
#' @param distrib A `PoissonDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter `mu`.
#' @param expected Logical; if `TRUE`, returns the expected fourth derivative.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list with the `mu_mu_mu_mu` component.
#' @seealso [poisson_distrib()]
S7::method(distrib_deriv4, PoissonDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) poisson_deriv4_expected_cpp(y, theta[[1]], threads)
  else poisson_deriv4_cpp(y, theta[[1]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Poisson Distribution Object
#'
#' @description
#' Creates a distribution object for the Poisson distribution parameterized by the mean parameter \eqn{\mu}.
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The Poisson distribution models counts \eqn{y \in \{0, 1, 2, \dots\}} with mean
#' (and variance) \eqn{\mu}.
#'
#' **Probability mass function:**
#' \deqn{P(Y=y; \mu) = \dfrac{\mu^y e^{-\mu}}{y!}}
#'
#' **Cumulative distribution function:**
#' \deqn{F(q; \mu) = \sum_{k=0}^{\lfloor q \rfloor} \dfrac{\mu^k e^{-\mu}}{k!}}
#'
#' **Quantile function:** the generalized inverse
#' \eqn{Q(p) = \min\{y \in \mathbb{N}_0 : F(y) \ge p\}}.
#'
#' **Score, observed and expected Hessian:**
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu}}
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{\mu}, skewness \eqn{1/\sqrt{\mu}},
#' excess kurtosis \eqn{1/\mu}.
#'
#' **Parameter domains:**
#'
#' - \eqn{\mu \in (0, +\infty)}
#'
#' Analytical third- and fourth-order derivatives ([distrib_deriv3()],
#' [distrib_deriv4()]) are also available.
#'
#' @seealso
#'
#' - [distrib_pdf.PoissonDistrib()] for the probability mass function.
#' - [distrib_cdf.PoissonDistrib()] for the cumulative distribution function.
#' - [distrib_quantile.PoissonDistrib()] for the quantile function.
#' - [distrib_rng.PoissonDistrib()] for random number generation.
#' - [distrib_gradient.PoissonDistrib()] for the analytical gradient.
#' - [distrib_hessian.PoissonDistrib()] for the analytical observed Hessian.
#' - [distrib_expected_hessian.PoissonDistrib()] for the analytical expected Hessian.
#'
#' @return An S7 object of class `PoissonDistrib` (inheriting from `discrete_distrib`) representing the Poisson distribution.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dpois ppois qpois rpois
#' @examples
#' d <- poisson_distrib()
#' d@params
#'
#' theta <- list(mu = 2)
#' distrib_pdf(d, 0:5, theta)
#' distrib_gradient(d, 0:5, theta)
#'
#' @export
poisson_distrib <- function(link_mu = log_link()) {
  PoissonDistrib(
    distrib_name = "poisson", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)), link_params = list(mu = link_mu)
  )
}
