#' @include distrib.R generics.R
NULL

#' @title S7 Class for Gaussian Distribution
#' @name Gaussian1Distrib
#' 
#' @description A subclass of `continuous_distrib` representing the Gaussian (Normal) distribution.
#' @inheritParams distrib
#' @return An object of class `Gaussian1Distrib`.
#' @seealso [gaussian1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.Gaussian1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gaussian1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gaussian1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gaussian1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gaussian1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gaussian1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gaussian1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gaussian1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gaussian1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gaussian1Distrib],
#'   [`distrib_rng()`][distrib_rng.Gaussian1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
Gaussian1Distrib <- S7::new_class("Gaussian1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Probability Density Function
#' @name distrib_pdf.Gaussian1Distrib
#' @description
#' Computes the probability density function for the Gaussian distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\}}
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_pdf, Gaussian1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(
    x = y,
    mean = theta[[1]],
    sd = theta[[2]],
    log = log
  )
}

#' @title Gaussian Cumulative Distribution Function
#' @name distrib_cdf.Gaussian1Distrib
#' @description
#' Computes the cumulative distribution function for the Gaussian distribution:
#' \deqn{F(q; \mu, \sigma) = \Phi\left(\dfrac{q-\mu}{\sigma}\right)}
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_cdf, Gaussian1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pnorm(
    q = q,
    mean = theta[[1]],
    sd = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gaussian Quantile Function
#' @name distrib_quantile.Gaussian1Distrib
#' @description
#' Computes the quantile function (inverse CDF) for the Gaussian distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \Phi^{-1}(p)}
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_quantile, Gaussian1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qnorm(
    p = p,
    mean = theta[[1]],
    sd = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gaussian Random Number Generator
#' @name distrib_rng.Gaussian1Distrib
#' @description
#' Generates random numbers from the Gaussian distribution.
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector of random draws.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_rng, Gaussian1Distrib) <- function(distrib, n, theta) {
  stats::rnorm(
    n = n,
    mean = theta[[1]],
    sd = theta[[2]]
  )
}

#' @title Gaussian Analytical Gradient
#' @name distrib_gradient.Gaussian1Distrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}}
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of first derivatives.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_gradient, Gaussian1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Observed Hessian
#' @name distrib_hessian.Gaussian1Distrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y - \mu)^2}{\sigma^4}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{2(y - \mu)}{\sigma^3}}
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of second derivatives.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_hessian, Gaussian1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Expected Hessian
#' @name distrib_expected_hessian.Gaussian1Distrib
#' @description
#' Computes the analytical expected Hessian of the Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#' 
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_expected_hessian, Gaussian1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  gaussian_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Third-Order Derivatives
#' @name distrib_deriv3.Gaussian1Distrib
#' @description
#' Computes the closed-form third-order partial derivatives of the Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\sigma} (observed, or expected when
#' `expected = TRUE`).
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param expected Logical; if `TRUE`, returns the expected third derivatives.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of third-derivative component vectors.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_deriv3, Gaussian1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.Gaussian1Distrib
#' @description
#' Computes the closed-form fourth-order partial derivatives of the Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\sigma} (observed, or expected when
#' `expected = TRUE`).
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param expected Logical; if `TRUE`, returns the expected fourth derivatives.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_deriv4, Gaussian1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian Response Derivatives
#' @name distrib_grad_y.Gaussian1Distrib
#' @description
#' Closed-form derivatives of the Gaussian log-density with respect to the response:
#' \eqn{\partial \ell / \partial y = -(y-\mu)/\sigma^2} and
#' \eqn{\partial^2 \ell / \partial y^2 = -1/\sigma^2}.
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_grad_y, Gaussian1Distrib) <- function(distrib, y, theta) {
  -(y - theta[[1]]) / theta[[2]]^2
}

#' @title Gaussian Response Second Derivative
#' @name distrib_hess_y.Gaussian1Distrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -1/\sigma^2}.
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector.
#' @seealso [gaussian1_distrib()]
S7::method(distrib_hess_y, Gaussian1Distrib) <- function(distrib, y, theta) {
  rep(-1 / theta[[2]]^2, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Gaussian Distribution Object (Standard Deviation Parameterization)
#'
#' @description
#' Creates a distribution object for the Gaussian distribution parameterized by mean (\eqn{\mu}) and standard deviation (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to [linkfunctions7::identity_link()].
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The Gaussian (Normal) distribution is parameterized by its mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}.
#'
#' **Probability density function:**
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\,\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\}}
#'
#' **Cumulative distribution function** (\eqn{\Phi} the standard normal CDF):
#' \deqn{F(q; \mu, \sigma) = \Phi\left(\dfrac{q-\mu}{\sigma}\right)}
#'
#' **Quantile function:**
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma\,\Phi^{-1}(p)}
#'
#' **Score** (gradient of the log-density):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}}
#'
#' **Observed Hessian:**
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y-\mu)^2}{\sigma^4}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma} = -\dfrac{2(y-\mu)}{\sigma^3}}
#'
#' **Expected Hessian** (negative Fisher information):
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0}
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{\sigma^2}, skewness 0, excess kurtosis 0.
#'
#' **Parameter domains:**
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives ([distrib_deriv3()],
#' [distrib_deriv4()]) and derivatives with respect to the response
#' ([distrib_grad_y()], [distrib_hess_y()]) are also available.
#'
#' @seealso
#' \itemize{
#'   \item [distrib_pdf.Gaussian1Distrib()] for the density function.
#'   \item [distrib_cdf.Gaussian1Distrib()] for the cumulative distribution function.
#'   \item [distrib_quantile.Gaussian1Distrib()] for the quantile function.
#'   \item [distrib_rng.Gaussian1Distrib()] for random number generation.
#'   \item [distrib_gradient.Gaussian1Distrib()] for the analytical gradient.
#'   \item [distrib_hessian.Gaussian1Distrib()] for the analytical observed Hessian.
#'   \item [distrib_expected_hessian.Gaussian1Distrib()] for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class `Gaussian1Distrib` (inheriting from `continuous_distrib`) representing the Gaussian distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @examples
#' d <- gaussian1_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' @export
gaussian1_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  
  Gaussian1Distrib(
    distrib_name = "gaussian1",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean", sigma = "standard deviation"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      sigma = link_sigma
    )
  )
  
}
