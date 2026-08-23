#' @include distrib.R generics.R
NULL

#' @title S7 Class for Inverse-Gaussian Distribution
#' @name InvGauss1Distrib
#' 
#' @description A subclass of `continuous_distrib` representing the Inverse-Gaussian distribution.
#' @inheritParams distrib
#' @return An object of class `InvGauss1Distrib`.
#' @seealso [invgauss1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.InvGauss1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.InvGauss1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.InvGauss1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.InvGauss1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.InvGauss1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.InvGauss1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.InvGauss1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.InvGauss1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.InvGauss1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.InvGauss1Distrib],
#'   [`distrib_rng()`][distrib_rng.InvGauss1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
InvGauss1Distrib <- S7::new_class("InvGauss1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Inverse-Gaussian Probability Density Function
#' @name distrib_pdf.InvGauss1Distrib
#' @description
#' Computes the probability density function for the Inverse-Gaussian distribution:
#' \deqn{f(y; \mu, \phi) = \sqrt{\dfrac{1}{2\pi\phi y^3}} \exp\left\{-\dfrac{(y-\mu)^2}{2\phi\mu^2 y}\right\}}
#' 
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations (\eqn{y > 0}).
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_pdf, InvGauss1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  statmod::dinvgauss(
    x = y,
    mean = theta[[1]],
    dispersion = theta[[2]],
    log = log
  )
}

#' @title Inverse-Gaussian Cumulative Distribution Function
#' @name distrib_cdf.InvGauss1Distrib
#' @description
#' Computes the cumulative distribution function for the Inverse-Gaussian distribution:
#' \deqn{F(q; \mu, \phi) = \Phi\!\left(\sqrt{\dfrac{1}{\phi q}}\left(\dfrac{q}{\mu}-1\right)\right) + e^{2/(\phi\mu)}\, \Phi\!\left(-\sqrt{\dfrac{1}{\phi q}}\left(\dfrac{q}{\mu}+1\right)\right)}
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' @param distrib An `InvGauss1Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_cdf, InvGauss1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  statmod::pinvgauss(
    q = q,
    mean = theta[[1]],
    dispersion = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Inverse-Gaussian Quantile Function
#' @name distrib_quantile.InvGauss1Distrib
#' @description
#' Computes the quantile function for the Inverse-Gaussian distribution as the
#' inverse of the CDF, \eqn{Q(p; \mu, \phi) = F^{-1}(p; \mu, \phi)}. There is no
#' closed form; it is obtained numerically via [statmod::qinvgauss()].
#'
#' @param distrib An `InvGauss1Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_quantile, InvGauss1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  statmod::qinvgauss(
    p = p,
    mean = theta[[1]],
    dispersion = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Inverse-Gaussian Random Number Generator
#' @name distrib_rng.InvGauss1Distrib
#' @description
#' Generates random numbers from the Inverse-Gaussian distribution.
#' 
#' @param distrib An `InvGauss1Distrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @return A numeric vector of random draws.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_rng, InvGauss1Distrib) <- function(distrib, n, theta) {
  statmod::rinvgauss(
    n = n,
    mean = theta[[1]],
    dispersion = theta[[2]]
  )
}

#' @title Inverse-Gaussian Analytical Gradient
#' @name distrib_gradient.InvGauss1Distrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Inverse-Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^3}}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2}}
#' 
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of first derivatives.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_gradient, InvGauss1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  invgauss_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse-Gaussian Analytical Observed Hessian
#' @name distrib_hessian.InvGauss1Distrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Inverse-Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{3y - 2\mu}{\phi\mu^4}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \phi^2} = \dfrac{\phi - 2(y-\mu)^2/(\mu^2 y)}{2\phi^3}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \phi} = -\dfrac{y - \mu}{\phi^2\mu^3}}
#' 
#' **Note:** The observed Hessian with respect to \eqn{\phi} is not guaranteed to be negative for all
#' observed values of \eqn{y}. Specifically, \eqn{\partial^2 \ell/\partial \phi^2 < 0} only when
#' \eqn{\phi < 2(y-\mu)^2/(\mu^2 y)}. This condition may be violated when observations are far from the mean
#' or when the dispersion parameter is large, potentially causing numerical instability in optimization
#' algorithms that rely on the observed Hessian (e.g., Newton-Raphson). In such cases, using the expected
#' Hessian (`distrib_expected_hessian`) is recommended for more stable convergence.
#' 
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of second derivatives.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_hessian, InvGauss1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  invgauss_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse-Gaussian Analytical Expected Hessian
#' @name distrib_expected_hessian.InvGauss1Distrib
#' @description
#' Computes the analytical expected Hessian of the Inverse-Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\phi\mu^3}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \phi^2}\right] = -\dfrac{1}{2\phi^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \phi}\right] = 0}
#' 
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_expected_hessian, InvGauss1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  invgauss_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse-Gaussian Analytical Third-Order Derivatives
#' @name distrib_deriv3.InvGauss1Distrib
#' @description Closed-form third-order derivatives of the Inverse-Gaussian log-density (observed, or expected when `expected = TRUE`).
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param expected Logical; if `TRUE`, returns the expected third derivatives.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of third-derivative component vectors.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_deriv3, InvGauss1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) invgauss_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else invgauss_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse-Gaussian Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.InvGauss1Distrib
#' @description Closed-form fourth-order derivatives of the Inverse-Gaussian log-density (observed, or expected when `expected = TRUE`).
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @param expected Logical; if `TRUE`, returns the expected fourth derivatives.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_deriv4, InvGauss1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) invgauss_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else invgauss_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse-Gaussian Response Derivatives
#' @name distrib_grad_y.InvGauss1Distrib
#' @description
#' Closed-form derivatives of the Inverse-Gaussian log-density with respect to the
#' response:
#' \eqn{\partial \ell / \partial y = -\dfrac{3}{2y} - \dfrac{y^2 - \mu^2}{2\phi\mu^2 y^2}} and
#' \eqn{\partial^2 \ell / \partial y^2 = \dfrac{3}{2y^2} - \dfrac{1}{\phi y^3}}.
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @return A numeric vector.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_grad_y, InvGauss1Distrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; phi <- theta[[2]]
  -1.5 / y - (y^2 - mu^2) / (2 * phi * mu^2 * y^2)
}

#' @title Inverse-Gaussian Response Second Derivative
#' @name distrib_hess_y.InvGauss1Distrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = 3/(2y^2) - 1/(\phi y^3)}.
#' @param distrib An `InvGauss1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `phi`.
#' @return A numeric vector.
#' @seealso [invgauss1_distrib()]
S7::method(distrib_hess_y, InvGauss1Distrib) <- function(distrib, y, theta) {
  phi <- theta[[2]]
  1.5 / y^2 - 1 / (phi * y^3)
}

# --- CONSTRUCTOR WRAPPER ---

#' Inverse-Gaussian Distribution Object (Mean-Dispersion Parameterization)
#'
#' @description
#' Creates a distribution object for the Inverse-Gaussian distribution parameterized by mean (\eqn{\mu}) and dispersion (\eqn{\phi}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to [linkfunctions7::log_link()] to ensure positivity.
#' @param link_phi A link function object for the dispersion parameter \eqn{\phi}.
#'   Defaults to [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The Inverse-Gaussian distribution is given a mean/dispersion parameterization,
#' with mean \eqn{\mu} and dispersion \eqn{\phi}.
#'
#' **Probability density function:**
#' \deqn{f(y; \mu, \phi) = \sqrt{\dfrac{1}{2\pi\phi y^3}} \exp\left\{-\dfrac{(y-\mu)^2}{2\phi\mu^2 y}\right\}, \quad y > 0}
#'
#' **Cumulative distribution function** (\eqn{\Phi} the standard normal CDF):
#' \deqn{F(q; \mu, \phi) = \Phi\!\left(\sqrt{\tfrac{1}{\phi q}}\left(\tfrac{q}{\mu}-1\right)\right) + e^{2/(\phi\mu)}\, \Phi\!\left(-\sqrt{\tfrac{1}{\phi q}}\left(\tfrac{q}{\mu}+1\right)\right)}
#'
#' **Quantile function:** no closed form; the numerical inverse of the CDF.
#'
#' **Score:**
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^3}, \qquad
#'       \dfrac{\partial \ell}{\partial \phi} = \dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2}}
#'
#' **Observed Hessian:**
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{3y - 2\mu}{\phi\mu^4}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \phi^2} = \dfrac{\phi - 2(y-\mu)^2/(\mu^2 y)}{2\phi^3}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \phi} = -\dfrac{y - \mu}{\phi^2\mu^3}}
#'
#' **Expected Hessian:**
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\phi\mu^3}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \phi^2}\right] = -\dfrac{1}{2\phi^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \phi}\right] = 0}
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{\phi\mu^3},
#' skewness \eqn{3\sqrt{\phi\mu}}, excess kurtosis \eqn{15\phi\mu}.
#'
#' **Parameter domains:**
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\phi \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives ([distrib_deriv3()],
#' [distrib_deriv4()]) and response derivatives ([distrib_grad_y()],
#' [distrib_hess_y()]) are also available.
#'
#' @seealso
#' \itemize{
#'   \item [distrib_pdf.InvGauss1Distrib()] for the probability density function.
#'   \item [distrib_cdf.InvGauss1Distrib()] for the cumulative distribution function.
#'   \item [distrib_quantile.InvGauss1Distrib()] for the quantile function.
#'   \item [distrib_rng.InvGauss1Distrib()] for random number generation.
#'   \item [distrib_gradient.InvGauss1Distrib()] for the analytical gradient.
#'   \item [distrib_hessian.InvGauss1Distrib()] for the analytical observed Hessian.
#'   \item [distrib_expected_hessian.InvGauss1Distrib()] for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class `InvGauss1Distrib` (inheriting from `continuous_distrib`) representing the Inverse-Gaussian distribution.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom statmod dinvgauss pinvgauss qinvgauss rinvgauss
#' @examples
#' d <- invgauss1_distrib()
#' d@params
#'
#' theta <- list(mu = 1, phi = 2)
#' distrib_pdf(d, c(0.5, 1, 2), theta)
#' distrib_gradient(d, c(0.5, 1, 2), theta)
#'
#' @export
invgauss1_distrib <- function(link_mu = log_link(), link_phi = log_link()) {
  InvGauss1Distrib(
    distrib_name = "invgauss1", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu", "phi"), params_interpretation = c(mu = "mean", phi = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), phi = c(0, Inf)),
    link_params = list(mu = link_mu, phi = link_phi)
  )
}
