#' @include distrib.R generics.R
NULL

#' @title S7 Class for Logistic Distribution
#' @name LogisticDistrib
#' 
#' @description A subclass of `continuous_distrib` representing the Logistic distribution.
#' @inheritParams distrib
#' @return An object of class `LogisticDistrib`.
#' @seealso [logistic_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.LogisticDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.LogisticDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.LogisticDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.LogisticDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.LogisticDistrib],
#'   [`distrib_gradient()`][distrib_gradient.LogisticDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.LogisticDistrib],
#'   [`distrib_hessian()`][distrib_hessian.LogisticDistrib],
#'   [`distrib_pdf()`][distrib_pdf.LogisticDistrib],
#'   [`distrib_quantile()`][distrib_quantile.LogisticDistrib],
#'   [`distrib_rng()`][distrib_rng.LogisticDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
LogisticDistrib <- S7::new_class("LogisticDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Logistic Probability Density Function
#' @name distrib_pdf.LogisticDistrib
#' @description
#' Computes the probability density function for the Logistic distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left[1 + \exp\left(-\dfrac{y-\mu}{\sigma}\right)\right]^2}}
#' 
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [logistic_distrib()]
S7::method(distrib_pdf, LogisticDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dlogis(
    x = y,
    location = theta[[1]],
    scale = theta[[2]],
    log = log
  )
}

#' @title Logistic Cumulative Distribution Function
#' @name distrib_cdf.LogisticDistrib
#' @description
#' Computes the cumulative distribution function for the Logistic distribution:
#' \deqn{F(q; \mu, \sigma) = \dfrac{1}{1 + \exp\left(-\dfrac{q-\mu}{\sigma}\right)}}
#'
#' @param distrib A `LogisticDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [logistic_distrib()]
S7::method(distrib_cdf, LogisticDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::plogis(
    q = q,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Logistic Quantile Function
#' @name distrib_quantile.LogisticDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Logistic distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \log\left(\dfrac{p}{1-p}\right)}
#'
#' @param distrib A `LogisticDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso [logistic_distrib()]
S7::method(distrib_quantile, LogisticDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qlogis(
    p = p,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Logistic Random Number Generator
#' @name distrib_rng.LogisticDistrib
#' @description
#' Generates random numbers from the Logistic distribution.
#' 
#' @param distrib A `LogisticDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector of random draws.
#' @seealso [logistic_distrib()]
S7::method(distrib_rng, LogisticDistrib) <- function(distrib, n, theta) {
  stats::rlogis(
    n = n,
    location = theta[[1]],
    scale = theta[[2]]
  )
}

#' @title Logistic Analytical Gradient
#' @name distrib_gradient.LogisticDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Logistic log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right)}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma} \left[ 1 - \dfrac{y-\mu}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' 
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A list containing the vectors of first derivatives.
#' @seealso [logistic_distrib()]
S7::method(distrib_gradient, LogisticDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  logistic_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Logistic Analytical Observed Hessian
#' @name distrib_hessian.LogisticDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Logistic log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{2\sigma^2} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right)}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1}{\sigma^2} \left[ 1 - \dfrac{(y-\mu)^2}{2\sigma^2} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) - \dfrac{2(y-\mu)}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{1}{\sigma^2} \left[ \tanh\left(\dfrac{y-\mu}{2\sigma}\right) + \dfrac{y-\mu}{2\sigma} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' 
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A list containing the vectors of second derivatives.
#' @seealso [logistic_distrib()]
S7::method(distrib_hessian, LogisticDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  logistic_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Logistic Analytical Expected Hessian
#' @name distrib_expected_hessian.LogisticDistrib
#' @description
#' Computes the analytical expected Hessian of the Logistic log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{3\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{3+\pi^2}{9\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#' 
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [logistic_distrib()]
S7::method(distrib_expected_hessian, LogisticDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  logistic_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Logistic Analytical Third-Order Derivatives
#' @name distrib_deriv3.LogisticDistrib
#' @description
#' Closed-form observed third-order derivatives of the Logistic log-density.
#'
#' Writing \eqn{z = (y-\mu)/\sigma}, \eqn{t = 1/(1+e^{-z})} and \eqn{u = 1-t}, the
#' log-density is \eqn{\ell = -\log\sigma + g(z)} with
#' \eqn{g(z) = -z - 2\log(1+e^{-z})}, whose derivatives in \eqn{z} are
#' \deqn{g_1 = 1-2t, \quad g_2 = -2tu, \quad g_3 = -2tu(1-2t), \quad g_4 = -2tu(1-6tu).}
#' Since \eqn{z} depends on both parameters, each component is a polynomial in
#' \eqn{z} with the \eqn{g_j} as coefficients:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = -\dfrac{g_3}{\sigma^3}}
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^2 \partial \sigma} = -\dfrac{2g_2 + z g_3}{\sigma^3}}
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu \partial \sigma^2} = -\dfrac{2g_1 + 4z g_2 + z^2 g_3}{\sigma^3}}
#' \deqn{\dfrac{\partial^3 \ell}{\partial \sigma^3} = -\dfrac{2 + 6z g_1 + 6z^2 g_2 + z^3 g_3}{\sigma^3}}
#'
#' The expected third derivatives are not available in closed form, so
#' `expected = TRUE` is routed to [expected_derivative_methods()].
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param expected Logical; if `TRUE`, returns the approximated expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso [logistic_distrib()]
S7::method(distrib_deriv3, LogisticDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    logistic_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Logistic Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.LogisticDistrib
#' @description
#' Closed-form observed fourth-order derivatives of the Logistic log-density, in the
#' notation of [distrib_deriv3.LogisticDistrib()]:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4} = \dfrac{g_4}{\sigma^4}}
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^3 \partial \sigma} = \dfrac{3g_3 + z g_4}{\sigma^4}}
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^2 \partial \sigma^2} = \dfrac{6g_2 + 6z g_3 + z^2 g_4}{\sigma^4}}
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu \partial \sigma^3} = \dfrac{6g_1 + 18z g_2 + 9z^2 g_3 + z^3 g_4}{\sigma^4}}
#' \deqn{\dfrac{\partial^4 \ell}{\partial \sigma^4} = \dfrac{6 + 24z g_1 + 36z^2 g_2 + 12z^3 g_3 + z^4 g_4}{\sigma^4}}
#'
#' The expected fourth derivatives are not available in closed form, so
#' `expected = TRUE` is routed to [expected_derivative_methods()].
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param expected Logical; if `TRUE`, returns the approximated expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [logistic_distrib()]
S7::method(distrib_deriv4, LogisticDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    logistic_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Logistic Response Derivatives
#' @name distrib_grad_y.LogisticDistrib
#' @description
#' Closed-form derivatives of the Logistic log-density with respect to the response.
#' With \eqn{z = (y-\mu)/\sigma}: \eqn{\partial \ell / \partial y = -\tanh(z/2)/\sigma}
#' and \eqn{\partial^2 \ell / \partial y^2 = -\mathrm{sech}^2(z/2)/(2\sigma^2)}.
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector.
#' @seealso [logistic_distrib()]
S7::method(distrib_grad_y, LogisticDistrib) <- function(distrib, y, theta) {
  s <- theta[[2]]
  -tanh(0.5 * (y - theta[[1]]) / s) / s
}

#' @title Logistic Response Second Derivative
#' @name distrib_hess_y.LogisticDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -\mathrm{sech}^2(z/2)/(2\sigma^2)}, \eqn{z = (y-\mu)/\sigma}.
#' @param distrib A `LogisticDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector.
#' @seealso [logistic_distrib()]
S7::method(distrib_hess_y, LogisticDistrib) <- function(distrib, y, theta) {
  s <- theta[[2]]
  th <- tanh(0.5 * (y - theta[[1]]) / s)
  -(1 - th^2) / (2 * s^2)
}

# --- CONSTRUCTOR WRAPPER ---

#' Logistic Distribution Object
#'
#' @description
#' Creates a distribution object for the Logistic distribution parameterized by location (\eqn{\mu})
#' and scale (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to [linkfunctions7::identity_link()].
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The Logistic distribution is a location-scale distribution with location
#' \eqn{\mu} and scale \eqn{\sigma}.
#'
#' **Probability density function:**
#' \deqn{f(y; \mu, \sigma) = \dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma\left[1 + \exp\left(-\dfrac{y-\mu}{\sigma}\right)\right]^2}}
#'
#' **Cumulative distribution function:**
#' \deqn{F(q; \mu, \sigma) = \dfrac{1}{1 + \exp\left(-\dfrac{q-\mu}{\sigma}\right)}}
#'
#' **Quantile function:**
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \log\left(\dfrac{p}{1-p}\right)}
#'
#' **Score:**
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right), \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma}\left[1 - \dfrac{y-\mu}{\sigma}\tanh\left(\dfrac{y-\mu}{2\sigma}\right)\right]}
#'
#' **Expected Hessian:**
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{3\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{3+\pi^2}{9\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0}
#' The observed Hessian is available via [distrib_hessian.LogisticDistrib()].
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{\pi^2\sigma^2/3}, skewness 0,
#' excess kurtosis \eqn{6/5}.
#'
#' **Parameter domains:**
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' Response derivatives ([distrib_grad_y()], [distrib_hess_y()])
#' and observed third- and fourth-order parameter derivatives
#' ([distrib_deriv3()], [distrib_deriv4()]) are available in
#' closed form. The corresponding expected derivatives are not: they are obtained
#' through [expected_derivative_methods()]. Seven of the nine components are
#' nonetheless known exactly, and the location-scale structure of the family forbids
#' any of them from depending on \eqn{\mu}:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^3 \ell}{\partial \mu^2 \partial \sigma}\right] = \dfrac{1}{2\sigma^3}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^3 \ell}{\partial \sigma^3}\right] = \dfrac{\pi^2+2}{2\sigma^3}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^4 \ell}{\partial \mu^4}\right] = \dfrac{1}{15\sigma^4},}
#' the remaining ones in that list being zero by symmetry.
#'
#' @seealso
#' \itemize{
#'   \item [distrib_pdf.LogisticDistrib()] for the probability density function.
#'   \item [distrib_cdf.LogisticDistrib()] for the cumulative distribution function.
#'   \item [distrib_quantile.LogisticDistrib()] for the quantile function.
#'   \item [distrib_rng.LogisticDistrib()] for random number generation.
#'   \item [distrib_gradient.LogisticDistrib()] for the analytical gradient.
#'   \item [distrib_hessian.LogisticDistrib()] for the analytical observed Hessian.
#'   \item [distrib_expected_hessian.LogisticDistrib()] for the analytical expected Hessian.
#'   \item [distrib_deriv3.LogisticDistrib()] and [distrib_deriv4.LogisticDistrib()] for the observed higher-order derivatives.
#' }
#'
#' @return An S7 object of class `LogisticDistrib` (inheriting from `continuous_distrib`) representing the Logistic distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dlogis plogis qlogis rlogis
#' @examples
#' d <- logistic_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' @export
logistic_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  LogisticDistrib(
    distrib_name = "logistic", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "mean", sigma = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma)
  )
}
