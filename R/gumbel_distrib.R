#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Gumbel Distribution
#' @name GumbelDistrib
#'
#' @description A subclass of `continuous_distrib` representing the Gumbel
#' (type I extreme value, maximum) distribution.
#' @inheritParams distrib
#' @return An object of class `GumbelDistrib`.
#' @seealso [gumbel_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.GumbelDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GumbelDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.GumbelDistrib],
#'   [`distrib_gradient()`][distrib_gradient.GumbelDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.GumbelDistrib],
#'   [`distrib_hessian()`][distrib_hessian.GumbelDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.GumbelDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.GumbelDistrib],
#'   [`distrib_pdf()`][distrib_pdf.GumbelDistrib],
#'   [`distrib_quantile()`][distrib_quantile.GumbelDistrib],
#'   [`distrib_rng()`][distrib_rng.GumbelDistrib],
#'   [`kurtosis()`][kurtosis],
#'   [`mean()`][mean.distrib],
#'   [`skewness()`][skewness],
#'   [`variance()`][variance]
#'
#' Everything else is inherited from [continuous_distrib()].
GumbelDistrib <- S7::new_class("GumbelDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gumbel Probability Density Function
#' @name distrib_pdf.GumbelDistrib
#' @description
#' Computes the probability density function for the Gumbel distribution, with
#' \eqn{z = (y - \mu)/\sigma}:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sigma}
#'   \exp\left\{-z - e^{-z}\right\}}
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [gumbel_distrib()]
S7::method(distrib_pdf, GumbelDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  z <- (y - mu) / sigma
  log_d <- -log(sigma) - z - exp(-z)
  if (log) log_d else exp(log_d)
}

#' @title Gumbel Cumulative Distribution Function
#' @name distrib_cdf.GumbelDistrib
#' @description
#' Computes the cumulative distribution function for the Gumbel distribution:
#' \deqn{F(q; \mu, \sigma) = \exp\left\{-e^{-z}\right\}, \qquad
#'       z = (q - \mu)/\sigma}
#' @details
#' The lower tail is exact on the log scale, since \eqn{\log F = -e^{-z}}, and
#' the upper tail uses `expm1`, so neither loses precision in its own tail.
#' @param distrib A `GumbelDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [gumbel_distrib()]
S7::method(distrib_cdf, GumbelDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  z <- (q - theta[[1]]) / theta[[2]]
  w <- exp(-z)
  if (lower.tail) {
    # log F = -w exactly, so the lower tail never underflows to log(0).
    if (log.p) -w else exp(-w)
  } else {
    if (log.p) log(-expm1(-w)) else -expm1(-w)
  }
}

#' @title Gumbel Quantile Function
#' @name distrib_quantile.GumbelDistrib
#' @description
#' Computes the quantile function for the Gumbel distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu - \sigma \log(-\log p)}
#' @param distrib A `GumbelDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso [gumbel_distrib()]
S7::method(distrib_quantile, GumbelDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  # log p is what the formula needs, so it is kept rather than exponentiated
  # and logged again; for the upper tail log(1 - p) comes from log1p.
  lp <- if (lower.tail) {
    if (log.p) p else log(p)
  } else {
    if (log.p) log(-expm1(p)) else log1p(-p)
  }
  theta[[1]] - theta[[2]] * log(-lp)
}

#' @title Gumbel Random Number Generator
#' @name distrib_rng.GumbelDistrib
#' @description
#' Generates random numbers as \eqn{\mu - \sigma \log E} with \eqn{E} standard
#' exponential, which is inverse transform written without a logarithm of a
#' uniform.
#' @param distrib A `GumbelDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector of random draws.
#' @seealso [gumbel_distrib()]
S7::method(distrib_rng, GumbelDistrib) <- function(distrib, n, theta) {
  theta[[1]] - theta[[2]] * log(stats::rexp(n))
}

#' @title Gumbel Analytical Gradient
#' @name distrib_gradient.GumbelDistrib
#' @description
#' Closed-form first derivatives of the Gumbel log-density, written in
#' \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}}:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1 - w}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{z(1 - w) - 1}{\sigma}}
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [gumbel_distrib()]
S7::method(distrib_gradient, GumbelDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  z <- (y - mu) / sigma
  w <- exp(-z)
  list(
    mu = (1 - w) / sigma,
    sigma = (z * (1 - w) - 1) / sigma
  )
}

#' @title Gumbel Analytical Observed Hessian
#' @name distrib_hessian.GumbelDistrib
#' @description
#' Closed-form second derivatives of the Gumbel log-density:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{w}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = -\dfrac{1 - w + z w}{\sigma^2},}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = \dfrac{1 - 2z + 2zw - z^2 w}{\sigma^2}.}
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [gumbel_distrib()]
S7::method(distrib_hessian, GumbelDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  z <- (y - mu) / sigma
  w <- exp(-z)
  s2 <- sigma^2
  list(
    mu_mu = -w / s2,
    sigma_sigma = (1 - 2 * z + 2 * z * w - z^2 * w) / s2,
    mu_sigma = -(1 - w + z * w) / s2
  )
}

#' @title Gumbel Analytical Expected Hessian
#' @name distrib_expected_hessian.GumbelDistrib
#' @description
#' Closed form. Under the model \eqn{w = e^{-Z}} is standard exponential, so
#' every expectation needed is a derivative of \eqn{\Gamma} at 2:
#' \eqn{E[w] = 1}, \eqn{E[w \log w] = 1 - \gamma} and
#' \eqn{E[w(\log w)^2] = (1-\gamma)^2 + \pi^2/6 - 1}, together with
#' \eqn{E[Z] = \gamma}. Hence
#' \deqn{E\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'         = -\dfrac{1}{\sigma^2},
#'       \qquad
#'       E\left[\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}\right]
#'         = \dfrac{1 - \gamma}{\sigma^2},}
#' \deqn{E\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right]
#'         = -\dfrac{(1-\gamma)^2 + \pi^2/6}{\sigma^2},}
#' with \eqn{\gamma} the Euler-Mascheroni constant. Because the closed form
#' exists, the `approx` argument is ignored.
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param approx Ignored; the expectation is exact.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso [gumbel_distrib()]
S7::method(distrib_expected_hessian, GumbelDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  sigma <- theta[[2]]
  n <- length(y)
  s2 <- sigma^2
  eg <- -digamma(1)   # the Euler-Mascheroni constant
  list(
    mu_mu = rep(-1 / s2, length.out = n),
    sigma_sigma = rep(-((1 - eg)^2 + pi^2 / 6) / s2, length.out = n),
    mu_sigma = rep((1 - eg) / s2, length.out = n)
  )
}

#' @title Gumbel Analytical Third-Order Derivatives
#' @name distrib_deriv3.GumbelDistrib
#' @description
#' Closed-form third-order derivatives of the Gumbel log-density (observed, or
#' expected when `expected = TRUE`). With \eqn{z = (y-\mu)/\sigma} and
#' \eqn{w = e^{-z}}, every derivative is a polynomial in \eqn{z} and
#' \eqn{z^j w}; the expected values use \eqn{E[z^k w] = (-1)^k \Gamma^{(k)}(2)}
#' and \eqn{E[z] = \gamma}, since \eqn{w} is standard exponential under the
#' model.
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param expected Logical; if `TRUE`, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso [gumbel_distrib()]
S7::method(distrib_deriv3, GumbelDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) gumbel_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gumbel_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gumbel Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.GumbelDistrib
#' @description
#' Closed-form fourth-order derivatives of the Gumbel log-density (observed,
#' or expected when `expected = TRUE`), in the notation of
#' [distrib_deriv3.GumbelDistrib()].
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @param expected Logical; if `TRUE`, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [gumbel_distrib()]
S7::method(distrib_deriv4, GumbelDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) gumbel_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gumbel_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gumbel Response Derivative
#' @name distrib_grad_y.GumbelDistrib
#' @description
#' Closed form: \eqn{\partial \ell / \partial y = (w - 1)/\sigma}, which is
#' minus the derivative in \eqn{\mu}, as it must be for a location family.
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector.
#' @seealso [gumbel_distrib()]
S7::method(distrib_grad_y, GumbelDistrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  (exp(-(y - theta[[1]]) / sigma) - 1) / sigma
}

#' @title Gumbel Response Second Derivative
#' @name distrib_hess_y.GumbelDistrib
#' @description
#' Closed form: \eqn{\partial^2 \ell / \partial y^2 = -w/\sigma^2}.
#' @param distrib A `GumbelDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `sigma`.
#' @return A numeric vector.
#' @seealso [gumbel_distrib()]
S7::method(distrib_hess_y, GumbelDistrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  -exp(-(y - theta[[1]]) / sigma) / sigma^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Gumbel Distribution Object
#'
#' @description
#' Creates a distribution object for the Gumbel (type I extreme value)
#' distribution, in the form for **maxima**, with location \eqn{\mu} and
#' scale \eqn{\sigma}.
#'
#' @param link_mu A link function object for the location \eqn{\mu}. Defaults to
#'   [linkfunctions7::identity_link()].
#' @param link_sigma A link function object for the scale \eqn{\sigma}. Defaults
#'   to [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The Gumbel distribution is the limit law of the maximum of a sample from a
#' light-tailed distribution, which is what it is used for. It is a location-scale
#' family on the whole real line, skewed to the right with a fixed shape:
#' unlike the gaussian, its skewness and kurtosis are constants and cannot be
#' fitted.
#'
#' **Probability density function**, with \eqn{z = (y-\mu)/\sigma}:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sigma}\exp\left\{-z - e^{-z}\right\}}
#'
#' **Cumulative distribution function:**
#' \deqn{F(q; \mu, \sigma) = \exp\left\{-e^{-z}\right\}}
#'
#' **Quantile function:**
#' \deqn{Q(p; \mu, \sigma) = \mu - \sigma\log(-\log p)}
#'
#' **Score**, with \eqn{w = e^{-z}}:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1 - w}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{z(1 - w) - 1}{\sigma}}
#'
#' **Observed Hessian:**
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{w}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = -\dfrac{1 - w + zw}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = \dfrac{1 - 2z + 2zw - z^2 w}{\sigma^2}}
#'
#' **Expected Hessian:** see
#' [distrib_expected_hessian.GumbelDistrib()]. Note that
#' \eqn{E[\partial^2 \ell / \partial \mu \, \partial \sigma}] does not vanish:
#' the location and the scale are **not** orthogonal here, unlike in a
#' symmetric location-scale family, because the density is skewed.
#'
#' **Moments:** mean \eqn{\mu + \gamma\sigma}, variance
#' \eqn{\pi^2\sigma^2/6}, skewness \eqn{12\sqrt{6}\,\zeta(3)/\pi^3 \approx 1.1395}
#' and excess kurtosis \eqn{12/5}, the last two free of both parameters.
#'
#' **Relation to the Weibull.** If \eqn{Y} is Gumbel then \eqn{e^{-Y}} is
#' Weibull, so [weibull1_distrib()] is this family on the log scale and
#' reversed; the two share the expectations that produce their information
#' matrices.
#'
#' **Minima.** The distribution of minima is the reflection: fit this
#' family to \eqn{-Y}, or use [transformation()] with
#' [`affine_transform(scale = -1)`][affine_transform].
#'
#' **Higher orders.** Third and fourth derivatives are closed form,
#' observed and expected: with \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}},
#' every derivative is a polynomial in \eqn{z} and \eqn{z^j w}, and every
#' expectation is a derivative of \eqn{\Gamma} at 2.
#'
#' **Parameter Domains:**
#'
#' - \eqn{\mu \in (-\infty, +\infty)}
#' - \eqn{\sigma \in (0, +\infty)}
#'
#' @return An S7 object of class [GumbelDistrib()] (inheriting from
#'   `continuous_distrib`).
#'
#' @references
#' Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
#' Values*, chapter 3. Springer.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats rexp
#'
#' @examples
#' d <- gumbel_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' # the mean is shifted by Euler's constant, and the shape is fixed
#' c(mean = mean(d, theta), skewness = skewness(d, theta))
#'
#' # exp(-Y) is Weibull: the two families are one on the log scale
#' stats::sd(exp(-distrib_rng(d, 1000, theta))) > 0
#'
#' @seealso [weibull1_distrib()], [gpd_distrib()]
#' @export
gumbel_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  GumbelDistrib(
    distrib_name = "gumbel",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "sigma"),
    params_interpretation = c(mu = "location", sigma = "scale"),
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
