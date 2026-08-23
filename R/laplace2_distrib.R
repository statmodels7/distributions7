#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Laplace Distribution in Location and Rate
#' @name Laplace2Distrib
#'
#' @description A subclass of `continuous_distrib` for the Laplace
#' (double-exponential) distribution written in its location \eqn{\mu} and its
#' **rate** \eqn{\lambda = 1/\sigma}. Like [LaplaceDistrib()],
#' its log-likelihood is not differentiable in \eqn{\mu}.
#' @inheritParams distrib
#' @return An object of class `Laplace2Distrib`.
#' @seealso [laplace2_distrib()], [laplace_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.Laplace2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Laplace2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Laplace2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Laplace2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Laplace2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Laplace2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Laplace2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Laplace2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Laplace2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Laplace2Distrib],
#'   [`distrib_rng()`][distrib_rng.Laplace2Distrib],
#'   [`kurtosis()`][kurtosis],
#'   [`mean()`][mean.distrib],
#'   [`skewness()`][skewness],
#'   [`variance()`][variance]
#'
#' Everything else is inherited from [continuous_distrib()].
Laplace2Distrib <- S7::new_class("Laplace2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Laplace Density in Location and Rate
#' @name distrib_pdf.Laplace2Distrib
#' @description
#' \deqn{f(y; \mu, \lambda) = \dfrac{\lambda}{2} \exp\left(-\lambda|y-\mu|\right)}
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [laplace2_distrib()]
S7::method(distrib_pdf, Laplace2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  log_d <- base::log(lam / 2) - lam * abs(y - mu)
  if (log) log_d else exp(log_d)
}

#' @title Laplace Distribution Function in Location and Rate
#' @name distrib_cdf.Laplace2Distrib
#' @description
#' \deqn{F(q; \mu, \lambda) = \begin{cases} \dfrac{1}{2}\exp\left(\lambda(q-\mu)\right) & q < \mu \\ 1 - \dfrac{1}{2}\exp\left(-\lambda(q-\mu)\right) & q \ge \mu \end{cases}}
#' @param distrib A `Laplace2Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [laplace2_distrib()]
S7::method(distrib_cdf, Laplace2Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  res <- ifelse(q < mu, 0.5 * exp(lam * (q - mu)), 1 - 0.5 * exp(-lam * (q - mu)))
  if (!lower.tail) res <- 1 - res
  if (log.p) base::log(res) else res
}

#' @title Laplace Quantile Function in Location and Rate
#' @name distrib_quantile.Laplace2Distrib
#' @description
#' \deqn{Q(p; \mu, \lambda) = \mu - \dfrac{1}{\lambda}\,\mathrm{sign}(p - \tfrac{1}{2})\,\log\left(1 - 2\left|p - \tfrac{1}{2}\right|\right)}
#' @param distrib A `Laplace2Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso [laplace2_distrib()]
S7::method(distrib_quantile, Laplace2Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  mu - sign(p - 0.5) * base::log(1 - 2 * abs(p - 0.5)) / lam
}

#' @title Laplace Random Generation in Location and Rate
#' @name distrib_rng.Laplace2Distrib
#' @description Generates draws by inverse-transform sampling.
#' @param distrib A `Laplace2Distrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @return A numeric vector of random draws.
#' @seealso [laplace2_distrib()]
S7::method(distrib_rng, Laplace2Distrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Laplace Analytical Gradient in Location and Rate
#' @name distrib_gradient.Laplace2Distrib
#' @description
#' The derivative with respect to \eqn{\mu} exists only almost everywhere
#' (there is a kink at \eqn{y = \mu}, a set of probability zero); the
#' subgradient value 0 is returned there.
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \lambda\,\mathrm{sign}(y-\mu), \qquad
#'       \dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{\lambda} - |y-\mu|}
#'
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @return A list containing the vectors of first derivatives.
#' @seealso [laplace2_distrib()]
S7::method(distrib_gradient, Laplace2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  r <- y - mu
  list(
    mu = lam * sign(r),
    lambda = 1 / lam - abs(r)
  )
}

#' @title Laplace Analytical Observed Hessian in Location and Rate
#' @name distrib_hessian.Laplace2Distrib
#' @description
#' Almost everywhere,
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \partial \lambda} = \mathrm{sign}(y-\mu), \qquad
#'       \dfrac{\partial^2 \ell}{\partial \lambda^2} = -\dfrac{1}{\lambda^2}}
#'
#' The log-likelihood is linear in \eqn{\lambda} apart from the
#' \eqn{\log\lambda} term, so every derivative in \eqn{\lambda} beyond the
#' first is free of the data. As for [laplace_distrib()], the zero
#' second derivative in \eqn{\mu} means Newton-Raphson cannot update
#' \eqn{\mu}; [distrib_expected_hessian()] supplies the Fisher
#' information \eqn{\lambda^2}.
#'
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @return A list containing the vectors of second derivatives.
#' @seealso [laplace2_distrib()]
S7::method(distrib_hessian, Laplace2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  n <- length(y)
  list(
    mu_mu = rep(0, n),
    lambda_lambda = rep(-1 / lam^2, length.out = n),
    mu_lambda = sign(y - mu) + rep(0, n)
  )
}

#' @title Laplace Analytical Expected Hessian in Location and Rate
#' @name distrib_expected_hessian.Laplace2Distrib
#' @description
#' As for [laplace_distrib()], the second Bartlett identity fails in
#' \eqn{\mu} and the expected Hessian is defined from the variance of the
#' score:
#'
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\lambda^2, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \lambda}\right] = 0, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \lambda^2}\right] = -\dfrac{1}{\lambda^2}}
#'
#' Because the closed form exists, the `approx` argument is ignored.
#'
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [laplace2_distrib()]
S7::method(distrib_expected_hessian, Laplace2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  lam <- theta[[2]]
  n <- length(y)
  list(
    mu_mu = rep(-lam^2, length.out = n),
    lambda_lambda = rep(-1 / lam^2, length.out = n),
    mu_lambda = rep(0, n)
  )
}

#' @title Laplace Third-Order Derivatives in Location and Rate
#' @name distrib_deriv3.Laplace2Distrib
#' @description
#' Closed form, almost everywhere. The only non-zero component is
#' \eqn{\ell^{(\lambda\lambda\lambda)} = 2/\lambda^3}, which is free of the
#' data, so the observed and the expected derivatives coincide.
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @param expected Logical; if `TRUE`, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso [laplace2_distrib()]
S7::method(distrib_deriv3, Laplace2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  lam <- theta[[2]]
  n <- length(y)
  zero <- rep(0, n)
  list(
    mu_mu_mu = zero,
    mu_mu_lambda = zero,
    mu_lambda_lambda = zero,
    lambda_lambda_lambda = rep(2 / lam^3, length.out = n)
  )
}

#' @title Laplace Fourth-Order Derivatives in Location and Rate
#' @name distrib_deriv4.Laplace2Distrib
#' @description
#' Closed form, almost everywhere. The only non-zero component is
#' \eqn{\ell^{(\lambda\lambda\lambda\lambda)} = -6/\lambda^4}; the observed and
#' the expected derivatives coincide.
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @param expected Logical; if `TRUE`, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [laplace2_distrib()]
S7::method(distrib_deriv4, Laplace2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  lam <- theta[[2]]
  n <- length(y)
  zero <- rep(0, n)
  list(
    mu_mu_mu_mu = zero,
    mu_mu_mu_lambda = zero,
    mu_mu_lambda_lambda = zero,
    mu_lambda_lambda_lambda = zero,
    lambda_lambda_lambda_lambda = rep(-6 / lam^4, length.out = n)
  )
}

#' @title Laplace Response Derivative in Location and Rate
#' @name distrib_grad_y.Laplace2Distrib
#' @description
#' \eqn{\partial \ell / \partial y = -\lambda\,\mathrm{sign}(y-\mu)}, almost
#' everywhere; the analytic form is provided because finite differences would
#' be inaccurate across the kink at \eqn{y = \mu}.
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @return A numeric vector.
#' @seealso [laplace2_distrib()]
S7::method(distrib_grad_y, Laplace2Distrib) <- function(distrib, y, theta) {
  -theta[[2]] * sign(y - theta[[1]])
}

#' @title Laplace Response Second Derivative in Location and Rate
#' @name distrib_hess_y.Laplace2Distrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = 0} (almost everywhere).
#' @param distrib A `Laplace2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters `mu` and `lambda`.
#' @return A numeric vector of zeros.
#' @seealso [laplace2_distrib()]
S7::method(distrib_hess_y, Laplace2Distrib) <- function(distrib, y, theta) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Laplace Distribution in Location and Rate
#'
#' @description
#' Creates a Laplace (double-exponential) distribution object parametrized by
#' location (\eqn{\mu}) and **rate** (\eqn{\lambda}).
#'
#' @details
#' This is the same law as [laplace_distrib()] in different
#' coordinates: \eqn{\lambda} here is \eqn{1/\sigma} there. The rate form is
#' the one a penalty consumes, because the negative log-density at fixed
#' \eqn{\mu = 0} is \eqn{\lambda \lvert y \rvert} up to a constant: the lasso
#' penalty is linear in \eqn{\lambda}, and its derivatives in \eqn{\lambda}
#' beyond the first carry no data at all.
#'
#' **Probability density function:**
#' \deqn{f(y; \mu, \lambda) = \dfrac{\lambda}{2} \exp\left(-\lambda|y-\mu|\right)}
#'
#' **Score** (defined almost everywhere):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \lambda\,\mathrm{sign}(y-\mu), \qquad
#'       \dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{\lambda} - |y-\mu|}
#'
#' **Observed Hessian** (almost everywhere):
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \lambda} = \mathrm{sign}(y-\mu), \quad
#'       \dfrac{\partial^2 \ell}{\partial \lambda^2} = -\dfrac{1}{\lambda^2}}
#'
#' **Expected Hessian** (Fisher information from the score variance):
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\lambda^2, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \lambda^2}\right] = -\dfrac{1}{\lambda^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \lambda}\right] = 0}
#'
#' **Moments:** mean \eqn{\mu}, variance \eqn{2/\lambda^2}, skewness 0,
#' excess kurtosis 3.
#'
#' The kink at \eqn{y = \mu} and its consequences are those of
#' [laplace_distrib()]: `params_smooth = c(mu = FALSE, lambda
#' = TRUE)`, the observed Hessian cannot update \eqn{\mu}, and the expected
#' Hessian comes from the variance of the score.
#'
#' **Parameter Domains:**
#'
#' - \eqn{\mu \in (-\infty, +\infty)}
#' - \eqn{\lambda \in (0, +\infty)}
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to [linkfunctions7::identity_link()].
#' @param link_lambda A link function object for the rate parameter \eqn{\lambda}.
#'   Defaults to [linkfunctions7::log_link()] to ensure positivity.
#'
#' @return An S7 object of class `Laplace2Distrib` (inheriting from `continuous_distrib`).
#'
#' @seealso [laplace_distrib()]
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats runif
#' @examples
#' d <- laplace2_distrib()
#' theta <- list(mu = 0, lambda = 2)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#'
#' # the same law as laplace_distrib with sigma = 1/2
#' distrib_pdf(laplace_distrib(), 0.7, list(mu = 0, sigma = 0.5))
#'
#' @export
laplace2_distrib <- function(link_mu = identity_link(), link_lambda = log_link()) {
  Laplace2Distrib(
    distrib_name = "laplace2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "lambda"),
    params_interpretation = c(mu = "location", lambda = "rate"),
    n_params = 2,

    params_bounds = list(
      mu = c(-Inf, Inf),
      lambda = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      lambda = link_lambda
    ),

    params_smooth = c(mu = FALSE, lambda = TRUE)
  )
}
