#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Skew Normal Distribution
#' @name SkewNormal1Distrib
#'
#' @description A subclass of `continuous_distrib` representing Azzalini's
#' skew normal distribution: a gaussian with a shape parameter controlling the
#' asymmetry.
#' @inheritParams distrib
#' @return An object of class `SkewNormal1Distrib`.
#' @seealso [skewnormal1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.SkewNormal1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.SkewNormal1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.SkewNormal1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewNormal1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.SkewNormal1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.SkewNormal1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.SkewNormal1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.SkewNormal1Distrib],
#'   [`distrib_rng()`][distrib_rng.SkewNormal1Distrib],
#'   [`kurtosis()`][kurtosis],
#'   [`mean()`][mean.distrib],
#'   [`skewness()`][skewness],
#'   [`variance()`][variance]
#'
#' Everything else is inherited from [continuous_distrib()], including
#' the quantile function, which is obtained by root finding on the distribution
#' function, and the expected information, which has no closed form here.
SkewNormal1Distrib <- S7::new_class("SkewNormal1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Skew Normal Probability Density Function
#' @name distrib_pdf.SkewNormal1Distrib
#' @description
#' Computes the probability density function, with \eqn{z = (y-\mu)/\sigma}:
#' \deqn{f(y; \mu, \sigma, \alpha) = \dfrac{2}{\sigma}\,\phi(z)\,\Phi(\alpha z)}
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_pdf, SkewNormal1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - mu) / sigma
  # The two logarithms are taken separately, so the far tail of the skewed side
  # gives a large negative number rather than log(0).
  log_d <- log(2) - log(sigma) + stats::dnorm(z, log = TRUE) +
    stats::pnorm(alpha * z, log.p = TRUE)
  if (log) log_d else exp(log_d)
}

#' @title Skew Normal Cumulative Distribution Function
#' @name distrib_cdf.SkewNormal1Distrib
#' @description
#' Computes the distribution function through Owen's T function:
#' \deqn{F(q; \mu, \sigma, \alpha) = \Phi(z) - 2\,T(z, \alpha),
#'       \qquad z = (q - \mu)/\sigma.}
#' @details
#' The identity is Azzalini's. Evaluating it costs one bounded one-dimensional
#' quadrature per observation, which is cheaper and more accurate than the base
#' class's route of integrating the density over a semi-infinite range.
#' @param distrib A `SkewNormal1Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if `TRUE`, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_cdf, SkewNormal1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  z <- (q - theta[[1]]) / theta[[2]]
  res <- stats::pnorm(z) - 2 * numericals7::owen_t(z, theta[[3]])
  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Skew Normal Random Number Generator
#' @name distrib_rng.SkewNormal1Distrib
#' @description
#' Generates draws exactly, from the stochastic representation
#' \eqn{Z = \delta|U_0| + \sqrt{1-\delta^2}\,U_1} with \eqn{U_0, U_1}
#' independent standard normal and \eqn{\delta = \alpha/\sqrt{1+\alpha^2}}.
#' No inversion or rejection is involved.
#' @param distrib A `SkewNormal1Distrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @return A numeric vector of random draws.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_rng, SkewNormal1Distrib) <- function(distrib, n, theta) {
  alpha <- theta[[3]]
  delta <- alpha / sqrt(1 + alpha^2)
  z <- delta * abs(stats::rnorm(n)) + sqrt(1 - delta^2) * stats::rnorm(n)
  theta[[1]] + theta[[2]] * z
}

#' @title Skew Normal Analytical Gradient
#' @name distrib_gradient.SkewNormal1Distrib
#' @description
#' Closed-form first derivatives, with \eqn{z = (y-\mu)/\sigma},
#' \eqn{t = \alpha z} and \eqn{R = \phi(t)/\Phi(t)}:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{z - \alpha R}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{z^2 - 1 - \alpha z R}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \alpha} = z R.}
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_gradient, SkewNormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - mu) / sigma
  m <- numericals7::mills_ratio(alpha * z)
  list(
    mu = (z - alpha * m$r) / sigma,
    sigma = (z^2 - 1 - alpha * z * m$r) / sigma,
    alpha = z * m$r
  )
}

#' @title Skew Normal Analytical Observed Hessian
#' @name distrib_hessian.SkewNormal1Distrib
#' @description
#' Closed-form second derivatives, with \eqn{R' = -R(t + R)}:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{\alpha^2 R' - 1}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{\alpha^2 z R' - 2z + \alpha R}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \alpha}
#'         = -\dfrac{R + \alpha z R'}{\sigma},}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = \dfrac{1 - 3z^2 + 2\alpha z R + \alpha^2 z^2 R'}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma \, \partial \alpha}
#'         = -\dfrac{z R + \alpha z^2 R'}{\sigma},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \alpha^2} = z^2 R'.}
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_hessian, SkewNormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - mu) / sigma
  m <- numericals7::mills_ratio(alpha * z)
  r <- m$r
  dr <- m$dr
  s2 <- sigma^2
  list(
    mu_mu = (alpha^2 * dr - 1) / s2,
    sigma_sigma = (1 - 3 * z^2 + 2 * alpha * z * r + alpha^2 * z^2 * dr) / s2,
    alpha_alpha = z^2 * dr,
    mu_sigma = (alpha^2 * z * dr - 2 * z + alpha * r) / s2,
    mu_alpha = -(r + alpha * z * dr) / sigma,
    sigma_alpha = -(z * r + alpha * z^2 * dr) / sigma
  )
}

#' @title Skew Normal Analytical Third-Order Derivatives
#' @name distrib_deriv3.SkewNormal1Distrib
#' @description
#' Closed-form third-order derivatives of the skew normal log-density. With
#' \eqn{t = \alpha z} and \eqn{R} the inverse Mills ratio, the derivatives of
#' \eqn{\log \Phi(t)} follow from \eqn{R' = -R(t+R)} and stay polynomials in
#' \eqn{t} and \eqn{R}, so every component is elementary. The expected
#' derivatives have no closed form (the same integrals as the expected
#' information) and come from `expected_derivative()`.
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @param expected Logical; if `TRUE`, the expectation is approximated
#'   numerically.
#' @param approx Strategy for the expectation; see [distrib_deriv3()].
#' @param nsim Monte Carlo sample size when `approx = "mc"`.
#' @return A named list of third-derivative component vectors.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_deriv3, SkewNormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    skewnormal_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

#' @title Skew Normal Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.SkewNormal1Distrib
#' @description
#' Closed-form fourth-order derivatives of the skew normal log-density, in the
#' notation of [distrib_deriv3.SkewNormal1Distrib()]. The expected
#' derivatives are approximated numerically, as at third order.
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @param expected Logical; if `TRUE`, the expectation is approximated
#'   numerically.
#' @param approx Strategy for the expectation; see [distrib_deriv4()].
#' @param nsim Monte Carlo sample size when `approx = "mc"`.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_deriv4, SkewNormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    skewnormal_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

#' @title Skew Normal Response Derivative
#' @name distrib_grad_y.SkewNormal1Distrib
#' @description
#' Closed form: \eqn{\partial \ell / \partial y = (\alpha R - z)/\sigma}, which
#' is minus the derivative in \eqn{\mu}, as it must be for a location family.
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @return A numeric vector.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_grad_y, SkewNormal1Distrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - theta[[1]]) / sigma
  (alpha * numericals7::mills_ratio(alpha * z)$r - z) / sigma
}

#' @title Skew Normal Response Second Derivative
#' @name distrib_hess_y.SkewNormal1Distrib
#' @description
#' Closed form: \eqn{\partial^2 \ell / \partial y^2 = (\alpha^2 R' - 1)/\sigma^2}.
#' @param distrib A `SkewNormal1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma` and `alpha`.
#' @return A numeric vector.
#' @seealso [skewnormal1_distrib()]
S7::method(distrib_hess_y, SkewNormal1Distrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - theta[[1]]) / sigma
  (alpha^2 * numericals7::mills_ratio(alpha * z)$dr - 1) / sigma^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Skew Normal Distribution Object
#'
#' @description
#' Creates a distribution object for Azzalini's skew normal distribution, with
#' location \eqn{\mu}, scale \eqn{\sigma} and shape \eqn{\alpha}. The gaussian
#' is the special case \eqn{\alpha = 0}.
#'
#' @param link_mu A link function object for the location \eqn{\mu}. Defaults to
#'   [linkfunctions7::identity_link()].
#' @param link_sigma A link function object for the scale \eqn{\sigma}. Defaults
#'   to [linkfunctions7::log_link()] to ensure positivity.
#' @param link_alpha A link function object for the shape \eqn{\alpha}, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()].
#'
#' @details
#' **Probability density function**, with \eqn{z = (y-\mu)/\sigma}:
#' \deqn{f(y; \mu, \sigma, \alpha) = \dfrac{2}{\sigma}\,\phi(z)\,\Phi(\alpha z)}
#' The factor \eqn{2\Phi(\alpha z)} tilts the gaussian: it is above one where
#' \eqn{\alpha z > 0} and below one where \eqn{\alpha z < 0}, so positive
#' \eqn{\alpha} skews the density to the right. At \eqn{\alpha = 0} it is
#' identically one and the family reduces to the gaussian.
#'
#' **Cumulative distribution function:**
#' \deqn{F(q; \mu, \sigma, \alpha) = \Phi(z) - 2\,T(z, \alpha)}
#' with \eqn{T} Owen's T function; see [numericals7::owen_t()]. The quantile
#' function has no closed form and comes from the base class by root finding.
#'
#' **Score and observed Hessian** are closed form, written in the inverse
#' Mills ratio \eqn{R(t) = \phi(t)/\Phi(t)} at \eqn{t = \alpha z} and its
#' derivative \eqn{R' = -R(t+R)}; see
#' [distrib_gradient.SkewNormal1Distrib()] and
#' [distrib_hessian.SkewNormal1Distrib()].
#'
#' **Expected information.** There is no elementary closed form, so none
#' is registered and [distrib_expected_hessian()] approximates it by
#' the strategy named in `approx`, the default being the score variance.
#'
#' **Singularity at the gaussian.** At \eqn{\alpha = 0} the expected
#' information of this parametrization is **singular**: the derivative in
#' \eqn{\alpha} becomes collinear with the derivative in \eqn{\mu}, so the two
#' cannot be separated there. This is a property of the family and not of the
#' implementation, and it is why the profile log-likelihood in \eqn{\alpha} is
#' flat at the origin. A fit whose true shape is near zero will report a large
#' standard error for \eqn{\alpha}; the centered parametrization of Azzalini and
#' Capitanio removes the singularity and is a different object, not a
#' reparametrization this class performs.
#'
#' **Moments.** With \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and
#' \eqn{b = \sqrt{2/\pi}}, the mean is \eqn{\mu + \sigma b \delta} and the
#' variance \eqn{\sigma^2(1 - b^2\delta^2)}. The skewness is bounded: it lies in
#' \eqn{(-0.9953, 0.9953)} whatever \eqn{\alpha} is, which is the limitation of
#' the family and the reason the skew \eqn{t} exists.
#'
#' **Higher orders.** The observed third and fourth derivatives are
#' closed form, every derivative of \eqn{\log\Phi(t)} being a polynomial in
#' \eqn{t} and the inverse Mills ratio through \eqn{R' = -R(t+R)}; their
#' expected values share the obstruction of the expected information and are
#' approximated numerically.
#'
#' **Parameter Domains:**
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\alpha \in (-\infty, +\infty)}
#' }
#'
#' @return An S7 object of class [SkewNormal1Distrib()] (inheriting from
#'   `continuous_distrib`).
#'
#' @references
#' Azzalini, A. (1985). A class of distributions which includes the normal ones.
#' *Scandinavian Journal of Statistics* 12, 171-178.
#'
#' Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
#' Families*. Cambridge University Press.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dnorm pnorm rnorm integrate
#'
#' @examples
#' d <- skewnormal1_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1, alpha = 3)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' # shape zero is the gaussian
#' max(abs(distrib_pdf(d, c(-1, 0, 1), list(mu = 0, sigma = 1, alpha = 0)) -
#'         stats::dnorm(c(-1, 0, 1))))
#'
#' # the skewness the family can reach is bounded
#' c(alpha_3 = skewness(d, theta),
#'   alpha_50 = skewness(d, list(mu = 0, sigma = 1, alpha = 50)))
#'
#' @seealso [skewnormal2_distrib()], [skewt_distrib()]
#' @export
skewnormal1_distrib <- function(link_mu = identity_link(),
                               link_sigma = log_link(),
                               link_alpha = identity_link()) {
  SkewNormal1Distrib(
    distrib_name = "skew normal1",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "sigma", "alpha"),
    params_interpretation = c(mu = "location", sigma = "scale", alpha = "shape"),
    n_params = 3,

    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf),
      alpha = c(-Inf, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma,
      alpha = link_alpha
    )
  )
}
