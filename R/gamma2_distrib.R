#' @include distrib.R generics.R
NULL

#' @title Gamma Distribution Class, Mean and Variance
#' @name Gamma2Distrib
#'
#' @description
#' The S7 class of the gamma family parametrized by its mean \eqn{\mu > 0} and
#' its variance \eqn{\sigma^2 > 0}, so that the shape is
#' \eqn{\alpha = \mu^2/\sigma^2} and the rate \eqn{\lambda = \mu/\sigma^2}. It
#' inherits from `continuous_distrib`, so it answers every generic of the
#' `distrib` contract; the eleven methods listed below are registered on it
#' directly and everything else comes from the parent.
#'
#' Build one with [gamma2_distrib()], which supplies the two link functions and
#' fills the properties in. This page documents the raw S7 constructor, which
#' takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Gamma2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [gamma2_distrib()] they hold `"gamma2"`, `"univariate"`,
#'   `c(0, Inf)`, `c("mu", "sigma2")`, the interpretations
#'   `c(mu = "mean", sigma2 = "variance")`, `2`, the domain \eqn{(0, \infty)}
#'   for both parameters, and the two links.
#'
#' @seealso [gamma2_distrib()] to build one;
#'   [gamma1_distrib()] for the same law in mean and dispersion, where the two
#'   parameters are orthogonal and these are not;
#'   [distrib_pdf.Gamma2Distrib()] and [distrib_gradient.Gamma2Distrib()] for
#'   the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Gamma2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gamma2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gamma2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gamma2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gamma2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gamma2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gamma2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gamma2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gamma2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gamma2Distrib],
#'   [`distrib_rng()`][distrib_rng.Gamma2Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- gamma2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # The second parameter is the variance itself, so it is what variance()
#' # returns and what the fitted standard error describes.
#' variance(d, list(mu = 3, sigma2 = 2))
Gamma2Distrib <- S7::new_class("Gamma2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gamma Probability Density Function in Mean and Variance
#' @name distrib_pdf.Gamma2Distrib
#' @description
#' Computes the gamma density
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{\lambda^{\alpha}}{\Gamma(\alpha)}
#'       y^{\alpha - 1} e^{-\lambda y}, \qquad
#'       \alpha = \dfrac{\mu^2}{\sigma^2}, \quad
#'       \lambda = \dfrac{\mu}{\sigma^2}, \quad y > 0,}
#' by calling [stats::dgamma()] at that shape and rate. With `log = TRUE` the
#' logarithm is formed inside `dgamma()` and stays finite where the density
#' itself underflows.
#'
#' The density is unbounded at the origin when \eqn{\sigma^2 > \mu^2}, where
#' the shape falls below one; it is flat there when \eqn{\sigma^2 = \mu^2},
#' the exponential case, and vanishes at the origin when
#' \eqn{\sigma^2 < \mu^2}.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of observations. The support is \eqn{(0, \infty)};
#'   a negative value gives 0.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_cdf.Gamma2Distrib()] for the distribution function,
#'   [distrib_gradient.Gamma2Distrib()] for the derivatives of the log-density,
#'   [distrib_pdf.Gamma1Distrib()] for the same density in the dispersion, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#'
#' # The method is stats::dgamma at shape mu^2/sigma2 and rate mu/sigma2.
#' all.equal(distrib_pdf(d, y, th),
#'           dgamma(y, shape = 9 / 2, rate = 3 / 2))
#'
#' # The same law as gamma1 at phi = sigma2/mu^2.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(gamma1_distrib(), y, list(mu = 3, phi = 2 / 9)))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(1, 3, 9), sigma2 = 2))
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e4, th)
#' distrib_pdf(d, 1e4, th, log = TRUE)
S7::method(distrib_pdf, Gamma2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dgamma(
    x = y,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    log = log
  )
}

#' @title Gamma Cumulative Distribution Function in Mean and Variance
#' @name distrib_cdf.Gamma2Distrib
#' @description
#' Computes the gamma distribution function, the regularized incomplete gamma
#' function
#' \deqn{F(q; \mu, \sigma^2) = \dfrac{\gamma(\alpha, \lambda q)}{\Gamma(\alpha)},
#'       \qquad \alpha = \dfrac{\mu^2}{\sigma^2}, \quad
#'       \lambda = \dfrac{\mu}{\sigma^2},}
#' with \eqn{\gamma} the lower incomplete gamma function, by calling
#' [stats::pgamma()] at that shape and rate. Both tails are available exactly:
#' `lower.tail = FALSE` evaluates \eqn{1 - F} without forming the difference,
#' and `log.p = TRUE` returns a logarithm that stays finite where the
#' probability itself underflows to zero.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param q A numeric vector of quantiles. A value at or below zero gives a
#'   lower-tail probability of 0.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma2))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Gamma2Distrib()] for the inverse,
#'   [distrib_pdf.Gamma2Distrib()] for the density, [distrib_grad_cdf()] for
#'   the derivatives of this function in the parameters, which the gamma takes
#'   by finite difference because the derivative of an incomplete gamma in its
#'   shape is hypergeometric, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' th <- list(mu = 3, sigma2 = 2)
#'
#' # The method is stats::pgamma at the implied shape and rate.
#' all.equal(distrib_cdf(d, c(1, 3, 5), th),
#'           pgamma(c(1, 3, 5), shape = 9 / 2, rate = 3 / 2))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 5, th) + distrib_cdf(d, 5, th, lower.tail = FALSE)
#'
#' # The gamma is right skewed, so less than half the mass lies below the mean.
#' distrib_cdf(d, 3, th)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 3000, th, lower.tail = FALSE)
#' distrib_cdf(d, 3000, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, Gamma2Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pgamma(
    q = q,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gamma Quantile Function in Mean and Variance
#' @name distrib_quantile.Gamma2Distrib
#' @description
#' Computes the gamma quantile function, the inverse of the regularized
#' incomplete gamma function in its argument, by calling [stats::qgamma()] at
#' shape \eqn{\alpha = \mu^2/\sigma^2} and rate \eqn{\lambda = \mu/\sigma^2}.
#' There is no elementary closed form; `qgamma()` inverts the distribution
#' function numerically. The gamma distribution function is strictly
#' increasing on \eqn{(0, \infty)}, so the round trip through
#' [distrib_cdf.Gamma2Distrib()] returns `p`.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities, which is how a quantile deep in a tail is
#'   requested without the probability underflowing. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty]}, of length
#'   `max(length(p), length(mu), length(sigma2))`.
#'
#' @seealso [distrib_cdf.Gamma2Distrib()], which this inverts;
#'   [distrib_rng.Gamma2Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' th <- list(mu = 3, sigma2 = 2)
#'
#' # A central 95 percent interval, visibly asymmetric about the mean of 3.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The median falls below the mean, the gamma being right skewed.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, Gamma2Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qgamma(
    p = p,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gamma Random Number Generator in Mean and Variance
#' @name distrib_rng.Gamma2Distrib
#' @description
#' Draws `n` independent gamma variates by calling [stats::rgamma()] at shape
#' \eqn{\alpha = \mu^2/\sigma^2} and rate \eqn{\lambda = \mu/\sigma^2}, so the
#' draws come from R's own gamma generator and depend on `.Random.seed` in the
#' usual way. The ratio-of-uniforms fallback the base class supplies is
#' bypassed.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. Both
#'   must be strictly positive.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_quantile.Gamma2Distrib()] for the inverse-transform route,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#'
#' # Same generator as stats::rgamma at the implied shape and rate.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 3, sigma2 = 2))
#' set.seed(2)
#' identical(a, rgamma(3, shape = 9 / 2, rate = 3 / 2))
#'
#' # The sample moments recover both parameters directly.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, sigma2 = 2))
#' c(mu = mean(z), sigma2 = var(z))
S7::method(distrib_rng, Gamma2Distrib) <- function(distrib, n, theta) {
  stats::rgamma(
    n = n,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]]
  )
}

#' @title Gamma Score in Mean and Variance
#' @name distrib_gradient.Gamma2Distrib
#' @description
#' Computes the first derivatives of the gamma log-density with respect to
#' \eqn{\mu} and \eqn{\sigma^2}, one value per observation, in closed form.
#' With \eqn{\alpha = \mu^2/\sigma^2} and \eqn{\lambda = \mu/\sigma^2},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} =
#'       \dfrac{-2\mu\psi(\alpha) + 2\mu\log\lambda + \mu + 2\mu\log y - y}
#'             {\sigma^2},}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} =
#'       -\dfrac{\mu\left\{-\mu\psi(\alpha) + \mu
#'         + \mu(\log\lambda + \log y) - y\right\}}{(\sigma^2)^2},}
#' with \eqn{\psi} the digamma function. Both components carry the digamma,
#' because both parameters move the shape; in [gamma1_distrib()], where the
#' second parameter is a dispersion, only the second component does.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma2`, each of
#'   length `max(length(y), length(mu), length(sigma2))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\sigma^2 > 0} the variance. \eqn{\psi} is the digamma function,
#' \eqn{\psi = (\log\Gamma)'}.
#'
#' @seealso [distrib_hessian.Gamma2Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Gamma2Distrib()] for their expectation,
#'   [distrib_gradient.Gamma1Distrib()] for the same score in the dispersion,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out with the digamma function.
#' al <- 9 / 2
#' lam <- 3 / 2
#' all.equal(g$mu,
#'           (-2 * 3 * digamma(al) + 2 * 3 * log(lam) + 3 +
#'              2 * 3 * log(y) - y) / 2)
#' all.equal(g$sigma2,
#'           -3 * (-3 * digamma(al) + 3 + 3 * (log(lam) + log(y)) - y) / 2^2)
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
S7::method(distrib_gradient, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  gamma_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Observed Hessian in Mean and Variance
#' @name distrib_hessian.Gamma2Distrib
#' @description
#' Computes the three distinct second derivatives of the gamma log-density
#' with respect to \eqn{\mu} and \eqn{\sigma^2}, one value per observation, in
#' closed form, by differentiating the score of
#' [distrib_gradient.Gamma2Distrib()] once more. Every component carries the
#' trigamma function of \eqn{\alpha = \mu^2/\sigma^2}, both parameters moving
#' the shape.
#'
#' The curvature in \eqn{\mu} is not negative at every data point, and neither
#' is the curvature in \eqn{\sigma^2}: the observed information of a gamma is
#' not positive definite everywhere. Its expectation is; see
#' [distrib_expected_hessian.Gamma2Distrib()].
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
#'   `mu_sigma2`, in that order, each of length
#'   `max(length(y), length(mu), length(sigma2))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives. \eqn{\psi}
#' and \eqn{\psi_1} are the digamma and trigamma functions.
#'
#' @seealso [distrib_gradient.Gamma2Distrib()] for the score,
#'   [distrib_expected_hessian.Gamma2Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Gamma2Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # The curvature in the variance is positive at y = 3, so the observed
#' # information is not positive definite at every observation.
#' h$sigma2_sigma2
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 3 + eps, sigma2 = 2))$mu
#' dn <- distrib_gradient(d, y, list(mu = 3 - eps, sigma2 = 2))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
S7::method(distrib_hessian, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  gamma_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Expected Hessian in Mean and Variance
#' @name distrib_expected_hessian.Gamma2Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. With \eqn{\alpha =
#' \mu^2/\sigma^2} and \eqn{\psi_1} the trigamma function,
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] =
#'         \dfrac{3\sigma^2 - 4\mu^2\psi_1(\alpha)}{(\sigma^2)^2}, \qquad
#'       \mathbb{E}\left[\ell^{(\sigma^2\sigma^2)}\right] =
#'         -\dfrac{\mu^2\left\{\mu^2\psi_1(\alpha) - \sigma^2\right\}}
#'                {(\sigma^2)^4},}
#' \deqn{\mathbb{E}\left[\ell^{(\mu\sigma^2)}\right] =
#'         \dfrac{2\mu\left\{\mu^2\psi_1(\alpha) - \sigma^2\right\}}
#'               {(\sigma^2)^3}.}
#' They follow from \eqn{\mathbb{E}[Y] = \mu} and
#' \eqn{\mathbb{E}[\log Y] = \psi(\alpha) - \log\lambda}.
#'
#' **The mixed entry does not vanish**, so the mean and the variance are not
#' orthogonal in this parametrization; their maximum likelihood estimates are
#' asymptotically correlated. [gamma1_distrib()] carries the same law with a
#' dispersion in place of the variance, and there the mixed entry is exactly
#' zero. The matrix is negative definite throughout, so the information it
#' negates is positive definite.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available. Accepted so that the
#'   signature matches the generic's, where it selects between the Bartlett,
#'   quadrature, Monte Carlo and outer-product routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
#'   `mu_sigma2`, in that order, each of length
#'   `max(length(y), length(mu), length(sigma2))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The gamma is a
#' regular family, so the second Bartlett identity holds and this equals the
#' variance of the score.
#'
#' @seealso [distrib_hessian.Gamma2Distrib()] for the observed quantity this is
#'   the expectation of, [distrib_expected_hessian.Gamma1Distrib()] for the
#'   orthogonal parametrization of the same law, [fisher_scoring()], which
#'   inverts it at each step, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' th <- list(mu = 3, sigma2 = 2)
#'
#' # The three constants, one value per observation.
#' eh <- lapply(distrib_expected_hessian(d, c(1, 3, 5), th), unique)
#' eh
#'
#' # Written out with the trigamma function.
#' al <- 9 / 2
#' c((3 * 2 - 4 * 9 * trigamma(al)) / 2^2,
#'   -9 * (9 * trigamma(al) - 2) / 2^4,
#'   2 * 3 * (9 * trigamma(al) - 2) / 2^3)
#'
#' # The mixed entry is not zero: the mean and the variance are correlated.
#' # In gamma1, the same law in the dispersion, it is exactly zero.
#' c(gamma2 = eh$mu_sigma2,
#'   gamma1 = distrib_expected_hessian(gamma1_distrib(), 0,
#'                                     list(mu = 3, phi = 2 / 9))$mu_phi)
#'
#' # Negative definite all the same, so the information is positive definite.
#' M <- matrix(c(eh$mu_mu, eh$mu_sigma2, eh$mu_sigma2, eh$sigma2_sigma2), 2)
#' eigen(M, only.values = TRUE)$values
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_expected_hessian, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  gamma_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Third-Order Derivatives in Mean and Variance
#' @name distrib_deriv3.Gamma2Distrib
#' @description
#' Computes the four distinct third derivatives of the gamma log-density with
#' respect to \eqn{\mu} and \eqn{\sigma^2}, in closed form. They come from
#' differentiating the log-density in the shape \eqn{\alpha = \mu^2/\sigma^2}
#' and the rate \eqn{\lambda = \mu/\sigma^2} and carrying the result across by
#' the chain rule, so each polygamma function of \eqn{\alpha} is evaluated
#' once. Every component of this order carries \eqn{\psi_2}, the second
#' derivative of the digamma function.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{Y} with \eqn{\mu} and \eqn{\log Y} with \eqn{\psi(\alpha) -
#' \log\lambda}. Both routes are closed form, so no quadrature is run and
#' `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available for both the observed
#'   and the expected values.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma2`,
#'   `mu_sigma2_sigma2` and `sigma2_sigma2_sigma2`, each of length
#'   `max(length(y), length(mu), length(sigma2))`. The names enumerate the
#'   distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\psi} is the digamma function and
#' \eqn{\psi_m} its \eqn{m}th derivative.
#'
#' @seealso [distrib_hessian.Gamma2Distrib()] for the order below and
#'   [distrib_deriv4.Gamma2Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic and for the numerical route a family
#'   without a closed form takes.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Every component varies with the observation.
#' d3$mu_mu_mu
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 3 + eps, sigma2 = 2))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 3 - eps, sigma2 = 2))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
S7::method(distrib_deriv3, Gamma2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) gamma_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gamma_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Fourth-Order Derivatives in Mean and Variance
#' @name distrib_deriv4.Gamma2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the gamma log-density with
#' respect to \eqn{\mu} and \eqn{\sigma^2}, in closed form, by the same chain
#' rule from the shape \eqn{\alpha = \mu^2/\sigma^2} and the rate
#' \eqn{\lambda = \mu/\sigma^2} that the lower orders use. Every component
#' carries \eqn{\psi_3}, the third derivative of the digamma function.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{Y} with \eqn{\mu} and \eqn{\log Y} with \eqn{\psi(\alpha) -
#' \log\lambda}. Both routes are closed form, so `approx` and `nsim` are
#' ignored.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available for both the observed
#'   and the expected values.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_sigma2`, `mu_mu_sigma2_sigma2`, `mu_sigma2_sigma2_sigma2` and
#'   `sigma2_sigma2_sigma2_sigma2`, each of length
#'   `max(length(y), length(mu), length(sigma2))`.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\psi} is the digamma
#' function and \eqn{\psi_m} its \eqn{m}th derivative.
#'
#' @seealso [distrib_deriv3.Gamma2Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.Gamma2Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 3 + eps, sigma2 = 2))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 3 - eps, sigma2 = 2))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
S7::method(distrib_deriv4, Gamma2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) gamma_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gamma_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma First Derivative in the Response, Mean and Variance
#' @name distrib_grad_y.Gamma2Distrib
#' @description
#' Computes the first derivative of the gamma log-density with respect to the
#' response, in closed form at the implied shape \eqn{\alpha = \mu^2/\sigma^2}
#' and rate \eqn{\lambda = \mu/\sigma^2}:
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\alpha - 1}{y} - \lambda.}
#' It changes sign at the mode \eqn{y = (\alpha-1)/\lambda}, so it is positive
#' below the mode and negative above it. When \eqn{\sigma^2 = \mu^2} the shape
#' is 1, the first term drops out and the derivative is the constant
#' \eqn{-1/\mu} of an exponential; when \eqn{\sigma^2 > \mu^2} the density has
#' no interior mode and the derivative is negative throughout.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of strictly positive observations. At `y = 0` the
#'   value is infinite unless the shape is exactly 1.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Gamma2Distrib()] for the second derivative in the
#'   response, [distrib_gradient.Gamma2Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#'
#' # Written out at the implied shape and rate.
#' all.equal(distrib_grad_y(d, y, th), (9 / 2 - 1) / y - 3 / 2)
#'
#' # Zero at the mode, positive below it and negative above.
#' mode <- (9 / 2 - 1) / (3 / 2)
#' c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#'
#' # At sigma2 = mu^2 the family is exponential and the derivative is constant.
#' distrib_grad_y(d, y, list(mu = 3, sigma2 = 9))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, Gamma2Distrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; s2 <- theta[[2]]
  (mu^2 / s2 - 1) / y - mu / s2
}

#' @title Gamma Second Derivative in the Response, Mean and Variance
#' @name distrib_hess_y.Gamma2Distrib
#' @description
#' Computes the second derivative of the gamma log-density with respect to the
#' response, in closed form at the implied shape \eqn{\alpha = \mu^2/\sigma^2}:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{\alpha - 1}{y^2}.}
#' The rate drops out, the log-density being linear in \eqn{y} apart from the
#' \eqn{(\alpha-1)\log y} term. The sign follows the shape: the log-density is
#' concave in the response when \eqn{\sigma^2 < \mu^2}, exactly flat when
#' \eqn{\sigma^2 = \mu^2}, where the family is exponential, and convex when
#' \eqn{\sigma^2 > \mu^2}.
#'
#' @param distrib A `Gamma2Distrib` object, from [gamma2_distrib()].
#' @param y A numeric vector of strictly positive observations. At `y = 0` the
#'   value is infinite unless the shape is exactly 1.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive; they enter only through the shape.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_grad_y.Gamma2Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Gamma2Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#'
#' all.equal(distrib_hess_y(d, y, th), -(9 / 2 - 1) / y^2)
#'
#' # Only the shape survives, so a setting with the same mu^2/sigma2 agrees.
#' all.equal(distrib_hess_y(d, y, th),
#'           distrib_hess_y(d, y, list(mu = 6, sigma2 = 8)))
#'
#' # Concave below sigma2 = mu^2, flat at it, convex above.
#' vapply(c(2, 9, 20),
#'        function(v) distrib_hess_y(d, 3, list(mu = 3, sigma2 = v)),
#'        numeric(1))
S7::method(distrib_hess_y, Gamma2Distrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; s2 <- theta[[2]]
  -(mu^2 / s2 - 1) / y^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Gamma Distribution, Mean and Variance
#'
#' @description
#' Builds the distribution object for the gamma family parametrized by its mean
#' \eqn{\mu > 0} and its variance \eqn{\sigma^2 > 0}, so that the shape is
#' \eqn{\alpha = \mu^2/\sigma^2} and the rate \eqn{\lambda = \mu/\sigma^2}. The
#' returned object carries closed-form derivatives of the log-density to fourth
#' order, in the parameters and in the response, and closed-form moments, so
#' every generic of the toolkit answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. Both default to the logarithm,
#' both parameters being positive.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted mean positive.
#' @param link_sigma2 A `link` object from `linkfunctions7` for the variance
#'   \eqn{\sigma^2}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, \infty)} is
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{\lambda^{\alpha}}{\Gamma(\alpha)}\,
#'       y^{\alpha-1} e^{-\lambda y}, \qquad
#'       \alpha = \dfrac{\mu^2}{\sigma^2}, \quad
#'       \lambda = \dfrac{\mu}{\sigma^2},}
#' the distribution function \eqn{F(q) = \gamma(\alpha, \lambda q)/\Gamma(\alpha)}
#' with \eqn{\gamma} the lower incomplete gamma function, and the quantile
#' function its numerical inverse. The mean is \eqn{\mu}, the variance
#' \eqn{\sigma^2}, the skewness \eqn{2\sqrt{\sigma^2}/\mu} and the excess
#' kurtosis \eqn{6\sigma^2/\mu^2}.
#'
#' This is the same law as [gamma1_distrib()], which carries the mean and a
#' *dispersion*, the two being related by \eqn{\sigma^2 = \phi\mu^2}. They are
#' separate families because the second parameter is a different quantity in
#' each, with its own interpretation, standard error and interval. The choice
#' between them is not only cosmetic: see the next section.
#'
#' # Derivatives, and orthogonality
#'
#' Writing \eqn{\psi} for the digamma function and \eqn{\psi_1} for the
#' trigamma, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} =
#'       \dfrac{-2\mu\psi(\alpha) + 2\mu\log\lambda + \mu + 2\mu\log y - y}
#'             {\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma^2} =
#'       -\dfrac{\mu\left\{-\mu\psi(\alpha) + \mu + \mu(\log\lambda + \log y)
#'         - y\right\}}{(\sigma^2)^2},}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] =
#'         \dfrac{3\sigma^2 - 4\mu^2\psi_1(\alpha)}{(\sigma^2)^2}, \qquad
#'       \mathbb{E}\left[\ell^{(\sigma^2\sigma^2)}\right] =
#'         -\dfrac{\mu^2\left\{\mu^2\psi_1(\alpha) - \sigma^2\right\}}
#'                {(\sigma^2)^4},}
#' \deqn{\mathbb{E}\left[\ell^{(\mu\sigma^2)}\right] =
#'         \dfrac{2\mu\left\{\mu^2\psi_1(\alpha) - \sigma^2\right\}}
#'               {(\sigma^2)^3}.}
#'
#' **The mixed entry does not vanish.** The mean and the variance are not
#' orthogonal here, so their estimates are asymptotically correlated and a
#' mean equation fitted with the variance held at a wrong value is biased. In
#' [gamma1_distrib()], where the second parameter is the dispersion
#' \eqn{\phi = \sigma^2/\mu^2}, the same entry is exactly zero. That is the
#' reason a generalized linear model uses the dispersion.
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.Gamma2Distrib()] and [distrib_deriv4.Gamma2Distrib()], as
#' are the derivatives in the response, [distrib_grad_y.Gamma2Distrib()] and
#' [distrib_hess_y.Gamma2Distrib()]. The derivatives of the *distribution*
#' function in the parameters have no elementary form, the derivative of an
#' incomplete gamma in its shape being hypergeometric, and are taken by finite
#' difference on the analytic cdf.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. The mean has
#' the closed-form estimate \eqn{\hat\mu = \bar y}; the variance has none and
#' is reached numerically, landing near but not at the sample variance. The
#' example below shows both.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\sigma^2 > 0} the variance. \eqn{\alpha} and \eqn{\lambda} are the
#' implied shape and rate. \eqn{\psi} is the digamma function and \eqn{\psi_m}
#' its \eqn{m}th derivative. \eqn{\eta} is a parameter on the unconstrained
#' scale of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Gamma2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"gamma2"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "sigma2")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [gamma1_distrib()] for the same law in the mean and the dispersion, which is
#' the orthogonal parametrization; [exponential_distrib()] for the case
#' \eqn{\sigma^2 = \mu^2}; [invgauss2_distrib()] and [lognormal2_distrib()] for
#' other positive families written in the mean and the variance;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [Gamma2Distrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 17.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dgamma pgamma qgamma rgamma
#'
#' @examples
#' d <- gamma2_distrib()
#' d
#'
#' # The density is stats::dgamma at shape mu^2/sigma2 and rate mu/sigma2.
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, sigma2 = 2)
#' all.equal(distrib_pdf(d, y, th), dgamma(y, shape = 9 / 2, rate = 3 / 2))
#'
#' # Moments: skewness 2 sqrt(sigma2)/mu, excess kurtosis 6 sigma2/mu^2.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(2 * sqrt(2) / 3, 6 * 2 / 9)
#'
#' # Fitting recovers the parameters. The mean is the sample mean exactly;
#' # the variance is close to the sample variance without equalling it.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted  = coef(fit),
#'       moments = c(mu = mean(z), sigma2 = var(z)))
#'
#' # The mean and the variance are correlated here, so the mixed entry of the
#' # expected Hessian is non-zero; in gamma1 it is exactly zero.
#' c(gamma2 = distrib_expected_hessian(d, 0, th)$mu_sigma2,
#'   gamma1 = distrib_expected_hessian(gamma1_distrib(), 0,
#'                                     list(mu = 3, phi = 2 / 9))$mu_phi)
#'
#' @export
gamma2_distrib <- function(link_mu = log_link(), link_sigma2 = log_link()) {
  
  Gamma2Distrib(
    distrib_name = "gamma2",
    dimension = "univariate",
    bounds = c(0, Inf),
    
    params = c("mu", "sigma2"),
    params_interpretation = c(mu = "mean", sigma2 = "variance"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(0, Inf),
      sigma2 = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      sigma2 = link_sigma2
    )
  )
  
}
