#' @include distrib.R generics.R
NULL

#' @title Exponential Distribution Class
#' @name ExponentialDistrib
#'
#' @description
#' The S7 class of the exponential family parametrized by its **mean**
#' \eqn{\mu > 0}, with density \eqn{f(y) = \mu^{-1}e^{-y/\mu}} on
#' \eqn{[0, \infty)}. It inherits from `continuous_distrib`, so it answers
#' every generic of the `distrib` contract; the eleven methods listed below are
#' registered on it directly.
#'
#' This is the only single-parameter continuous family in the package, so every
#' derivative array it returns has one component per order and its information
#' is a number.
#'
#' The parametrization is by the mean, not by the rate: R's own `dexp` takes
#' `rate`, and the methods pass `rate = 1/mu`. The mean is the scale of the
#' distribution, so \eqn{\mu} is also its standard deviation.
#'
#' Build one with [exponential_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `ExponentialDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [exponential_distrib()] they hold `"exponential"`,
#'   `"univariate"`, `c(0, Inf)`, `"mu"`, `c(mu = "mean")`, `1`, the domain
#'   \eqn{(0, \infty)}, and the one link.
#'
#' @seealso [exponential_distrib()] to build one;
#'   [gamma1_distrib()] and [weibull1_distrib()], both of which contain this
#'   family at a unit shape; [geometric_distrib()] for its discrete analogue.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.ExponentialDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.ExponentialDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.ExponentialDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ExponentialDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.ExponentialDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ExponentialDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.ExponentialDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ExponentialDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ExponentialDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ExponentialDistrib],
#'   [`distrib_rng()`][distrib_rng.ExponentialDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- exponential_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # One parameter, the mean, on the positive half line.
#' d@params
#' d@n_params
#' d@bounds
#'
#' # The mean is also the standard deviation, and the shape is fixed: the
#' # skewness is 2 and the excess kurtosis 6 at every mu.
#' th <- list(mu = 2)
#' c(mean = mean(d, th), sd = std_dev(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
ExponentialDistrib <- S7::new_class("ExponentialDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Exponential Probability Density Function
#' @name distrib_pdf.ExponentialDistrib
#' @description
#' Computes the exponential density
#' \deqn{f(y; \mu) = \dfrac{1}{\mu}\exp\left(-\dfrac{y}{\mu}\right), \qquad y \ge 0,}
#' by calling [stats::dexp()] at `rate = 1/mu`. The parametrization here is by
#' the mean, so the reciprocal is taken inside the method; a reader passing
#' `rate` where `mu` is expected gets the reciprocal distribution.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations. The support is \eqn{[0, \infty)};
#'   [stats::dexp()] returns 0 for a negative value.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. A value of length 1 is recycled.
#'   `mu` must be strictly positive; a non-positive value gives `NaN` with a
#'   warning from [stats::dexp()].
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu))`, one value per observation.
#'
#' @seealso [distrib_cdf.ExponentialDistrib()] for the distribution function,
#'   [distrib_pdf.Gamma1Distrib()], which contains this at a unit shape, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' y <- c(0.3, 1.1, 4.0)
#'
#' # The method is stats::dexp at rate = 1/mu.
#' all.equal(distrib_pdf(d, y, list(mu = 2)), dexp(y, rate = 1 / 2))
#'
#' # The log-density is exactly linear in y, with slope -1/mu.
#' diff(distrib_pdf(d, c(1, 2, 3), list(mu = 2), log = TRUE))
#' -1 / 2
#'
#' # The maximum is at the origin and equals 1/mu.
#' distrib_pdf(d, 0, list(mu = 2))
S7::method(distrib_pdf, ExponentialDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dexp(y, rate = 1 / theta[[1]], log = log)
}

#' @title Exponential Cumulative Distribution Function
#' @name distrib_cdf.ExponentialDistrib
#' @description
#' Computes the exponential distribution function
#' \deqn{F(q; \mu) = 1 - \exp\left(-\dfrac{q}{\mu}\right), \qquad q \ge 0,}
#' by calling [stats::pexp()] at `rate = 1/mu`. With `lower.tail = FALSE` the
#' survival function \eqn{\exp(-q/\mu)} is returned directly, without forming
#' the difference, so it stays exact far into the tail; combined with
#' `log.p = TRUE` it is simply \eqn{-q/\mu}.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param q A numeric vector of quantiles. A negative value gives 0.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `q`. `mu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are the survival
#'   function \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu))`.
#'
#' @seealso [distrib_quantile.ExponentialDistrib()] for the inverse,
#'   [distrib_pdf.ExponentialDistrib()] for the density, and [distrib_cdf()]
#'   for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' th <- list(mu = 2)
#'
#' # The method is stats::pexp at rate = 1/mu.
#' all.equal(distrib_cdf(d, c(0.3, 1.1, 4.0), th),
#'           pexp(c(0.3, 1.1, 4.0), rate = 1 / 2))
#'
#' # The survival function on the log scale is exactly -q/mu.
#' distrib_cdf(d, c(10, 100, 1000), th, lower.tail = FALSE, log.p = TRUE)
#' -c(10, 100, 1000) / 2
#'
#' # Memoryless: the chance of surviving another unit does not depend on how
#' # long the wait has already been.
#' c(distrib_cdf(d, 4, th, lower.tail = FALSE) /
#'     distrib_cdf(d, 3, th, lower.tail = FALSE),
#'   distrib_cdf(d, 1, th, lower.tail = FALSE))
S7::method(distrib_cdf, ExponentialDistrib) <- function(distrib, q, theta,
                                                        lower.tail = TRUE,
                                                        log.p = FALSE) {
  stats::pexp(q, rate = 1 / theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Exponential Quantile Function
#' @name distrib_quantile.ExponentialDistrib
#' @description
#' Computes the exponential quantile function
#' \deqn{Q(p; \mu) = -\mu \log(1 - p)}
#' by calling [stats::qexp()] at `rate = 1/mu`. The median is
#' \eqn{\mu \log 2 \approx 0.693\mu}, below the mean, the distribution being
#' right skewed at every parameter value.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. `p = 1` gives `Inf`; a value outside the
#'   range gives `NaN` with a warning.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `p`. `mu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}, and the quantile is
#'   then \eqn{-\mu\log p}, which is exact deep in the tail.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty)}, of length
#'   `max(length(p), length(mu))`.
#'
#' @seealso [distrib_cdf.ExponentialDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' th <- list(mu = 2)
#'
#' # The median is mu log 2, below the mean.
#' c(median = distrib_quantile(d, 0.5, th), mu_log2 = 2 * log(2), mean = 2)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # Asked from the upper tail the quantile is -mu log p, exact however small.
#' distrib_quantile(d, 1e-300, th, lower.tail = FALSE)
#' -2 * log(1e-300)
S7::method(distrib_quantile, ExponentialDistrib) <- function(distrib, p, theta,
                                                             lower.tail = TRUE,
                                                             log.p = FALSE) {
  stats::qexp(p, rate = 1 / theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Exponential Random Number Generator
#' @name distrib_rng.ExponentialDistrib
#' @description
#' Draws `n` independent exponential variates by calling [stats::rexp()] at
#' `rate = 1/mu`, so the draws come from R's own generator and depend on
#' `.Random.seed` in the usual way. The generalized ratio-of-uniforms fallback
#' the base class supplies is bypassed.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of length `n`. A value of length 1 is recycled, so a vector
#'   of length `n` draws one variate per mean. `mu` must be strictly positive.
#'
#' @return A numeric vector of `n` non-negative draws.
#'
#' @seealso [distrib_quantile.ExponentialDistrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the mean back, and [distrib_rng()] for
#'   the generic.
#'
#' @examples
#' d <- exponential_distrib()
#'
#' # Same generator as stats::rexp, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 2))
#' set.seed(2)
#' identical(a, rexp(3, rate = 1 / 2))
#'
#' # The sample mean and standard deviation both estimate mu.
#' set.seed(21)
#' z <- distrib_rng(d, 2e4, list(mu = 3))
#' c(mean = mean(z), sd = sd(z))
S7::method(distrib_rng, ExponentialDistrib) <- function(distrib, n, theta) {
  stats::rexp(n, rate = 1 / theta[[1]])
}

#' @title Exponential Score
#' @name distrib_gradient.ExponentialDistrib
#' @description
#' Computes the first derivative of the exponential log-density with respect to
#' the mean, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}.}
#' The log-density is \eqn{-\log\mu - y/\mu}, so the score is proportional to
#' the residual and its sum vanishes exactly at \eqn{\hat\mu = \bar y}. The
#' family has one parameter, so the returned list has one component.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{\mu > 0} the mean,
#' which is also the standard deviation of this family.
#'
#' @seealso [distrib_hessian.ExponentialDistrib()] for the second derivative,
#'   [distrib_expected_hessian.ExponentialDistrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' y <- c(0.3, 1.1, 4.0)
#' th <- list(mu = 2)
#'
#' # The closed form, written out.
#' all.equal(distrib_gradient(d, y, th)$mu, (y - 2) / 2^2)
#'
#' # One component, the family having one parameter.
#' names(distrib_gradient(d, y, th))
#'
#' # The summed score vanishes at the sample mean, which is the estimate.
#' set.seed(21)
#' z <- distrib_rng(d, 1000, list(mu = 3))
#' sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
#'
#' # On the link scale the component is multiplied by h' = mu, the derivative
#' # of the inverse log link.
#' distrib_gradient(d, y, th, scale = "link")$mu / distrib_gradient(d, y, th)$mu
S7::method(distrib_gradient, ExponentialDistrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ..., threads = 1L) {
  exponential_gradient_cpp(y, theta[[1]], threads)
}

#' @title Exponential Observed Hessian
#' @name distrib_hessian.ExponentialDistrib
#' @description
#' Computes the second derivative of the exponential log-density with respect
#' to the mean, one value per observation, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} - \dfrac{2y}{\mu^3}
#'       = \dfrac{\mu - 2y}{\mu^3}.}
#' It is positive wherever \eqn{y < \mu/2}, so a single observation below half
#' the mean contributes convexity; summed over a sample the curvature is
#' negative at the estimate, where \eqn{\bar y = \mu}. The expectation is
#' [distrib_expected_hessian.ExponentialDistrib()].
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @seealso [distrib_gradient.ExponentialDistrib()] for the score,
#'   [distrib_expected_hessian.ExponentialDistrib()] for the expectation of
#'   this quantity, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' y <- c(0.3, 1.1, 4.0)
#' th <- list(mu = 2)
#' h <- distrib_hessian(d, y, th)
#'
#' # The closed form, written out.
#' all.equal(h$mu_mu, (2 - 2 * y) / 2^3)
#'
#' # Positive below mu/2, negative above it.
#' data.frame(y = y, mu_mu = h$mu_mu, below_half = y < 2 / 2)
#'
#' # A central difference of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 2 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 2 - eps))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
S7::method(distrib_hessian, ExponentialDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  exponential_hessian_cpp(y, theta[[1]], threads)
}

#' @title Exponential Expected Hessian
#' @name distrib_expected_hessian.ExponentialDistrib
#' @description
#' Returns the expectation of the observed second derivative under the model,
#' in closed form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu^2},}
#' which follows from \eqn{\mathbb{E}[Y] = \mu} substituted into
#' \eqn{(\mu - 2y)/\mu^3}. The Fisher information for one observation is
#' \eqn{1/\mu^2}, so the asymptotic standard error of \eqn{\hat\mu} is
#' \eqn{\mu/\sqrt n}, which is the standard error of a sample mean of variance
#' \eqn{\mu^2}.
#'
#' Because the value does not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu))` and constant at \eqn{-1/\mu^2}.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' **observed information** under the model. The exponential is a regular
#' family, so the second Bartlett identity holds and this equals the variance
#' of the score.
#'
#' @seealso [distrib_hessian.ExponentialDistrib()] for the observed quantity
#'   this is the expectation of, [fisher_scoring()], which inverts it, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' th <- list(mu = 2)
#'
#' # A single number, -1/mu^2.
#' unique(distrib_expected_hessian(d, c(0.3, 1.1, 4.0), th)$mu_mu)
#' -1 / 2^2
#'
#' # The observed value averages onto it over a large sample.
#' set.seed(21)
#' z <- distrib_rng(d, 2e5, th)
#' c(observed = mean(distrib_hessian(d, z, th)$mu_mu),
#'   expected = distrib_expected_hessian(d, 0, th)$mu_mu)
#'
#' # It is the variance of the score, this family being regular.
#' c(var_of_score = mean(distrib_gradient(d, z, th)$mu^2),
#'   information = -distrib_expected_hessian(d, 0, th)$mu_mu)
S7::method(distrib_expected_hessian, ExponentialDistrib) <- function(distrib, y, theta,
                                                                     scale = c("parameter", "link"),
                                                                     approx = c("bartlett", "integrate", "mc", "opg"),
                                                                     nsim = 10000, ..., threads = 1L) {
  exponential_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Exponential Third-Order Derivative
#' @name distrib_deriv3.ExponentialDistrib
#' @description
#' Computes the third derivative of the exponential log-density with respect to
#' the mean, in closed form:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = -\dfrac{2}{\mu^3} + \dfrac{6y}{\mu^4}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}, which gives \eqn{4/\mu^3}. Both routes are closed
#' form, so no quadrature is run and `approx` and `nsim` are ignored.
#'
#' The family has one parameter, so there is one component, where a
#' two-parameter family carries four at this order.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the observed value. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available for both the observed
#'   and the expected value.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu_mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell^{(\mu\mu\mu)}} is the third derivative of the log-density with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives; a
#' subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.ExponentialDistrib()] for the order below and
#'   [distrib_deriv4.ExponentialDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' y <- c(0.3, 1.1, 4.0)
#' th <- list(mu = 2)
#'
#' # The closed form, written out.
#' all.equal(distrib_deriv3(d, y, th)$mu_mu_mu, -2 / 2^3 + 6 * y / 2^4)
#'
#' # The expected value is 4/mu^3.
#' unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#' 4 / 2^3
#'
#' # A central difference of the Hessian reproduces the observed value.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 2 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 2 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv3, ExponentialDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ..., threads = 1L) {
  if (expected) exponential_deriv3_expected_cpp(y, theta[[1]], threads)
  else exponential_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Exponential Fourth-Order Derivative
#' @name distrib_deriv4.ExponentialDistrib
#' @description
#' Computes the fourth derivative of the exponential log-density with respect
#' to the mean, in closed form:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4} = \dfrac{6}{\mu^4} - \dfrac{24y}{\mu^5}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}, which gives \eqn{-18/\mu^4}. Both routes are
#' closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the observed value. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available for both the observed
#'   and the expected value.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu_mu_mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell^{(\mu\mu\mu\mu)}} is the fourth derivative of the log-density with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.ExponentialDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' y <- c(0.3, 1.1, 4.0)
#' th <- list(mu = 2)
#'
#' # The closed form, written out.
#' all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu, 6 / 2^4 - 24 * y / 2^5)
#'
#' # The expected value is -18/mu^4.
#' unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
#' -18 / 2^4
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 2 + eps))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 2 - eps))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, ExponentialDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ..., threads = 1L) {
  if (expected) exponential_deriv4_expected_cpp(y, theta[[1]], threads)
  else exponential_deriv4_cpp(y, theta[[1]], threads)
}

#' @title Exponential First Derivative in the Response
#' @name distrib_grad_y.ExponentialDistrib
#' @description
#' Returns \eqn{-1/\mu} for every observation. The log-density is
#' \eqn{-\log\mu - y/\mu}, which is **linear** in \eqn{y}, so its derivative in
#' the response is the constant slope \eqn{-1/\mu} and carries no information
#' about where the observation fell.
#'
#' That linearity is the memorylessness of the family stated on the log scale:
#' the log survival function is \eqn{-q/\mu}, a straight line, so the hazard is
#' constant.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `length(y)`, every entry \eqn{-1/\mu}.
#'
#' @seealso [distrib_hess_y.ExponentialDistrib()] for the second derivative,
#'   which is zero; [distrib_gradient.ExponentialDistrib()] for the score in
#'   the mean; [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' y <- c(0.3, 1.1, 4.0)
#' th <- list(mu = 2)
#'
#' distrib_grad_y(d, y, th)
#' -1 / 2
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#'
#' # The same slope appears in the log survival function: a constant hazard.
#' diff(distrib_cdf(d, c(1, 2, 3), th, lower.tail = FALSE, log.p = TRUE))
S7::method(distrib_grad_y, ExponentialDistrib) <- function(distrib, y, theta, ...) {
  rep(-1 / theta[[1]], length.out = length(y))
}

#' @title Exponential Second Derivative in the Response
#' @name distrib_hess_y.ExponentialDistrib
#' @description
#' Returns zero for every observation. The exponential log-density is linear in
#' \eqn{y} on the whole support, so its second derivative in the response
#' vanishes everywhere; there is no kink to qualify the statement, unlike the
#' Laplace, whose zero holds only away from a point.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with the single component `mu`, which is not read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length `length(y)`.
#'
#' @seealso [distrib_grad_y.ExponentialDistrib()] for the constant first
#'   derivative; [distrib_hess_y.LaplaceDistrib()], which is zero for a
#'   different reason; [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- exponential_distrib()
#' th <- list(mu = 2)
#'
#' distrib_hess_y(d, c(0.3, 1.1, 4.0), th)
#'
#' # The first derivative is constant, so its difference is zero everywhere,
#' # with no exceptional point.
#' eps <- 1e-6
#' y <- c(1e-8, 0.5, 5, 50)
#' (distrib_grad_y(d, y + eps, th) - distrib_grad_y(d, y - eps, th)) / (2 * eps)
S7::method(distrib_hess_y, ExponentialDistrib) <- function(distrib, y, theta, ...) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Exponential Distribution
#'
#' @description
#' Builds the distribution object for the exponential family parametrized by
#' its **mean** \eqn{\mu > 0}, with density \eqn{\mu^{-1}e^{-y/\mu}} on
#' \eqn{[0, \infty)}. The returned object carries closed-form derivatives of
#' the log-density to fourth order and closed-form moments.
#'
#' The family has one parameter and a fixed shape: the standard deviation
#' equals the mean, the skewness is 2 and the excess kurtosis is 6 whatever
#' \eqn{\mu} is. It is the only continuous law with a constant hazard, which is
#' the memorylessness the details set out.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in [0, \infty)} is
#' \deqn{f(y; \mu) = \dfrac{1}{\mu}\exp\left(-\dfrac{y}{\mu}\right),}
#' with \eqn{\mu \in (0, \infty)}. The distribution function is
#' \eqn{F(q) = 1 - e^{-q/\mu}} and the quantile function
#' \eqn{Q(p) = -\mu\log(1-p)}.
#'
#' The mean and the standard deviation are both \eqn{\mu}, so the coefficient
#' of variation is 1; the skewness is 2 and the excess kurtosis 6, neither
#' depending on the parameter. The median is \eqn{\mu\log 2}, below the mean.
#'
#' R parametrizes its own `dexp` by the **rate**, \eqn{1/\mu}, and the methods
#' convert. [laplace2_distrib()] is the corresponding choice for the Laplace,
#' where the rate is what the package carries.
#'
#' # Memorylessness
#'
#' The survival function is \eqn{P(Y > q) = e^{-q/\mu}}, so
#' \eqn{P(Y > s + t \mid Y > s) = P(Y > t)} for every \eqn{s, t \ge 0}: a wait
#' already endured tells nothing about the wait remaining. Equivalently the
#' hazard \eqn{f/(1-F)} is the constant \eqn{1/\mu}. On the log scale this is
#' the statement that both \eqn{\log f} and \eqn{\log(1-F)} are straight lines
#' in \eqn{y}, which is why [distrib_grad_y.ExponentialDistrib()] is a constant
#' and [distrib_hess_y.ExponentialDistrib()] is zero.
#'
#' The exponential is the only continuous law with this property, so a fitted
#' exponential is also the statement that the process has no memory. When it
#' does, [weibull1_distrib()] and [gamma1_distrib()] both contain this family
#' at a unit shape and let the hazard rise or fall.
#'
#' # Derivatives
#'
#' Writing \eqn{\ell = -\log\mu - y/\mu}, the four orders are
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\mu^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{\mu - 2y}{\mu^3}, \qquad
#'       \dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{6y - 2\mu}{\mu^4}, \qquad
#'       \dfrac{\partial^4 \ell}{\partial \mu^4} = \dfrac{6\mu - 24y}{\mu^5},}
#' and their expectations follow by substituting \eqn{\mathbb{E}[Y] = \mu}:
#' \eqn{0}, \eqn{-1/\mu^2}, \eqn{4/\mu^3} and \eqn{-18/\mu^4}. The information
#' is \eqn{1/\mu^2}, so the asymptotic standard error of \eqn{\hat\mu} is
#' \eqn{\mu/\sqrt n}.
#'
#' # Estimation
#'
#' The maximum likelihood estimate is the sample mean, \eqn{\hat\mu = \bar y},
#' in closed form. [fit_distrib()] reaches it on the link scale; the example
#' below checks it.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{\mu > 0} the mean,
#' which is also the standard deviation and the reciprocal of the rate. The
#' **hazard** is \eqn{f(y)/(1 - F(y))}, the instantaneous failure rate given
#' survival to \eqn{y}.
#'
#' @return An S7 object of class `ExponentialDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"exponential"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`,
#'   `params_interpretation` `c(mu = "mean")`, `n_params` `1`, `params_bounds`
#'   the list of \eqn{(0, \infty)}, and `link_params` the one link given here.
#'
#' @seealso
#' [gamma1_distrib()] and [weibull1_distrib()], which contain this family at a
#' unit shape and let the hazard vary; [geometric_distrib()], its
#' memoryless discrete counterpart; [gpd_distrib()], which contains it at a
#' zero shape; [fixed()], which produces it from a larger family by holding a
#' shape; [fit_distrib()] to estimate the mean; [ExponentialDistrib] for the
#' class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 19.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dexp pexp qexp rexp
#'
#' @examples
#' d <- exponential_distrib()
#' d
#'
#' # The density is R's own at rate = 1/mu.
#' y <- c(0.3, 1.1, 4.0)
#' all.equal(distrib_pdf(d, y, list(mu = 2)), dexp(y, rate = 1 / 2))
#'
#' # A fixed shape: the standard deviation is the mean, and the skewness and
#' # excess kurtosis do not move with mu.
#' vapply(c(0.5, 2, 10), function(m) {
#'   th <- list(mu = m)
#'   c(sd = std_dev(d, th), skew = skewness(d, th), kurt = kurtosis(d, th))
#' }, numeric(3))
#'
#' # Memoryless: the chance of surviving one more unit is the same at every
#' # elapsed time.
#' vapply(c(0, 1, 5, 20), function(s)
#'   distrib_cdf(d, s + 1, list(mu = 2), lower.tail = FALSE) /
#'     distrib_cdf(d, s, list(mu = 2), lower.tail = FALSE), numeric(1))
#'
#' # The estimate is the sample mean, in closed form.
#' set.seed(21)
#' z <- distrib_rng(d, 1000, list(mu = 3))
#' c(fitted = unname(coef(fit_distrib(d, z))), sample_mean = mean(z))
#'
#' # It is a Weibull of unit shape, which fixed() produces from the larger
#' # family; the two densities agree exactly.
#' all.equal(distrib_pdf(d, y, list(mu = 2)),
#'           distrib_pdf(fixed(weibull1_distrib(), sigma = 1), y,
#'                       list(mu = 2)))
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
