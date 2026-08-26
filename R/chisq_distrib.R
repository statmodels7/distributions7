#' @include distrib.R generics.R
NULL

#' @title Chi-Squared Distribution Class
#' @name ChisqDistrib
#'
#' @description
#' The S7 class of the chi-squared family on \eqn{(0, \infty)}, parametrized by
#' its mean \eqn{\mu > 0}, which is the degrees of freedom. They are treated as
#' a continuous positive parameter, so the one-parameter family is estimable;
#' the variance is then \eqn{2\mu}. It
#' inherits from `continuous_distrib`, so it answers every generic of the
#' `distrib` contract; the eleven methods listed below are registered on it
#' directly and everything else comes from the parent.
#'
#' Build one with [chisq_distrib()], which supplies the link function and fills
#' the properties in. This page documents the raw S7 constructor, which takes
#' the parent's properties and validates none of the relationships between
#' them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `ChisqDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [chisq_distrib()] they hold `"chisq"`, `"univariate"`,
#'   `c(0, Inf)`, `"mu"`, the interpretation `c(mu = "mean")`, `1`, the domain
#'   \eqn{(0, \infty)}, and the one link.
#'
#' @seealso [chisq_distrib()] to build one;
#'   [gamma2_distrib()], which contains this family at
#'   \eqn{\sigma^2 = 2\mu}; [exponential_distrib()] for the case
#'   \eqn{\mu = 2}; [distrib_pdf.ChisqDistrib()] and
#'   [distrib_gradient.ChisqDistrib()] for the closed forms this class
#'   supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.ChisqDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.ChisqDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.ChisqDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ChisqDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ChisqDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.ChisqDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ChisqDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.ChisqDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ChisqDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ChisqDistrib],
#'   [`distrib_rng()`][distrib_rng.ChisqDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- chisq_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # One parameter, and it is the mean.
#' d@params
#' d@n_params
#' d@params_interpretation
#'
#' # The variance is tied to it at 2 mu, so the family has no free spread.
#' vapply(c(1, 4, 20), function(m) variance(d, list(mu = m)), numeric(1))
ChisqDistrib <- S7::new_class("ChisqDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Chi-Squared Probability Density Function
#' @name distrib_pdf.ChisqDistrib
#' @description
#' Computes the chi-squared density
#' \deqn{f(y; \mu) = \dfrac{y^{\mu/2 - 1} e^{-y/2}}{2^{\mu/2}\,\Gamma(\mu/2)},
#'       \qquad y > 0,}
#' by calling [stats::dchisq()] at `df = mu`. With `log = TRUE` the logarithm
#' is formed inside `dchisq()` and stays finite where the density itself
#' underflows.
#'
#' The density is unbounded at the origin for \eqn{\mu < 2}, flat there at
#' \eqn{\mu = 2}, where the family is the exponential with mean 2, and vanishes
#' at the origin for \eqn{\mu > 2}.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of observations. The support is
#'   \eqn{(0, \infty)}; a value at or below zero gives 0, `Inf` or a finite
#'   value according to \eqn{\mu}.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive and need not be a whole number.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu))`, one value per observation.
#'
#' @seealso [distrib_cdf.ChisqDistrib()] for the distribution function,
#'   [distrib_gradient.ChisqDistrib()] for the derivative of the log-density,
#'   [distrib_pdf.Gamma2Distrib()] for the family this sits inside, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#'
#' # The method is stats::dchisq at df = mu.
#' all.equal(distrib_pdf(d, y, th), dchisq(y, df = 4))
#'
#' # It is the gamma with sigma2 = 2 mu, and the exponential at mu = 2.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(gamma2_distrib(), y, list(mu = 4, sigma2 = 8)))
#' all.equal(distrib_pdf(d, y, list(mu = 2)),
#'           distrib_pdf(exponential_distrib(), y, list(mu = 2)))
#'
#' # The degrees of freedom need not be a whole number.
#' distrib_pdf(d, y, list(mu = 3.7))
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 2000, th)
#' distrib_pdf(d, 2000, th, log = TRUE)
S7::method(distrib_pdf, ChisqDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dchisq(y, df = theta[[1]], log = log)
}

#' @title Chi-Squared Cumulative Distribution Function
#' @name distrib_cdf.ChisqDistrib
#' @description
#' Computes the chi-squared distribution function, the regularized incomplete
#' gamma function at shape \eqn{\mu/2} and argument \eqn{q/2}, by calling
#' [stats::pchisq()] at `df = mu`. Both tails are available exactly:
#' `lower.tail = FALSE` evaluates \eqn{1 - F} without forming the difference,
#' and `log.p = TRUE` returns a logarithm that stays finite where the
#' probability itself underflows to zero.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param q A numeric vector of quantiles. A value at or below zero gives a
#'   lower-tail probability of 0.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `q`, recycled if of length 1. It must be
#'   strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu))`. With `log.p = TRUE` the values are
#'   logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.ChisqDistrib()] for the inverse,
#'   [distrib_pdf.ChisqDistrib()] for the density, [distrib_grad_cdf()] for the
#'   derivatives of this function in the parameter, which the chi-squared takes
#'   by finite difference because the derivative of an incomplete gamma in its
#'   shape is hypergeometric, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' th <- list(mu = 4)
#'
#' # The method is stats::pchisq at df = mu.
#' all.equal(distrib_cdf(d, c(1, 4, 9), th), pchisq(c(1, 4, 9), df = 4))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 9, th) + distrib_cdf(d, 9, th, lower.tail = FALSE)
#'
#' # Right skewed, so less than half the mass lies below the mean.
#' distrib_cdf(d, 4, th)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 2000, th, lower.tail = FALSE)
#' distrib_cdf(d, 2000, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, ChisqDistrib) <- function(distrib, q, theta,
                                                  lower.tail = TRUE,
                                                  log.p = FALSE, ...) {
  stats::pchisq(q, df = theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Chi-Squared Quantile Function
#' @name distrib_quantile.ChisqDistrib
#' @description
#' Computes the chi-squared quantile function, the inverse of the regularized
#' incomplete gamma function in its argument, by calling [stats::qchisq()] at
#' `df = mu`. There is no elementary closed form; `qchisq()` inverts the
#' distribution function numerically. The distribution function is strictly
#' increasing on \eqn{(0, \infty)}, so the round trip through
#' [distrib_cdf.ChisqDistrib()] returns `p`.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `p`, recycled if of length 1. It must be
#'   strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty]}, of length
#'   `max(length(p), length(mu))`.
#'
#' @seealso [distrib_cdf.ChisqDistrib()], which this inverts;
#'   [distrib_rng.ChisqDistrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' th <- list(mu = 4)
#'
#' # The critical values a test reads, here at four degrees of freedom.
#' distrib_quantile(d, c(0.9, 0.95, 0.99), th)
#'
#' # A central 95 percent interval, asymmetric about the mean of 4.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The median falls below the mean, the family being right skewed.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, ChisqDistrib) <- function(distrib, p, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  stats::qchisq(p, df = theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Chi-Squared Random Number Generator
#' @name distrib_rng.ChisqDistrib
#' @description
#' Draws `n` independent chi-squared variates by calling [stats::rchisq()] at
#' `df = mu`, so the draws come from R's own generator and depend on
#' `.Random.seed` in the usual way. The ratio-of-uniforms fallback the base
#' class supplies is bypassed.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of length `n`. A component of length 1 is recycled, so a
#'   vector of length `n` draws one variate per parameter setting. It must be
#'   strictly positive and need not be a whole number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_quantile.ChisqDistrib()] for the inverse-transform route,
#'   [fit_distrib()] to estimate the parameter back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#'
#' # Same generator as stats::rchisq at df = mu.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 4))
#' set.seed(2)
#' identical(a, rchisq(3, df = 4))
#'
#' # The sample mean recovers the degrees of freedom, and the variance is
#' # about twice it, the two being tied.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 4))
#' c(mean = mean(z), var = var(z))
S7::method(distrib_rng, ChisqDistrib) <- function(distrib, n, theta, ...) {
  stats::rchisq(n, df = theta[[1]])
}

#' @title Chi-Squared Score
#' @name distrib_gradient.ChisqDistrib
#' @description
#' Computes the derivative of the chi-squared log-density with respect to
#' \eqn{\mu}, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} =
#'       \dfrac{\log y - \log 2 - \psi(\mu/2)}{2},}
#' with \eqn{\psi} the digamma function. The family is a one-parameter
#' exponential family in \eqn{\log y}, so this is the sufficient statistic
#' minus its expectation, and **it is the only order that involves the
#' response at all**: every derivative from the second up is a polygamma
#' function of \eqn{\mu/2} and nothing else.
#'
#' The score has mean zero because
#' \eqn{\mathbb{E}[\log Y] = \psi(\mu/2) + \log 2}.
#'
#' With `scale = "link"` the generic applies the chain rule for the link the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of strictly positive observations. A value at
#'   zero makes the logarithm infinite and the score non-finite.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list with one numeric vector, `mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{\mu > 0} the mean,
#' which is also the degrees of freedom. \eqn{\psi} is the digamma function,
#' \eqn{\psi = (\log\Gamma)'}.
#'
#' @seealso [distrib_hessian.ChisqDistrib()] for the second derivative,
#'   [distrib_expected_hessian.ChisqDistrib()], which returns the same number,
#'   [distrib_grad_y.ChisqDistrib()] for the derivative in the response, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#' g <- distrib_gradient(d, y, th)
#'
#' # The closed form, written out with the digamma function.
#' all.equal(g$mu, (log(y) - log(2) - digamma(2)) / 2)
#'
#' # It is a sufficient statistic minus its expectation, so the sample mean of
#' # log y matches psi(mu/2) + log 2.
#' set.seed(9)
#' z <- distrib_rng(d, 2e5, th)
#' c(sample = mean(log(z)), theory = digamma(2) + log(2))
#'
#' # The score vanishes where log y equals that expectation.
#' distrib_gradient(d, exp(digamma(2) + log(2)), th)$mu
#'
#' # Summed over a fitted sample it is at the optimizer's tolerance.
#' set.seed(7)
#' zz <- distrib_rng(d, 2000, th)
#' sum(distrib_gradient(d, zz, as.list(coef(fit_distrib(d, zz))))$mu)
S7::method(distrib_gradient, ChisqDistrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"), ..., threads = 1L) {
  chisq_gradient_cpp(y, theta[[1]], threads)
}

#' @title Chi-Squared Observed Hessian
#' @name distrib_hessian.ChisqDistrib
#' @description
#' Computes the second derivative of the chi-squared log-density with respect
#' to \eqn{\mu}, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\psi'(\mu/2)}{4},}
#' with \eqn{\psi'} the trigamma function. **It does not involve the
#' response.** The family is a one-parameter exponential family in
#' \eqn{\log y}, so the data reach the log-density only through a term linear
#' in \eqn{\mu} and the second derivative kills it. The value is constant
#' within a parameter setting, is recycled to the length of `y`, and equals
#' [distrib_expected_hessian.ChisqDistrib()] exactly.
#'
#' The coincidence is on the **parameter** scale. On the link scale the two
#' differ, the second-order chain rule adding a term
#' \eqn{h''(\eta)\,\partial\ell/\partial\mu} that the expected version drops
#' and a sample does not; the example below measures it.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list with one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu))` and constant within itself when the parameter
#'   is.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#' \eqn{\psi'} is the trigamma function and \eqn{h(\eta)} the inverse link.
#'
#' @seealso [distrib_gradient.ChisqDistrib()] for the score, which is the one
#'   quantity here that reads the data;
#'   [distrib_expected_hessian.ChisqDistrib()], which returns the same number
#'   on the parameter scale; [distrib_deriv3.ChisqDistrib()] for the order
#'   above; and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#'
#' # A constant across the observations, written out with the trigamma.
#' distrib_hessian(d, y, th)
#' -trigamma(2) / 4
#'
#' # Free of the response, so it equals its own expectation to the bit.
#' identical(distrib_hessian(d, y, th), distrib_expected_hessian(d, y, th))
#'
#' # On the link scale they differ by h''(eta) times the score, which is what
#' # makes Fisher scoring and Newton's method take different steps.
#' obs <- distrib_hessian(d, y, th, scale = "link")$mu_mu
#' exp_ <- distrib_expected_hessian(d, y, th, scale = "link")$mu_mu
#' rbind(observed = obs, expected = exp_,
#'       difference = obs - exp_,
#'       h2_times_score = 4 * distrib_gradient(d, y, th)$mu)
S7::method(distrib_hessian, ChisqDistrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ..., threads = 1L) {
  chisq_hessian_cpp(y, theta[[1]], threads)
}

#' @title Chi-Squared Expected Hessian
#' @name distrib_expected_hessian.ChisqDistrib
#' @description
#' Returns the same number as [distrib_hessian.ChisqDistrib()],
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'       = -\dfrac{\psi'(\mu/2)}{4},}
#' the observed second derivative being free of the response and so equal to
#' its own expectation. Nothing is averaged, integrated or simulated, so
#' `approx` and `nsim` are ignored and `y` is read only for its length. The
#' value is negative at every \eqn{\mu}, the trigamma function being positive,
#' so the information is positive throughout.
#'
#' The identity holds on the **parameter** scale. On the link scale the
#' second-order chain rule adds \eqn{h''(\eta)\,\partial\ell/\partial\mu} to
#' the observed Hessian; the expected version drops that term because the score
#' has mean zero, and a finite sample does not. Fisher scoring and Newton's
#' method therefore take different steps here and agree at the optimum, where
#' the summed score vanishes.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, the expectation being exact. Accepted so that the
#'   signature matches the generic's, where it selects between the Bartlett,
#'   quadrature, Monte Carlo and outer-product routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list with one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu))` and constant within itself when the parameter
#'   is.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' **observed information** under the model. The chi-squared is a regular
#' family, so the second Bartlett identity holds and this equals the variance
#' of the score. \eqn{\psi'} is the trigamma function.
#'
#' @seealso [distrib_hessian.ChisqDistrib()], which returns the same number;
#'   [fisher_scoring()], which inverts it at each step; and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' th <- list(mu = 4)
#'
#' # A constant, and identical to the observed Hessian.
#' unique(distrib_expected_hessian(d, c(1, 4, 9), th)$mu_mu)
#' identical(distrib_expected_hessian(d, c(1, 4, 9), th),
#'           distrib_hessian(d, c(1, 4, 9), th))
#'
#' # Negative at every mu, the trigamma function being positive.
#' vapply(c(0.5, 4, 40),
#'        function(m) distrib_expected_hessian(d, 0, list(mu = m))$mu_mu,
#'        numeric(1))
#'
#' # The two fitting methods differ on the link scale and land together.
#' set.seed(7)
#' z <- distrib_rng(d, 2000, th)
#' c(newton = coef(fit_distrib(d, z, method = "newton")),
#'   fisher = coef(fit_distrib(d, z, method = "fisher")))
S7::method(distrib_expected_hessian, ChisqDistrib) <- function(distrib, y, theta,
                                                               scale = c("parameter", "link"),
                                                               approx = c("bartlett", "integrate", "mc", "opg"),
                                                               nsim = 10000, ..., threads = 1L) {
  chisq_hessian_cpp(y, theta[[1]], threads)
}

#' @title Chi-Squared Third-Order Derivative
#' @name distrib_deriv3.ChisqDistrib
#' @description
#' Computes the single third derivative of the chi-squared log-density with
#' respect to \eqn{\mu}, in closed form:
#' \deqn{\ell^{(\mu\mu\mu)} = -\dfrac{\psi''(\mu/2)}{8},}
#' with \eqn{\psi''} the second derivative of the digamma function. It is one
#' case of the general pattern
#' \eqn{\ell^{(k)} = -\psi^{(k-2)}(\mu/2)/2^{k}} for \eqn{k \ge 2}, which holds
#' because the family is a one-parameter exponential family in \eqn{\log y} and
#' the response leaves the derivatives after the first order.
#'
#' The value is free of the response, so `expected` selects nothing: the same
#' kernel runs either way and the two results are identical to the bit.
#' `approx` and `nsim` are ignored for the same reason.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected third derivatives being the same number. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list with one numeric vector, `mu_mu_mu`, of length
#'   `max(length(y), length(mu))` and constant within itself when the parameter
#'   is.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\psi} is the digamma function and
#' \eqn{\psi^{(m)}} its \eqn{m}th derivative.
#'
#' @seealso [distrib_hessian.ChisqDistrib()] for the order below and
#'   [distrib_deriv4.ChisqDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#'
#' # A constant, written out with the second derivative of the digamma.
#' unique(distrib_deriv3(d, y, th)$mu_mu_mu)
#' -psigamma(2, deriv = 2) / 8
#'
#' # Free of the response, so asking for the expectation changes nothing.
#' identical(distrib_deriv3(d, y, th), distrib_deriv3(d, y, th, expected = TRUE))
#'
#' # A central difference of the Hessian reproduces it.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 4 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 4 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
#'           tolerance = 1e-6)
S7::method(distrib_deriv3, ChisqDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ..., threads = 1L) {
  chisq_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Chi-Squared Fourth-Order Derivative
#' @name distrib_deriv4.ChisqDistrib
#' @description
#' Computes the single fourth derivative of the chi-squared log-density with
#' respect to \eqn{\mu}, in closed form:
#' \deqn{\ell^{(\mu\mu\mu\mu)} = -\dfrac{\psi'''(\mu/2)}{16},}
#' with \eqn{\psi'''} the third derivative of the digamma function, one case of
#' \eqn{\ell^{(k)} = -\psi^{(k-2)}(\mu/2)/2^{k}}.
#'
#' The value is free of the response, so `expected`, `approx` and `nsim` are
#' all without effect and the result is identical to the bit whichever is asked
#' for.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected fourth derivatives being the same number. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list with one numeric vector, `mu_mu_mu_mu`, of length
#'   `max(length(y), length(mu))` and constant within itself when the parameter
#'   is.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\psi} is the digamma
#' function and \eqn{\psi^{(m)}} its \eqn{m}th derivative.
#'
#' @seealso [distrib_deriv3.ChisqDistrib()] for the order below and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#'
#' # A constant, written out with the third derivative of the digamma.
#' unique(distrib_deriv4(d, y, th)$mu_mu_mu_mu)
#' -psigamma(2, deriv = 3) / 16
#'
#' # Free of the response, so asking for the expectation changes nothing.
#' identical(distrib_deriv4(d, y, th),
#'           distrib_deriv4(d, y, th, expected = TRUE))
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, y, list(mu = 4 + eps))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 4 - eps))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, ChisqDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ..., threads = 1L) {
  chisq_deriv4_cpp(y, theta[[1]], threads)
}

#' @title Chi-Squared First Derivative in the Response
#' @name distrib_grad_y.ChisqDistrib
#' @description
#' Computes the first derivative of the chi-squared log-density with respect to
#' the response, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\mu/2 - 1}{y}
#'       - \dfrac{1}{2}.}
#' It changes sign at \eqn{y = \mu - 2}, which is the mode of the density for
#' \eqn{\mu > 2}. At \eqn{\mu = 2} the first term drops out and the derivative
#' is the constant \eqn{-1/2} of an exponential with mean 2; below \eqn{\mu =
#' 2} the density has no interior mode and the derivative is negative
#' throughout.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of strictly positive observations. At `y = 0` the
#'   value is infinite unless \eqn{\mu} is exactly 2.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `max(length(y), length(mu))`, one value
#'   per observation.
#'
#' @seealso [distrib_hess_y.ChisqDistrib()] for the second derivative in the
#'   response, [distrib_gradient.ChisqDistrib()] for the derivative in the
#'   parameter, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#'
#' # Written out.
#' all.equal(distrib_grad_y(d, y, th), (4 / 2 - 1) / y - 0.5)
#'
#' # Zero at the mode mu - 2, positive below it and negative above.
#' c(mode = 4 - 2, at_mode = distrib_grad_y(d, 2, th))
#' distrib_grad_y(d, c(1, 9), th)
#'
#' # At mu = 2 the family is exponential and the derivative is constant.
#' distrib_grad_y(d, y, list(mu = 2))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, ChisqDistrib) <- function(distrib, y, theta, ...) {
  (theta[[1]] / 2 - 1) / y - 0.5
}

#' @title Chi-Squared Second Derivative in the Response
#' @name distrib_hess_y.ChisqDistrib
#' @description
#' Computes the second derivative of the chi-squared log-density with respect
#' to the response, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{\mu/2 - 1}{y^2}.}
#' The exponential term drops out, the log-density being linear in \eqn{y}
#' apart from the \eqn{(\mu/2 - 1)\log y} term. The sign follows \eqn{\mu}: the
#' log-density is concave in the response for \eqn{\mu > 2}, exactly flat at
#' \eqn{\mu = 2}, where the family is the exponential with mean 2, and convex
#' for \eqn{\mu < 2}.
#'
#' @param distrib A `ChisqDistrib` object, from [chisq_distrib()].
#' @param y A numeric vector of strictly positive observations. At `y = 0` the
#'   value is infinite unless \eqn{\mu} is exactly 2.
#' @param theta A named list with one component `mu`, a numeric vector of
#'   length 1 or of the length of `y`, recycled if of length 1. It must be
#'   strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `max(length(y), length(mu))`, one value
#'   per observation.
#'
#' @seealso [distrib_grad_y.ChisqDistrib()] for the first derivative in the
#'   response, [distrib_hessian.ChisqDistrib()] for the curvature in the
#'   parameter, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- chisq_distrib()
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#'
#' all.equal(distrib_hess_y(d, y, th), -(4 / 2 - 1) / y^2)
#'
#' # Concave above mu = 2, flat at it, convex below.
#' vapply(c(4, 2, 1), function(m) distrib_hess_y(d, 4, list(mu = m)),
#'        numeric(1))
#'
#' # Negative everywhere at four degrees of freedom.
#' all(distrib_hess_y(d, y, th) < 0)
S7::method(distrib_hess_y, ChisqDistrib) <- function(distrib, y, theta, ...) {
  -(theta[[1]] / 2 - 1) / (y * y)
}

# --- CONSTRUCTOR WRAPPER ---

#' Chi-Squared Distribution
#'
#' @description
#' Builds the distribution object for the chi-squared family on
#' \eqn{(0, \infty)}, parametrized by its mean \eqn{\mu > 0}, which is the
#' degrees of freedom. The returned object carries closed-form derivatives of
#' the log-density to fourth order, in the parameter and in the response, and
#' closed-form moments, so every generic of the toolkit answers without a
#' numerical fallback.
#'
#' The one argument chooses the link that carries the parameter to the
#' unconstrained scale an optimizer works on, and defaults to the logarithm.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The degrees of freedom are treated as a continuous positive parameter, so
#' the one-parameter family is estimable. The density on
#' \eqn{y \in (0, \infty)} is
#' \deqn{f(y; \mu) = \dfrac{y^{\mu/2 - 1} e^{-y/2}}{2^{\mu/2}\,\Gamma(\mu/2)},}
#' the mean is \eqn{\mu}, the variance \eqn{2\mu}, the skewness
#' \eqn{2\sqrt{2/\mu}} and the excess kurtosis \eqn{12/\mu}. The spread is
#' tied to the mean, so the family has no free scale.
#'
#' The law is a gamma with shape \eqn{\mu/2} and scale 2, and it is **not** a
#' gamma with a parameter held: this package writes the gamma in
#' \eqn{(\mu, \sigma^2)}, and a scale of 2 is the relation
#' \eqn{\sigma^2 = 2\mu} between two parameters rather than a value one of them
#' can be fixed at. At \eqn{\mu = 2} it is the exponential with mean 2.
#'
#' # The response leaves after the first derivative
#'
#' The family is a one-parameter exponential family in \eqn{\log y}, so
#' \deqn{\dfrac{\partial \ell}{\partial \mu}
#'         = \dfrac{\log y - \log 2 - \psi(\mu/2)}{2}, \qquad
#'       \ell^{(k)} = -\dfrac{\psi^{(k-2)}(\mu/2)}{2^{k}}, \quad k \ge 2,}
#' with \eqn{\psi} the digamma function. From the second order on nothing
#' involves the response at all. Two things follow. On the parameter scale the
#' observed information is exactly the expected information, and the same holds
#' at third and fourth order, so asking any of those methods for
#' `expected = TRUE` returns the same numbers. And
#' \eqn{\mathbb{E}[\log Y] = \psi(\mu/2) + \log 2} is what gives the score mean
#' zero.
#'
#' That coincidence does not carry to the scale a fit optimizes on. The
#' second-order chain rule adds a term
#' \eqn{h''(\eta)\,\partial\ell/\partial\mu} to the link-scale Hessian; the
#' expected version drops it because the score has mean zero, and a finite
#' sample does not. Fisher scoring and Newton's method therefore take different
#' steps here and agree at the optimum, where the summed score vanishes. The
#' example on [distrib_hessian.ChisqDistrib()] measures that difference.
#'
#' The derivatives of the *distribution* function in the parameter have no
#' elementary form, the derivative of an incomplete gamma in its shape being
#' hypergeometric, and are taken by finite difference on the analytic cdf.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. The estimate
#' solves \eqn{\psi(\hat\mu/2) = \overline{\log y} - \log 2} and has no closed
#' form; the sample mean is the method of moments starting value and lands
#' beside it.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{\mu > 0} the mean,
#' which is also the degrees of freedom. \eqn{\psi} is the digamma function and
#' \eqn{\psi^{(m)}} its \eqn{m}th derivative. \eqn{\eta} is the parameter on
#' the unconstrained scale of its link, with \eqn{\mu = h(\eta)}.
#'
#' @return An S7 object of class `ChisqDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"chisq"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`, `n_params` `1`,
#'   `params_bounds` the domain \eqn{(0, \infty)}, and `link_params` the link
#'   given here.
#'
#' @seealso
#' [gamma2_distrib()] for the two-parameter family this sits inside at
#' \eqn{\sigma^2 = 2\mu}, and [gamma1_distrib()] for its dispersion form;
#' [exponential_distrib()] for the case \eqn{\mu = 2};
#' [gengamma1_distrib()] for a wider family still; [fit_distrib()] to estimate
#' the parameter; [check_distrib()] to validate a family of your own against
#' the same battery this one passes; [ChisqDistrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 18.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dchisq pchisq qchisq rchisq
#'
#' @examples
#' d <- chisq_distrib()
#' d
#'
#' # The density is stats::dchisq at df = mu.
#' y <- c(1, 4, 9)
#' th <- list(mu = 4)
#' all.equal(distrib_pdf(d, y, th), dchisq(y, df = 4))
#'
#' # Moments: the variance is 2 mu, so the spread is not free.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(2 * sqrt(2 / 4), 12 / 4)
#'
#' # It is the gamma at sigma2 = 2 mu, and the exponential at mu = 2.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(gamma2_distrib(), y, list(mu = 4, sigma2 = 8)))
#' all.equal(distrib_pdf(d, y, list(mu = 2)),
#'           distrib_pdf(exponential_distrib(), y, list(mu = 2)))
#'
#' # From the second order on nothing involves the response, so the observed
#' # and expected information coincide on the parameter scale.
#' identical(distrib_hessian(d, y, th), distrib_expected_hessian(d, y, th))
#'
#' # Fitting recovers the degrees of freedom; the sample mean starts it off.
#' set.seed(7)
#' z <- distrib_rng(d, 2000, th)
#' rbind(fitted = coef(fit_distrib(d, z)), moment = c(mu = mean(z)))
#'
#' @export
chisq_distrib <- function(link_mu = log_link()) {
  ChisqDistrib(
    distrib_name = "chisq", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = link_mu)
  )
}
