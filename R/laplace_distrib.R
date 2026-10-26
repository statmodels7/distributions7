#' @include distrib.R generics.R
NULL

#' @title Laplace Distribution Class, Location and Scale
#' @name LaplaceDistrib
#'
#' @description
#' This is the S7 class of the Laplace (double exponential) family with location
#' \eqn{\mu} and scale \eqn{\sigma > 0}, with density \eqn{f(y) =
#' (2\sigma)^{-1}\exp(-|y-\mu|/\sigma)} on the whole real line. It inherits from
#' `continuous_distrib`, so every generic of the `distrib` class applies to it;
#' the methods listed below are registered on it directly and everything else
#' comes from the parent.
#'
#' The density has a kink at \eqn{y = \mu}, where \eqn{|y - \mu|} is not
#' differentiable. This makes the family non-regular in its location: the
#' observed second derivative in \eqn{\mu} is 0 almost everywhere while the
#' information is \eqn{1/\sigma^2}. The class records this by setting
#' `params_smooth` to `c(mu = FALSE, sigma = TRUE)`, a property that
#' [check_distrib()] and the finite-difference guards read.
#'
#' Build one with [laplace_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties; the validator of
#' `distrib` checks them, and the constructor adds no check of its own.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `LaplaceDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [laplace_distrib()] they hold `"laplace"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma")`, the interpretations
#'   `c(mu = "location", sigma = "scale")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, the two links, and
#'   `c(mu = FALSE, sigma = TRUE)` for `params_smooth`.
#'
#' @seealso [laplace_distrib()] to build one;
#'   [laplace2_distrib()] for the same law written by its rate
#'   \eqn{\lambda = 1/\sigma};
#'   [distrib_expected_hessian.LaplaceDistrib()] for the information at the
#'   kink; [pseudohuber_distrib()] for a smooth relative.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.LaplaceDistrib],
#'   [`distrib_cdf()`][distrib_cdf.LaplaceDistrib],
#'   [`distrib_quantile()`][distrib_quantile.LaplaceDistrib] and
#'   [`distrib_rng()`][distrib_rng.LaplaceDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.LaplaceDistrib],
#'   [`distrib_hessian()`][distrib_hessian.LaplaceDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.LaplaceDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.LaplaceDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.LaplaceDistrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.LaplaceDistrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.LaplaceDistrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.LaplaceDistrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.LaplaceDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.LaplaceDistrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.LaplaceDistrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.LaplaceDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.LaplaceDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.LaplaceDistrib],
#'   [`distrib_cross3_y()`][distrib_cross3_y.LaplaceDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.LaplaceDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.LaplaceDistrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.LaplaceDistrib],
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.LaplaceDistrib],
#'   [`distrib_deriv3_cdf()`][distrib_deriv3_cdf.LaplaceDistrib] and
#'   [`distrib_deriv4_cdf()`][distrib_deriv4_cdf.LaplaceDistrib];
#' - the moments: [`mean()`][mean.LaplaceDistrib],
#'   [`variance()`][variance.LaplaceDistrib],
#'   [`skewness()`][skewness.LaplaceDistrib] and
#'   [`kurtosis()`][kurtosis.LaplaceDistrib];
#' - the decomposition at the kink:
#'   [`kink_decomposition()`][kink_decomposition.LaplaceDistrib].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The class declares that the log-likelihood has a kink in mu.
#' d@params_smooth
#'
#' # sigma is a scale, not a standard deviation: the variance is 2 sigma^2.
#' th <- list(mu = 0.4, sigma = 1.5)
#' c(variance = variance(d, th), two_sigma_sq = 2 * 1.5^2)
#'
#' # The tails are heavier than a Gaussian's: the excess kurtosis is 3.
#' kurtosis(d, th)
LaplaceDistrib <- S7::new_class("LaplaceDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Laplace Probability Density Function
#' @name distrib_pdf.LaplaceDistrib
#' @description
#' Computes the Laplace density
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{2\sigma} \exp\left(-\dfrac{|y - \mu|}{\sigma}\right)}
#' from the log-density \eqn{-\log 2 - \log\sigma - |y-\mu|/\sigma}, which is
#' formed first and exponentiated only when `log = FALSE`. The two logarithms
#' are taken separately, so the value stays finite where \eqn{2\sigma}
#' overflows. The density is continuous
#' everywhere and has a corner at \eqn{y = \mu}, where the two exponential arms
#' meet. That corner makes the family non-regular in its location.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so every value is accepted.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive; the generic signals an
#'   error for a zero or negative value before the method is called.
#' @param log Logical of length 1. When `TRUE` the log-density is returned,
#'   which is exact and finite at every real `y`. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.LaplaceDistrib()] for the distribution function,
#'   [distrib_grad_y.LaplaceDistrib()] for the kink at \eqn{y = \mu},
#'   [laplace2_distrib()] for the rate parametrization, and [distrib_pdf()]
#'   for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The density is compared with its formula.
#' all.equal(distrib_pdf(d, y, th), exp(-abs(y - 0.4) / 1.5) / (2 * 1.5))
#'
#' # The log-density is exactly linear in |y - mu|, so it is a pair of
#' # straight lines meeting at mu.
#' round(distrib_pdf(d, 0.4 + c(-3, -1.5, 0, 1.5, 3), th, log = TRUE), 6)
#'
#' # The density is symmetric about mu and has its maximum there.
#' distrib_pdf(d, 0.4 + c(-2, 2), th)
#' c(at_mu = distrib_pdf(d, 0.4, th), one_over_2sigma = 1 / (2 * 1.5))
S7::method(distrib_pdf, LaplaceDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  b <- theta[[2]]
  log_d <- -log(2) - log(b) - abs(y - mu) / b
  if (log) log_d else exp(log_d)
}

#' @title Laplace Cumulative Distribution Function
#' @name distrib_cdf.LaplaceDistrib
#' @description
#' Computes the Laplace distribution function, which is one exponential on each
#' side of the location:
#' \deqn{F(q; \mu, \sigma) = \begin{cases}
#'   \tfrac{1}{2}\exp\!\left(\dfrac{q-\mu}{\sigma}\right), & q < \mu,\\[4pt]
#'   1 - \tfrac{1}{2}\exp\!\left(-\dfrac{q-\mu}{\sigma}\right), & q \ge \mu.
#' \end{cases}}
#' Both arms are continuous at \eqn{q = \mu}, where the value is \eqn{1/2}.
#' The method computes the logarithm of each tail at the standardized
#' residual \eqn{r = (q-\mu)/\sigma}: the tail on the far side of the
#' location is \eqn{\log(1/2) - |r|}, and the tail on the near side is
#' \eqn{\log(1 - e^{-|r|}/2)}, computed with `log1p()`. The requested tail is
#' returned on the log scale or exponentiated, so the method does not subtract
#' from one when `lower.tail = FALSE` or `log.p = TRUE`, and a tail
#' probability far below machine epsilon, or its logarithm, is returned
#' accurately. The helper is
#' [laplace_cdf_tails()].
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are
#'   \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.LaplaceDistrib()] for the inverse,
#'   [distrib_pdf.LaplaceDistrib()] for the density, and [distrib_cdf()] for
#'   the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The value is one half at the location, where the two arms meet.
#' distrib_cdf(d, 0.4, th)
#'
#' # Each arm is an exponential, written out.
#' q <- c(-2, 3)
#' all.equal(distrib_cdf(d, q, th),
#'           ifelse(q < 0.4, 0.5 * exp((q - 0.4) / 1.5),
#'                  1 - 0.5 * exp(-(q - 0.4) / 1.5)))
#'
#' # The distribution is symmetric: F(mu - a) + F(mu + a) = 1.
#' distrib_cdf(d, 0.4 - 2, th) + distrib_cdf(d, 0.4 + 2, th)
#'
#' # Fifty scale units out, the upper tail is exp(-50)/2, and its logarithm
#' # is exact as well.
#' distrib_cdf(d, 0.4 + 50 * 1.5, th, lower.tail = FALSE)
#' distrib_cdf(d, 0.4 + 50 * 1.5, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, LaplaceDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  mu <- theta[[1]]
  b <- theta[[2]]
  laplace_cdf_tails((q - mu) / b, lower.tail, log.p)
}

#' @title The Laplace Tails on the Log Scale
#' @name laplace_cdf_tails
#' @description
#' These are internal helpers of the two Laplace families. `laplace_cdf_tails()`
#' returns the probability that `lower.tail` and `log.p` select at the
#' standardized residual `r`, and computes each tail from its own arm:
#' \eqn{\log(1/2) - |r|} for the far tail and `log1p(-exp(-|r|)/2)` for the
#' near one.
#' `laplace_std_quantile()` returns the standard Laplace quantile from the
#' logarithm of the smaller tail, \eqn{\pm(\log P + \log 2)}.
#' @param r A numeric vector of standardized residuals.
#' @param p A numeric vector of probabilities, or their logarithms.
#' @param lower.tail,log.p Logical of length 1.
#' @return A numeric vector.
#' @examples
#' distributions7:::laplace_cdf_tails(c(-50, 0, 50), TRUE, TRUE)
#' distributions7:::laplace_std_quantile(1e-20, TRUE, FALSE)
#' @keywords internal
laplace_cdf_tails <- function(r, lower.tail, log.p) {
  far <- -log(2) - abs(r)
  near <- log1p(-0.5 * exp(-abs(r)))
  cdf_from_log_tails(ifelse(r < 0, far, near), ifelse(r < 0, near, far),
                     lower.tail, log.p)
}

#' @rdname laplace_cdf_tails
laplace_std_quantile <- function(p, lower.tail, log.p) {
  lt <- quantile_log_tails(p, lower.tail, log.p)
  # by indexing, so that a NaN tail gives NaN and not NA
  out <- -(lt$upper + log(2))
  lo <- !is.na(lt$lower) & lt$lower <= lt$upper
  out[lo] <- lt$lower[lo] + log(2)
  out
}

#' @title Laplace Quantile Function
#' @name distrib_quantile.LaplaceDistrib
#' @description
#' Computes the Laplace quantile function
#' \deqn{Q(p; \mu, \sigma) = \mu - \sigma\,\mathrm{sign}(p - \tfrac{1}{2})\,\log\left(1 - 2\left|p - \tfrac{1}{2}\right|\right).}
#' The median is \eqn{\mu}, and the two quartiles are
#' \eqn{\mu \pm \sigma \log 2}. The method works from the logarithm of the
#' smaller tail, \eqn{\log P}: the quantile is \eqn{\mu + \sigma(\log P + \log 2)}
#' when \eqn{P} is the lower tail and \eqn{\mu - \sigma(\log P + \log 2)} when
#' it is the upper one. The other tail's logarithm is taken as
#' \eqn{\log(1 - e^{\log P})} with `expm1()` or `log1p()`, so a probability
#' far below machine epsilon, given directly, through `lower.tail = FALSE` or
#' through `log.p = TRUE`, gives an accurate quantile. The helper is
#' [laplace_std_quantile()].
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. The endpoints give `-Inf` and `Inf`; a
#'   value outside the range gives `NaN` with a warning.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of length
#'   `max(length(p), length(mu), length(sigma))`.
#'
#' @seealso [distrib_cdf.LaplaceDistrib()], which this inverts;
#'   [distrib_rng.LaplaceDistrib()], which draws through it; and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The median is mu and the quartiles are mu -/+ sigma log 2.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#' 0.4 + c(-1, 0, 1) * 1.5 * log(2)
#'
#' # The quantile is an exact inverse: the round trip returns the
#' # probabilities that it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # A far tail given as a logarithm is compared with mu + sigma (log p + log 2).
#' distrib_quantile(d, -800, th, log.p = TRUE)
#' 0.4 + 1.5 * (-800 + log(2))
S7::method(distrib_quantile, LaplaceDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  mu <- theta[[1]]
  b <- theta[[2]]
  mu + b * laplace_std_quantile(p, lower.tail, log.p)
}

#' @title Laplace Random Number Generator
#' @name distrib_rng.LaplaceDistrib
#' @description
#' Draws `n` independent Laplace variates by inverse transform: `n` uniform
#' variates from [stats::runif()] are passed through
#' [distrib_quantile.LaplaceDistrib()]. The base class draws the same way for
#' a family with its own quantile method, so the method gives the draws that
#' the base class would give.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param n A single positive integer, the number of draws. One uniform variate
#'   is consumed per draw.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled.
#'   `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.LaplaceDistrib()], through which the draws are
#'   made; [fit_distrib()] to estimate the parameters back; and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#'
#' # Inverse transform, so the draws are the quantiles of the uniforms drawn.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.5))
#' set.seed(2)
#' identical(a, distrib_quantile(d, runif(3), list(mu = 0.4, sigma = 1.5)))
#'
#' # The sample variance recovers 2 sigma^2, not sigma^2.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
#' c(mean = mean(z), var = var(z), two_sigma_sq = 2 * 4)
S7::method(distrib_rng, LaplaceDistrib) <- function(distrib, n, theta, ...) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Laplace Score
#' @name distrib_gradient.LaplaceDistrib
#' @description
#' Computes the first derivatives of the Laplace log-density with respect to
#' \eqn{\mu} and \eqn{\sigma}, one value per observation, in closed form.
#' Writing \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\mathrm{sign}(r)}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}\left(\dfrac{|r|}{\sigma} - 1\right).}
#'
#' The score in \eqn{\mu} takes only the three values \eqn{-1/\sigma}, 0 and
#' \eqn{+1/\sigma}: it carries the sign of the residual but not its
#' size, which is why the maximum likelihood estimate of \eqn{\mu} is the
#' sample median. At \eqn{r = 0} exactly, `sign(0)` is 0 and the method returns
#' 0, the midpoint of the subdifferential \eqn{[-1/\sigma, 1/\sigma]}; the
#' derivative does not exist there.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count is. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each of
#'   length `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\sigma > 0} the scale, with variance \eqn{2\sigma^2}. \eqn{r = y - \mu}
#' is the residual.
#'
#' @seealso [distrib_hessian.LaplaceDistrib()] for the second derivatives,
#'   which vanish in \eqn{\mu};
#'   [distrib_expected_hessian.LaplaceDistrib()] for the information, which
#'   does not; and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms are written out and compared.
#' r <- y - 0.4
#' all.equal(g$mu, sign(r) / 1.5)
#' all.equal(g$sigma, (abs(r) / 1.5 - 1) / 1.5)
#'
#' # The score in mu is a sign: three values, whatever the residual.
#' unique(distrib_gradient(d, 0.4 + c(-100, -1, 0, 1, 100), th)$mu)
#'
#' # It vanishes summed at the sample median, which is the estimate of mu.
#' set.seed(12)
#' z <- distrib_rng(d, 401, list(mu = 3, sigma = 2))
#' sum(distrib_gradient(d, z, list(mu = median(z), sigma = 2))$mu)
S7::method(distrib_gradient, LaplaceDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  laplace_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Laplace Observed Hessian
#' @name distrib_hessian.LaplaceDistrib
#' @description
#' Computes the three distinct second derivatives of the Laplace log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, one value per observation, in
#' closed form. Writing \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = 0,
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma - 2|r|}{\sigma^3},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma} = -\dfrac{\mathrm{sign}(r)}{\sigma^2}.}
#'
#' The first entry is exactly zero for every observation, and that is the
#' correct value: the log-density is piecewise linear in \eqn{\mu}, so away
#' from \eqn{r = 0} its second derivative vanishes, and at \eqn{r = 0} the
#' derivative does not exist. The curvature of the log-likelihood in the
#' location is concentrated at the single point \eqn{r = 0}, so the
#' expectation of this second derivative is 0 although the information is
#' positive.
#'
#' Hence the observed information here is not the
#' information. [distrib_expected_hessian.LaplaceDistrib()] returns
#' \eqn{-1/\sigma^2}, obtained from the variance of the score, and that page
#' explains which identity holds and which fails.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count is. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
#'   `mu_mu` is a vector of zeros.
#'
#' @seealso [distrib_expected_hessian.LaplaceDistrib()], which returns the
#'   information and is not the expectation of this;
#'   [distrib_gradient.LaplaceDistrib()] for the score;
#'   [fisher_scoring()], the estimation route that this family needs; and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' h <- distrib_hessian(d, y, th)
#'
#' # The curvature in the location is exactly zero at every observation.
#' h$mu_mu
#'
#' # The other two components are compared with their formulas.
#' r <- y - 0.4
#' all.equal(h$sigma_sigma, (1.5 - 2 * abs(r)) / 1.5^3)
#' all.equal(h$mu_sigma, -sign(r) / 1.5^2)
#'
#' # Averaging the observed curvature in mu gives 0, and the information is
#' # 1/sigma^2. The two disagree because this family is not regular in mu.
#' set.seed(12)
#' z <- distrib_rng(d, 1e5, th)
#' c(observed_mean = mean(distrib_hessian(d, z, th)$mu_mu),
#'   expected = distrib_expected_hessian(d, 0, th)$mu_mu)
S7::method(distrib_hessian, LaplaceDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  laplace_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Laplace Expected Hessian
#' @name distrib_expected_hessian.LaplaceDistrib
#' @description
#' Returns the negative of the Fisher information, in closed form and with no
#' quadrature or simulation:
#' \deqn{-I(\mu) = -\dfrac{1}{\sigma^2}, \qquad
#'       -I(\sigma) = -\dfrac{1}{\sigma^2}, \qquad
#'       -I(\mu, \sigma) = 0.}
#'
#' For \eqn{\mu} this is not the expectation of
#' [distrib_hessian.LaplaceDistrib()], which is identically zero. It is the
#' variance of the score: \eqn{\partial\ell/\partial\mu = \mathrm{sign}(r)/\sigma}
#' has mean 0 and square \eqn{1/\sigma^2} almost surely, so
#' \eqn{\mathrm{Var} = 1/\sigma^2}. The mixed entry vanishes because the family
#' is symmetric about \eqn{\mu}, so \eqn{\mathbb{E}[\mathrm{sign}(r)] = 0}.
#'
#' Because the values do not depend on the data, `y` is read only for its
#' length.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL`: the generic applies an object before dispatching on
#'   the family, so this method receives `NULL` and returns the closed form.
#'   See [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count is. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length `length(y)` and constant within itself when the
#'   parameters are.
#'
#' @section Notation:
#' The observed information is \eqn{-\partial^2\ell/\partial\theta\,\partial\theta^\top}
#' at the data. The expected information is usually its expectation, and
#' for a regular family the two agree with the variance of the score by the
#' second Bartlett identity. The Laplace is not regular in \eqn{\mu}: the
#' first identity \eqn{\mathbb{E}[\partial\ell/\partial\mu] = 0} still holds,
#' the second does not, and the information is defined as the variance of
#' the score, which this method returns.
#'
#' @seealso [distrib_hessian.LaplaceDistrib()], whose `mu_mu` is zero;
#' [distrib_gradient.LaplaceDistrib()], whose variance this is;
#' [fisher_scoring()], which inverts this matrix and so fits this family where
#' a Newton step could not; and [distrib_expected_hessian()] for the generic
#' and the methods that it offers for a family with no closed form.
#'
#' @examples
#' d <- laplace_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # Both diagonal entries are -1/sigma^2; the mixed entry is 0.
#' lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#' -1 / 1.5^2
#'
#' # It is the variance of the score, which the sample confirms; the mean of
#' # the observed second derivative is 0 and does not.
#' set.seed(12)
#' z <- distrib_rng(d, 1e5, th)
#' s <- distrib_gradient(d, z, th)$mu
#' c(var_of_score = mean(s^2),
#'   information = -distrib_expected_hessian(d, 0, th)$mu_mu,
#'   mean_observed = mean(distrib_hessian(d, z, th)$mu_mu))
#'
#' # Fisher scoring fits the family because this matrix is nonsingular; the
#' # observed Hessian has a zero mu_mu entry.
#' coef(fit_distrib(d, z))
S7::method(distrib_expected_hessian, LaplaceDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), method = NULL, ..., threads = 1L) {
  laplace_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Laplace Third-Order Derivatives
#' @name distrib_deriv3.LaplaceDistrib
#' @description
#' Computes the four distinct third derivatives of the Laplace log-density with
#' respect to \eqn{\mu} and \eqn{\sigma}, in closed form. The log-density is
#' \eqn{-\log\sigma - \log 2 - |r|/\sigma} with \eqn{r = y - \mu}, so it is
#' linear in \eqn{\mu} on each side of the kink and every component that
#' differentiates twice or more in \eqn{\mu} is zero. The remaining components
#' come from \eqn{-\log\sigma} and from \eqn{|r|/\sigma}.
#'
#' With `expected = TRUE` the expectations under the model are returned, also
#' in closed form, using \eqn{\mathbb{E}|r| = \sigma} and
#' \eqn{\mathbb{E}[\mathrm{sign}(r)] = 0}. Both routes are closed form and the
#' method runs no quadrature.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL`: the generic applies an object before dispatching on
#'   the family, so this method receives `NULL` and returns the closed form.
#'   See [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does. These are the derivatives
#' that exist away from the kink at \eqn{y = \mu}.
#'
#' @seealso [distrib_hessian.LaplaceDistrib()] for the order below and
#'   [distrib_deriv4.LaplaceDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Linear in mu on each side of the kink, so every component that
#' # differentiates twice in mu is zero.
#' c(d3$mu_mu_mu[1], d3$mu_mu_sigma[1])
#'
#' # The expected values are also in closed form.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the observed component,
#' # away from the kink.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 + eps))$sigma_sigma
#' dn <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 - eps))$sigma_sigma
#' all.equal((up - dn) / (2 * eps), d3$sigma_sigma_sigma, tolerance = 1e-6)
S7::method(distrib_deriv3, LaplaceDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), method = NULL, ...) {
  if (expected) laplace_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  else laplace_deriv3_cpp(y, theta[[1]], theta[[2]])
}

#' @title Laplace Fourth-Order Derivatives
#' @name distrib_deriv4.LaplaceDistrib
#' @description
#' Computes the five distinct fourth derivatives of the Laplace log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, in closed form. As at third
#' order, the log-density is linear in \eqn{\mu} on each side of the kink, so
#' every component differentiating twice or more in \eqn{\mu} is zero and the
#' remaining components come from the \eqn{-\log\sigma} term and from
#' \eqn{|r|/\sigma}.
#'
#' With `expected = TRUE` the expectations under the model are returned, also
#' in closed form. Both routes are closed form.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL`: the generic applies an object before dispatching on
#'   the family, so this method receives `NULL` and returns the closed form.
#'   See [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_sigma`, `mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
#'   `sigma_sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k l)}} is the fourth derivative of the log-density with
#' respect to parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.LaplaceDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' names(distrib_deriv4(d, y, th))
#'
#' # The pure-sigma component is 6/sigma^4 from the -log(sigma) term, less the
#' # contribution of |r|/sigma.
#' distrib_deriv4(d, y, th)$sigma_sigma_sigma_sigma
#'
#' # The expected values are also in closed form.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, y, list(mu = 0.4, sigma = 1.5 + eps))$sigma_sigma_sigma
#' dn <- distrib_deriv3(d, y, list(mu = 0.4, sigma = 1.5 - eps))$sigma_sigma_sigma
#' all.equal((up - dn) / (2 * eps),
#'           distrib_deriv4(d, y, th)$sigma_sigma_sigma_sigma, tolerance = 1e-5)
S7::method(distrib_deriv4, LaplaceDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), method = NULL, ...) {
  if (expected) laplace_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  else laplace_deriv4_cpp(y, theta[[1]], theta[[2]])
}

#' @title Laplace Fifth-Order Derivatives
#' @name distrib_deriv5.LaplaceDistrib
#' @description
#' Computes the six distinct fifth derivatives of the Laplace log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, in closed form from a
#' compiled kernel. Away from the kink at \eqn{y = \mu}, with
#' \eqn{r = y - \mu}, the four components with two or more derivatives in
#' \eqn{\mu} are zero and
#' \deqn{\ell^{(\mu\sigma\sigma\sigma\sigma)} = \dfrac{24\,\mathrm{sign}(r)}{\sigma^5},
#'       \qquad \ell^{(\sigma\sigma\sigma\sigma\sigma)} =
#'       \dfrac{120 |r|}{\sigma^6} - \dfrac{24}{\sigma^5}.}
#' The kernel forms the powers of \eqn{1/\sigma}, so that a small scale does
#' not underflow an intermediate power where the value is finite.
#' With `scale = "link"` the components are composed with the link
#' derivatives by the routine that every family with a fifth-order kernel
#' uses.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric vector
#'   of length 1 or of the length of `y`.
#' @param scale One of `"parameter"` (the default) or `"link"`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of six numeric vectors, `mu_mu_mu_mu_mu`,
#'   `mu_mu_mu_mu_sigma`, `mu_mu_mu_sigma_sigma`, `mu_mu_sigma_sigma_sigma`,
#'   `mu_sigma_sigma_sigma_sigma` and `sigma_sigma_sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @seealso [distrib_deriv4()] for the order below and [distrib_deriv5()]
#'   for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' d5 <- distrib_deriv5(d, y, list(mu = 0.4, sigma = 1.5))
#' all.equal(d5$sigma_sigma_sigma_sigma_sigma,
#'           120 * abs(y - 0.4) / 1.5^6 - 24 / 1.5^5)
S7::method(distrib_deriv5, LaplaceDistrib) <- function(distrib, y, theta,
    scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, laplace_deriv5_cpp(y, theta[[1]], theta[[2]], threads),
               match.arg(scale))
}

#' @title Laplace First Derivative in the Response
#' @name distrib_grad_y.LaplaceDistrib
#' @description
#' Computes the first derivative of the Laplace log-density with respect to the
#' response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{\mathrm{sign}(y - \mu)}{\sigma},}
#' in closed form. The Laplace is a location family in \eqn{\mu}, so this is
#' the negative of the score in \eqn{\mu}. It takes only three values, and at
#' \eqn{y = \mu} exactly the method returns 0, the midpoint of the
#' subdifferential; the derivative does not exist there, because the
#' log-density has a corner.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, taking the values
#'   \eqn{-1/\sigma}, 0 and \eqn{1/\sigma} only.
#'
#' @seealso [distrib_hess_y.LaplaceDistrib()] for the second derivative, which
#'   is zero; [distrib_gradient.LaplaceDistrib()] for the score in the
#'   parameters; and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' all.equal(distrib_grad_y(d, y, th), -sign(y - 0.4) / 1.5)
#'
#' # This is a location family, so the derivative is minus the score in the
#' # location.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # The derivative takes three values only, and it is 0 at the kink itself.
#' distrib_grad_y(d, 0.4 + c(-100, -1e-9, 0, 1e-9, 100), th)
S7::method(distrib_grad_y, LaplaceDistrib) <- function(distrib, y, theta, ...) {
  -sign(y - theta[[1]]) / theta[[2]]
}

#' @title Laplace Second Derivative in the Response
#' @name distrib_hess_y.LaplaceDistrib
#' @description
#' Returns zero for every observation. The Laplace log-density is
#' \eqn{-\log(2\sigma) - |y-\mu|/\sigma}, which is linear in \eqn{y} on
#' each side of the location, so its second derivative in the response vanishes
#' wherever it exists. At \eqn{y = \mu} the first derivative drops from
#' \eqn{1/\sigma} to \eqn{-1/\sigma} and the second derivative does not exist;
#' the returned zero is the value away from that single point.
#'
#' The same fact appears in the parameters as the zero `mu_mu` component of
#' [distrib_hessian.LaplaceDistrib()], because the family is a location
#' family.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma`. The method does not read either component.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length `length(y)`.
#'
#' @seealso [distrib_grad_y.LaplaceDistrib()] for the first derivative, which
#'   jumps at the location; [distrib_hessian.LaplaceDistrib()] for the same
#'   vanishing curvature in the parameters; [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- laplace_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' distrib_hess_y(d, y, th)
#'
#' # This is a location family, so the second derivative equals the
#' # curvature in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # The first derivative drops by 2/sigma across the location, which is the
#' # curvature that the zero above does not capture.
#' diff(distrib_grad_y(d, 0.4 + c(-1e-9, 1e-9), th))
#' -2 / 1.5
S7::method(distrib_hess_y, LaplaceDistrib) <- function(distrib, y, theta, ...) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Laplace Distribution, Location and Scale
#'
#' @description
#' Builds the distribution object for the Laplace (double exponential) family
#' with location \eqn{\mu} and scale \eqn{\sigma > 0}. The returned object
#' carries closed-form derivatives of the log-density to fifth order in the
#' parameters and to fourth order in the response, and closed-form
#' moments.
#'
#' The family is symmetric about \eqn{\mu} with variance \eqn{2\sigma^2} and
#' excess kurtosis 3, so it is heavier tailed than a Gaussian. Its maximum
#' likelihood estimates are the sample median and the mean absolute deviation
#' about it, and both are available in closed form.
#'
#' This family is not regular in \eqn{\mu}: the density has a corner at
#' \eqn{y = \mu} and the observed second derivative there is zero almost
#' everywhere, while the information is \eqn{1/\sigma^2}. The object records
#' this in `params_smooth`, and the details below describe the consequences.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the location
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], because the
#'   location already ranges over the whole line.
#' @param link_sigma A `link` object from `linkfunctions7` for the scale
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{2\sigma}\exp\left(-\dfrac{|y-\mu|}{\sigma}\right),}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\sigma \in (0, \infty)}. The
#' distribution function is one exponential on each side of \eqn{\mu} and the
#' quantile function inverts it in closed form.
#'
#' The mean and the median are \eqn{\mu}, the variance is \eqn{2\sigma^2}, the
#' skewness is 0 and the excess kurtosis is 3. [laplace2_distrib()] is the same
#' law written by its rate \eqn{\lambda = 1/\sigma}, which is the form that a
#' lasso penalty uses.
#'
#' # The kink, and what it costs
#'
#' \eqn{|y - \mu|} is not differentiable at \eqn{y = \mu}, so the log-density
#' is piecewise linear in \eqn{\mu} with a corner. Three consequences follow,
#' and all of them appear in the methods:
#'
#' - the score in \eqn{\mu} is \eqn{\mathrm{sign}(r)/\sigma}, carrying only
#'   the sign of the residual;
#' - the observed second derivative in \eqn{\mu} is 0 wherever it exists, so
#'   [distrib_hessian.LaplaceDistrib()] returns a vector of zeros;
#' - the second Bartlett identity fails, and the information is defined as
#'   the variance of the score, \eqn{1/\sigma^2}, which
#'   [distrib_expected_hessian.LaplaceDistrib()] returns.
#'
#' `params_smooth` is `c(mu = FALSE, sigma = TRUE)`. When an entry is `FALSE`,
#' [check_distrib()] recomputes each finite-difference reference with half
#' the step and drops the observations at which the two disagree, because a
#' central difference that straddles the corner is not a derivative there;
#' the comparison still covers every parameter, including \eqn{\mu}.
#'
#' # Estimation
#'
#' Both estimates are closed form: \eqn{\hat\mu} is the sample median and
#' \eqn{\hat\sigma} the mean absolute deviation about it,
#' \eqn{n^{-1}\sum|y_i - \hat\mu|}. [fit_distrib()] reaches them by Fisher
#' scoring, which inverts the information above; a Newton step could not,
#' because the observed Hessian is singular in \eqn{\mu}. The example below checks
#' the fitted values against both closed forms.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\sigma > 0} the scale, with variance \eqn{2\sigma^2}. \eqn{r = y - \mu}
#' is the residual. The observed information is
#' \eqn{-\partial^2\ell/\partial\theta\,\partial\theta^\top} at the data; the
#' expected information is the variance of the score, which for a regular
#' family equals the expectation of the observed one and here does not.
#'
#' @return An S7 object of class `LaplaceDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"laplace"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, `link_params` the two links given here, and
#'   `params_smooth` `c(mu = FALSE, sigma = TRUE)`.
#'
#' @seealso
#' [laplace2_distrib()] for the rate parametrization;
#' [gaussian1_distrib()] for the light-tailed comparison;
#' [cauchy_distrib()] and [pseudohuber_distrib()] for the other robust location
#' families, the second of which smooths this one's corner away;
#' [enet_distrib()], whose density tends to that of [laplace2_distrib()] as
#' its mixing weight \eqn{\alpha} tends to 1;
#' [fit_distrib()] and [fisher_scoring()] to estimate the parameters;
#' [LaplaceDistrib] for the class.
#'
#' @references
#' Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001).
#' *The Laplace Distribution and Generalizations*. Birkhauser, Boston.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats runif
#'
#' @examples
#' d <- laplace_distrib()
#' d
#'
#' # The density is compared with its formula.
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' all.equal(distrib_pdf(d, y, th), exp(-abs(y - 0.4) / 1.5) / (2 * 1.5))
#'
#' # sigma is a scale: the variance is 2 sigma^2 and the excess kurtosis is 3.
#' c(variance = variance(d, th), kurtosis = kurtosis(d, th))
#'
#' # Both estimates are closed form: the median and the mean absolute
#' # deviation about it.
#' set.seed(12)
#' z <- distrib_rng(d, 2000, list(mu = 3, sigma = 2))
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = median(z), sigma = mean(abs(z - median(z)))))
#'
#' # The family is not regular in mu: the observed curvature there is 0 and
#' # the information is 1/sigma^2.
#' c(observed = mean(distrib_hessian(d, z, th)$mu_mu),
#'   expected = distrib_expected_hessian(d, 0, th)$mu_mu)
#'
#' @export
laplace_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  LaplaceDistrib(
    distrib_name = "laplace",
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
    ),

    params_smooth = c(mu = FALSE, sigma = TRUE)
  )
}

#' @name kink_decomposition.LaplaceDistrib
#' @title The Laplace's Kink
#'
#' @description
#' The log-density is \eqn{-\log(2\sigma) - \lvert y - \mu \rvert / \sigma}, so
#' the composition is the absolute value at \eqn{v = y - \mu} with
#' \eqn{c(\theta) = -1/\sigma}.
#'
#' @details
#' Only the location moves the kink: \eqn{\partial v/\partial \mu = -1} and
#' \eqn{\partial v/\partial \sigma = 0}, so [params_order()] reports
#' \eqn{m_\mu = 0} and \eqn{m_\sigma = \infty}, which agrees with
#' the value `c(mu = FALSE, sigma = TRUE)` of `params_smooth`.
#'
#' The scale is smooth although it multiplies the absolute value: the
#' derivative in \eqn{\sigma} of \eqn{c(\theta)\phi(v)} is
#' \eqn{c'(\theta)\phi(v)}, which exists at every \eqn{y}. A derivative fails
#' only where the argument \eqn{v} crosses the origin.
#'
#' @param distrib A `LaplaceDistrib` object.
#'
#' @return A [kink_spec()] with `phi = "abs"`.
#'
#' @examples
#' kink_decomposition(laplace_distrib())
#' params_order(laplace_distrib())
#' check_kink(laplace_distrib())
#'
#' @seealso [kink_decomposition()] for the generic, [laplace_distrib()] for the
#'   family.
#' @keywords internal
S7::method(kink_decomposition, LaplaceDistrib) <- function(distrib) {
  kink_spec(
    phi  = "abs",
    v    = function(y, theta) y - theta$mu,
    dv   = function(y, theta) list(mu = -1, sigma = 0),
    coef = function(theta) -1 / theta$sigma)
}
