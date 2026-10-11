#' @include distrib.R generics.R
NULL

#' @title Logistic Distribution Class
#' @name LogisticDistrib
#'
#' @description
#' The S7 class of the logistic family with mean \eqn{\mu} and scale
#' \eqn{\sigma > 0}, whose distribution function is the logistic sigmoid
#' \eqn{F(q) = [1 + e^{-(q-\mu)/\sigma}]^{-1}} on the whole real line. It
#' inherits from `continuous_distrib`, so every generic of the `distrib` class
#' applies to it; the methods listed below are registered on it directly and
#' everything else comes from the parent.
#'
#' The family is symmetric about \eqn{\mu}, which is therefore both the mean
#' and the median. Its variance is \eqn{\pi^2\sigma^2/3}, so \eqn{\sigma} is
#' **not** the standard deviation, and its excess kurtosis is \eqn{6/5},
#' slightly heavier than a Gaussian's.
#'
#' Build one with [logistic_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties; the validator of
#' `distrib` checks them, and the constructor adds no check of its own.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `LogisticDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [logistic_distrib()] they hold `"logistic"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma")`, the interpretations
#'   `c(mu = "mean", sigma = "scale")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [logistic_distrib()] to build one;
#'   [gaussian1_distrib()], which it resembles with slightly heavier tails;
#'   [distrib_expected_hessian.LogisticDistrib()] for the information.
#'
#' @section Methods:
#' Registered on this class, each linked to the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.LogisticDistrib],
#'   [`distrib_cdf()`][distrib_cdf.LogisticDistrib],
#'   [`distrib_quantile()`][distrib_quantile.LogisticDistrib] and
#'   [`distrib_rng()`][distrib_rng.LogisticDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.LogisticDistrib],
#'   [`distrib_hessian()`][distrib_hessian.LogisticDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.LogisticDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.LogisticDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.LogisticDistrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.LogisticDistrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.LogisticDistrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.LogisticDistrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.LogisticDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.LogisticDistrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.LogisticDistrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.LogisticDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.LogisticDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.LogisticDistrib],
#'   [`distrib_cross3_y()`][distrib_cross3_y.LogisticDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.LogisticDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.LogisticDistrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.LogisticDistrib],
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.LogisticDistrib],
#'   [`distrib_deriv3_cdf()`][distrib_deriv3_cdf.LogisticDistrib] and
#'   [`distrib_deriv4_cdf()`][distrib_deriv4_cdf.LogisticDistrib];
#' - the moments: [`mean()`][mean.LogisticDistrib],
#'   [`variance()`][variance.LogisticDistrib],
#'   [`skewness()`][skewness.LogisticDistrib] and
#'   [`kurtosis()`][kurtosis.LogisticDistrib].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params_interpretation
#'
#' # sigma is a scale, not a standard deviation: the variance is pi^2 sigma^2/3.
#' th <- list(mu = 0.4, sigma = 1.5)
#' c(variance = variance(d, th), pi_sq_over_3 = pi^2 * 1.5^2 / 3)
#'
#' # Symmetric, so the mean is the median; the excess kurtosis is 6/5.
#' c(mean = mean(d, th), median = distrib_quantile(d, 0.5, th),
#'   excess_kurtosis = kurtosis(d, th))
LogisticDistrib <- S7::new_class("LogisticDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Logistic Probability Density Function
#' @name distrib_pdf.LogisticDistrib
#' @description
#' Computes the logistic density
#' \deqn{f(y; \mu, \sigma) = \dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left[1 + \exp\left(-\dfrac{y-\mu}{\sigma}\right)\right]^2}}
#' by calling [stats::dlogis()] at `location = mu` and `scale = sigma`. The
#' density is symmetric about \eqn{\mu} and its tails decay exponentially,
#' like \eqn{e^{-|y-\mu|/\sigma}}, so they are heavier than a Gaussian's.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive; the generic signals an
#'   error for a zero or negative value before the method is called.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.LogisticDistrib()] for the distribution function,
#'   [distrib_gradient.LogisticDistrib()] for the derivatives of the
#'   log-density, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#'
#' # The method is stats::dlogis at this parametrization.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)),
#'           dlogis(y, location = 0.4, scale = 1.5))
#'
#' # Symmetric about mu: equal densities at equal distances either side.
#' distrib_pdf(d, 0.4 + c(-2, 2), list(mu = 0.4, sigma = 1.5))
#'
#' # The density and the distribution function are linked by
#' # f = F (1 - F) / sigma, the logistic sigmoid's own derivative.
#' F <- distrib_cdf(d, y, list(mu = 0.4, sigma = 1.5))
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)), F * (1 - F) / 1.5)
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
#' Computes the logistic distribution function
#' \deqn{F(q; \mu, \sigma) = \dfrac{1}{1 + \exp\left(-\dfrac{q-\mu}{\sigma}\right)}}
#' by calling [stats::plogis()]. This is the logistic sigmoid, the inverse of
#' the logit link, so the same curve appears here as a distribution function
#' and in `linkfunctions7` as an inverse link.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned, which stays finite far into either tail.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.LogisticDistrib()] for the inverse,
#'   [distrib_pdf.LogisticDistrib()] for the density,
#'   [linkfunctions7::logit_link()], whose inverse is this curve, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The method is stats::plogis at this parametrization.
#' all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
#'           plogis(c(-1.2, 0.3, 2.5), location = 0.4, scale = 1.5))
#'
#' # Symmetric: F(mu - a) + F(mu + a) = 1.
#' distrib_cdf(d, 0.4 - 2, th) + distrib_cdf(d, 0.4 + 2, th)
#'
#' # It is the inverse logit link, so the two agree exactly.
#' all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), list(mu = 0, sigma = 1)),
#'           linkfunctions7::linkinv(linkfunctions7::logit_link(),
#'                                   c(-1.2, 0.3, 2.5)))
S7::method(distrib_cdf, LogisticDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
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
#' Computes the logistic quantile function
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \log\left(\dfrac{p}{1-p}\right)}
#' by calling [stats::qlogis()]. The argument of the logarithm is the odds, so
#' this function is the logit link applied to `p` and then carried onto the
#' location and scale.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; the endpoints give `-Inf` and `Inf`.
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
#' @seealso [distrib_cdf.LogisticDistrib()], which this inverts;
#'   [linkfunctions7::logit_link()], which it applies to `p`; and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The median is mu, the family being symmetric.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # It is mu + sigma times the log odds.
#' all.equal(distrib_quantile(d, p, th), 0.4 + 1.5 * log(p / (1 - p)))
S7::method(distrib_quantile, LogisticDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
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
#' Draws `n` independent logistic variates by calling [stats::rlogis()], so the
#' draws come from R's own generator and depend on `.Random.seed` in the usual
#' way. The base class would draw by inversion through
#' [distrib_quantile.LogisticDistrib()], which gives a different sequence
#' for the same seed.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled.
#'   `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.LogisticDistrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#'
#' # Same generator as stats::rlogis, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.5))
#' set.seed(2)
#' identical(a, rlogis(3, location = 0.4, scale = 1.5))
#'
#' # The sample variance recovers pi^2 sigma^2 / 3, not sigma^2.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
#' c(mean = mean(z), var = var(z), pi_sq_over_3 = pi^2 * 4 / 3)
S7::method(distrib_rng, LogisticDistrib) <- function(distrib, n, theta, ...) {
  stats::rlogis(
    n = n,
    location = theta[[1]],
    scale = theta[[2]]
  )
}

#' @title Logistic Score
#' @name distrib_gradient.LogisticDistrib
#' @description
#' Computes the first derivatives of the logistic log-density with respect to
#' \eqn{\mu} and \eqn{\sigma}, one value per observation, in closed form.
#' Writing \eqn{z = (y - \mu)/\sigma},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma} \tanh\left(\dfrac{z}{2}\right),
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma}\left[1 - z \tanh\left(\dfrac{z}{2}\right)\right].}
#'
#' The score in \eqn{\mu} is bounded by \eqn{1/\sigma} and tends to
#' \eqn{\pm 1/\sigma} as \eqn{|z|} grows, so a distant observation contributes
#' a fixed amount. Under a Gaussian the contribution grows without bound, and
#' under a Cauchy it tends to zero.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each of
#'   length `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean,
#' \eqn{\sigma > 0} the scale and \eqn{z = (y-\mu)/\sigma} the standardized
#' residual. \eqn{\sigma} is not the standard deviation, which is
#' \eqn{\pi\sigma/\sqrt{3}}.
#'
#' @seealso [distrib_hessian.LogisticDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.LogisticDistrib()] for their expectation,
#'   [distrib_gradient.CauchyDistrib()] for a redescending score and
#'   [distrib_gradient.Gaussian1Distrib()] for an unbounded one, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' z <- (y - 0.4) / 1.5
#' all.equal(g$mu, tanh(z / 2) / 1.5)
#' all.equal(g$sigma, -(1 - z * tanh(z / 2)) / 1.5)
#'
#' # The score in mu saturates at 1/sigma instead of growing.
#' round(distrib_gradient(d, 0.4 + c(0, 1.5, 3, 6, 12, 60), th)$mu, 4)
#' 1 / 1.5
#'
#' # The summed score is close to zero at the maximum likelihood estimate,
#' # within the tolerance at which the fit stops.
#' set.seed(9)
#' zz <- distrib_rng(d, 3000, list(mu = 2, sigma = 1))
#' fit <- fit_distrib(d, zz)
#' round(vapply(distrib_gradient(d, zz, as.list(coef(fit))), sum, numeric(1)), 8)
S7::method(distrib_gradient, LogisticDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  logistic_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Logistic Observed Hessian
#' @name distrib_hessian.LogisticDistrib
#' @description
#' Computes the three distinct second derivatives of the logistic log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, one value per observation, in
#' closed form. Writing \eqn{z = (y - \mu)/\sigma},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{2\sigma^2} \mathrm{sech}^2\left(\dfrac{z}{2}\right),}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1}{\sigma^2}\left[1 - \dfrac{z^2}{2}\,\mathrm{sech}^2\left(\dfrac{z}{2}\right) - 2 z \tanh\left(\dfrac{z}{2}\right)\right],}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma} = -\dfrac{1}{\sigma^2}\left[\tanh\left(\dfrac{z}{2}\right) + \dfrac{z}{2}\,\mathrm{sech}^2\left(\dfrac{z}{2}\right)\right].}
#'
#' The curvature in \eqn{\mu} is negative everywhere and vanishes as
#' \eqn{|z|} grows, so a distant observation carries almost no information
#' about the location. The kernel computes
#' \eqn{\mathrm{sech}^2(z/2)} as \eqn{4e/(1+e)^2} with \eqn{e = e^{-|z|}},
#' so the value keeps its relative accuracy far into the tails, where
#' \eqn{1 - \tanh^2(z/2)} would cancel to zero. The expectations are
#' [distrib_expected_hessian.LogisticDistrib()].
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
#'   The three name the distinct entries of a symmetric \eqn{2 \times 2} matrix
#'   per observation.
#'
#' @seealso [distrib_gradient.LogisticDistrib()] for the score,
#'   [distrib_expected_hessian.LogisticDistrib()] for the expectation of this
#'   quantity, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' h <- distrib_hessian(d, y, th)
#'
#' # The closed form for the location, written out.
#' z <- (y - 0.4) / 1.5
#' all.equal(h$mu_mu, -(1 - tanh(z / 2)^2) / (2 * 1.5^2))
#'
#' # Concave in mu everywhere, and flattening as the residual grows.
#' round(distrib_hessian(d, 0.4 + c(0, 1.5, 3, 6, 12), th)$mu_mu, 6)
#'
#' # A central difference of the score reproduces it.
#' eps <- 1e-5
#' up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
S7::method(distrib_hessian, LogisticDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  logistic_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Logistic Expected Hessian
#' @name distrib_expected_hessian.LogisticDistrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{3\sigma^2},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{3 + \pi^2}{9\sigma^2},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}\right] = 0.}
#' The mixed entry vanishes because the family is symmetric about \eqn{\mu},
#' so the location and the scale are orthogonal and their estimates are
#' asymptotically independent. The information in the location is
#' \eqn{1/(3\sigma^2)}, a third of a Gaussian's \eqn{1/\sigma^2} at the same
#' \eqn{\sigma}. Compared at equal variance the two are much closer, the
#' logistic variance being \eqn{\pi^2\sigma^2/3}.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The logistic is
#' a regular family, so the second Bartlett identity holds and this equals the
#' variance of the score.
#'
#' @seealso [distrib_hessian.LogisticDistrib()] for the observed quantity this
#'   is the expectation of, [fisher_scoring()], which inverts it at each step,
#'   and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The three closed forms.
#' lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#' c(-1 / (3 * 1.5^2), -(3 + pi^2) / (9 * 1.5^2))
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(8)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_expected_hessian, LogisticDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ..., threads = 1L) {
  logistic_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Logistic Third-Order Derivatives
#' @name distrib_deriv3.LogisticDistrib
#' @description
#' Computes the four distinct third derivatives of the logistic log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}. The observed values are closed
#' form. Writing \eqn{z = (y-\mu)/\sigma}, \eqn{t = 1/(1+e^{-z})} and
#' \eqn{u = 1-t}, the log-density is \eqn{\ell = -\log\sigma + g(z)} with
#' \eqn{g(z) = -z - 2\log(1+e^{-z})}, whose derivatives in \eqn{z} are
#' \deqn{g_1 = 1-2t, \quad g_2 = -2tu, \quad g_3 = -2tu(1-2t), \quad g_4 = -2tu(1-6tu).}
#' Both parameters enter only through \eqn{z}, so each component is a
#' polynomial in \eqn{z} with the \eqn{g_j} as coefficients:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = -\dfrac{g_3}{\sigma^3}, \qquad
#'       \dfrac{\partial^3 \ell}{\partial \mu^2 \partial \sigma} = -\dfrac{2g_2 + z g_3}{\sigma^3},}
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu \partial \sigma^2} = -\dfrac{2g_1 + 4z g_2 + z^2 g_3}{\sigma^3}, \qquad
#'       \dfrac{\partial^3 \ell}{\partial \sigma^3} = -\dfrac{2 + 6z g_1 + 6z^2 g_2 + z^3 g_3}{\sigma^3}.}
#'
#' With `expected = TRUE` the method calls [expected_derivative()], so the
#' values are computed numerically according to `approx` (by default
#' `"integrate"`, one quadrature per component), and `approx` and `nsim`
#' are read. The expectations at this order are elementary
#' (\eqn{1/(2\sigma^3)} for `mu_mu_sigma`, \eqn{(2+\pi^2)/(2\sigma^3)} for
#' `sigma_sigma_sigma`, and zero for the two components odd in \eqn{z}), but
#' the method does not use these closed forms.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned, by the numerical route described above. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default), `"bartlett"`, `"mc"` or
#'   `"opg"`, matched by [base::match.arg()]. Read only when
#'   `expected = TRUE`, where it selects how the expectation is approximated.
#'   At this order `"opg"` is evaluated as `"bartlett"`, the outer product
#'   of scores being an identity of order two.
#' @param nsim A single positive integer, the sample size used when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use for the observed values. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.LogisticDistrib()] for the order below and
#'   [distrib_deriv4.LogisticDistrib()] for the order above;
#'   [expected_derivative()] for the numerical route the expected values take;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the observed component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#'
#' # The expected values are numerical here: the pure-mu component is zero by
#' # symmetry and comes back at the quadrature's own accuracy, not at zero.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
S7::method(distrib_deriv3, LogisticDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    logistic_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Logistic Fourth-Order Derivatives
#' @name distrib_deriv4.LogisticDistrib
#' @description
#' Computes the five distinct fourth derivatives of the logistic log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}. The observed values are closed
#' form, as polynomials in \eqn{z = (y-\mu)/\sigma} with the \eqn{g_j} of
#' [distrib_deriv3.LogisticDistrib()] as coefficients:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4} = \dfrac{g_4}{\sigma^4}, \qquad
#'       \dfrac{\partial^4 \ell}{\partial \mu^3 \partial \sigma} = \dfrac{3g_3 + z g_4}{\sigma^4}, \qquad
#'       \dfrac{\partial^4 \ell}{\partial \mu^2 \partial \sigma^2} = \dfrac{6g_2 + 6z g_3 + z^2 g_4}{\sigma^4},}
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu \partial \sigma^3} = \dfrac{6g_1 + 18z g_2 + 9z^2 g_3 + z^3 g_4}{\sigma^4}, \qquad
#'       \dfrac{\partial^4 \ell}{\partial \sigma^4} = \dfrac{6 + 24z g_1 + 36z^2 g_2 + 12z^3 g_3 + z^4 g_4}{\sigma^4}.}
#'
#' With `expected = TRUE` the method calls [expected_derivative()], as at
#' third order, so the values are computed numerically according to
#' `approx`, and `approx` and `nsim` are read. The expectations are
#' elementary here too (\eqn{1/(15\sigma^4)} for `mu_mu_mu_mu`,
#' \eqn{(\pi^2-60)/(45\sigma^4)} for `mu_mu_sigma_sigma`,
#' \eqn{(7\pi^4-600\pi^2-810)/(225\sigma^4)} for `sigma_sigma_sigma_sigma`,
#' and zero for the two components odd in \eqn{z}), but the method does not
#' use these closed forms.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned, by the numerical route. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default), `"bartlett"`, `"mc"` or
#'   `"opg"`, matched by [base::match.arg()]. Read only when `expected = TRUE`;
#'   `"opg"` is evaluated as `"bartlett"` at this order.
#' @param nsim A single positive integer, the sample size used when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use for the observed values. Defaults to `1L`.
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
#' @seealso [distrib_deriv3.LogisticDistrib()] for the order below and for the
#'   \eqn{g_j} used here; [expected_derivative()] for the numerical route;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' names(distrib_deriv4(d, y, th))
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, LogisticDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    logistic_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Logistic Fifth-Order Derivatives
#' @name distrib_deriv5.LogisticDistrib
#' @description
#' Computes the six distinct fifth derivatives of the logistic log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, in closed form from a
#' compiled kernel.
#' Writing the log-density as \eqn{\ell = -\log\sigma + g(z)} with
#' \eqn{z = (y - \mu)/\sigma}, each component is \eqn{\sigma^{-5}} times a
#' sum of terms \eqn{c\, z^p g^{(i)}(z)} with integer coefficients, plus
#' \eqn{-24\sigma^{-5}} in the component that is pure in \eqn{\sigma}. The
#' kernel evaluates the sum from the derivatives of
#' \eqn{g(z) = -z - 2\log(1 + e^{-z})}, written in the sigmoid
#' \eqn{t = 1/(1 + e^{-z})} and \eqn{u = 1 - t}, so that every factor
#' \eqn{g^{(i)}} lies in \eqn{[-1, 1]}.
#' With `scale = "link"` the components are composed with the link
#' derivatives by the routine that every family with a fifth-order kernel
#' uses.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
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
#' d <- logistic_distrib()
#' d5 <- distrib_deriv5(d, c(-1.2, 0.3, 2.5), list(mu = 0.4, sigma = 1.5))
#' names(d5)
S7::method(distrib_deriv5, LogisticDistrib) <- function(distrib, y, theta,
    scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, logistic_deriv5_cpp(y, theta[[1]], theta[[2]], threads),
               match.arg(scale))
}

#' @title Logistic First Derivative in the Response
#' @name distrib_grad_y.LogisticDistrib
#' @description
#' Computes the first derivative of the logistic log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{1}{\sigma}\tanh\left(\dfrac{z}{2}\right),
#'       \qquad z = \dfrac{y-\mu}{\sigma},}
#' in closed form. The logistic is a location family in \eqn{\mu}, so the
#' response enters the log-density only through \eqn{z} and this derivative is
#' the negative of the score in \eqn{\mu}. It is bounded by \eqn{1/\sigma}.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_hess_y.LogisticDistrib()] for the second derivative in the
#'   response, [distrib_gradient.LogisticDistrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' all.equal(distrib_grad_y(d, y, th), -tanh((y - 0.4) / (2 * 1.5)) / 1.5)
#'
#' # A location family: the derivative in the response is minus the score in
#' # the location.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Bounded by 1/sigma however far out the observation is.
#' max(abs(distrib_grad_y(d, seq(-1e3, 1e3, length.out = 1e4), th)))
#' 1 / 1.5
S7::method(distrib_grad_y, LogisticDistrib) <- function(distrib, y, theta, ...) {
  s <- theta[[2]]
  -tanh(0.5 * (y - theta[[1]]) / s) / s
}

#' @title Logistic Second Derivative in the Response
#' @name distrib_hess_y.LogisticDistrib
#' @description
#' Computes the second derivative of the logistic log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{1}{2\sigma^2}\,\mathrm{sech}^2\left(\dfrac{z}{2}\right),
#'       \qquad z = \dfrac{y-\mu}{\sigma},}
#' in closed form. It is negative everywhere, so the log-density is concave in
#' the response throughout, and it decays to zero in either tail. Being a
#' location family, the logistic has the same curvature in the response as in
#' its location, so this equals the `mu_mu` component of
#' [distrib_hessian.LogisticDistrib()].
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_grad_y.LogisticDistrib()] for the first derivative in the
#'   response, [distrib_hessian.LogisticDistrib()] for the curvature in the
#'   parameters, [distrib_hess_y.CauchyDistrib()] for a family that is convex
#'   in its tails, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- logistic_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same curvature as in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Concave everywhere, unlike the Cauchy, which turns convex in its tails.
#' all(distrib_hess_y(d, seq(-50, 50, length.out = 1e3), th) < 0)
#' any(distrib_hess_y(cauchy_distrib(), seq(-50, 50, length.out = 1e3), th) > 0)
S7::method(distrib_hess_y, LogisticDistrib) <- function(distrib, y, theta, ...) {
  s <- theta[[2]]
  # sech^2(z/2) = 4e/(1 + e)^2 with e = exp(-|z|), without 1 - tanh^2
  e <- exp(-abs(y - theta[[1]]) / s)
  -2 * e / ((1 + e)^2 * s^2)
}

# --- CONSTRUCTOR WRAPPER ---

#' Logistic Distribution
#'
#' @description
#' Builds the distribution object for the logistic family with mean \eqn{\mu}
#' and scale \eqn{\sigma > 0}. The returned object carries closed-form
#' derivatives of the log-density to fifth order in the parameters and to
#' fourth order in the response, and closed-form moments.
#'
#' The family is symmetric about \eqn{\mu} with variance
#' \eqn{\pi^2\sigma^2/3}, so \eqn{\sigma} is a scale, not a standard
#' deviation. Its distribution function is the logistic sigmoid, the same curve
#' `linkfunctions7` uses as the inverse logit link.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the mean ranging
#'   over the whole line already.
#' @param link_sigma A `link` object from `linkfunctions7` for the scale
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma) = \dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left[1 + \exp\left(-\dfrac{y-\mu}{\sigma}\right)\right]^2},}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\sigma \in (0, \infty)}. The
#' distribution function is \eqn{F(q) = [1 + e^{-(q-\mu)/\sigma}]^{-1}} and the
#' quantile function \eqn{Q(p) = \mu + \sigma \log(p/(1-p))}.
#'
#' The mean and the median are \eqn{\mu}, the variance is
#' \eqn{\pi^2 \sigma^2/3}, the skewness is 0 and the excess kurtosis is
#' \eqn{6/5}. A logistic and a Gaussian matched on their variance are close in
#' the body and differ in the tails: the logistic density decays like
#' \eqn{e^{-|z|}} and the Gaussian density like \eqn{e^{-z^2/2}}.
#'
#' # Derivatives
#'
#' With \eqn{z = (y-\mu)/\sigma} the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma}\tanh\left(\dfrac{z}{2}\right), \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma}\left[1 - z\tanh\left(\dfrac{z}{2}\right)\right],}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{3\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{3+\pi^2}{9\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0.}
#'
#' The observed derivatives of orders three to five are closed form, in
#' [distrib_deriv3.LogisticDistrib()], [distrib_deriv4.LogisticDistrib()] and
#' [distrib_deriv5.LogisticDistrib()]. The `expected = TRUE` branch of the
#' third and fourth orders calls [expected_derivative()] and computes the
#' expectations numerically, although they are elementary.
#'
#' # Estimation
#'
#' There is no closed-form estimate; [fit_distrib()] maximizes the
#' log-likelihood on the link scale. The log-likelihood is concave in
#' \eqn{\mu} for a fixed \eqn{\sigma}, so the location is well determined.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\sigma > 0} the scale, with standard deviation \eqn{\pi\sigma/\sqrt 3}.
#' \eqn{z = (y-\mu)/\sigma} is the standardized residual. \eqn{\eta} is a
#' parameter on the unconstrained scale of its link.
#'
#' @return An S7 object of class `LogisticDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"logistic"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [gaussian1_distrib()] for the light-tailed comparison;
#' [gumbel_distrib()], the asymmetric extreme-value relative;
#' [linkfunctions7::logit_link()], whose inverse is this distribution function;
#' [fit_distrib()] to estimate the parameters; [LogisticDistrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995).
#' *Continuous Univariate Distributions*, Volume 2, 2nd edition, Chapter 23.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dlogis plogis qlogis rlogis
#'
#' @examples
#' d <- logistic_distrib()
#' d
#'
#' # The density is R's own.
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' all.equal(distrib_pdf(d, y, th), dlogis(y, 0.4, 1.5))
#'
#' # sigma is a scale: the standard deviation is pi sigma / sqrt(3).
#' c(sd_from_variance = sqrt(variance(d, th)), pi * 1.5 / sqrt(3))
#'
#' # Fitting recovers the parameters.
#' set.seed(9)
#' z <- distrib_rng(d, 3000, list(mu = 2, sigma = 1))
#' coef(fit_distrib(d, z))
#'
#' # Matched on variance, a logistic sits close to a Gaussian in the body and
#' # puts more mass in the tails.
#' g <- gaussian1_distrib()
#' s <- pi * 1.5 / sqrt(3)
#' rbind(logistic = distrib_cdf(d, 0.4 + c(1, 2, 3, 4) * s, th,
#'                              lower.tail = FALSE),
#'       gaussian = distrib_cdf(g, 0.4 + c(1, 2, 3, 4) * s,
#'                              list(mu = 0.4, sigma = s), lower.tail = FALSE))
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
