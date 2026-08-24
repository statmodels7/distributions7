#' @include distrib.R generics.R
NULL

#' @title Cauchy Distribution Class
#' @name CauchyDistrib
#'
#' @description
#' The S7 class of the Cauchy family with location \eqn{\mu} and scale
#' \eqn{\sigma > 0}, the Student's t on one degree of freedom, with density
#' \eqn{f(y) = [\pi\sigma(1 + ((y-\mu)/\sigma)^2)]^{-1}} on the whole real
#' line. It inherits from `continuous_distrib`, so it answers every generic of
#' the `distrib` contract; the eleven methods listed below are registered on it
#' directly and everything else comes from the parent.
#'
#' No moment of this family exists, so `mean()`, `variance()`, `skewness()` and
#' `kurtosis()` return `NaN`. Its location and scale are the median and the
#' half-interquartile range, both available from [distrib_quantile()].
#'
#' Build one with [cauchy_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `CauchyDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [cauchy_distrib()] they hold `"cauchy"`, `"univariate"`,
#'   `c(-Inf, Inf)`, `c("mu", "sigma")`, the interpretations
#'   `c(mu = "location", sigma = "scale")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [cauchy_distrib()] to build one;
#'   [student_t1_distrib()], of which this is the case \eqn{\nu = 1};
#'   [skewness.CauchyDistrib()] for why the moments are `NaN`;
#'   [distrib_gradient.CauchyDistrib()] for the redescending score.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.CauchyDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.CauchyDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.CauchyDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.CauchyDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.CauchyDistrib],
#'   [`distrib_gradient()`][distrib_gradient.CauchyDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.CauchyDistrib],
#'   [`distrib_hessian()`][distrib_hessian.CauchyDistrib],
#'   [`distrib_pdf()`][distrib_pdf.CauchyDistrib],
#'   [`distrib_quantile()`][distrib_quantile.CauchyDistrib],
#'   [`distrib_rng()`][distrib_rng.CauchyDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- cauchy_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params
#' d@params_interpretation
#'
#' # The interpretations say "location" and "scale", not "mean" and "standard
#' # deviation": no moment of this family exists.
#' th <- list(mu = 0.4, sigma = 1.5)
#' c(mean(d, th), variance(d, th))
#'
#' # What does exist is the median and the half-interquartile range.
#' q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#' c(median = q[2], half_iqr = (q[3] - q[1]) / 2)
CauchyDistrib <- S7::new_class("CauchyDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Cauchy Probability Density Function
#' @name distrib_pdf.CauchyDistrib
#' @description
#' Computes the Cauchy density
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right]}}
#' by calling [stats::dcauchy()] at `location = mu` and `scale = sigma`. The
#' density decays like \eqn{y^{-2}}, so it stays representable far out where a
#' Gaussian of the same scale has underflowed to zero.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive; a zero or negative value
#'   gives `NaN` with a warning from [stats::dcauchy()].
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.CauchyDistrib()] for the distribution function,
#'   [distrib_gradient.CauchyDistrib()] for the derivatives of the log-density,
#'   [distrib_pdf.Gaussian1Distrib()] for the light-tailed comparison, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#'
#' # The method is stats::dcauchy at this parametrization.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)),
#'           dcauchy(y, location = 0.4, scale = 1.5))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0, 1, 2), sigma = c(1, 1.5, 2)))
#'
#' # Forty scale units out the Cauchy density is still an ordinary number,
#' # where the Gaussian of the same scale has underflowed.
#' distrib_pdf(d, 40, list(mu = 0, sigma = 1))
#' distrib_pdf(gaussian1_distrib(), 40, list(mu = 0, sigma = 1))
S7::method(distrib_pdf, CauchyDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dcauchy(
    x = y,
    location = theta[[1]],
    scale = theta[[2]],
    log = log
  )
}

#' @title Cauchy Cumulative Distribution Function
#' @name distrib_cdf.CauchyDistrib
#' @description
#' Computes the Cauchy distribution function
#' \deqn{F(q; \mu, \sigma) = \dfrac{1}{2} + \dfrac{1}{\pi}\arctan\left(\dfrac{q-\mu}{\sigma}\right)}
#' by calling [stats::pcauchy()]. Far out in either tail the probability
#' behaves like \eqn{1/(\pi |z|)} with \eqn{z = (q-\mu)/\sigma}, so it decays
#' at a polynomial rate and `lower.tail = FALSE` returns an ordinary number
#' where a light-tailed family returns zero.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.CauchyDistrib()] for the inverse,
#'   [distrib_pdf.CauchyDistrib()] for the density, and [distrib_cdf()] for
#'   the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The method is stats::pcauchy at this parametrization.
#' all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
#'           pcauchy(c(-1.2, 0.3, 2.5), location = 0.4, scale = 1.5))
#'
#' # A hundred scale units out the tail is 1/(100 pi), not zero.
#' distrib_cdf(d, 0.4 + 100 * 1.5, th, lower.tail = FALSE)
#' 1 / (100 * pi)
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
S7::method(distrib_cdf, CauchyDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pcauchy(
    q = q,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Cauchy Quantile Function
#' @name distrib_quantile.CauchyDistrib
#' @description
#' Computes the Cauchy quantile function
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \tan\left(\pi\left(p - \dfrac{1}{2}\right)\right)}
#' by calling [stats::qcauchy()]. The median is \eqn{\mu} and the quartiles are
#' \eqn{\mu \pm \sigma}, so the two parameters are read off this function
#' directly: they are the median and the half-interquartile range, the
#' location and spread this family has in place of a mean and a standard
#' deviation.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
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
#'
#' @return A numeric vector of quantiles, of length
#'   `max(length(p), length(mu), length(sigma))`.
#'
#' @seealso [distrib_cdf.CauchyDistrib()], which this inverts;
#'   [skewness.CauchyDistrib()] for why the median and the half-interquartile
#'   range are what this family offers; [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The median is mu and the quartiles are mu -/+ sigma exactly.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.01, 0.5, 0.99)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The tails are heavy, so the extreme quantiles are far out: the 99.9th
#' # percentile sits at over 300 scale units.
#' (distrib_quantile(d, 0.999, th) - 0.4) / 1.5
S7::method(distrib_quantile, CauchyDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qcauchy(
    p = p,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Cauchy Random Number Generator
#' @name distrib_rng.CauchyDistrib
#' @description
#' Draws `n` independent Cauchy variates by calling [stats::rcauchy()], so the
#' draws come from R's own generator and depend on `.Random.seed` in the usual
#' way. The generalized ratio-of-uniforms fallback the base class supplies is
#' bypassed.
#'
#' A sample from this family carries extreme values at a rate that does not
#' diminish with `n`, and its running mean does not settle: the sample mean of
#' `n` Cauchy draws is itself Cauchy with the same scale, whatever `n` is.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled.
#'   `sigma` must be strictly positive.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.CauchyDistrib()] for the inverse-transform route,
#'   [skewness.CauchyDistrib()] for the non-convergence of the sample mean, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#'
#' # Same generator as stats::rcauchy, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.5))
#' set.seed(2)
#' identical(a, rcauchy(3, location = 0.4, scale = 1.5))
#'
#' # The sample median converges on mu; the sample mean does not converge at
#' # all, being itself Cauchy at every sample size.
#' set.seed(4)
#' z <- distrib_rng(d, 1e5, list(mu = 3, sigma = 2))
#' n <- c(1e2, 1e3, 1e4, 1e5)
#' rbind(median = vapply(n, function(k) median(z[1:k]), numeric(1)),
#'       mean   = vapply(n, function(k) mean(z[1:k]), numeric(1)))
S7::method(distrib_rng, CauchyDistrib) <- function(distrib, n, theta) {
  stats::rcauchy(
    n = n,
    location = theta[[1]],
    scale = theta[[2]]
  )
}

#' @title Cauchy Score
#' @name distrib_gradient.CauchyDistrib
#' @description
#' Computes the first derivatives of the Cauchy log-density with respect to
#' \eqn{\mu} and \eqn{\sigma}, one value per observation, in closed form.
#' Writing \eqn{r = y - \mu} and \eqn{d = \sigma^2 + r^2},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{2r}{d},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{r^2 - \sigma^2}{\sigma d}.}
#'
#' The score in \eqn{\mu} is **redescending**: it rises to \eqn{1/\sigma} at
#' \eqn{r = \sigma} and falls back towards zero as \eqn{|r|} grows, so an
#' observation far from the location contributes almost nothing to the
#' estimating equation. That is the mechanism behind the robustness of a Cauchy
#' likelihood, and it is also why the log-likelihood can have several local
#' maxima.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
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
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\sigma > 0} the scale. Neither is a moment: no moment of this family
#' exists.
#'
#' @seealso [distrib_hessian.CauchyDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.CauchyDistrib()] for their expectation,
#'   [distrib_gradient.Gaussian1Distrib()] for the unbounded score of a
#'   light-tailed family, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' r <- y - 0.4; dd <- 1.5^2 + r^2
#' all.equal(g$mu, 2 * r / dd)
#' all.equal(g$sigma, (r^2 - 1.5^2) / (1.5 * dd))
#'
#' # The score redescends: it peaks at r = sigma, where it is 1/sigma, and
#' # decays afterwards. An outlier is discounted, not chased.
#' r <- c(0, 0.75, 1.5, 3, 6, 12)
#' round(distrib_gradient(d, 0.4 + r, th)$mu, 4)
#' 1 / 1.5
#'
#' # A Gaussian score at the same residuals grows without bound instead.
#' round(distrib_gradient(gaussian1_distrib(), 0.4 + r, th)$mu, 4)
S7::method(distrib_gradient, CauchyDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  cauchy_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Observed Hessian
#' @name distrib_hessian.CauchyDistrib
#' @description
#' Computes the three distinct second derivatives of the Cauchy log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, one value per observation, in
#' closed form. Writing \eqn{r = y - \mu} and \eqn{d = \sigma^2 + r^2},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(r^2 - \sigma^2)}{d^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 r^2 - r^4}{\sigma^2 d^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma} = -\dfrac{4\sigma r}{d^2}.}
#'
#' The curvature in \eqn{\mu} turns **positive** wherever \eqn{|r| > \sigma}, so
#' a single observation beyond one scale unit contributes convexity and the
#' observed information can fail to be positive definite. The expected values
#' are well behaved, which is why Fisher scoring is the steadier route on this
#' family; see [distrib_expected_hessian.CauchyDistrib()].
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
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
#' @seealso [distrib_gradient.CauchyDistrib()] for the score,
#'   [distrib_expected_hessian.CauchyDistrib()] for the expectation of this
#'   quantity, [fisher_scoring()] for the estimation route that uses it, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' h <- distrib_hessian(d, y, th)
#'
#' # The closed forms, written out.
#' r <- y - 0.4; dd <- 1.5^2 + r^2
#' all.equal(h$mu_mu, 2 * (r^2 - 1.5^2) / dd^2)
#' all.equal(h$mu_sigma, -4 * 1.5 * r / dd^2)
#'
#' # The curvature in mu is positive wherever |r| exceeds sigma.
#' data.frame(r = r, mu_mu = h$mu_mu, beyond_sigma = abs(r) > 1.5)
#'
#' # A central difference of the score reproduces it.
#' eps <- 1e-5
#' up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
S7::method(distrib_hessian, CauchyDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  cauchy_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Expected Hessian
#' @name distrib_expected_hessian.CauchyDistrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] =
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{2\sigma^2},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}\right] = 0.}
#' The information in the location is \eqn{1/(2\sigma^2)}, half of a Gaussian's
#' \eqn{1/\sigma^2} at the same scale, which is the price the heavy tails
#' charge. The zero off-diagonal says the location and the scale are
#' orthogonal. The expectations exist even though no moment of the family does:
#' they are expectations of bounded functions of \eqn{y}, not of \eqn{y}.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
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
#' expectation of the **observed information** under the model. The Cauchy is a
#' regular family in both parameters, so the second Bartlett identity holds and
#' this equals the variance of the score.
#'
#' @seealso [distrib_hessian.CauchyDistrib()] for the observed quantity this is
#'   the expectation of, [distrib_expected_hessian.Gaussian1Distrib()] for the
#'   light-tailed comparison, [fisher_scoring()], which inverts it at each step,
#'   and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # Both diagonal entries are -1/(2 sigma^2); the mixed entry is 0.
#' lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#' -1 / (2 * 1.5^2)
#'
#' # Half the information a Gaussian of the same scale carries in its location.
#' distrib_expected_hessian(gaussian1_distrib(), 0, th)$mu_mu
#'
#' # The observed Hessian averages onto it over a large sample, even though
#' # the sample mean of the same draws does not converge at all.
#' set.seed(3)
#' z <- distrib_rng(d, 2e5, th)
#' vapply(distrib_hessian(d, z, th), mean, numeric(1))
S7::method(distrib_expected_hessian, CauchyDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  cauchy_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Third-Order Derivatives
#' @name distrib_deriv3.CauchyDistrib
#' @description
#' Computes the four distinct third derivatives of the Cauchy log-density with
#' respect to \eqn{\mu} and \eqn{\sigma}, in closed form. Every component is a
#' rational function of \eqn{r = y - \mu} and \eqn{\sigma} over a power of
#' \eqn{d = \sigma^2 + r^2}, so all of them are bounded in \eqn{y} and decay as
#' \eqn{|r|} grows.
#'
#' With `expected = TRUE` the expectations under the model are returned, also
#' in closed form: the two components odd in \eqn{r} vanish by symmetry, and
#' the other two are \eqn{1/(2\sigma^3)} and \eqn{3/(2\sigma^3)}. Both routes
#' are closed form, so no quadrature is run and `approx` and `nsim` are ignored.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
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
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.CauchyDistrib()] for the order below and
#'   [distrib_deriv4.CauchyDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Expected values: the two components odd in the residual vanish.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#' c(1 / (2 * 1.5^3), 3 / (2 * 1.5^3))
#'
#' # A central difference of the Hessian reproduces the observed component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
S7::method(distrib_deriv3, CauchyDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) cauchy_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else cauchy_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Fourth-Order Derivatives
#' @name distrib_deriv4.CauchyDistrib
#' @description
#' Computes the five distinct fourth derivatives of the Cauchy log-density with
#' respect to \eqn{\mu} and \eqn{\sigma}, in closed form. As at third order,
#' every component is a rational function of \eqn{r = y - \mu} and \eqn{\sigma}
#' over a power of \eqn{d = \sigma^2 + r^2} and so is bounded in \eqn{y}.
#'
#' With `expected = TRUE` the expectations under the model are returned, also
#' in closed form; the two components odd in \eqn{r} vanish by symmetry. Both
#' routes are closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
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
#'   `mu_mu_mu_sigma`, `mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
#'   `sigma_sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k l)}} is the fourth derivative of the log-density with
#' respect to parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.CauchyDistrib()] for the order below,
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' names(distrib_deriv4(d, y, th))
#'
#' # Expected values: the two components odd in the residual vanish.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, CauchyDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) cauchy_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else cauchy_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy First Derivative in the Response
#' @name distrib_grad_y.CauchyDistrib
#' @description
#' Computes the first derivative of the Cauchy log-density with respect to the
#' response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{2r}{\sigma^2 + r^2},
#'       \qquad r = y - \mu,}
#' in closed form. The Cauchy is a location family in \eqn{\mu}, so the response
#' enters the log-density only through \eqn{r} and this derivative is the
#' negative of the score in \eqn{\mu}. It is bounded by \eqn{1/\sigma}, unlike
#' the corresponding derivative of a Gaussian, which grows without bound.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_hess_y.CauchyDistrib()] for the second derivative in the
#'   response, [distrib_gradient.CauchyDistrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' r <- y - 0.4
#' all.equal(distrib_grad_y(d, y, th), -2 * r / (1.5^2 + r^2))
#'
#' # A location family: the derivative in the response is minus the score in
#' # the location.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Bounded by 1/sigma however far out the observation is.
#' max(abs(distrib_grad_y(d, seq(-1e3, 1e3, length.out = 1e4), th)))
#' 1 / 1.5
S7::method(distrib_grad_y, CauchyDistrib) <- function(distrib, y, theta) {
  r <- y - theta[[1]]
  -2 * r / (theta[[2]]^2 + r^2)
}

#' @title Cauchy Second Derivative in the Response
#' @name distrib_hess_y.CauchyDistrib
#' @description
#' Computes the second derivative of the Cauchy log-density with respect to the
#' response,
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = \dfrac{2(r^2 - \sigma^2)}{(\sigma^2 + r^2)^2},
#'       \qquad r = y - \mu,}
#' in closed form. Being a location family, the Cauchy has the same curvature in
#' the response as in its location, so this equals the `mu_mu` component of
#' [distrib_hessian.CauchyDistrib()]. It is negative only for \eqn{|r| <
#' \sigma}: the log-density is concave near the location and convex in the
#' tails, which is how heavy tails appear on the log scale.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_grad_y.CauchyDistrib()] for the first derivative in the
#'   response, [distrib_hessian.CauchyDistrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- cauchy_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same curvature as in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Concave within one scale unit of the location, convex beyond it.
#' z <- c(0, 1, 1.5, 2, 5)
#' data.frame(r = z * 1.5,
#'            hess_y = distrib_hess_y(d, 0.4 + z * 1.5, th),
#'            concave = distrib_hess_y(d, 0.4 + z * 1.5, th) < 0)
S7::method(distrib_hess_y, CauchyDistrib) <- function(distrib, y, theta) {
  r <- y - theta[[1]]
  s2 <- theta[[2]]^2
  2 * (r^2 - s2) / (s2 + r^2)^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Cauchy Distribution
#'
#' @description
#' Builds the distribution object for the Cauchy family with location
#' \eqn{\mu} and scale \eqn{\sigma > 0}, the Student's t on one degree of
#' freedom. The returned object carries closed-form derivatives of the
#' log-density to fourth order, in the parameters and in the response, so every
#' generic of the toolkit answers without a numerical fallback.
#'
#' No moment of this family exists, so the two parameters are the median and
#' the half-interquartile range. Its tails decay like \eqn{y^{-2}}, and the
#' resulting score is bounded, which makes a Cauchy likelihood a standard
#' choice when the data carry gross outliers.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the location
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the location
#'   ranging over the whole line already.
#' @param link_sigma A `link` object from `linkfunctions7` for the scale
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right]},}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\sigma \in (0, \infty)}. The
#' distribution function is
#' \eqn{F(q) = 1/2 + \pi^{-1}\arctan((q-\mu)/\sigma)} and the quantile function
#' \eqn{Q(p) = \mu + \sigma\tan(\pi(p - 1/2))}.
#'
#' # No moments
#'
#' The density decays like \eqn{y^{-2}}, so \eqn{\int |y|^p f(y)\,dy} diverges
#' for every \eqn{p \ge 1} and the family has no mean, no variance and no
#' higher moment. `mean()`, `variance()`, `skewness()` and `kurtosis()` return
#' `NaN` directly, without attempting a quadrature over a divergent integral;
#' [skewness.CauchyDistrib()] gives the argument in full. What the family does
#' have is a median, equal to \eqn{\mu}, and quartiles at \eqn{\mu \pm \sigma},
#' so \eqn{\sigma} is the half-interquartile range.
#'
#' # Derivatives
#'
#' With \eqn{r = y - \mu} and \eqn{d = \sigma^2 + r^2} the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{2r}{d}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{r^2 - \sigma^2}{\sigma d},}
#' the observed Hessian
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(r^2 - \sigma^2)}{d^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 r^2 - r^4}{\sigma^2 d^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma} = -\dfrac{4\sigma r}{d^2},}
#' and its expectation
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] =
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{2\sigma^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0.}
#' The information in the location is half a Gaussian's at the same scale. Every
#' derivative is bounded in \eqn{y}, so an expectation of any of them exists
#' even though no moment of \eqn{Y} does.
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.CauchyDistrib()] and [distrib_deriv4.CauchyDistrib()], as are
#' the derivatives in the response.
#'
#' # Estimation
#'
#' There is no closed-form estimate: the likelihood equations are polynomial in
#' \eqn{\mu} of degree \eqn{2n - 1} and are solved numerically. Because the
#' score redescends, the log-likelihood can carry several local maxima, and
#' [fit_distrib()] starts from a data-based value and takes Fisher scoring
#' steps, whose expected information is positive definite where the observed
#' one need not be.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\sigma > 0} the scale. Neither is a moment. \eqn{\eta} is a parameter on
#' the unconstrained scale of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `CauchyDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"cauchy"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [student_t1_distrib()], of which this is the case \eqn{\nu = 1}, and which
#' estimates the degrees of freedom instead of fixing them;
#' [laplace_distrib()] and [pseudohuber_distrib()] for the other robust
#' location families; [gaussian1_distrib()] for the light-tailed comparison;
#' [fit_distrib()] to estimate the parameters; [CauchyDistrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 16.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dcauchy pcauchy qcauchy rcauchy
#'
#' @examples
#' d <- cauchy_distrib()
#' d
#'
#' # The density and the distribution function are R's own.
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' all.equal(distrib_pdf(d, y, th), dcauchy(y, 0.4, 1.5))
#'
#' # No moment exists; the median and the half-interquartile range do.
#' c(mean = mean(d, th), variance = variance(d, th))
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#'
#' # Fitting a contaminated sample: the Cauchy location tracks the bulk where
#' # the sample mean is dragged away by the outliers.
#' set.seed(5)
#' z <- c(rnorm(200, mean = 3, sd = 1), rnorm(10, mean = 60, sd = 1))
#' fit <- fit_distrib(d, z)
#' c(cauchy_mu = unname(coef(fit)["mu"]),
#'   median = median(z), mean = mean(z))
#'
#' @export
cauchy_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  
  CauchyDistrib(
    distrib_name = "cauchy",
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
