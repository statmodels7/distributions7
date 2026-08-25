#' @include distrib.R generics.R
NULL

#' @title Gumbel Distribution Class
#' @name GumbelDistrib
#'
#' @description
#' The S7 class of the Gumbel (type I extreme value) family in the form for
#' **maxima**, a location-scale family on the whole real line with location
#' \eqn{\mu} and scale \eqn{\sigma > 0} and density
#' \eqn{\sigma^{-1}\exp\{-z - e^{-z}\}} at \eqn{z = (y-\mu)/\sigma}. Its shape
#' is fixed: the skewness and the excess kurtosis are constants and cannot be
#' fitted. It inherits from `continuous_distrib`, so it answers every generic
#' of the `distrib` contract; the eleven methods listed below are registered on
#' it directly and everything else comes from the parent.
#'
#' Build one with [gumbel_distrib()], which supplies the two link functions and
#' fills the properties in. This page documents the raw S7 constructor, which
#' takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `GumbelDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [gumbel_distrib()] they hold `"gumbel"`, `"univariate"`,
#'   `c(-Inf, Inf)`, `c("mu", "sigma")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [gumbel_distrib()] to build one;
#'   [weibull1_distrib()], which is this family under \eqn{e^{-Y}};
#'   [gpd_distrib()] for the other limit law of extremes;
#'   [distrib_pdf.GumbelDistrib()] and [distrib_gradient.GumbelDistrib()] for
#'   the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.GumbelDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.GumbelDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.GumbelDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GumbelDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.GumbelDistrib],
#'   [`distrib_gradient()`][distrib_gradient.GumbelDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.GumbelDistrib],
#'   [`distrib_hessian()`][distrib_hessian.GumbelDistrib],
#'   [`distrib_pdf()`][distrib_pdf.GumbelDistrib],
#'   [`distrib_quantile()`][distrib_quantile.GumbelDistrib],
#'   [`distrib_rng()`][distrib_rng.GumbelDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_bounds
#' d@bounds
#'
#' # The shape is fixed: skewness and excess kurtosis do not move with the
#' # parameters.
#' rbind(a = c(skew = skewness(d, list(mu = 0, sigma = 1)),
#'             kurt = kurtosis(d, list(mu = 0, sigma = 1))),
#'       b = c(skewness(d, list(mu = 5, sigma = 3)),
#'             kurtosis(d, list(mu = 5, sigma = 3))))
GumbelDistrib <- S7::new_class("GumbelDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gumbel Probability Density Function
#' @name distrib_pdf.GumbelDistrib
#' @description
#' Computes the Gumbel density, with \eqn{z = (y - \mu)/\sigma},
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sigma}
#'       \exp\left\{-z - e^{-z}\right\},}
#' written out rather than delegated, base R carrying no Gumbel. The two tails
#' are of very different weight: the density falls like \eqn{e^{-z}} to the
#' right and like \eqn{e^{-e^{-z}}} to the left, so the left one is doubly
#' exponential and dies far faster.
#'
#' With `log = TRUE` the exponent is returned directly and stays finite where
#' the density itself underflows.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.GumbelDistrib()] for the distribution function,
#'   [distrib_gradient.GumbelDistrib()] for the derivatives of the
#'   log-density, [distrib_pdf.Weibull1Distrib()] for the family this becomes
#'   under \eqn{e^{-Y}}, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#'
#' # The closed form, written out.
#' z <- (y - 0) / 1
#' all.equal(distrib_pdf(d, y, th), exp(-z - exp(-z)))
#'
#' # The mode is at mu, where the density peaks at 1/(sigma e).
#' c(at_mode = distrib_pdf(d, 0, th), one_over_e = 1 / exp(1))
#'
#' # The two tails are of very different weight.
#' distrib_pdf(d, c(-5, 5), th)
#'
#' # Far to the left the density underflows and its logarithm does not.
#' distrib_pdf(d, -6, th)
#' distrib_pdf(d, -6, th, log = TRUE)
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
#' Computes the Gumbel distribution function, with \eqn{z = (q - \mu)/\sigma},
#' \deqn{F(q; \mu, \sigma) = \exp\left\{-e^{-z}\right\}.}
#' This is the defining property of the family: a maximum of \eqn{n}
#' independent light-tailed variables has, after centering and scaling, a
#' distribution function that is the \eqn{n}th power of one distribution
#' function, and \eqn{\exp\{-e^{-z}\}} is the fixed point of that operation.
#'
#' @details
#' Both tails are computed on the scale that keeps them accurate. The lower
#' tail is exact on the log scale, \eqn{\log F = -e^{-z}}, and the upper tail
#' goes through [base::expm1()], so neither loses precision where it is small.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
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
#' @seealso [distrib_quantile.GumbelDistrib()] for the inverse,
#'   [distrib_pdf.GumbelDistrib()] for the density, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' th <- list(mu = 0, sigma = 1)
#' q <- c(-1, 0, 1)
#'
#' # The closed form, written out.
#' all.equal(distrib_cdf(d, q, th), exp(-exp(-q)))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 1, th) + distrib_cdf(d, 1, th, lower.tail = FALSE)
#'
#' # Max-stability: the cdf raised to the nth power is the same law shifted
#' # by sigma log n.
#' n <- 10
#' all.equal(distrib_cdf(d, q, th)^n,
#'           distrib_cdf(d, q, list(mu = log(n), sigma = 1)))
#'
#' # Deep in the lower tail the probability underflows and its log does not.
#' distrib_cdf(d, -700, th)
#' distrib_cdf(d, -700, th, log.p = TRUE)
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
#' Computes the Gumbel quantile function in closed form,
#' \deqn{Q(p; \mu, \sigma) = \mu - \sigma \log(-\log p).}
#' The distribution function is strictly increasing on the whole line, so this
#' is its exact inverse and the round trip through
#' [distrib_cdf.GumbelDistrib()] returns `p`. The median is
#' \eqn{\mu - \sigma\log\log 2}, about \eqn{\mu + 0.3665\sigma}, and lies
#' between the mode \eqn{\mu} and the mean \eqn{\mu + \gamma\sigma}.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`;
#'   the endpoints give `-Inf` and `Inf`.
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
#' @seealso [distrib_cdf.GumbelDistrib()], which this inverts;
#'   [distrib_rng.GumbelDistrib()], which uses the same inverse-transform
#'   identity; [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' th <- list(mu = 0, sigma = 1)
#'
#' # The closed form, written out.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_quantile(d, p, th), -log(-log(p)))
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # Mode, median and mean, in that order.
#' c(mode = 0, median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#'
#' # A return level: the value exceeded once in 100 periods.
#' distrib_quantile(d, 1 - 1 / 100, th)
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
#' Draws `n` independent Gumbel variates as \eqn{\mu - \sigma\log E} with
#' \eqn{E} standard exponential, using [stats::rexp()]. That is the
#' inverse-transform identity written without a logarithm of a uniform, which
#' keeps the left tail accurate: a uniform near zero has few significant digits
#' left after one logarithm and none after two.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting.
#'   `sigma` must be strictly positive.
#'
#' @return A numeric vector of `n` draws on the whole real line.
#'
#' @seealso [distrib_quantile.GumbelDistrib()] for the inverse it rests on,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#'
#' # The identity the method uses: mu - sigma log(E), E standard exponential.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0, sigma = 1))
#' set.seed(2)
#' identical(a, 0 - 1 * log(rexp(3)))
#'
#' # The sample moments recover the parameters through the fixed shape: the
#' # scale is sd sqrt(6)/pi and the location the mean less gamma times it.
#' set.seed(8)
#' z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
#' s <- sd(z) * sqrt(6) / pi
#' c(mu = mean(z) + digamma(1) * s, sigma = s)
S7::method(distrib_rng, GumbelDistrib) <- function(distrib, n, theta) {
  theta[[1]] - theta[[2]] * log(stats::rexp(n))
}

#' @title Gumbel Score
#' @name distrib_gradient.GumbelDistrib
#' @description
#' Computes the first derivatives of the Gumbel log-density with respect to
#' \eqn{\mu} and \eqn{\sigma}, one value per observation, in closed form. With
#' \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1 - w}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} =
#'         \dfrac{z(1 - w) - 1}{\sigma}.}
#'
#' The whole family is written in \eqn{w}, and under the model \eqn{w} is
#' **standard exponential** whatever the parameters. Every expectation this
#' family needs is therefore a derivative of \eqn{\Gamma} at 2, and the same
#' fact bounds the score above in \eqn{\mu} and leaves it unbounded below: a
#' value far to the left makes \eqn{w} enormous.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each of
#'   length `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\sigma > 0} the scale. \eqn{z = (y-\mu)/\sigma} is the standardized
#' value and \eqn{w = e^{-z}}, which is standard exponential under the model.
#'
#' @seealso [distrib_hessian.GumbelDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.GumbelDistrib()] for their expectation,
#'   [distrib_grad_y.GumbelDistrib()] for the derivative in the response, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' z <- (y - 0) / 1
#' w <- exp(-z)
#' all.equal(g$mu, (1 - w) / 1)
#' all.equal(g$sigma, (z * (1 - w) - 1) / 1)
#'
#' # The location component vanishes at y = mu, where w = 1.
#' distrib_gradient(d, 0, th)$mu
#'
#' # It is bounded above by 1/sigma and unbounded below.
#' distrib_gradient(d, c(-3, 10), th)$mu
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' set.seed(8)
#' s <- distrib_rng(d, 2000, list(mu = 3, sigma = 2))
#' fit <- fit_distrib(d, s)
#' vapply(distrib_gradient(d, s, as.list(coef(fit))), sum, numeric(1))
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

#' @title Gumbel Observed Hessian
#' @name distrib_hessian.GumbelDistrib
#' @description
#' Computes the three distinct second derivatives of the Gumbel log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, one value per observation, in
#' closed form. With \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}},
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{w}{\sigma^2}, \qquad
#'       \ell^{(\mu\sigma)} = -\dfrac{1 - w + zw}{\sigma^2}, \qquad
#'       \ell^{(\sigma\sigma)} =
#'         \dfrac{1 - 2z + 2zw - z^2 w}{\sigma^2}.}
#' The curvature in the location is negative at every observation, \eqn{w}
#' being positive, so the log-likelihood is concave in \eqn{\mu} at any data
#' set. The curvature in the scale is not: it is positive wherever
#' \eqn{1 - 2z + 2zw - z^2w > 0}, which includes \eqn{y = \mu}.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, in that order, each of length
#'   `max(length(y), length(mu), length(sigma))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#' \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}}.
#'
#' @seealso [distrib_gradient.GumbelDistrib()] for the score,
#'   [distrib_expected_hessian.GumbelDistrib()] for the expectation of this
#'   quantity, [distrib_deriv3.GumbelDistrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # The three closed forms, written out.
#' z <- (y - 0) / 1
#' w <- exp(-z)
#' all.equal(h$mu_mu, -w)
#' all.equal(h$mu_sigma, -(1 - w + z * w))
#' all.equal(h$sigma_sigma, 1 - 2 * z + 2 * z * w - z^2 * w)
#'
#' # Concave in mu everywhere, and positive in sigma at y = mu.
#' c(all_negative_in_mu = all(h$mu_mu < 0),
#'   sigma_at_mu = distrib_hessian(d, 0, th)$sigma_sigma)
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 0 + eps, sigma = 1))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0 - eps, sigma = 1))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
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

#' @title Gumbel Expected Hessian
#' @name distrib_expected_hessian.GumbelDistrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. With \eqn{\gamma} the
#' Euler-Mascheroni constant,
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\sigma^2},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\sigma)}\right] =
#'         \dfrac{1 - \gamma}{\sigma^2}, \qquad
#'       \mathbb{E}\left[\ell^{(\sigma\sigma)}\right] =
#'         -\dfrac{(1-\gamma)^2 + \pi^2/6}{\sigma^2}.}
#'
#' They are closed form for one reason, and it is the same reason the Weibull's
#' are: under the model \eqn{w = e^{-Z}} is standard exponential whatever the
#' parameters, so every expectation the family needs is a derivative of
#' \eqn{\Gamma} at 2. Here that is \eqn{\mathbb{E}[w] = 1},
#' \eqn{\mathbb{E}[w\log w] = 1 - \gamma} and
#' \eqn{\mathbb{E}[w(\log w)^2] = (1-\gamma)^2 + \pi^2/6 - 1}, with
#' \eqn{\mathbb{E}[Z] = \gamma}.
#'
#' **The mixed entry does not vanish.** The location and the scale are not
#' orthogonal here, where in a symmetric location-scale family such as the
#' Gaussian they are, and the reason is that this density is skewed. Their
#' maximum likelihood estimates are asymptotically correlated.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `sigma`
#'   must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available. Accepted so that the
#'   signature matches the generic's, where it selects between the Bartlett,
#'   quadrature, Monte Carlo and outer-product routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, in that order, each of length `length(y)` and constant within
#'   itself when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The Gumbel is a
#' regular family, so the second Bartlett identity holds and this equals the
#' variance of the score. \eqn{\gamma} is the Euler-Mascheroni constant,
#' \eqn{-\psi(1) \approx 0.5772}.
#'
#' @seealso [distrib_hessian.GumbelDistrib()] for the observed quantity this is
#'   the expectation of, [distrib_expected_hessian.Weibull1Distrib()] for the
#'   family that shares these expectations, [fisher_scoring()], which inverts
#'   it at each step, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' th <- list(mu = 0, sigma = 1)
#' e <- distrib_expected_hessian(d, c(-1, 0, 1), th)
#' lapply(e, unique)
#'
#' # Written out with the Euler-Mascheroni constant.
#' g <- -digamma(1)
#' c(-1, -((1 - g)^2 + pi^2 / 6), 1 - g)
#'
#' # The mixed entry is not zero, so the location and the scale are correlated.
#' e$mu_sigma[1]
#'
#' # Negative definite all the same.
#' M <- matrix(c(e$mu_mu[1], e$mu_sigma[1], e$mu_sigma[1], e$sigma_sigma[1]), 2)
#' eigen(M, only.values = TRUE)$values
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
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

#' @title Gumbel Third-Order Derivatives
#' @name distrib_deriv3.GumbelDistrib
#' @description
#' Computes the four distinct third derivatives of the Gumbel log-density with
#' respect to \eqn{\mu} and \eqn{\sigma}, in closed form. With
#' \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}}, every component is a
#' polynomial in \eqn{z} and in \eqn{z^j w} divided by a power of \eqn{\sigma}.
#'
#' With `expected = TRUE` the expectations are returned, also in closed form.
#' They rest on \eqn{w} being standard exponential under the model, which makes
#' \eqn{\mathbb{E}[z^k w] = (-1)^k \Gamma^{(k)}(2)} and
#' \eqn{\mathbb{E}[z] = \gamma}: every expectation the family needs is a
#' derivative of \eqn{\Gamma} at 2, assembled from polygamma functions there.
#' Since both routes are closed form, `approx` and `nsim` are ignored.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
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
#'   `max(length(y), length(mu), length(sigma))`. The names enumerate the
#'   distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{z = (y-\mu)/\sigma}, \eqn{w = e^{-z}} and
#' \eqn{\gamma} is the Euler-Mascheroni constant.
#'
#' @seealso [distrib_hessian.GumbelDistrib()] for the order below and
#'   [distrib_deriv4.GumbelDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # The pure location component is -w/sigma^3, which is the curvature in mu
#' # divided by sigma. Shown at a scale where the two differ.
#' th2 <- list(mu = 0, sigma = 2)
#' all.equal(distrib_deriv3(d, y, th2)$mu_mu_mu,
#'           distrib_hessian(d, y, th2)$mu_mu / 2)
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 0 + eps, sigma = 1))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0 - eps, sigma = 1))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-5)
S7::method(distrib_deriv3, GumbelDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) gumbel_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gumbel_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gumbel Fourth-Order Derivatives
#' @name distrib_deriv4.GumbelDistrib
#' @description
#' Computes the five distinct fourth derivatives of the Gumbel log-density with
#' respect to \eqn{\mu} and \eqn{\sigma}, in closed form, in the notation of
#' [distrib_deriv3.GumbelDistrib()]: with \eqn{z = (y-\mu)/\sigma} and
#' \eqn{w = e^{-z}}, every component is a polynomial in \eqn{z} and in
#' \eqn{z^j w} divided by a power of \eqn{\sigma}.
#'
#' With `expected = TRUE` the expectations are returned, also closed form,
#' resting on the same fact: \eqn{w} is standard exponential under the model,
#' so \eqn{\mathbb{E}[z^k w] = (-1)^k \Gamma^{(k)}(2)}. `approx` and `nsim` are
#' ignored either way.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
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
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}.
#' \eqn{z = (y-\mu)/\sigma} and \eqn{w = e^{-z}}.
#'
#' @seealso [distrib_deriv3.GumbelDistrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.GumbelDistrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 0 + eps, sigma = 1))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0 - eps, sigma = 1))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-4)
S7::method(distrib_deriv4, GumbelDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) gumbel_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gumbel_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gumbel First Derivative in the Response
#' @name distrib_grad_y.GumbelDistrib
#' @description
#' Computes the first derivative of the Gumbel log-density with respect to the
#' response, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{w - 1}{\sigma},
#'       \qquad w = e^{-(y-\mu)/\sigma}.}
#' The Gumbel is a location family in \eqn{\mu}, so the response enters the
#' log-density only through \eqn{y - \mu} and this derivative is the negative
#' of the score in \eqn{\mu}. It vanishes at \eqn{y = \mu}, which is the mode.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_hess_y.GumbelDistrib()] for the second derivative in the
#'   response, [distrib_gradient.GumbelDistrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#'
#' all.equal(distrib_grad_y(d, y, th), (exp(-y) - 1) / 1)
#'
#' # A location family: the derivative in the response is minus the score in
#' # the location.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Zero at the mode, which is mu itself.
#' distrib_grad_y(d, 0, th)
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, GumbelDistrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  (exp(-(y - theta[[1]]) / sigma) - 1) / sigma
}

#' @title Gumbel Second Derivative in the Response
#' @name distrib_hess_y.GumbelDistrib
#' @description
#' Computes the second derivative of the Gumbel log-density with respect to the
#' response, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{w}{\sigma^2},
#'       \qquad w = e^{-(y-\mu)/\sigma}.}
#' It is negative at every observation, \eqn{w} being positive, so the
#' log-density is concave in the response over the whole line. Being a location
#' family, the Gumbel has the same curvature in the response as in its
#' location, and this equals the `mu_mu` component of
#' [distrib_hessian.GumbelDistrib()].
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation,
#'   all strictly negative.
#'
#' @seealso [distrib_grad_y.GumbelDistrib()] for the first derivative in the
#'   response, [distrib_hessian.GumbelDistrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- gumbel_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same curvature as in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Negative everywhere, so the log-density is concave in the response.
#' all(distrib_hess_y(d, c(-20, 0, 20), th) < 0)
S7::method(distrib_hess_y, GumbelDistrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  -exp(-(y - theta[[1]]) / sigma) / sigma^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Gumbel Distribution
#'
#' @description
#' Builds the distribution object for the Gumbel (type I extreme value) family
#' in the form for **maxima**: a location-scale family on the whole real line
#' with location \eqn{\mu} and scale \eqn{\sigma > 0}. The returned object
#' carries closed-form derivatives of the log-density to fourth order, in the
#' parameters and in the response, observed and expected, and closed-form
#' moments, so every generic of the toolkit answers without a numerical
#' fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. The defaults are the identity for
#' the location, which ranges over the whole line, and the logarithm for the
#' scale, which keeps it positive at every predictor.
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
#' With \eqn{z = (y-\mu)/\sigma} the density on the whole real line is
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sigma}
#'       \exp\left\{-z - e^{-z}\right\},}
#' the distribution function \eqn{F(q) = \exp\{-e^{-z}\}} and the quantile
#' function \eqn{Q(p) = \mu - \sigma\log(-\log p)}, both closed form. The mean
#' is \eqn{\mu + \gamma\sigma} with \eqn{\gamma} the Euler-Mascheroni
#' constant, the variance \eqn{\pi^2\sigma^2/6}, the skewness
#' \eqn{12\sqrt{6}\,\zeta(3)/\pi^3 \approx 1.1395} and the excess kurtosis
#' \eqn{12/5}. The last two are constants: **this family has a fixed shape and
#' only its location and spread can be fitted.** The mode is \eqn{\mu}, the
#' median \eqn{\mu - \sigma\log\log 2}, and the mean is above both.
#'
#' The law is the limit of the maximum of a sample from a light-tailed
#' distribution, once centered and scaled, and that is what it is fitted to. It
#' is max-stable: raising its distribution function to the \eqn{n}th power
#' gives the same law with the location shifted by \eqn{\sigma\log n}.
#'
#' # Everything rests on w being standard exponential
#'
#' Write \eqn{w = e^{-z}}. Under the model \eqn{w} is standard exponential
#' whatever \eqn{\mu} and \eqn{\sigma} are, and the whole family is written in
#' it. The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1 - w}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} =
#'         \dfrac{z(1 - w) - 1}{\sigma},}
#' the observed Hessian is
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{w}{\sigma^2}, \quad
#'       \ell^{(\mu\sigma)} = -\dfrac{1 - w + zw}{\sigma^2}, \quad
#'       \ell^{(\sigma\sigma)} = \dfrac{1 - 2z + 2zw - z^2 w}{\sigma^2},}
#' and every expectation the family needs is a derivative of \eqn{\Gamma} at 2:
#' \eqn{\mathbb{E}[w] = 1}, \eqn{\mathbb{E}[w\log w] = 1-\gamma},
#' \eqn{\mathbb{E}[w(\log w)^2] = (1-\gamma)^2 + \pi^2/6 - 1}. The expected
#' Hessian follows in one line, and so do the third and fourth orders, through
#' \eqn{\mathbb{E}[z^k w] = (-1)^k\Gamma^{(k)}(2)}.
#'
#' **The location and the scale are not orthogonal.** The mixed entry of the
#' expected information is \eqn{(1-\gamma)/\sigma^2}, which is about
#' \eqn{0.4228/\sigma^2} and never zero. In a symmetric location-scale family
#' such as the Gaussian it vanishes; here the density is skewed and it does
#' not, so the two estimates are asymptotically correlated.
#'
#' # Relatives
#'
#' If \eqn{Y} is Gumbel then \eqn{e^{-Y}} is Weibull, so
#' [weibull1_distrib()] is this family on the log scale and reversed, and the
#' two share the expectations that produce their information matrices. For
#' **minima** the law is the reflection: fit this family to \eqn{-Y}, or use
#' [transformation()] with
#' [`affine_transform(scale = -1)`][affine_transform]. [gpd_distrib()] is the
#' other limit law of extremes, for exceedances over a threshold in place of
#' block maxima.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' estimate is closed form. The fixed shape supplies the method of moments
#' starting values: the scale is \eqn{s\sqrt{6}/\pi} and the location
#' \eqn{\bar y - \gamma s\sqrt{6}/\pi}, with \eqn{s} the sample standard
#' deviation. The example below shows them landing beside the estimates.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\sigma > 0} the scale. \eqn{z = (y-\mu)/\sigma} is the standardized
#' value and \eqn{w = e^{-z}}, standard exponential under the model.
#' \eqn{\gamma} is the Euler-Mascheroni constant, \eqn{-\psi(1) \approx
#' 0.5772}. \eqn{\eta} is a parameter on the unconstrained scale of its link,
#' with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `GumbelDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"gumbel"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [weibull1_distrib()] for the family this becomes under \eqn{e^{-Y}};
#' [gpd_distrib()] for exceedances over a threshold; [transformation()] for
#' the reflection that gives minima; [gaussian1_distrib()] for the symmetric
#' location-scale family whose parameters are orthogonal; [fit_distrib()] to
#' estimate the parameters; [check_distrib()] to validate a family of your own
#' against the same battery this one passes; [GumbelDistrib] for the class.
#'
#' @references
#' Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
#' Values*, Chapter 3. Springer, London.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats rexp
#'
#' @examples
#' d <- gumbel_distrib()
#' d
#'
#' # The density, written out.
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#' all.equal(distrib_pdf(d, y, th), exp(-y - exp(-y)))
#'
#' # The mean is shifted by Euler's constant and the shape is fixed: the
#' # skewness and kurtosis are the same at any parameter setting.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(skewness(d, list(mu = 5, sigma = 3)), kurtosis(d, list(mu = 5, sigma = 3)))
#'
#' # Max-stable: the cdf to the nth power is the same law shifted by
#' # sigma log n.
#' all.equal(distrib_cdf(d, y, th)^10,
#'           distrib_cdf(d, y, list(mu = log(10), sigma = 1)))
#'
#' # exp(-Y) is Weibull with scale exp(-mu) and shape 1/sigma; the densities
#' # agree once the Jacobian is applied.
#' set.seed(3)
#' yy <- distrib_rng(d, 5, list(mu = 0, sigma = 1 / 2))
#' all.equal(distrib_pdf(d, yy, list(mu = 0, sigma = 1 / 2)),
#'           distrib_pdf(weibull1_distrib(), exp(-yy),
#'                       list(mu = 1, sigma = 2)) * exp(-yy))
#'
#' # The location and the scale are not orthogonal: the mixed entry of the
#' # expected information is (1 - gamma)/sigma^2.
#' c(mixed = distrib_expected_hessian(d, 0, th)$mu_sigma,
#'   one_minus_gamma = 1 + digamma(1))
#'
#' # Fitting recovers the parameters; the moment estimates start it off.
#' set.seed(8)
#' s <- distrib_rng(d, 2000, list(mu = 3, sigma = 2))
#' sc <- sd(s) * sqrt(6) / pi
#' rbind(fitted  = coef(fit_distrib(d, s)),
#'       moments = c(mu = mean(s) + digamma(1) * sc, sigma = sc))
#'
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
