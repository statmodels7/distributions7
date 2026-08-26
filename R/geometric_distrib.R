#' @include distrib.R generics.R
NULL

#' @title Geometric Distribution Class
#' @name GeometricDistrib
#'
#' @description
#' The S7 class of the geometric family parametrized by its **mean**
#' \eqn{\mu > 0}: the number of failures before the first success in
#' independent trials, each succeeding with probability
#' \eqn{p = 1/(1+\mu)}. It inherits from `discrete_distrib`, so its support is
#' a set of isolated points, so expectations are sums and no derivative with
#' respect to the response is defined.
#'
#' The parametrization is by the mean, the quantity the modeling layer above
#' models: a link carries \eqn{\mu} to the unconstrained scale and the
#' probability follows. R's own `dgeom` takes
#' `prob`, and the methods convert through [geom_prob()].
#'
#' The variance is \eqn{\mu(1+\mu)}, always above the mean, so this family is
#' overdispersed relative to a Poisson of the same mean. It is the negative
#' binomial at a dispersion of one, and it is the only discrete law that is
#' memoryless.
#'
#' Build one with [geometric_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `GeometricDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [geometric_distrib()] they hold `"geometric"`,
#'   `"univariate"`, `c(0, Inf)`, `"mu"`, `c(mu = "mean")`, `1`, the domain
#'   \eqn{(0, \infty)}, and the one link.
#'
#' @seealso [geometric_distrib()] to build one;
#'   [negbin2_distrib()], of which this is the case \eqn{\theta = 1};
#'   [exponential_distrib()], its memoryless continuous counterpart;
#'   [poisson_distrib()] for the equidispersed alternative.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.GeometricDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.GeometricDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.GeometricDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GeometricDistrib],
#'   [`distrib_gradient()`][distrib_gradient.GeometricDistrib],
#'   [`distrib_hessian()`][distrib_hessian.GeometricDistrib],
#'   [`distrib_pdf()`][distrib_pdf.GeometricDistrib],
#'   [`distrib_quantile()`][distrib_quantile.GeometricDistrib],
#'   [`distrib_rng()`][distrib_rng.GeometricDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
#'
#' @examples
#' d <- geometric_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#' d@params
#' d@params_interpretation
#'
#' # Overdispersed by construction: the variance is mu(1+mu), above the mean.
#' vapply(c(0.5, 3, 20), function(m) {
#'   th <- list(mu = m)
#'   c(mean = mean(d, th), var = variance(d, th))
#' }, numeric(2))
#'
#' # The mass is R's own at prob = 1/(1+mu).
#' all.equal(distrib_pdf(d, c(0, 2, 7), list(mu = 3)),
#'           dgeom(c(0, 2, 7), prob = 1 / 4))
GeometricDistrib <- S7::new_class("GeometricDistrib", parent = discrete_distrib)

#' The Success Probability Behind a Geometric Mean
#'
#' @description
#' Converts a geometric mean into the success probability that the base R
#' functions [stats::dgeom()], [stats::pgeom()], [stats::qgeom()] and
#' [stats::rgeom()] take as their `prob` argument.
#'
#' @details
#' The mean number of failures before the first success is \eqn{(1-p)/p}, so
#' \eqn{p = 1/(1+\mu)}. Writing it once keeps the four methods that call
#' \pkg{stats} from each repeating the algebra, and keeps them from drifting
#' apart if the parametrization ever changes.
#'
#' The map is decreasing and carries \eqn{(0, \infty)} onto \eqn{(0, 1)}: a
#' large mean is a small success probability. No validation is performed, so a
#' non-positive `mu` returns a value outside \eqn{(0, 1)} and the caller sees
#' the failure at the base R function.
#'
#' @param mu The mean, a positive numeric vector of any length.
#'
#' @return A numeric vector of probabilities in \eqn{(0, 1)}, of the length of
#'   `mu`.
#'
#' @seealso [geometric_distrib()] for the family, and
#'   [distrib_pdf.GeometricDistrib()] for the first of the four callers.
#'
#' @keywords internal
geom_prob <- function(mu) 1 / (1 + mu)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Geometric Probability Mass Function
#' @name distrib_pdf.GeometricDistrib
#' @description
#' Computes the geometric probability mass
#' \deqn{P(Y = y; \mu) = \dfrac{1}{1+\mu}\left(\dfrac{\mu}{1+\mu}\right)^{y},
#'       \qquad y = 0, 1, 2, \dots}
#' by calling [stats::dgeom()] at `prob = 1/(1+mu)`, which [geom_prob()]
#' supplies. The mass falls geometrically, by the constant factor
#' \eqn{\mu/(1+\mu)} at every step, so its maximum is always at zero however
#' large the mean is.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param y A numeric vector of counts. A non-integer or negative value gives 0
#'   with a warning from [stats::dgeom()].
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. A value of length 1 is recycled.
#'   `mu` must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned, which
#'   stays finite for a count far above the mean. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(y), length(mu))`, one value per observation.
#'
#' @seealso [geom_prob()] for the conversion, [distrib_cdf.GeometricDistrib()]
#'   for the distribution function, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' y <- c(0, 2, 7)
#'
#' # The method is stats::dgeom at prob = 1/(1+mu).
#' all.equal(distrib_pdf(d, y, list(mu = 3)), dgeom(y, prob = 1 / (1 + 3)))
#'
#' # A probability mass: it sums to one over the support.
#' sum(distrib_pdf(d, 0:400, list(mu = 3)))
#'
#' # The ratio of consecutive masses is the constant mu/(1+mu), so the mode is
#' # always at zero.
#' m <- distrib_pdf(d, 0:5, list(mu = 3))
#' c(ratios = unique(round(m[-1] / m[-length(m)], 12)), mu_over_1_plus_mu = 3 / 4)
S7::method(distrib_pdf, GeometricDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dgeom(y, prob = geom_prob(theta[[1]]), log = log)
}

#' @title Geometric Cumulative Distribution Function
#' @name distrib_cdf.GeometricDistrib
#' @description
#' Computes the geometric distribution function
#' \deqn{F(q; \mu) = 1 - \left(\dfrac{\mu}{1+\mu}\right)^{\lfloor q \rfloor + 1}}
#' by calling [stats::pgeom()] at `prob = 1/(1+mu)`. The survival function is
#' exactly \eqn{(\mu/(1+\mu))^{q+1}}, a geometric decay, and
#' `lower.tail = FALSE` returns it without forming the difference.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param q A numeric vector of quantiles. A non-integer is floored, and a
#'   negative value gives 0.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `q`. `mu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are the survival
#'   function \eqn{P(Y > q)}, exact far into the tail.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu))`.
#'
#' @seealso [distrib_quantile.GeometricDistrib()] for the generalized inverse,
#'   [distrib_pdf.GeometricDistrib()] for the mass, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- geometric_distrib()
#' th <- list(mu = 3)
#'
#' # The method is stats::pgeom at prob = 1/(1+mu).
#' all.equal(distrib_cdf(d, c(0, 2, 7), th), pgeom(c(0, 2, 7), prob = 1 / 4))
#'
#' # The survival function is a geometric decay, written out.
#' all.equal(distrib_cdf(d, c(0, 2, 7), th, lower.tail = FALSE),
#'           (3 / 4)^(c(0, 2, 7) + 1))
#'
#' # Memoryless: the chance of at least one more failure does not depend on
#' # how many have already been seen.
#' vapply(c(0, 2, 10, 50), function(s)
#'   distrib_cdf(d, s + 1, th, lower.tail = FALSE) /
#'     distrib_cdf(d, s, th, lower.tail = FALSE), numeric(1))
S7::method(distrib_cdf, GeometricDistrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  stats::pgeom(q, prob = geom_prob(theta[[1]]), lower.tail = lower.tail,
               log.p = log.p)
}

#' @title Geometric Quantile Function
#' @name distrib_quantile.GeometricDistrib
#' @description
#' Computes the generalized inverse of the geometric distribution function,
#' \deqn{Q(p; \mu) = \min\{k \in \{0, 1, 2, \dots\} : F(k; \mu) \ge p\},}
#' by calling [stats::qgeom()] at `prob = 1/(1+mu)`. The distribution function
#' is a step function, so the round trip through
#' [distrib_cdf.GeometricDistrib()] returns a probability **at least** `p`.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. `p = 1` gives `Inf`.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `p`. `mu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of non-negative integers, of length
#'   `max(length(p), length(mu))`.
#'
#' @seealso [distrib_cdf.GeometricDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' th <- list(mu = 3)
#'
#' # Integers, and the median well below the mean, the mass being at zero.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # The round trip overshoots, the support being a lattice.
#' p <- c(0.025, 0.5, 0.975)
#' rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
S7::method(distrib_quantile, GeometricDistrib) <- function(distrib, p, theta,
                                                           lower.tail = TRUE,
                                                           log.p = FALSE, ...) {
  stats::qgeom(p, prob = geom_prob(theta[[1]]), lower.tail = lower.tail,
               log.p = log.p)
}

#' @title Geometric Random Number Generator
#' @name distrib_rng.GeometricDistrib
#' @description
#' Draws `n` independent geometric counts by calling [stats::rgeom()] at
#' `prob = 1/(1+mu)`, so the draws come from R's own generator and depend on
#' `.Random.seed` in the usual way.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of length `n`. A value of length 1 is recycled, so a vector
#'   of length `n` draws one count per mean. `mu` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An integer vector of `n` non-negative counts.
#'
#' @seealso [fit_distrib()] to estimate the mean back, and [distrib_rng()] for
#'   the generic.
#'
#' @examples
#' d <- geometric_distrib()
#'
#' # Same generator as stats::rgeom, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 5, list(mu = 3))
#' set.seed(2)
#' identical(a, rgeom(5, prob = 1 / 4))
#'
#' # The sample mean estimates mu and the variance mu(1+mu).
#' set.seed(5)
#' z <- distrib_rng(d, 2e4, list(mu = 2.5))
#' c(mean = mean(z), var = var(z), mu_1_plus_mu = 2.5 * 3.5)
S7::method(distrib_rng, GeometricDistrib) <- function(distrib, n, theta, ...) {
  stats::rgeom(n, prob = geom_prob(theta[[1]]))
}

#' @title Geometric Score
#' @name distrib_gradient.GeometricDistrib
#' @description
#' Computes the first derivative of the geometric log-mass with respect to the
#' mean, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1+\mu)}.}
#' The residual is divided by the variance \eqn{\mu(1+\mu)}, and the sum
#' vanishes exactly at \eqn{\hat\mu = \bar y}.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param y A numeric vector of counts.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation and \eqn{\mu > 0} the mean,
#' with variance \eqn{\mu(1+\mu)}. The success probability is
#' \eqn{p = 1/(1+\mu)}.
#'
#' @seealso [distrib_hessian.GeometricDistrib()] for the second derivative,
#'   [distrib_expected_hessian.GeometricDistrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # The closed form, written out: the residual over the variance.
#' all.equal(distrib_gradient(d, y, th)$mu, (y - 3) / (3 * (1 + 3)))
#'
#' # The summed score vanishes at the sample mean, which is the estimate.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, list(mu = 2.5))
#' sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
S7::method(distrib_gradient, GeometricDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  geometric_gradient_cpp(y, theta[[1]], threads)
}

#' @title Geometric Observed Hessian
#' @name distrib_hessian.GeometricDistrib
#' @description
#' Computes the second derivative of the geometric log-mass with respect to the
#' mean, one value per observation, in closed form. Differentiating
#' \eqn{(y-\mu)/(\mu(1+\mu))} gives
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2}
#'       = \dfrac{-\mu(1+\mu) - (y-\mu)(1+2\mu)}{\mu^2(1+\mu)^2},}
#' which depends on the data through \eqn{y} alone. Its expectation is
#' [distrib_expected_hessian.GeometricDistrib()].
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param y A numeric vector of counts.
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
#' @seealso [distrib_gradient.GeometricDistrib()] for the score,
#'   [distrib_expected_hessian.GeometricDistrib()] for the expectation of this
#'   quantity, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # A central difference of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 3 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 3 - eps))$mu
#' all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
#'           tolerance = 1e-6)
#'
#' # It varies with the count, unlike the expected value.
#' rbind(observed = distrib_hessian(d, y, th)$mu_mu,
#'       expected = distrib_expected_hessian(d, y, th)$mu_mu)
S7::method(distrib_hessian, GeometricDistrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ..., threads = 1L) {
  geometric_hessian_cpp(y, theta[[1]], threads)
}

#' @title Geometric Expected Hessian
#' @name distrib_expected_hessian.GeometricDistrib
#' @description
#' Returns the expectation of the observed second derivative under the model,
#' in closed form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu(1+\mu)},}
#' which follows from \eqn{\mathbb{E}[Y] = \mu} killing the data term. The
#' Fisher information is \eqn{1/(\mu(1+\mu))}, the reciprocal of the variance,
#' so it is smaller than a Poisson's \eqn{1/\mu} at the same mean: the extra
#' dispersion is paid for in precision.
#'
#' Because the value does not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param y A numeric vector of counts. Only its length is used.
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
#'   `max(length(y), length(mu))` and constant at \eqn{-1/(\mu(1+\mu))}.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' **observed information** under the model. The geometric is a regular family,
#' so the second Bartlett identity holds and this equals the variance of the
#' score.
#'
#' @seealso [distrib_hessian.GeometricDistrib()] for the observed quantity this
#'   is the expectation of, [distrib_expected_hessian.PoissonDistrib()] for the
#'   equidispersed comparison, [fisher_scoring()], which inverts it, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' th <- list(mu = 3)
#'
#' # A single number, the reciprocal of the variance.
#' unique(distrib_expected_hessian(d, c(0, 2, 7), th)$mu_mu)
#' -1 / (3 * (1 + 3))
#'
#' # Less information than a Poisson of the same mean, the variance being
#' # larger.
#' c(geometric = -distrib_expected_hessian(d, 0, th)$mu_mu,
#'   poisson = -distrib_expected_hessian(poisson_distrib(), 0, th)$mu_mu)
#'
#' # The observed value averages onto it over a large sample.
#' set.seed(5)
#' z <- distrib_rng(d, 2e5, th)
#' c(observed = mean(distrib_hessian(d, z, th)$mu_mu),
#'   expected = distrib_expected_hessian(d, 0, th)$mu_mu)
S7::method(distrib_expected_hessian, GeometricDistrib) <- function(distrib, y, theta,
                                                                   scale = c("parameter", "link"),
                                                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                                                   nsim = 10000, ..., threads = 1L) {
  geometric_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Geometric Third-Order Derivative
#' @name distrib_deriv3.GeometricDistrib
#' @description
#' Computes the third derivative of the geometric log-mass with respect to the
#' mean, in closed form. The log-mass is
#' \eqn{y\log\mu - (y+1)\log(1+\mu)}, so every derivative is a combination of
#' \eqn{\mu^{-k}} and \eqn{(1+\mu)^{-k}} with the data entering linearly
#' through \eqn{y}.
#'
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}. Both routes are closed form, so no quadrature is
#' run and `approx` and `nsim` are ignored.
#'
#' The family has one parameter, so there is one component, where a
#' two-parameter family carries four at this order.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the observed value. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available for both routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu_mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell^{(\mu\mu\mu)}} is the third derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_hessian.GeometricDistrib()] for the order below and
#'   [distrib_deriv4.GeometricDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # A central difference of the Hessian reproduces the observed value.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 3 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 3 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
#'           tolerance = 1e-5)
#'
#' # The expected value is one number, the data term having been averaged out.
#' unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
S7::method(distrib_deriv3, GeometricDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ..., threads = 1L) {
  if (expected) geometric_deriv3_expected_cpp(y, theta[[1]], threads)
  else geometric_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Geometric Fourth-Order Derivative
#' @name distrib_deriv4.GeometricDistrib
#' @description
#' Computes the fourth derivative of the geometric log-mass with respect to the
#' mean, in closed form. As at third order it is a combination of
#' \eqn{\mu^{-k}} and \eqn{(1+\mu)^{-k}} with the data entering linearly
#' through \eqn{y}.
#'
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}. Both routes are closed form, so `approx` and
#' `nsim` are ignored.
#'
#' @param distrib A `GeometricDistrib` object, from [geometric_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the observed value. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available for both routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu_mu_mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @section Notation:
#' \eqn{\ell^{(\mu\mu\mu\mu)}} is the fourth derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.GeometricDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- geometric_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 3 + eps))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 3 - eps))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
#'
#' # The expected value is one number.
#' unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
S7::method(distrib_deriv4, GeometricDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ..., threads = 1L) {
  if (expected) geometric_deriv4_expected_cpp(y, theta[[1]], threads)
  else geometric_deriv4_cpp(y, theta[[1]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Geometric Distribution
#'
#' @description
#' Builds the distribution object for the geometric family parametrized by its
#' **mean** \eqn{\mu > 0}: the number of failures before the first success in
#' independent trials of success probability \eqn{p = 1/(1+\mu)}. The returned
#' object carries closed-form derivatives of the log-mass to fourth order,
#' observed and expected, and closed-form moments.
#'
#' The variance is \eqn{\mu(1+\mu)}, always above the mean, so a geometric fit
#' is an overdispersed alternative to a Poisson with no extra parameter to
#' estimate. It is the negative binomial at a dispersion of exactly one, and it
#' is the only discrete law that is memoryless.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The mass on \eqn{y \in \{0, 1, 2, \dots\}} is
#' \deqn{P(Y = y; \mu) = \dfrac{1}{1+\mu}\left(\dfrac{\mu}{1+\mu}\right)^{y},}
#' with \eqn{\mu \in (0, \infty)}. The family is parametrized by the mean and
#' not by the success probability, so that a link can carry it to the
#' unconstrained scale; [geom_prob()] converts, giving \eqn{p = 1/(1+\mu)} for
#' the base R functions.
#'
#' The mean is \eqn{\mu}, the variance \eqn{\mu(1+\mu)}, and the mode is always
#' 0: consecutive masses differ by the constant factor \eqn{\mu/(1+\mu)},
#' whatever the mean is. The distribution is right skewed at every parameter
#' value.
#'
#' # Overdispersion and memorylessness
#'
#' The variance exceeds the mean by \eqn{\mu^2}, so this family sits above a
#' Poisson of the same mean and is a one-parameter answer to overdispersion.
#' The price is precision: the information is \eqn{1/(\mu(1+\mu))} against a
#' Poisson's \eqn{1/\mu}. Where the amount of overdispersion should itself be
#' estimated, [negbin2_distrib()] adds the parameter and contains this family
#' at \eqn{\theta = 1}.
#'
#' The survival function is \eqn{P(Y > q) = (\mu/(1+\mu))^{q+1}}, so
#' \eqn{P(Y > s + t \mid Y > s) = P(Y > t)}: the count already accumulated
#' tells nothing about the count remaining. The geometric is the only discrete
#' law with this property, as the exponential is the only continuous one, and
#' the two are counterparts.
#'
#' # Estimation
#'
#' The maximum likelihood estimate is the sample mean, \eqn{\hat\mu = \bar y},
#' in closed form; the score is the residual divided by the variance and its
#' sum vanishes there.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation and \eqn{\mu > 0} the mean,
#' with variance \eqn{\mu(1+\mu)}. The success probability is
#' \eqn{p = 1/(1+\mu)}. Counting starts at zero: \eqn{Y} is the number of
#' failures **before** the first success, not the number of trials.
#'
#' @return An S7 object of class `GeometricDistrib`, inheriting from
#'   `discrete_distrib`, with `distrib_name` `"geometric"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`,
#'   `params_interpretation` `c(mu = "mean")`, `n_params` `1`, `params_bounds`
#'   the list of \eqn{(0, \infty)}, and `link_params` the one link given here.
#'
#' @seealso
#' [negbin2_distrib()], of which this is the case \eqn{\theta = 1};
#' [poisson_distrib()] for the equidispersed alternative;
#' [exponential_distrib()], the memoryless continuous counterpart;
#' [geom_prob()] for the conversion to the success probability;
#' [fit_distrib()] to estimate the mean; [GeometricDistrib] for the class.
#'
#' @references
#' Johnson, N. L., Kemp, A. W. and Kotz, S. (2005).
#' *Univariate Discrete Distributions*, 3rd edition, Chapter 5.
#' Wiley, Hoboken.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dgeom pgeom qgeom rgeom
#'
#' @examples
#' d <- geometric_distrib()
#' d
#'
#' # The mass is R's own at prob = 1/(1+mu), and sums to one.
#' all.equal(distrib_pdf(d, c(0, 2, 7), list(mu = 3)),
#'           dgeom(c(0, 2, 7), prob = 1 / 4))
#' sum(distrib_pdf(d, 0:400, list(mu = 3)))
#'
#' # Overdispersed by construction: the variance is mu(1+mu).
#' vapply(c(0.5, 3, 20), function(m) {
#'   th <- list(mu = m)
#'   c(mean = mean(d, th), var = variance(d, th), skew = skewness(d, th))
#' }, numeric(3))
#'
#' # Memoryless, as the exponential is on the continuous side.
#' vapply(c(0, 2, 10, 50), function(s)
#'   distrib_cdf(d, s + 1, list(mu = 3), lower.tail = FALSE) /
#'     distrib_cdf(d, s, list(mu = 3), lower.tail = FALSE), numeric(1))
#'
#' # It is a negative binomial at a dispersion of one.
#' all.equal(distrib_pdf(d, c(0, 2, 7), list(mu = 3)),
#'           distrib_pdf(fixed(negbin2_distrib(), theta = 1), c(0, 2, 7),
#'                       list(mu = 3)))
#'
#' # The estimate is the sample mean, in closed form.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, list(mu = 2.5))
#' c(fitted = unname(coef(fit_distrib(d, z))), sample_mean = mean(z))
#'
#' @export
geometric_distrib <- function(link_mu = log_link()) {
  GeometricDistrib(
    distrib_name = "geometric", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = link_mu)
  )
}
