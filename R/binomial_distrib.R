#' @include distrib.R generics.R
NULL

#' @title Binomial Distribution Class
#' @name BinomialDistrib
#'
#' @description
#' The S7 class of the binomial family: the number of successes in `size`
#' independent trials, each succeeding with probability \eqn{\mu \in (0, 1)}.
#' It inherits from `discrete_distrib`, so its support is counted rather than
#' integrated and no derivative with respect to the response is defined.
#'
#' The class adds one property to the parent's, `size`, and that is what makes
#' it worth documenting separately from [bernoulli_distrib()]: **`size` is a
#' constant of the object, not a parameter**. It is not estimated, it carries
#' no link and no bound, and it does not appear in `params`. It may be a single
#' number or one value per observation, which is how grouped binary data with
#' unequal group sizes is described.
#'
#' The default link is the logit, the **canonical** link here: on its scale the
#' observed and the expected information coincide.
#'
#' Build one with [binomial_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties plus `size` and validates
#' none of the relationships between them.
#'
#' @inheritParams distrib
#' @param size A numeric vector, the number of trials. A single value applies
#'   to every observation; a vector of the length of the response gives one
#'   count of trials per observation. Not a parameter: it is fixed data, and
#'   the object stores it in the `size` property.
#'
#' @return An S7 object of class `BinomialDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's —
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params`, `params_smooth` — plus
#'   `size`. For an object built by [binomial_distrib()] they hold
#'   `"binomial"`, `"univariate"`, `c(0, max(size))`, `"mu"`,
#'   `c(mu = "probability")`, `1`, the domain \eqn{(0, 1)}, the one link, and
#'   the `size` given.
#'
#' @seealso [binomial_distrib()] to build one;
#'   [bernoulli_distrib()], the case `size = 1`;
#'   [betabinom1_distrib()] and [betabinom2_distrib()] when the probability
#'   varies between groups; [poisson_distrib()], the limit of many trials at a
#'   small probability.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.BinomialDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.BinomialDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.BinomialDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BinomialDistrib],
#'   [`distrib_gradient()`][distrib_gradient.BinomialDistrib],
#'   [`distrib_hessian()`][distrib_hessian.BinomialDistrib],
#'   [`distrib_pdf()`][distrib_pdf.BinomialDistrib],
#'   [`distrib_quantile()`][distrib_quantile.BinomialDistrib],
#'   [`distrib_rng()`][distrib_rng.BinomialDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # size is a property of the object; only mu is a parameter.
#' d@size
#' d@params
#' d@bounds
#'
#' # The mean is n p and the variance n p (1-p).
#' th <- list(mu = 0.3)
#' c(mean = mean(d, th), var = variance(d, th), n_p = 10 * 0.3)
#'
#' # size = 1 is the Bernoulli, exactly.
#' all.equal(distrib_pdf(binomial_distrib(size = 1), c(0, 1), th),
#'           distrib_pdf(bernoulli_distrib(), c(0, 1), th))
BinomialDistrib <- S7::new_class("BinomialDistrib", 
  parent = discrete_distrib,
  properties = list(
    size = S7::class_numeric
  )
)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Binomial Probability Mass Function
#' @name distrib_pdf.BinomialDistrib
#' @description
#' Computes the binomial probability mass
#' \deqn{P(Y = y; \mu) = \binom{n}{y}\mu^{y}(1-\mu)^{n-y}, \qquad y = 0, 1, \dots, n,}
#' by calling [stats::dbinom()] at `size = distrib@size`. The number of trials
#' comes from the object rather than from `theta`, `size` being fixed data and
#' not a parameter.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes, integers between 0 and
#'   `size`. Any other value gives 0 with a warning from [stats::dbinom()].
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. A value of length 1 is recycled.
#'   `mu` must lie in \eqn{(0, 1)}.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned, which
#'   stays finite where the mass underflows. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(y), length(mu), length(distrib@size))`.
#'
#' @seealso [distrib_pdf.BernoulliDistrib()] for one trial,
#'   [distrib_cdf.BinomialDistrib()] for the distribution function, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#'
#' # The method is stats::dbinom at the object's own size.
#' all.equal(distrib_pdf(d, c(0, 4, 10), list(mu = 0.3)),
#'           dbinom(c(0, 4, 10), size = 10, prob = 0.3))
#'
#' # A probability mass: it sums to one over 0:n.
#' sum(distrib_pdf(d, 0:10, list(mu = 0.3)))
#'
#' # size may be one value per observation, which is what grouped binary data
#' # with unequal group sizes needs.
#' g <- binomial_distrib(size = c(5, 10, 20))
#' distrib_pdf(g, c(1, 4, 8), list(mu = 0.3))
S7::method(distrib_pdf, BinomialDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbinom(
    x = y,
    size = distrib@size,
    prob = theta[[1]],
    log = log
  )
}

#' @title Binomial Cumulative Distribution Function
#' @name distrib_cdf.BinomialDistrib
#' @description
#' Computes the binomial distribution function
#' \deqn{F(q; \mu) = \sum_{k = 0}^{\lfloor q \rfloor} \binom{n}{k}\mu^{k}(1-\mu)^{n-k}}
#' by calling [stats::pbinom()] at `size = distrib@size`, which evaluates it
#' through the incomplete beta function rather than by summing. The function is
#' a step function, constant between integers, and reaches 1 at `size`.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param q A numeric vector of quantiles. A non-integer is floored, a negative
#'   value gives 0, and a value at or above `size` gives 1.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `q`. `mu` must lie in \eqn{(0, 1)}.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   computed directly and so exact in the upper tail.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(distrib@size))`.
#'
#' @seealso [distrib_quantile.BinomialDistrib()] for the generalized inverse,
#'   [distrib_pdf.BinomialDistrib()] for the mass, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' th <- list(mu = 0.3)
#'
#' # The method is stats::pbinom at the object's own size.
#' all.equal(distrib_cdf(d, c(0, 4, 10), th),
#'           pbinom(c(0, 4, 10), size = 10, prob = 0.3))
#'
#' # A step function reaching one at size.
#' distrib_cdf(d, c(3, 3.5, 3.9, 10, 20), th)
#'
#' # The jump at an integer is exactly the mass there.
#' c(jump = distrib_cdf(d, 4, th) - distrib_cdf(d, 3, th),
#'   mass = distrib_pdf(d, 4, th))
S7::method(distrib_cdf, BinomialDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbinom(
    q = q,
    size = distrib@size,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Binomial Quantile Function
#' @name distrib_quantile.BinomialDistrib
#' @description
#' Computes the generalized inverse of the binomial distribution function,
#' \deqn{Q(p; \mu) = \min\{k \in \{0, \dots, n\} : F(k; \mu) \ge p\},}
#' by calling [stats::qbinom()] at `size = distrib@size`. The distribution
#' function is a step function, so the round trip through
#' [distrib_cdf.BinomialDistrib()] returns a probability **at least** `p`.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `p`. `mu` must lie in \eqn{(0, 1)}.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#'
#' @return A numeric vector of integers between 0 and `size`, of length
#'   `max(length(p), length(mu), length(distrib@size))`.
#'
#' @seealso [distrib_cdf.BinomialDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' th <- list(mu = 0.3)
#'
#' # Integers, bounded by size.
#' distrib_quantile(d, c(0.025, 0.5, 0.975, 1), th)
#'
#' # The round trip overshoots, the support being a lattice.
#' p <- c(0.025, 0.5, 0.975)
#' rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
S7::method(distrib_quantile, BinomialDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbinom(
    p = p,
    size = distrib@size,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Binomial Random Number Generator
#' @name distrib_rng.BinomialDistrib
#' @description
#' Draws `n` independent counts of successes by calling [stats::rbinom()] at
#' `size = distrib@size`, so the draws come from R's own generator and depend
#' on `.Random.seed` in the usual way.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials, and may be one value per
#'   draw.
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of length `n`. A value of length 1 is recycled. `mu` must
#'   lie in \eqn{(0, 1)}.
#'
#' @return An integer vector of `n` counts, each between 0 and the
#'   corresponding `size`.
#'
#' @seealso [fit_distrib()] to estimate the probability back, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#'
#' # Same generator as stats::rbinom, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 5, list(mu = 0.3))
#' set.seed(2)
#' identical(a, rbinom(5, size = 10, prob = 0.3))
#'
#' # The sample mean estimates n p and the variance n p (1-p).
#' set.seed(4)
#' z <- distrib_rng(d, 2e4, list(mu = 0.3))
#' c(mean = mean(z), n_p = 10 * 0.3,
#'   var = var(z), n_p_q = 10 * 0.3 * 0.7)
#'
#' # Unequal group sizes, one per draw.
#' set.seed(5)
#' g <- binomial_distrib(size = c(5, 10, 20, 50))
#' rbind(size = c(5, 10, 20, 50), draw = distrib_rng(g, 4, list(mu = 0.3)))
S7::method(distrib_rng, BinomialDistrib) <- function(distrib, n, theta) {
  stats::rbinom(
    n = n,
    size = distrib@size,
    prob = theta[[1]]
  )
}

#' @title Binomial Score
#' @name distrib_gradient.BinomialDistrib
#' @description
#' Computes the first derivative of the binomial log-mass with respect to the
#' probability, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - \dfrac{n - y}{1-\mu}
#'       = \dfrac{y - n\mu}{\mu(1-\mu)}.}
#' The residual is measured against \eqn{n\mu}, the expected count, and divided
#' by \eqn{\mu(1-\mu)}. Its sum vanishes at \eqn{\hat\mu = \sum y_i / \sum n_i},
#' the pooled proportion.
#'
#' On the **link** scale with the default logit the generic's chain rule gives
#' \eqn{\partial\ell/\partial\eta = y - n\mu}, the raw residual: the logit is
#' the canonical link of this family.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{n} the number of trials
#' and \eqn{\mu \in (0,1)} the success probability, with mean \eqn{n\mu} and
#' variance \eqn{n\mu(1-\mu)}. \eqn{\eta = \log(\mu/(1-\mu))} is the log odds.
#'
#' @seealso [distrib_hessian.BinomialDistrib()] for the second derivative,
#'   [distrib_expected_hessian.BinomialDistrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' y <- c(0, 4, 10)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out.
#' all.equal(distrib_gradient(d, y, th)$mu, (y - 10 * 0.3) / (0.3 * 0.7))
#'
#' # On the canonical logit link the score is the raw residual.
#' distrib_gradient(d, y, th, scale = "link")$mu
#' y - 10 * 0.3
#'
#' # The summed score vanishes at the pooled proportion.
#' set.seed(4)
#' z <- distrib_rng(d, 500, list(mu = 0.3))
#' sum(distrib_gradient(d, z, list(mu = sum(z) / (500 * 10)))$mu)
S7::method(distrib_gradient, BinomialDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  binomial_gradient_cpp(y, theta[[1]], distrib@size, threads)
}

#' @title Binomial Observed Hessian
#' @name distrib_hessian.BinomialDistrib
#' @description
#' Computes the second derivative of the binomial log-mass with respect to the
#' probability, one value per observation, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{n-y}{(1-\mu)^2}.}
#' It is negative for every admissible \eqn{\mu} and depends on the data
#' through \eqn{y} alone.
#'
#' On the **link** scale with the default logit the chain rule gives
#' \eqn{\partial^2\ell/\partial\eta^2 = -n\mu(1-\mu)}, which carries **no data
#' at all**: the defining property of a canonical link, and what makes the
#' observed and the expected information the same matrix there. See
#' [distrib_expected_hessian.BinomialDistrib()].
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`.
#'
#' @seealso [distrib_gradient.BinomialDistrib()] for the score,
#'   [distrib_expected_hessian.BinomialDistrib()] for the expectation of this
#'   quantity, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' y <- c(0, 4, 10)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out.
#' all.equal(distrib_hessian(d, y, th)$mu_mu,
#'           -y / 0.3^2 - (10 - y) / 0.7^2)
#'
#' # On the canonical logit link the curvature carries no data: one value,
#' # repeated, equal to -n mu (1-mu).
#' distrib_hessian(d, y, th, scale = "link")$mu_mu
#' -10 * 0.3 * 0.7
#'
#' # A central difference of the score reproduces the parameter-scale value.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 0.3 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.3 - eps))$mu
#' all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
#'           tolerance = 1e-6)
S7::method(distrib_hessian, BinomialDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  binomial_hessian_cpp(y, theta[[1]], distrib@size, threads)
}

#' @title Binomial Expected Hessian
#' @name distrib_expected_hessian.BinomialDistrib
#' @description
#' Returns the expectation of the observed second derivative under the model,
#' in closed form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{n}{\mu(1-\mu)},}
#' which follows from \eqn{\mathbb{E}[Y] = n\mu} in
#' \eqn{-y/\mu^2 - (n-y)/(1-\mu)^2}. The Fisher information is
#' \eqn{n/(\mu(1-\mu))}, so it grows in proportion to the number of trials and
#' is smallest at \eqn{\mu = 1/2}.
#'
#' On the **link** scale with the default logit the value is
#' \eqn{-n\mu(1-\mu)}, and so is the observed Hessian: the logit is the
#' canonical link and the two coincide exactly, which is what makes iteratively
#' reweighted least squares Fisher scoring and Newton's method at once.
#'
#' Because the value does not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
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
#'   `max(length(y), length(mu), length(distrib@size))`, constant at
#'   \eqn{-n/(\mu(1-\mu))} where `size` is.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' **observed information** under the model. The binomial is a regular family,
#' so the second Bartlett identity holds and this equals the variance of the
#' score.
#'
#' @seealso [distrib_hessian.BinomialDistrib()] for the observed quantity this
#'   is the expectation of, [fisher_scoring()], which inverts it, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' th <- list(mu = 0.3)
#'
#' # A single number, -n divided by the Bernoulli variance.
#' unique(distrib_expected_hessian(d, c(0, 4, 10), th)$mu_mu)
#' -10 / (0.3 * 0.7)
#'
#' # It grows in proportion to the number of trials.
#' vapply(c(1, 10, 100), function(n)
#'   -distrib_expected_hessian(binomial_distrib(size = n), 0, th)$mu_mu,
#'   numeric(1))
#'
#' # On the canonical logit link the observed and expected values agree
#' # exactly, at every observation.
#' rbind(observed = distrib_hessian(d, c(0, 4, 10), th, scale = "link")$mu_mu,
#'       expected = distrib_expected_hessian(d, c(0, 4, 10), th,
#'                                           scale = "link")$mu_mu)
S7::method(distrib_expected_hessian, BinomialDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  binomial_expected_hessian_cpp(y, theta[[1]], distrib@size, threads)
}

#' @title Binomial Third-Order Derivative
#' @name distrib_deriv3.BinomialDistrib
#' @description
#' Computes the third derivative of the binomial log-mass with respect to the
#' probability, in closed form:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{2y}{\mu^3} - \dfrac{2(n-y)}{(1-\mu)^3}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = n\mu}, which gives \eqn{2n/\mu^2 - 2n/(1-\mu)^2}. Both
#' routes are closed form, so no quadrature is run and `approx` and `nsim` are
#' ignored.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes. With `expected = TRUE`
#'   only its length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
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
#'   `max(length(y), length(mu), length(distrib@size))`.
#'
#' @section Notation:
#' \eqn{\ell^{(\mu\mu\mu)}} is the third derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_hessian.BinomialDistrib()] for the order below and
#'   [distrib_deriv4.BinomialDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' y <- c(0, 4, 10)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out.
#' all.equal(distrib_deriv3(d, y, th)$mu_mu_mu,
#'           2 * y / 0.3^3 - 2 * (10 - y) / 0.7^3)
#'
#' # The expected value, and zero at mu = 1/2 by symmetry.
#' unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#' 2 * 10 / 0.3^2 - 2 * 10 / 0.7^2
#' distrib_deriv3(d, 0, list(mu = 0.5), expected = TRUE)$mu_mu_mu
#'
#' # A central difference of the Hessian reproduces the observed value.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 0.3 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.3 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv3, BinomialDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) binomial_deriv3_expected_cpp(y, theta[[1]], distrib@size, threads)
  else binomial_deriv3_cpp(y, theta[[1]], distrib@size, threads)
}

#' @title Binomial Fourth-Order Derivative
#' @name distrib_deriv4.BinomialDistrib
#' @description
#' Computes the fourth derivative of the binomial log-mass with respect to the
#' probability, in closed form:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4} = -\dfrac{6y}{\mu^4} - \dfrac{6(n-y)}{(1-\mu)^4}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = n\mu}, which gives \eqn{-6n/\mu^3 - 6n/(1-\mu)^3}. Both
#' routes are closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes. With `expected = TRUE`
#'   only its length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
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
#'   `max(length(y), length(mu), length(distrib@size))`.
#'
#' @section Notation:
#' \eqn{\ell^{(\mu\mu\mu\mu)}} is the fourth derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.BinomialDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' y <- c(0, 4, 10)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out; negative at every admissible mu.
#' all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           -6 * y / 0.3^4 - 6 * (10 - y) / 0.7^4)
#'
#' # The expected value.
#' unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
#' -6 * 10 / 0.3^3 - 6 * 10 / 0.7^3
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 0.3 + eps))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.3 - eps))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, BinomialDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) binomial_deriv4_expected_cpp(y, theta[[1]], distrib@size, threads)
  else binomial_deriv4_cpp(y, theta[[1]], distrib@size, threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Binomial Distribution
#'
#' @description
#' Builds the distribution object for the binomial family: the number of
#' successes in `size` independent trials, each succeeding with probability
#' \eqn{\mu \in (0, 1)}. The returned object carries closed-form derivatives of
#' the log-mass to fourth order, observed and expected, and closed-form moments.
#'
#' **`size` is fixed data and not a parameter.** It is carried on the object,
#' it is not estimated, it has no link and no bound, and it does not appear in
#' `params`. Giving it one value per observation is how grouped binary data
#' with unequal group sizes is described.
#'
#' The default link is the logit, the canonical link here, so the observed and
#' the expected information coincide on its scale.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the probability
#'   \eqn{\mu}. Defaults to [linkfunctions7::logit_link()], the canonical link,
#'   which maps \eqn{(0, 1)} onto the whole line.
#'   [linkfunctions7::probit_link()] and [linkfunctions7::cloglog_link()] are
#'   the usual alternatives; under either the observed and expected information
#'   differ.
#' @param size A numeric vector of trial counts, positive integers. A single
#'   value applies to every observation; a vector of the length of the response
#'   gives one count per observation. Defaults to `1`, which makes the object a
#'   Bernoulli; pass the group sizes for grouped data. The value is stored in
#'   the `size` property and sets `bounds` to `c(0, max(size))`.
#'
#' @details
#' # The parametrization
#'
#' The mass on \eqn{y \in \{0, 1, \dots, n\}} is
#' \deqn{P(Y = y; \mu) = \binom{n}{y}\mu^{y}(1-\mu)^{n-y},}
#' with \eqn{n} the number of trials and \eqn{\mu \in (0, 1)}. The mean is
#' \eqn{n\mu}, the variance \eqn{n\mu(1-\mu)}, the skewness
#' \eqn{(1-2\mu)/\sqrt{n\mu(1-\mu)}} and the excess kurtosis
#' \eqn{(1-6\mu(1-\mu))/(n\mu(1-\mu))}. Both shape measures shrink like
#' \eqn{1/n}, which is the central limit theorem visible in the moments.
#'
#' As \eqn{n \to \infty} with \eqn{n\mu} held, the family tends to a Poisson of
#' that mean; [poisson_distrib()] is the limit and is the family to use when
#' the number of trials is large and unrecorded.
#'
#' # size is data
#'
#' The number of trials sits on the object rather than in `theta`, so it takes
#' no link, is not counted in `n_params`, and cannot be estimated. Two
#' consequences worth knowing: the same object cannot be reused across data
#' sets with different group sizes, and `bounds` records `c(0, max(size))`, so
#' a vector `size` gives a bound that is correct for the largest group and
#' loose for the rest.
#'
#' # The canonical link
#'
#' The log-mass is linear in \eqn{y} once \eqn{\eta = \log(\mu/(1-\mu))} is the
#' parameter, so the logit is canonical and on its scale
#' \deqn{\dfrac{\partial \ell}{\partial \eta} = y - n\mu, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \eta^2} = -n\mu(1-\mu),}
#' the second carrying no data. Under a probit or a complementary log-log link
#' the observed and expected information differ.
#'
#' # Estimation
#'
#' The maximum likelihood estimate is the pooled proportion,
#' \eqn{\hat\mu = \sum y_i / \sum n_i}, in closed form. If the observed counts
#' are more variable than \eqn{n\mu(1-\mu)} allows, the binomial is the wrong
#' family and [betabinom1_distrib()] adds the dispersion parameter it lacks.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{n} the number of trials
#' and \eqn{\mu \in (0,1)} the success probability.
#' \eqn{\eta = \log(\mu/(1-\mu))} is the log odds. The **canonical** link of an
#' exponential family makes the log-mass linear in the sufficient statistic.
#'
#' @return An S7 object of class `BinomialDistrib`, inheriting from
#'   `discrete_distrib`, with `distrib_name` `"binomial"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, max(size))`, `size` as given, `params`
#'   `"mu"`, `params_interpretation` `c(mu = "probability")`, `n_params` `1`,
#'   `params_bounds` the list of \eqn{(0, 1)}, and `link_params` the one link.
#'
#' @seealso
#' [bernoulli_distrib()], the case `size = 1`;
#' [betabinom1_distrib()] and [betabinom2_distrib()] for overdispersed grouped
#' binary data; [poisson_distrib()], the many-trials limit;
#' [multinomial_distrib()] for more than two categories;
#' [fit_distrib()] to estimate the probability; [BinomialDistrib] for the class.
#'
#' @references
#' Johnson, N. L., Kemp, A. W. and Kotz, S. (2005).
#' *Univariate Discrete Distributions*, 3rd edition, Chapter 3.
#' Wiley, Hoboken.
#'
#' @importFrom linkfunctions7 logit_link
#' @importFrom stats dbinom pbinom qbinom rbinom
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' d
#'
#' # The mass is R's own and sums to one over 0:n.
#' all.equal(distrib_pdf(d, c(0, 4, 10), list(mu = 0.3)),
#'           dbinom(c(0, 4, 10), size = 10, prob = 0.3))
#' sum(distrib_pdf(d, 0:10, list(mu = 0.3)))
#'
#' # size = 1 is the Bernoulli, exactly.
#' all.equal(distrib_pdf(binomial_distrib(size = 1), c(0, 1), list(mu = 0.3)),
#'           distrib_pdf(bernoulli_distrib(), c(0, 1), list(mu = 0.3)))
#'
#' # Skewness shrinks like 1/sqrt(n): the central limit theorem in the moments.
#' vapply(c(1, 10, 100), function(n)
#'   skewness(binomial_distrib(size = n), list(mu = 0.3)), numeric(1))
#'
#' # Many trials at a small probability is a Poisson of the same mean.
#' rbind(binomial = distrib_pdf(binomial_distrib(size = 1000),
#'                              0:5, list(mu = 3 / 1000)),
#'       poisson = distrib_pdf(poisson_distrib(), 0:5, list(mu = 3)))
#'
#' # Grouped data with unequal group sizes: one size per observation.
#' set.seed(5)
#' g <- binomial_distrib(size = c(5, 10, 20, 50))
#' rbind(size = c(5, 10, 20, 50), draw = distrib_rng(g, 4, list(mu = 0.3)))
#'
#' @export
binomial_distrib <- function(link_mu = logit_link(), size = 1) {
  
  BinomialDistrib(
    distrib_name = "binomial",
    dimension = "univariate",
    bounds = c(0, max(size)),
    size = size,
    
    params = c("mu"),
    params_interpretation = c(mu = "probability"),
    n_params = 1,
    
    params_bounds = list(
      mu = c(0, 1)
    ),
    
    link_params = list(
      mu = link_mu
    )
  )
  
}
