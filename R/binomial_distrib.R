#' @include distrib.R generics.R
NULL

#' @title Binomial Distribution Class
#' @name BinomialDistrib
#'
#' @description
#' This is the S7 class of the binomial family: the number of successes in
#' `size` independent trials, where each trial succeeds with probability
#' \eqn{\mu \in (0, 1)}. It inherits from `discrete_distrib`, so expectations
#' over its support are sums and the class defines no derivative with respect to
#' the response.
#'
#' The class adds one property to the parent's, `size`, and that property is
#' the reason that the class is documented separately from [bernoulli_distrib()]:
#' `size` is a constant of the object, not a parameter. It is not estimated,
#' it carries no link and no bound, and it does not appear in `params`. It may
#' be a single number or one value per observation, which describes grouped
#' binary data with unequal group sizes.
#'
#' The default link is the logit, the canonical link here: on its scale the
#' observed and the expected information coincide.
#'
#' Build one with [binomial_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties plus `size`; the
#' validator of `distrib` checks the parent's properties, and `size` is
#' checked by [binomial_distrib()], not by this constructor.
#'
#' @inheritParams distrib
#' @param size A numeric vector, the number of trials. A single value applies
#'   to every observation; a vector of the length of the response gives one
#'   count of trials per observation. It is fixed data and not a parameter,
#'   and the object stores it in the `size` property.
#'
#' @return An S7 object of class `BinomialDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's
#'   (`distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params`, `params_smooth`) plus
#'   `size`. For an object built by [binomial_distrib()] they hold
#'   `"binomial"`, `"univariate"`, `c(0, max(size))`, `"mu"`,
#'   `c(mu = "probability")`, `1`, the domain \eqn{(0, 1)}, the one link, and
#'   the supplied `size`.
#'
#' @seealso [binomial_distrib()] to build one;
#'   [bernoulli_distrib()], the case `size = 1`;
#'   [betabinom1_distrib()] and [betabinom2_distrib()] when the probability
#'   varies between groups; [poisson_distrib()], the limit of many trials at a
#'   small probability.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.BinomialDistrib],
#'   [`distrib_cdf()`][distrib_cdf.BinomialDistrib],
#'   [`distrib_quantile()`][distrib_quantile.BinomialDistrib] and
#'   [`distrib_rng()`][distrib_rng.BinomialDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.BinomialDistrib],
#'   [`distrib_hessian()`][distrib_hessian.BinomialDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BinomialDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.BinomialDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.BinomialDistrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.BinomialDistrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.BinomialDistrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.BinomialDistrib];
#' - the moments: [`mean()`][mean.BinomialDistrib],
#'   [`variance()`][variance.BinomialDistrib],
#'   [`skewness()`][skewness.BinomialDistrib] and
#'   [`kurtosis()`][kurtosis.BinomialDistrib];
#' - the fitting helpers: [distrib_scalar_route()].
#'
#' Every other generic is inherited from [discrete_distrib()].
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
#' comes from the object and not from `theta`, because `size` is fixed data and
#' not a parameter.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes, integers between 0 and
#'   `size`. Any other value gives 0, and a non-integer value also gives a
#'   warning from [stats::dbinom()]. A missing value gives `NA`.
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
#' # size may be one value per observation, as grouped binary data with
#' # unequal group sizes require.
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
#' through the incomplete beta function. The function is a step function,
#' constant between integers, and reaches 1 at `size`.
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
#' @param ... Unused, and accepted so that the signature matches the generic's.
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
#' # The distribution function is a step function that reaches one at size.
#' distrib_cdf(d, c(3, 3.5, 3.9, 10, 20), th)
#'
#' # The jump at an integer is exactly the mass there.
#' c(jump = distrib_cdf(d, 4, th) - distrib_cdf(d, 3, th),
#'   mass = distrib_pdf(d, 4, th))
S7::method(distrib_cdf, BinomialDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
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
#' [distrib_cdf.BinomialDistrib()] returns a probability at least `p`.
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
#' @param ... Unused, and accepted so that the signature matches the generic's.
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
#' # The quantiles are integers bounded by size.
#' distrib_quantile(d, c(0.025, 0.5, 0.975, 1), th)
#'
#' # The round trip overshoots, because the support is discrete.
#' p <- c(0.025, 0.5, 0.975)
#' rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
S7::method(distrib_quantile, BinomialDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
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
#' @param ... Unused, and accepted so that the signature matches the generic's.
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
#' # Unequal group sizes give one size per draw.
#' set.seed(5)
#' g <- binomial_distrib(size = c(5, 10, 20, 50))
#' rbind(size = c(5, 10, 20, 50), draw = distrib_rng(g, 4, list(mu = 0.3)))
S7::method(distrib_rng, BinomialDistrib) <- function(distrib, n, theta, ...) {
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
#' On the link scale with the default logit the generic's chain rule gives
#' \eqn{\partial\ell/\partial\eta = y - n\mu}, the raw residual: the logit is
#' the canonical link of this family.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`. Each of `y`, `mu` and
#'   `size` must have length 1 or that length; any other length signals an
#'   error.
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
#' # The closed form is written out.
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
  a <- binomial_recycle(y, theta[[1]], distrib@size)
  binomial_gradient_cpp(a$y, a$mu, a$size, threads)
}

#' @title Binomial Observed Hessian
#' @name distrib_hessian.BinomialDistrib
#' @description
#' Computes the second derivative of the binomial log-mass with respect to the
#' probability, one value per observation, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{n-y}{(1-\mu)^2}.}
#' It is negative for every admissible \eqn{\mu} and depends on the data
#' only through \eqn{y}. Each term is formed by dividing its coefficient by
#' the base once per power, and a term with a zero coefficient contributes
#' zero, so at \eqn{y = 0} the value is \eqn{-n/(1-\mu)^2} however small
#' \eqn{\mu} is.
#'
#' On the link scale with the default logit the chain rule gives
#' \eqn{\partial^2\ell/\partial\eta^2 = -n\mu(1-\mu)}, which carries no
#' data; this is the defining property of a canonical link, and in exact
#' arithmetic it makes the observed and the expected information the same
#' matrix there (see [distrib_expected_hessian.BinomialDistrib()]). The
#' generic forms the link-scale value from the derivatives on the parameter
#' scale, and for an observation on the far side of the probability
#' (\eqn{y = 0} with \eqn{\mu} near 1, or \eqn{y = n} with \eqn{\mu} near 0)
#' its two terms cancel to a value much smaller than either. The relative
#' error of the observed value then grows with \eqn{|\eta|} (about
#' \eqn{10^{-3}} at \eqn{\eta = 27.6}), while the expected value is computed
#' directly and stays exact.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()]. Its
#'   `size` property supplies the number of trials.
#' @param y A numeric vector of counts of successes.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`. Each of `y`, `mu` and
#'   `size` must have length 1 or that length; any other length signals an
#'   error.
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
#' # The closed form is written out.
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
  a <- binomial_recycle(y, theta[[1]], distrib@size)
  binomial_hessian_cpp(a$y, a$mu, a$size, threads)
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
#' On the link scale with the default logit the value is
#' \eqn{-n\mu(1-\mu)}, and so is the observed Hessian in exact arithmetic: the
#' logit is the canonical link, so iteratively reweighted least squares is
#' Fisher scoring and Newton's method at once. The generic composes the
#' observed value, which loses relative accuracy where \eqn{|\eta|} is large
#' (see [distrib_hessian.BinomialDistrib()]); the value returned here does not.
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
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param approx Ignored, because a closed form is available. It is accepted so
#'   that the signature matches the generic's.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`, equal to
#'   \eqn{-n/(\mu(1-\mu))} with the `size` of each observation. Each of `y`,
#'   `mu` and `size` must have length 1 or that length; any other length
#'   signals an error.
#'
#' @section Notation:
#' The expected information is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' observed information under the model. The binomial is a regular family,
#' so the second Bartlett identity holds and this equals the variance of the
#' score.
#'
#' @seealso [distrib_hessian.BinomialDistrib()] for the observed quantity whose
#'   expectation this is, [fisher_scoring()], which inverts it, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' th <- list(mu = 0.3)
#'
#' # The result is a single number, -n divided by the Bernoulli variance.
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
S7::method(distrib_expected_hessian, BinomialDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ..., threads = 1L) {
  a <- binomial_recycle(y, theta[[1]], distrib@size)
  binomial_expected_hessian_cpp(a$y, a$mu, a$size, threads)
}

#' @title Binomial Third-Order Derivative
#' @name distrib_deriv3.BinomialDistrib
#' @description
#' Computes the third derivative of the binomial log-mass with respect to the
#' probability, in closed form:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{2y}{\mu^3} - \dfrac{2(n-y)}{(1-\mu)^3}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = n\mu}, which gives \eqn{2n/\mu^2 - 2n/(1-\mu)^2}. Both
#' routes are closed form, so the method runs no quadrature and ignores
#' `approx` and `nsim`.
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
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param approx Ignored, because a closed form is available for both routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu_mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`. Each of `y`, `mu` and
#'   `size` must have length 1 or that length; any other length signals an
#'   error.
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
#' # The closed form is written out.
#' all.equal(distrib_deriv3(d, y, th)$mu_mu_mu,
#'           2 * y / 0.3^3 - 2 * (10 - y) / 0.7^3)
#'
#' # The expected value is checked against its closed form, and it is zero at
#' # mu = 1/2 by symmetry.
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
  a <- binomial_recycle(y, theta[[1]], distrib@size)
  if (expected) binomial_deriv3_expected_cpp(a$y, a$mu, a$size, threads)
  else binomial_deriv3_cpp(a$y, a$mu, a$size, threads)
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
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param approx Ignored, because a closed form is available for both routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu_mu_mu`, of length
#'   `max(length(y), length(mu), length(distrib@size))`. Each of `y`, `mu` and
#'   `size` must have length 1 or that length; any other length signals an
#'   error.
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
#' # The closed form is written out, and it is negative at every admissible mu.
#' all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           -6 * y / 0.3^4 - 6 * (10 - y) / 0.7^4)
#'
#' # The expected value is checked against its closed form.
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
  a <- binomial_recycle(y, theta[[1]], distrib@size)
  if (expected) binomial_deriv4_expected_cpp(a$y, a$mu, a$size, threads)
  else binomial_deriv4_cpp(a$y, a$mu, a$size, threads)
}

#' Recycle the Arguments of the Binomial Kernels
#'
#' @description
#' Recycles the response, the success probability and the number of trials
#' to their common length before they reach the compiled kernels, which read
#' each argument either as a scalar or at the observation's index. Each
#' argument must have length 1 or the common length.
#'
#' @param y The response, a numeric vector.
#' @param mu The success probability, a numeric vector.
#' @param size The number of trials, a numeric vector.
#'
#' @return A list with components `y`, `mu` and `size`, each of the common
#'   length.
#'
#' @keywords internal
binomial_recycle <- function(y, mu, size) {
  n <- max(length(y), length(mu), length(size))
  lens <- c(y = length(y), mu = length(mu), size = length(size))
  bad <- !(lens %in% c(1L, n))
  if (any(bad)) {
    stop(sprintf("'%s' has length %d; the response, 'mu' and 'size' must have length 1 or %d.",
                 names(lens)[bad][1], lens[bad][1], n), call. = FALSE)
  }
  list(y = rep_len(y, n), mu = rep_len(mu, n), size = rep_len(size, n))
}

# --- CONSTRUCTOR WRAPPER ---

#' Binomial Distribution
#'
#' @description
#' Builds the distribution object for the binomial family: the number of
#' successes in `size` independent trials, where each trial succeeds with
#' probability \eqn{\mu \in (0, 1)}. The returned object carries closed-form
#' derivatives of the log-mass to fourth order, observed and expected, and
#' closed-form moments.
#'
#' `size` is fixed data, not a parameter. It is carried on the object, it
#' is not estimated, it has no link and no bound, and it does not appear in
#' `params`. It may be given one value per observation, which describes
#' grouped binary data with unequal group sizes.
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
#'   the `size` property and sets `bounds` to `c(0, max(size))`. A value that
#'   is missing, infinite, not an integer or below 1 signals an error. The
#'   derivative methods signal an error when the length of `size` is neither
#'   1 nor the common length of the response and `mu`; the mass, the
#'   distribution function, the quantile and the draws recycle it as
#'   [stats::dbinom()] and its companions do.
#'
#' @details
#' # The parametrization
#'
#' The mass on \eqn{y \in \{0, 1, \dots, n\}} is
#' \deqn{P(Y = y; \mu) = \binom{n}{y}\mu^{y}(1-\mu)^{n-y},}
#' with \eqn{n} the number of trials and \eqn{\mu \in (0, 1)}. The mean is
#' \eqn{n\mu}, the variance \eqn{n\mu(1-\mu)}, the skewness
#' \eqn{(1-2\mu)/\sqrt{n\mu(1-\mu)}} and the excess kurtosis
#' \eqn{(1-6\mu(1-\mu))/(n\mu(1-\mu))}. The skewness shrinks like
#' \eqn{1/\sqrt{n}} and the excess kurtosis like \eqn{1/n}, as the central
#' limit theorem requires.
#'
#' As \eqn{n \to \infty} while \eqn{n\mu} stays fixed, the family tends to a
#' Poisson of that mean; [poisson_distrib()] is the limit and is the family to use when
#' the number of trials is large and unrecorded.
#'
#' # size is data
#'
#' The number of trials sits on the object and not in `theta`, so it takes no
#' link, is not counted in `n_params`, and cannot be estimated. The same
#' object therefore cannot be reused across data sets with different group
#' sizes. The property `bounds` records `c(0, max(size))`, which for a vector
#' `size` is the bound of the largest group and is loose for the others.
#'
#' # The canonical link
#'
#' The log-mass is linear in \eqn{y} once \eqn{\eta = \log(\mu/(1-\mu))} is the
#' parameter, so the logit is canonical and on its scale
#' \deqn{\dfrac{\partial \ell}{\partial \eta} = y - n\mu, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \eta^2} = -n\mu(1-\mu),}
#' where the second derivative does not depend on the data. Under a probit or a complementary log-log link
#' the observed and expected information differ.
#'
#' # Estimation
#'
#' The maximum likelihood estimate is the pooled proportion,
#' \eqn{\hat\mu = \sum y_i / \sum n_i}, in closed form. If the observed counts
#' are more variable than \eqn{n\mu(1-\mu)} allows, the binomial is the wrong
#' family and [betabinom1_distrib()] adds the dispersion parameter that the
#' binomial lacks.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{n} the number of trials
#' and \eqn{\mu \in (0,1)} the success probability.
#' \eqn{\eta = \log(\mu/(1-\mu))} is the log odds. The canonical link of an
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
#' # Grouped data with unequal group sizes take one size per observation.
#' set.seed(5)
#' g <- binomial_distrib(size = c(5, 10, 20, 50))
#' rbind(size = c(5, 10, 20, 50), draw = distrib_rng(g, 4, list(mu = 0.3)))
#'
#' @export
binomial_distrib <- function(link_mu = logit_link(), size = 1) {
  if (!is.numeric(size) || length(size) == 0L || anyNA(size) ||
      any(!is.finite(size)) || any(size < 1) || any(size != floor(size))) {
    stop("'size' must be a vector of positive integers.", call. = FALSE)
  }
  
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
