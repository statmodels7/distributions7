#' @include distrib.R generics.R
NULL

#' @title Bernoulli Distribution Class
#' @name BernoulliDistrib
#'
#' @description
#' The S7 class of the Bernoulli family parametrized by its success probability
#' \eqn{\mu \in (0, 1)}, with mass \eqn{P(Y = 1) = \mu} and
#' \eqn{P(Y = 0) = 1 - \mu}. It inherits from `discrete_distrib`, so its
#' support is counted rather than integrated and no derivative with respect to
#' the response is defined.
#'
#' The mean is \eqn{\mu} and the variance \eqn{\mu(1-\mu)}, so the two are tied
#' and the family has no dispersion parameter of its own. The default link is
#' the logit, which is the **canonical** link here: on its scale the observed
#' and the expected information coincide.
#'
#' A Bernoulli is a binomial on one trial. [binomial_distrib()] with
#' `size = 1` gives the same law, and the two are worth keeping apart only
#' because the binomial carries a `size` property this class does not.
#'
#' Build one with [bernoulli_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `BernoulliDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [bernoulli_distrib()] they hold `"bernoulli"`,
#'   `"univariate"`, `c(0, 1)`, `"mu"`, `c(mu = "probability")`, `1`, the
#'   domain \eqn{(0, 1)}, and the one link.
#'
#' @seealso [bernoulli_distrib()] to build one;
#'   [binomial_distrib()], of which this is the one-trial case;
#'   [betabinom1_distrib()] for a binomial whose probability varies;
#'   [linkfunctions7::logit_link()], the canonical link.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.BernoulliDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.BernoulliDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.BernoulliDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BernoulliDistrib],
#'   [`distrib_gradient()`][distrib_gradient.BernoulliDistrib],
#'   [`distrib_hessian()`][distrib_hessian.BernoulliDistrib],
#'   [`distrib_pdf()`][distrib_pdf.BernoulliDistrib],
#'   [`distrib_quantile()`][distrib_quantile.BernoulliDistrib],
#'   [`distrib_rng()`][distrib_rng.BernoulliDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
#'
#' @examples
#' d <- bernoulli_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#' d@params_interpretation
#' d@params_bounds
#'
#' # The default link is the logit, the canonical one.
#' d@link_params$mu@link_name
#'
#' # Mean and variance are tied: p and p(1-p), maximal at p = 1/2.
#' vapply(c(0.1, 0.5, 0.9), function(p) {
#'   th <- list(mu = p)
#'   c(mean = mean(d, th), var = variance(d, th))
#' }, numeric(2))
#'
#' # The two masses sum to one.
#' sum(distrib_pdf(d, c(0, 1), list(mu = 0.3)))
BernoulliDistrib <- S7::new_class("BernoulliDistrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Bernoulli Probability Mass Function
#' @name distrib_pdf.BernoulliDistrib
#' @description
#' Computes the Bernoulli probability mass
#' \deqn{P(Y = y; \mu) = \mu^{y}(1-\mu)^{1-y}, \qquad y \in \{0, 1\},}
#' by calling [stats::dbinom()] at `size = 1`. The `pdf` in the generic's name
#' is the density with respect to counting measure, so what this returns is a
#' probability.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param y A numeric vector of zeros and ones. Any other value gives 0 with a
#'   warning from [stats::dbinom()].
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. A value of length 1 is recycled.
#'   `mu` must lie in \eqn{(0, 1)}; a value outside \eqn{[0, 1]} gives `NaN`
#'   with a warning.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned, which
#'   stays finite for a probability so small that the mass underflows.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(y), length(mu))`, one value per observation.
#'
#' @seealso [distrib_pdf.BinomialDistrib()] for several trials,
#'   [distrib_gradient.BernoulliDistrib()] for the derivatives of the log-mass,
#'   and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # The method is stats::dbinom at size 1.
#' all.equal(distrib_pdf(d, c(0, 1, 1), list(mu = 0.3)),
#'           dbinom(c(0, 1, 1), size = 1, prob = 0.3))
#'
#' # The two masses are 1 - p and p, and sum to one.
#' distrib_pdf(d, c(0, 1), list(mu = 0.3))
#'
#' # A probability may vary by observation, one value each.
#' distrib_pdf(d, c(0, 1, 1), list(mu = c(0.1, 0.5, 0.9)))
S7::method(distrib_pdf, BernoulliDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbinom(
    x = y,
    size = 1,
    prob = theta[[1]],
    log = log
  )
}

#' @title Bernoulli Cumulative Distribution Function
#' @name distrib_cdf.BernoulliDistrib
#' @description
#' Computes the Bernoulli distribution function by calling [stats::pbinom()] at
#' `size = 1`. It takes three values only:
#' \deqn{F(q; \mu) = \begin{cases}
#'   0, & q < 0,\\
#'   1 - \mu, & 0 \le q < 1,\\
#'   1, & q \ge 1.
#' \end{cases}}
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param q A numeric vector of quantiles. A non-integer is floored.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `q`. `mu` must lie in \eqn{(0, 1)}.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu))`.
#'
#' @seealso [distrib_quantile.BernoulliDistrib()] for the generalized inverse,
#'   [distrib_pdf.BernoulliDistrib()] for the mass, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' th <- list(mu = 0.3)
#'
#' # Three values, whatever q is.
#' distrib_cdf(d, c(-1, 0, 0.5, 1, 2), th)
#'
#' # The jump at 1 is exactly the mass there.
#' c(jump = distrib_cdf(d, 1, th) - distrib_cdf(d, 0, th),
#'   mass = distrib_pdf(d, 1, th))
#'
#' # The method is stats::pbinom at size 1.
#' all.equal(distrib_cdf(d, c(0, 1), th), pbinom(c(0, 1), size = 1, prob = 0.3))
S7::method(distrib_cdf, BernoulliDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbinom(
    q = q,
    size = 1,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Bernoulli Quantile Function
#' @name distrib_quantile.BernoulliDistrib
#' @description
#' Computes the generalized inverse of the Bernoulli distribution function by
#' calling [stats::qbinom()] at `size = 1`. It returns 0 while
#' \eqn{p \le 1 - \mu} and 1 afterwards, so the whole function is one step at
#' \eqn{1 - \mu}.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `p`. `mu` must lie in \eqn{(0, 1)}.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#'
#' @return A numeric vector of zeros and ones, of length
#'   `max(length(p), length(mu))`.
#'
#' @seealso [distrib_cdf.BernoulliDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' th <- list(mu = 0.3)
#'
#' # One step, at p = 1 - mu.
#' distrib_quantile(d, c(0.1, 0.5, 0.69, 0.71, 0.99), th)
#'
#' # The median is 0 whenever mu < 1/2 and 1 whenever mu > 1/2.
#' vapply(c(0.3, 0.7), function(p) distrib_quantile(d, 0.5, list(mu = p)),
#'        numeric(1))
S7::method(distrib_quantile, BernoulliDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbinom(
    p = p,
    size = 1,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Bernoulli Random Number Generator
#' @name distrib_rng.BernoulliDistrib
#' @description
#' Draws `n` independent zero-one variates by calling [stats::rbinom()] at
#' `size = 1`, so the draws come from R's own generator and depend on
#' `.Random.seed` in the usual way.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of length `n`. A value of length 1 is recycled, so a vector
#'   of length `n` draws one variate per probability, which is what a
#'   regression on a Bernoulli response needs. `mu` must lie in \eqn{(0, 1)}.
#'
#' @return An integer vector of `n` zeros and ones.
#'
#' @seealso [fit_distrib()] to estimate the probability back, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # Same generator as stats::rbinom, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 5, list(mu = 0.3))
#' set.seed(2)
#' identical(a, rbinom(5, size = 1, prob = 0.3))
#'
#' # The sample proportion estimates mu.
#' set.seed(4)
#' z <- distrib_rng(d, 2e4, list(mu = 0.3))
#' c(proportion = mean(z), variance = var(z), p_times_q = 0.3 * 0.7)
#'
#' # A probability per observation, which is what a regression supplies.
#' set.seed(5)
#' p <- plogis(seq(-2, 2, length.out = 6))
#' rbind(p = round(p, 3), draw = distrib_rng(d, 6, list(mu = p)))
S7::method(distrib_rng, BernoulliDistrib) <- function(distrib, n, theta) {
  stats::rbinom(
    n = n,
    size = 1,
    prob = theta[[1]]
  )
}

#' @title Bernoulli Score
#' @name distrib_gradient.BernoulliDistrib
#' @description
#' Computes the first derivative of the Bernoulli log-mass with respect to the
#' probability, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - \dfrac{1-y}{1-\mu}
#'       = \dfrac{y - \mu}{\mu(1-\mu)}.}
#' The residual is divided by the variance, so an observation near a
#' probability of 0 or 1 that goes the other way contributes an arbitrarily
#' large score.
#'
#' On the **link** scale with the default logit the generic's chain rule gives
#' \eqn{\partial\ell/\partial\eta = y - \mu}, the raw residual: the logit is
#' the canonical link of this family.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param y A numeric vector of zeros and ones.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
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
#' \eqn{\ell} is the log-mass of one observation and \eqn{\mu \in (0,1)} the
#' success probability, with variance \eqn{\mu(1-\mu)}.
#' \eqn{\eta = \log(\mu/(1-\mu))} is the parameter on the link scale, the log
#' odds. A **canonical** link makes the log-mass linear in the sufficient
#' statistic times \eqn{\eta}.
#'
#' @seealso [distrib_hessian.BernoulliDistrib()] for the second derivative,
#'   [distrib_expected_hessian.BernoulliDistrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' y <- c(0, 1, 1)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out.
#' all.equal(distrib_gradient(d, y, th)$mu, (y - 0.3) / (0.3 * 0.7))
#'
#' # On the canonical logit link the score is the raw residual.
#' distrib_gradient(d, y, th, scale = "link")$mu
#' y - 0.3
#'
#' # The summed score vanishes at the sample proportion, which is the estimate.
#' set.seed(4)
#' z <- distrib_rng(d, 1000, list(mu = 0.3))
#' sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
S7::method(distrib_gradient, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  bernoulli_gradient_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Observed Hessian
#' @name distrib_hessian.BernoulliDistrib
#' @description
#' Computes the second derivative of the Bernoulli log-mass with respect to the
#' probability, one value per observation, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{1-y}{(1-\mu)^2}.}
#' It is negative for every admissible \eqn{\mu} and takes one of two values
#' according to whether the observation is 0 or 1.
#'
#' On the **link** scale with the default logit the chain rule gives
#' \eqn{\partial^2\ell/\partial\eta^2 = -\mu(1-\mu)}, which carries **no data
#' at all**. That is the defining property of a canonical link and makes the
#' observed and the expected information the same matrix there; see
#' [distrib_expected_hessian.BernoulliDistrib()].
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param y A numeric vector of zeros and ones.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list of one numeric vector, `mu_mu`, of length
#'   `max(length(y), length(mu))`.
#'
#' @seealso [distrib_gradient.BernoulliDistrib()] for the score,
#'   [distrib_expected_hessian.BernoulliDistrib()] for the expectation of this
#'   quantity, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' y <- c(0, 1, 1)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out; two values only.
#' all.equal(distrib_hessian(d, y, th)$mu_mu,
#'           -y / 0.3^2 - (1 - y) / 0.7^2)
#'
#' # On the canonical logit link the curvature carries no data: one value,
#' # repeated, equal to -mu(1-mu).
#' distrib_hessian(d, y, th, scale = "link")$mu_mu
#' -0.3 * 0.7
#'
#' # A central difference of the score reproduces the parameter-scale value.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 0.3 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.3 - eps))$mu
#' all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
#'           tolerance = 1e-6)
S7::method(distrib_hessian, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  bernoulli_hessian_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Expected Hessian
#' @name distrib_expected_hessian.BernoulliDistrib
#' @description
#' Returns the expectation of the observed second derivative under the model,
#' in closed form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu(1-\mu)},}
#' which follows from \eqn{\mathbb{E}[Y] = \mu} in
#' \eqn{-y/\mu^2 - (1-y)/(1-\mu)^2}. The Fisher information for one observation
#' is \eqn{1/(\mu(1-\mu))}, the reciprocal of the variance, so it is smallest
#' at \eqn{\mu = 1/2} and grows without bound as the probability approaches
#' either endpoint.
#'
#' On the **link** scale with the default logit the value is
#' \eqn{-\mu(1-\mu)}, and so is the observed Hessian: the logit is the
#' canonical link, the observed curvature there carries no data, and the two
#' coincide exactly. Fisher scoring and Newton's method therefore take the same
#' step on a logistic regression, which is what makes iteratively reweighted
#' least squares both at once.
#'
#' Because the value does not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
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
#'   `max(length(y), length(mu))` and constant at \eqn{-1/(\mu(1-\mu))}.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' **observed information** under the model. The Bernoulli is a regular family,
#' so the second Bartlett identity holds and this equals the variance of the
#' score.
#'
#' @seealso [distrib_hessian.BernoulliDistrib()] for the observed quantity this
#'   is the expectation of, [fisher_scoring()], which inverts it, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' th <- list(mu = 0.3)
#'
#' # A single number, the reciprocal of the variance.
#' unique(distrib_expected_hessian(d, c(0, 1, 1), th)$mu_mu)
#' -1 / (0.3 * 0.7)
#'
#' # Least information where the outcome is most uncertain.
#' vapply(c(0.5, 0.1, 0.01), function(p)
#'   -distrib_expected_hessian(d, 0, list(mu = p))$mu_mu, numeric(1))
#'
#' # On the canonical logit link the observed and the expected values agree
#' # exactly, at every observation.
#' rbind(observed = distrib_hessian(d, c(0, 1, 1), th, scale = "link")$mu_mu,
#'       expected = distrib_expected_hessian(d, c(0, 1, 1), th,
#'                                           scale = "link")$mu_mu)
S7::method(distrib_expected_hessian, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  bernoulli_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Third-Order Derivative
#' @name distrib_deriv3.BernoulliDistrib
#' @description
#' Computes the third derivative of the Bernoulli log-mass with respect to the
#' probability, in closed form:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{2y}{\mu^3} - \dfrac{2(1-y)}{(1-\mu)^3}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}, which gives
#' \eqn{2/\mu^2 - 2/(1-\mu)^2}. Both routes are closed form, so no quadrature
#' is run and `approx` and `nsim` are ignored.
#'
#' The family has one parameter, so there is one component rather than the four
#' a two-parameter family carries at this order.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param y A numeric vector of zeros and ones. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
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
#' \eqn{\ell^{(\mu\mu\mu)}} is the third derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_hessian.BernoulliDistrib()] for the order below and
#'   [distrib_deriv4.BernoulliDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' y <- c(0, 1, 1)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out.
#' all.equal(distrib_deriv3(d, y, th)$mu_mu_mu,
#'           2 * y / 0.3^3 - 2 * (1 - y) / 0.7^3)
#'
#' # The expected value, and zero at mu = 1/2 by symmetry.
#' unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#' 2 / 0.3^2 - 2 / 0.7^2
#' distrib_deriv3(d, 0, list(mu = 0.5), expected = TRUE)$mu_mu_mu
#'
#' # A central difference of the Hessian reproduces the observed value.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 0.3 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.3 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv3, BernoulliDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) bernoulli_deriv3_expected_cpp(y, theta[[1]], threads)
  else bernoulli_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Fourth-Order Derivative
#' @name distrib_deriv4.BernoulliDistrib
#' @description
#' Computes the fourth derivative of the Bernoulli log-mass with respect to the
#' probability, in closed form:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4} = -\dfrac{6y}{\mu^4} - \dfrac{6(1-y)}{(1-\mu)^4}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}, which gives \eqn{-6/\mu^3 - 6/(1-\mu)^3}. Both
#' routes are closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param y A numeric vector of zeros and ones. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. `mu` must lie in \eqn{(0, 1)}.
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
#' \eqn{\ell^{(\mu\mu\mu\mu)}} is the fourth derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.BernoulliDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- bernoulli_distrib()
#' y <- c(0, 1, 1)
#' th <- list(mu = 0.3)
#'
#' # The closed form, written out; negative at every admissible mu.
#' all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           -6 * y / 0.3^4 - 6 * (1 - y) / 0.7^4)
#'
#' # The expected value.
#' unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
#' -6 / 0.3^3 - 6 / 0.7^3
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 0.3 + eps))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.3 - eps))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, BernoulliDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) bernoulli_deriv4_expected_cpp(y, theta[[1]], threads)
  else bernoulli_deriv4_cpp(y, theta[[1]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Bernoulli Distribution
#'
#' @description
#' Builds the distribution object for the Bernoulli family parametrized by its
#' success probability \eqn{\mu \in (0, 1)}. The returned object carries
#' closed-form derivatives of the log-mass to fourth order, observed and
#' expected, and closed-form moments.
#'
#' This is the response distribution of logistic regression, and the default
#' link is the logit that gives that model its name. The logit is the canonical
#' link here, so the observed and the expected information coincide on its
#' scale and iteratively reweighted least squares is Fisher scoring and
#' Newton's method at once.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the probability
#'   \eqn{\mu}. Defaults to [linkfunctions7::logit_link()], the canonical link,
#'   which maps \eqn{(0, 1)} onto the whole line.
#'   [linkfunctions7::probit_link()] and [linkfunctions7::cloglog_link()] are
#'   the usual alternatives; under either the observed and expected information
#'   differ.
#'
#' @details
#' # The parametrization
#'
#' The mass on \eqn{y \in \{0, 1\}} is
#' \deqn{P(Y = y; \mu) = \mu^{y}(1-\mu)^{1-y},}
#' with \eqn{\mu \in (0, 1)}. The mean is \eqn{\mu}, the variance
#' \eqn{\mu(1-\mu)}, the skewness \eqn{(1-2\mu)/\sqrt{\mu(1-\mu)}} and the
#' excess kurtosis \eqn{(1 - 6\mu(1-\mu))/(\mu(1-\mu))}. The family has no
#' dispersion parameter: fixing the mean fixes the variance, which is why
#' overdispersion in binary data has to come from elsewhere, most often from a
#' random effect or from grouping.
#'
#' # The canonical link
#'
#' The log-mass is \eqn{y\log(\mu/(1-\mu)) + \log(1-\mu)}, linear in \eqn{y}
#' once \eqn{\eta = \log(\mu/(1-\mu))} is the parameter. The logit is therefore
#' the canonical link, and on its scale
#' \deqn{\dfrac{\partial \ell}{\partial \eta} = y - \mu, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \eta^2} = -\mu(1-\mu).}
#' The second carries no data, so the observed and the expected information
#' coincide. Under a probit or a complementary log-log link they do not, and
#' the choice between Fisher scoring and Newton's method becomes a real one.
#'
#' # Derivatives
#'
#' On the parameter scale,
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1-\mu)}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{1-y}{(1-\mu)^2},}
#' with expectations \eqn{0} and \eqn{-1/(\mu(1-\mu))}. The information is the
#' reciprocal of the variance, so it is smallest at \eqn{\mu = 1/2} and grows
#' without bound towards either endpoint: a near-certain outcome is very
#' informative about its own probability.
#'
#' # Estimation
#'
#' The maximum likelihood estimate is the sample proportion,
#' \eqn{\hat\mu = \bar y}, in closed form. A sample of all zeros or all ones
#' drives the estimate to the boundary, where the logit is infinite; that is
#' separation, and it is a property of the data.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation and \eqn{\mu \in (0,1)} the
#' success probability. \eqn{\eta = \log(\mu/(1-\mu))} is the log odds, the
#' parameter on the link scale. The **canonical** link of an exponential family
#' is the one that makes the log-mass linear in the sufficient statistic.
#'
#' @return An S7 object of class `BernoulliDistrib`, inheriting from
#'   `discrete_distrib`, with `distrib_name` `"bernoulli"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, 1)`, `params` `"mu"`,
#'   `params_interpretation` `c(mu = "probability")`, `n_params` `1`,
#'   `params_bounds` the list of \eqn{(0, 1)}, and `link_params` the one link
#'   given here.
#'
#' @seealso
#' [binomial_distrib()], of which this is the one-trial case;
#' [betabinom1_distrib()] and [betabinom2_distrib()] for grouped binary data
#' that is overdispersed; [linkfunctions7::logit_link()],
#' [linkfunctions7::probit_link()] and [linkfunctions7::cloglog_link()] for the
#' usual links; [fit_distrib()] to estimate the probability;
#' [BernoulliDistrib] for the class.
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
#' d <- bernoulli_distrib()
#' d
#'
#' # The two masses, summing to one.
#' distrib_pdf(d, c(0, 1), list(mu = 0.3))
#'
#' # Mean and variance are tied, the variance being maximal at 1/2.
#' vapply(c(0.1, 0.5, 0.9), function(p) {
#'   th <- list(mu = p)
#'   c(mean = mean(d, th), var = variance(d, th), skew = skewness(d, th))
#' }, numeric(3))
#'
#' # The estimate is the sample proportion, in closed form.
#' set.seed(4)
#' z <- distrib_rng(d, 1000, list(mu = 0.3))
#' c(fitted = unname(coef(fit_distrib(d, z))), proportion = mean(z))
#'
#' # On the canonical logit link the observed and expected information agree;
#' # under a probit they do not.
#' th <- list(mu = 0.3)
#' probit <- bernoulli_distrib(link_mu = linkfunctions7::probit_link())
#' rbind(logit = c(obs = distrib_hessian(d, 1, th, scale = "link")$mu_mu,
#'                 exp = distrib_expected_hessian(d, 1, th,
#'                                                scale = "link")$mu_mu),
#'       probit = c(obs = distrib_hessian(probit, 1, th, scale = "link")$mu_mu,
#'                  exp = distrib_expected_hessian(probit, 1, th,
#'                                                 scale = "link")$mu_mu))
#'
#' @export
bernoulli_distrib <- function(link_mu = logit_link()) {
  
  BernoulliDistrib(
    distrib_name = "bernoulli",
    dimension = "univariate",
    bounds = c(0, 1),
    
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
