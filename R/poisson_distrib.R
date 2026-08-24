#' @include distrib.R generics.R
NULL

#' @title Poisson Distribution Class
#' @name PoissonDistrib
#'
#' @description
#' The S7 class of the Poisson family parametrized by its mean \eqn{\mu > 0},
#' with mass \eqn{P(Y = y) = e^{-\mu}\mu^y/y!} on the non-negative integers. It
#' inherits from `discrete_distrib`, so its support is counted rather than
#' integrated: expectations are exact sums and no derivative with respect to
#' the response is defined.
#'
#' The mean and the variance are both \eqn{\mu}, which is the constraint a
#' count model most often has to relax; [negbin2_distrib()] and
#' [negbin1_distrib()] do that, and both contain this family in the limit.
#'
#' The default link is the logarithm, which is the **canonical** link here. On
#' its scale the observed and the expected information coincide, so Fisher
#' scoring and Newton's method take the same step.
#'
#' Build one with [poisson_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `PoissonDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [poisson_distrib()] they hold `"poisson"`,
#'   `"univariate"`, `c(0, Inf)`, `"mu"`, `c(mu = "mean")`, `1`, the domain
#'   \eqn{(0, \infty)}, and the one link.
#'
#' @seealso [poisson_distrib()] to build one;
#'   [negbin2_distrib()] and [negbin1_distrib()] for overdispersed counts;
#'   [zero_inflated()] and [zero_adjusted()] for excess zeros;
#'   [distrib_expected_hessian.PoissonDistrib()] for the information.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.PoissonDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.PoissonDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.PoissonDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.PoissonDistrib],
#'   [`distrib_gradient()`][distrib_gradient.PoissonDistrib],
#'   [`distrib_hessian()`][distrib_hessian.PoissonDistrib],
#'   [`distrib_pdf()`][distrib_pdf.PoissonDistrib],
#'   [`distrib_quantile()`][distrib_quantile.PoissonDistrib],
#'   [`distrib_rng()`][distrib_rng.PoissonDistrib]
#'
#' Everything else is inherited from [discrete_distrib()]. Derivatives with
#' respect to the response are refused by that parent: the support is a
#' lattice, so there is nothing to differentiate along.
#'
#' @examples
#' d <- poisson_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#' d@params
#' d@bounds
#'
#' # Equidispersion: the mean is the variance at every parameter value.
#' vapply(c(0.5, 3, 20), function(m) {
#'   th <- list(mu = m)
#'   c(mean = mean(d, th), var = variance(d, th))
#' }, numeric(2))
#'
#' # The log link is canonical, so on its scale the observed and the expected
#' # information are the same number.
#' th <- list(mu = 3)
#' rbind(observed = distrib_hessian(d, c(0, 2, 7), th, scale = "link")$mu_mu,
#'       expected = distrib_expected_hessian(d, c(0, 2, 7), th,
#'                                           scale = "link")$mu_mu)
PoissonDistrib <- S7::new_class("PoissonDistrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Poisson Probability Mass Function
#' @name distrib_pdf.PoissonDistrib
#' @description
#' Computes the Poisson probability mass
#' \deqn{P(Y = y; \mu) = \dfrac{e^{-\mu}\mu^{y}}{y!}, \qquad y = 0, 1, 2, \dots}
#' by calling [stats::dpois()]. The `pdf` in the generic's name is the density
#' with respect to counting measure, so what this returns is a probability and
#' is at most 1.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param y A numeric vector of counts. A non-integer or negative value gives 0
#'   with a warning from [stats::dpois()].
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `y`. A value of length 1 is recycled.
#'   `mu` must be strictly positive; a negative value gives `NaN` with a
#'   warning.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned, which
#'   stays finite for a count far above the mean where the mass itself
#'   underflows. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(y), length(mu))`, one value per observation.
#'
#' @seealso [distrib_cdf.PoissonDistrib()] for the distribution function,
#'   [distrib_gradient.PoissonDistrib()] for the derivatives of the log-mass,
#'   and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' y <- c(0, 2, 7)
#'
#' # The method is stats::dpois.
#' all.equal(distrib_pdf(d, y, list(mu = 3)), dpois(y, lambda = 3))
#'
#' # A probability mass: it sums to one over the support.
#' sum(distrib_pdf(d, 0:200, list(mu = 3)))
#'
#' # Far above the mean the mass underflows and its logarithm does not.
#' distrib_pdf(d, 400, list(mu = 3))
#' distrib_pdf(d, 400, list(mu = 3), log = TRUE)
S7::method(distrib_pdf, PoissonDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dpois(
    x = y,
    lambda = theta[[1]],
    log = log
  )
}

#' @title Poisson Cumulative Distribution Function
#' @name distrib_cdf.PoissonDistrib
#' @description
#' Computes the Poisson distribution function
#' \deqn{F(q; \mu) = P(Y \le q) = \sum_{k = 0}^{\lfloor q \rfloor} \dfrac{e^{-\mu}\mu^{k}}{k!}}
#' by calling [stats::ppois()], which evaluates it through the incomplete gamma
#' function rather than by summing. The function is a step function, constant
#' between integers, so `F(2)` and `F(2.9)` are the same number.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param q A numeric vector of quantiles. A non-integer is floored, and a
#'   negative value gives 0.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `q`. `mu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   computed directly and so exact far into the upper tail.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu))`.
#'
#' @seealso [distrib_quantile.PoissonDistrib()] for the generalized inverse,
#'   [distrib_pdf.PoissonDistrib()] for the mass, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- poisson_distrib()
#' th <- list(mu = 3)
#'
#' # The method is stats::ppois.
#' all.equal(distrib_cdf(d, c(0, 2, 7), th), ppois(c(0, 2, 7), lambda = 3))
#'
#' # A step function: constant between integers.
#' distrib_cdf(d, c(2, 2.5, 2.9), th)
#'
#' # The jump at an integer is exactly the mass there.
#' c(jump = distrib_cdf(d, 2, th) - distrib_cdf(d, 1, th),
#'   mass = distrib_pdf(d, 2, th))
#'
#' # The upper tail is exact where its complement would have rounded to one.
#' distrib_cdf(d, 60, th, lower.tail = FALSE)
S7::method(distrib_cdf, PoissonDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::ppois(
    q = q,
    lambda = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Poisson Quantile Function
#' @name distrib_quantile.PoissonDistrib
#' @description
#' Computes the generalized inverse of the Poisson distribution function,
#' \deqn{Q(p; \mu) = \min\{k \in \{0, 1, 2, \dots\} : F(k; \mu) \ge p\},}
#' by calling [stats::qpois()]. The distribution function is a step function,
#' so its inverse is a step function too and the round trip through
#' [distrib_cdf.PoissonDistrib()] returns a probability **at least** `p`, not
#' `p` itself.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. `p = 1` gives `Inf`.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of the length of `p`. `mu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#'
#' @return A numeric vector of non-negative integers, of length
#'   `max(length(p), length(mu))`.
#'
#' @seealso [distrib_cdf.PoissonDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' th <- list(mu = 3)
#'
#' # Integers, being the smallest count whose cumulative mass reaches p.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # The round trip overshoots, the support being a lattice.
#' p <- c(0.025, 0.5, 0.975)
#' rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
#'
#' # It is the smallest k with F(k) >= p, which the definition checks.
#' k <- distrib_quantile(d, 0.9, th)
#' c(k = k, F_k = distrib_cdf(d, k, th), F_k_minus_1 = distrib_cdf(d, k - 1, th))
S7::method(distrib_quantile, PoissonDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qpois(
    p = p,
    lambda = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Poisson Random Number Generator
#' @name distrib_rng.PoissonDistrib
#' @description
#' Draws `n` independent Poisson counts by calling [stats::rpois()], so the
#' draws come from R's own generator and depend on `.Random.seed` in the usual
#' way. The cumulative-table fallback the discrete base class supplies is
#' bypassed.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with the single component `mu`, a numeric vector
#'   of length 1 or of length `n`. A value of length 1 is recycled, so a vector
#'   of length `n` draws one count per mean. `mu` must be strictly positive.
#'
#' @return An integer vector of `n` non-negative counts.
#'
#' @seealso [fit_distrib()] to estimate the mean back, and [distrib_rng()] for
#'   the generic.
#'
#' @examples
#' d <- poisson_distrib()
#'
#' # Same generator as stats::rpois, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 5, list(mu = 3))
#' set.seed(2)
#' identical(a, rpois(5, lambda = 3))
#'
#' # The sample mean and variance both estimate mu, this family being
#' # equidispersed.
#' set.seed(4)
#' z <- distrib_rng(d, 2e4, list(mu = 4.2))
#' c(mean = mean(z), var = var(z))
S7::method(distrib_rng, PoissonDistrib) <- function(distrib, n, theta) {
  stats::rpois(
    n = n,
    lambda = theta[[1]]
  )
}

#' @title Poisson Score
#' @name distrib_gradient.PoissonDistrib
#' @description
#' Computes the first derivative of the Poisson log-mass with respect to the
#' mean, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - 1 = \dfrac{y - \mu}{\mu}.}
#' The log-mass is \eqn{y\log\mu - \mu - \log y!}, so the score is the residual
#' divided by the mean and its sum vanishes exactly at \eqn{\hat\mu = \bar y}.
#'
#' On the **link** scale with the default logarithm the generic's chain rule
#' gives \eqn{\partial\ell/\partial\eta = y - \mu}, the raw residual: the log
#' is the canonical link of this family and the score there is the sufficient
#' statistic minus its mean.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param y A numeric vector of counts.
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
#' \eqn{\ell} is the log-mass of one observation and \eqn{\mu > 0} the mean,
#' which is also the variance. \eqn{\eta = \log\mu} is the parameter on the
#' link scale. A **canonical** link is the one for which the log-mass is linear
#' in the sufficient statistic times \eqn{\eta}.
#'
#' @seealso [distrib_hessian.PoissonDistrib()] for the second derivative,
#'   [distrib_expected_hessian.PoissonDistrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # The closed form, written out.
#' all.equal(distrib_gradient(d, y, th)$mu, (y - 3) / 3)
#'
#' # On the canonical log link the score is the raw residual.
#' distrib_gradient(d, y, th, scale = "link")$mu
#' y - 3
#'
#' # The summed score vanishes at the sample mean, which is the estimate.
#' set.seed(4)
#' z <- distrib_rng(d, 1000, list(mu = 4.2))
#' sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
S7::method(distrib_gradient, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  poisson_gradient_cpp(y, theta[[1]], threads)
}

#' @title Poisson Observed Hessian
#' @name distrib_hessian.PoissonDistrib
#' @description
#' Computes the second derivative of the Poisson log-mass with respect to the
#' mean, one value per observation, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}.}
#' It depends on the data only through \eqn{y}, and is exactly zero at an
#' observed count of zero; the curvature at \eqn{y = 0} comes entirely from the
#' \eqn{-\mu} term of the log-mass, which is linear.
#'
#' On the **link** scale with the default logarithm the chain rule gives
#' \eqn{\partial^2\ell/\partial\eta^2 = -\mu}, which carries **no data at
#' all**. That is the defining property of a canonical link, and it makes the
#' observed and the expected information the same matrix there; see
#' [distrib_expected_hessian.PoissonDistrib()].
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
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
#' @seealso [distrib_gradient.PoissonDistrib()] for the score,
#'   [distrib_expected_hessian.PoissonDistrib()] for the expectation of this
#'   quantity, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # The closed form, written out; zero at an observed zero.
#' all.equal(distrib_hessian(d, y, th)$mu_mu, -y / 3^2)
#'
#' # On the canonical log link the curvature carries no data: one value,
#' # repeated, equal to -mu.
#' distrib_hessian(d, y, th, scale = "link")$mu_mu
#'
#' # A central difference of the score reproduces the parameter-scale value.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 3 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 3 - eps))$mu
#' all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
#'           tolerance = 1e-6)
S7::method(distrib_hessian, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  poisson_hessian_cpp(y, theta[[1]], threads)
}

#' @title Poisson Expected Hessian
#' @name distrib_expected_hessian.PoissonDistrib
#' @description
#' Returns the expectation of the observed second derivative under the model,
#' in closed form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu},}
#' which follows from \eqn{\mathbb{E}[Y] = \mu} in \eqn{-y/\mu^2}. The Fisher
#' information for one observation is \eqn{1/\mu}, so the asymptotic standard
#' error of \eqn{\hat\mu} is \eqn{\sqrt{\mu/n}}, the standard error of a
#' sample mean whose variance is \eqn{\mu}.
#'
#' On the **link** scale with the default logarithm the value is \eqn{-\mu},
#' and so is the observed Hessian: the log is the canonical link, the observed
#' curvature there carries no data, and the two coincide exactly. Fisher
#' scoring and Newton's method therefore take the same step on a Poisson model
#' with a log link, which is why the two are not distinguished in the classical
#' generalized-linear-model literature for this case.
#'
#' Because the value does not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
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
#'   `max(length(y), length(mu))` and constant at \eqn{-1/\mu}.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\mu^2]}, the expectation of the
#' **observed information** under the model. The Poisson is a regular family,
#' so the second Bartlett identity holds and this equals the variance of the
#' score.
#'
#' @seealso [distrib_hessian.PoissonDistrib()] for the observed quantity this
#'   is the expectation of, [fisher_scoring()], which inverts it, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' th <- list(mu = 3)
#'
#' # A single number, -1/mu.
#' unique(distrib_expected_hessian(d, c(0, 2, 7), th)$mu_mu)
#' -1 / 3
#'
#' # On the canonical log link the observed and the expected values agree
#' # exactly, at every observation.
#' rbind(observed = distrib_hessian(d, c(0, 2, 7), th, scale = "link")$mu_mu,
#'       expected = distrib_expected_hessian(d, c(0, 2, 7), th,
#'                                           scale = "link")$mu_mu)
#'
#' # It is the variance of the score, this family being regular.
#' set.seed(4)
#' z <- distrib_rng(d, 2e5, th)
#' c(var_of_score = mean(distrib_gradient(d, z, th)$mu^2),
#'   information = -distrib_expected_hessian(d, 0, th)$mu_mu)
S7::method(distrib_expected_hessian, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  poisson_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Poisson Third-Order Derivative
#' @name distrib_deriv3.PoissonDistrib
#' @description
#' Computes the third derivative of the Poisson log-mass with respect to the
#' mean, in closed form:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{2y}{\mu^3}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}, which gives \eqn{2/\mu^2}. Both routes are closed
#' form, so no quadrature is run and `approx` and `nsim` are ignored.
#'
#' The family has one parameter, so there is one component rather than the four
#' a two-parameter family carries at this order.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is used.
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
#' \eqn{\ell^{(\mu\mu\mu)}} is the third derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives; a
#' subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.PoissonDistrib()] for the order below and
#'   [distrib_deriv4.PoissonDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # The closed form, written out; zero at an observed zero.
#' all.equal(distrib_deriv3(d, y, th)$mu_mu_mu, 2 * y / 3^3)
#'
#' # The expected value is 2/mu^2.
#' unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#' 2 / 3^2
#'
#' # A central difference of the Hessian reproduces the observed value.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 3 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 3 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv3, PoissonDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) poisson_deriv3_expected_cpp(y, theta[[1]], threads)
  else poisson_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Poisson Fourth-Order Derivative
#' @name distrib_deriv4.PoissonDistrib
#' @description
#' Computes the fourth derivative of the Poisson log-mass with respect to the
#' mean, in closed form:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \mu^4} = -\dfrac{6y}{\mu^4}.}
#' With `expected = TRUE` the expectation is returned, obtained by substituting
#' \eqn{\mathbb{E}[Y] = \mu}, which gives \eqn{-6/\mu^3}. Both routes are
#' closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is used.
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
#' \eqn{\ell^{(\mu\mu\mu\mu)}} is the fourth derivative of the log-mass with
#' respect to \eqn{\mu}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.PoissonDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- poisson_distrib()
#' y <- c(0, 2, 7)
#' th <- list(mu = 3)
#'
#' # The closed form, written out.
#' all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu, -6 * y / 3^4)
#'
#' # The expected value is -6/mu^3.
#' unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
#' -6 / 3^3
#'
#' # A central difference of the third order reproduces the fourth.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 3 + eps))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 3 - eps))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
#'           tolerance = 1e-5)
S7::method(distrib_deriv4, PoissonDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) poisson_deriv4_expected_cpp(y, theta[[1]], threads)
  else poisson_deriv4_cpp(y, theta[[1]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Poisson Distribution
#'
#' @description
#' Builds the distribution object for the Poisson family parametrized by its
#' mean \eqn{\mu > 0}, with mass \eqn{e^{-\mu}\mu^y/y!} on the non-negative
#' integers. The returned object carries closed-form derivatives of the
#' log-mass to fourth order, observed and expected, and closed-form moments.
#'
#' The family is **equidispersed**: its variance equals its mean at every
#' parameter value. Real counts are usually more variable than that, so a
#' Poisson fit is often the null a count analysis argues against;
#' [negbin2_distrib()] and [negbin1_distrib()] add a dispersion parameter and
#' both contain this family in a limit.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and is the **canonical** link of this
#'   family.
#'
#' @details
#' # The parametrization
#'
#' The mass on \eqn{y \in \{0, 1, 2, \dots\}} is
#' \deqn{P(Y = y; \mu) = \dfrac{e^{-\mu}\mu^{y}}{y!},}
#' with \eqn{\mu \in (0, \infty)}. The mean and the variance are both
#' \eqn{\mu}, the skewness is \eqn{\mu^{-1/2}} and the excess kurtosis
#' \eqn{\mu^{-1}}, so the distribution is right skewed at a small mean and
#' close to Gaussian at a large one.
#'
#' # The canonical link, and what it buys
#'
#' The log-mass is \eqn{y\log\mu - \mu - \log y!}, which is linear in \eqn{y}
#' once \eqn{\eta = \log\mu} is the parameter. The logarithm is therefore the
#' canonical link, and on its scale
#' \deqn{\dfrac{\partial \ell}{\partial \eta} = y - \mu, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \eta^2} = -\mu.}
#' The second of these carries no data, so the observed and the expected
#' information coincide: Fisher scoring and Newton's method take the same step,
#' and the iteratively reweighted least squares of a Poisson log-linear model
#' is both at once. Under any other link the two differ and the choice matters.
#'
#' # Derivatives
#'
#' On the parameter scale, with \eqn{\ell = y\log\mu - \mu - \log y!},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\mu}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}, \qquad
#'       \dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{2y}{\mu^3}, \qquad
#'       \dfrac{\partial^4 \ell}{\partial \mu^4} = -\dfrac{6y}{\mu^4},}
#' and their expectations follow by substituting \eqn{\mathbb{E}[Y] = \mu}:
#' \eqn{0}, \eqn{-1/\mu}, \eqn{2/\mu^2} and \eqn{-6/\mu^3}.
#'
#' # Estimation
#'
#' The maximum likelihood estimate is the sample mean, \eqn{\hat\mu = \bar y},
#' in closed form. A quick check of the equidispersion assumption is to compare
#' the sample variance with it; the example below does that on data drawn from
#' the family and on data that is not.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation and \eqn{\mu > 0} the mean,
#' which is also the variance. \eqn{\eta = \log\mu} is the parameter on the
#' link scale. The **canonical** link of an exponential family is the one that
#' makes the log-mass linear in the sufficient statistic.
#'
#' @return An S7 object of class `PoissonDistrib`, inheriting from
#'   `discrete_distrib`, with `distrib_name` `"poisson"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`,
#'   `params_interpretation` `c(mu = "mean")`, `n_params` `1`, `params_bounds`
#'   the list of \eqn{(0, \infty)}, and `link_params` the one link given here.
#'
#' @seealso
#' [negbin2_distrib()] and [negbin1_distrib()] for overdispersed counts, both
#' of which contain this family in a limit; [binomial_distrib()], of which it
#' is the many-trials limit; [zero_inflated()] and [zero_adjusted()] for excess
#' zeros; [geometric_distrib()] for a memoryless count law;
#' [fit_distrib()] to estimate the mean; [PoissonDistrib] for the class.
#'
#' @references
#' Johnson, N. L., Kemp, A. W. and Kotz, S. (2005).
#' *Univariate Discrete Distributions*, 3rd edition, Chapter 4.
#' Wiley, Hoboken.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dpois ppois qpois rpois
#'
#' @examples
#' d <- poisson_distrib()
#' d
#'
#' # The mass is R's own and sums to one over the support.
#' y <- c(0, 2, 7)
#' all.equal(distrib_pdf(d, y, list(mu = 3)), dpois(y, 3))
#' sum(distrib_pdf(d, 0:200, list(mu = 3)))
#'
#' # Equidispersed, and closer to symmetric as the mean grows.
#' vapply(c(0.5, 3, 20), function(m) {
#'   th <- list(mu = m)
#'   c(mean = mean(d, th), var = variance(d, th), skew = skewness(d, th))
#' }, numeric(3))
#'
#' # The estimate is the sample mean, in closed form.
#' set.seed(4)
#' z <- distrib_rng(d, 1000, list(mu = 4.2))
#' c(fitted = unname(coef(fit_distrib(d, z))), sample_mean = mean(z))
#'
#' # Equidispersion is a testable assumption. Drawn from the family the
#' # sample variance matches the mean; drawn with extra variability it does
#' # not, and a negative binomial is the family to reach for.
#' set.seed(6)
#' over <- rpois(1000, lambda = rgamma(1000, shape = 2, rate = 2 / 4.2))
#' rbind(poisson = c(mean = mean(z), var = var(z)),
#'       overdispersed = c(mean = mean(over), var = var(over)))
#'
#' @export
poisson_distrib <- function(link_mu = log_link()) {
  PoissonDistrib(
    distrib_name = "poisson", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)), link_params = list(mu = link_mu)
  )
}
