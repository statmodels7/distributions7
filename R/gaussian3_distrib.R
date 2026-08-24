#' @include distrib.R generics.R
NULL

#' @title Gaussian Distribution Class, Mean and Precision
#' @name Gaussian3Distrib
#'
#' @description
#' The S7 class of the Gaussian (normal) family parametrized by its mean
#' \eqn{\mu} and its precision \eqn{\tau = 1/\sigma^2 > 0}, with density
#' \eqn{f(y) = (\tau/2\pi)^{1/2}\exp\{-\tau(y-\mu)^2/2\}} on the whole real
#' line. It inherits from `continuous_distrib`, so it answers every generic of
#' the `distrib` contract; the eleven methods listed below are registered on it
#' directly and everything else comes from the parent.
#'
#' Build one with [gaussian3_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Gaussian3Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [gaussian3_distrib()] they hold `"gaussian3"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "tau")`, the interpretations
#'   `c(mu = "mean", tau = "precision")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [gaussian3_distrib()] to build one;
#'   [gaussian1_distrib()] for the same law in mean and standard deviation and
#'   [gaussian2_distrib()] for mean and variance;
#'   [distrib_pdf.Gaussian3Distrib()] and [distrib_gradient.Gaussian3Distrib()]
#'   for the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Gaussian3Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gaussian3Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gaussian3Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gaussian3Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gaussian3Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gaussian3Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gaussian3Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gaussian3Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gaussian3Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gaussian3Distrib],
#'   [`distrib_rng()`][distrib_rng.Gaussian3Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- gaussian3_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@params_bounds
#'
#' # The precision is the reciprocal of the variance, so a large tau is a
#' # tight distribution.
#' variance(d, list(mu = 0, tau = 100))
Gaussian3Distrib <- S7::new_class("Gaussian3Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Probability Density Function in Mean and Precision
#' @name distrib_pdf.Gaussian3Distrib
#' @description
#' Computes the Gaussian density
#' \deqn{f(y; \mu, \tau) = \sqrt{\dfrac{\tau}{2\pi}}
#'       \exp\left\{-\dfrac{\tau(y-\mu)^2}{2}\right\}}
#' by calling [stats::dnorm()] at `mean = mu` and `sd = 1/sqrt(tau)`, so the
#' accuracy and the underflow behavior are R's own. With `log = TRUE` the
#' logarithm is formed inside `dnorm()` and stays finite far into the tails,
#' where the density itself underflows to zero.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `tau` must be strictly positive; a zero gives `sd = Inf` and a
#'   density of 0, and a negative value gives `NaN` with a warning from
#'   [stats::dnorm()].
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(tau))`, one value per observation.
#'
#' @seealso [distrib_cdf.Gaussian3Distrib()] for the distribution function,
#'   [distrib_gradient.Gaussian3Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf.Gaussian1Distrib()] for the same density in the
#'   standard deviation, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#'
#' # The method is stats::dnorm at the reciprocal square root of tau.
#' all.equal(distrib_pdf(d, y, list(mu = 1, tau = 0.25)),
#'           dnorm(y, mean = 1, sd = 2))
#'
#' # Same law as gaussian1 at sigma = 1/sqrt(tau), to the last bit.
#' all.equal(distrib_pdf(d, y, list(mu = 1, tau = 0.25)),
#'           distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0, 1, 2), tau = c(1, 0.25, 0.0625)))
#'
#' # In the far tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 40, list(mu = 0, tau = 1))
#' distrib_pdf(d, 40, list(mu = 0, tau = 1), log = TRUE)
S7::method(distrib_pdf, Gaussian3Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(y, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]), log = log)
}

#' @title Gaussian Cumulative Distribution Function in Mean and Precision
#' @name distrib_cdf.Gaussian3Distrib
#' @description
#' Computes the Gaussian distribution function
#' \deqn{F(q; \mu, \tau) = \Phi\left(\sqrt{\tau}\,(q-\mu)\right)}
#' with \eqn{\Phi} the standard normal distribution function, by calling
#' [stats::pnorm()] at `sd = 1/sqrt(tau)`. Both tails are available exactly:
#' `lower.tail = FALSE` evaluates \eqn{1 - F} without forming the difference,
#' and `log.p = TRUE` returns a logarithm that stays finite where the
#' probability itself underflows to zero.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `tau` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(tau))`. With `log.p = TRUE` the values
#'   are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Gaussian3Distrib()] for the inverse,
#'   [distrib_pdf.Gaussian3Distrib()] for the density,
#'   [distrib_grad_cdf()] for the derivatives of this function in the
#'   parameters, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' th <- list(mu = 1, tau = 0.25)
#'
#' # The method is stats::pnorm at sd = 1/sqrt(tau).
#' all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
#'           pnorm(c(-1.2, 0.3, 2.5), mean = 1, sd = 2))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
#'
#' # Forty standard deviations out the upper tail underflows; its log does not.
#' distrib_cdf(d, 40, list(mu = 0, tau = 1), lower.tail = FALSE)
#' distrib_cdf(d, 40, list(mu = 0, tau = 1), lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, Gaussian3Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  stats::pnorm(q, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Quantile Function in Mean and Precision
#' @name distrib_quantile.Gaussian3Distrib
#' @description
#' Computes the Gaussian quantile function
#' \deqn{Q(p; \mu, \tau) = \mu + \dfrac{\Phi^{-1}(p)}{\sqrt{\tau}}}
#' by calling [stats::qnorm()] at `sd = 1/sqrt(tau)`. The Gaussian distribution
#' function is strictly increasing on the whole line, so \eqn{Q} is its exact
#' inverse and the round trip through [distrib_cdf.Gaussian3Distrib()] returns
#' `p`.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; the endpoints give `-Inf` and `Inf`.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `tau` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities, which is how a quantile deep in a tail is
#'   requested without the probability underflowing. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of length
#'   `max(length(p), length(mu), length(tau))`.
#'
#' @seealso [distrib_cdf.Gaussian3Distrib()], which this inverts;
#'   [distrib_rng.Gaussian3Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' th <- list(mu = 1, tau = 0.25)
#'
#' # The median is mu and the quartiles are symmetric about it.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # Four times the precision halves the distance from the median.
#' distrib_quantile(d, 0.975, list(mu = 0, tau = 1)) -
#'   distrib_quantile(d, 0.975, list(mu = 0, tau = 4))
S7::method(distrib_quantile, Gaussian3Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  stats::qnorm(p, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Random Number Generator in Mean and Precision
#' @name distrib_rng.Gaussian3Distrib
#' @description
#' Draws `n` independent Gaussian variates by calling [stats::rnorm()] at
#' `sd = 1/sqrt(tau)`, so the draws come from R's own normal generator and
#' depend on `.Random.seed` in the usual way. The inverse-transform fallback
#' the base class supplies is bypassed.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. `tau`
#'   must be strictly positive.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.Gaussian3Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#'
#' # Same generator as stats::rnorm at sd = 1/sqrt(tau).
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = c(0, 1, 2), tau = c(1, 0.25, 0.0625)))
#' set.seed(2)
#' identical(a, rnorm(3, mean = c(0, 1, 2), sd = c(1, 2, 4)))
#'
#' # The sample moments recover the parameters, the precision as 1/var.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, tau = 0.25))
#' c(mean = mean(z), tau = 1 / var(z))
S7::method(distrib_rng, Gaussian3Distrib) <- function(distrib, n, theta) {
  stats::rnorm(n, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]))
}

#' @title Gaussian Score in Mean and Precision
#' @name distrib_gradient.Gaussian3Distrib
#' @description
#' Computes the first derivatives of the Gaussian log-density with respect to
#' \eqn{\mu} and \eqn{\tau}, one value per observation, in closed form. With
#' \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \tau r,
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \tau} = \dfrac{1}{2\tau} - \dfrac{r^2}{2}.}
#' The precision component is linear in \eqn{r^2}. That linearity is why every
#' third and fourth derivative of this parametrization is free of the response.
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning, giving \eqn{\partial \ell / \partial \eta_j
#' = h_j'(\eta_j)\, \partial \ell / \partial \theta_j}. This method always
#' returns the parameter scale; the transformation happens in the generic.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `tau` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `tau`, each of length
#'   `max(length(y), length(mu), length(tau))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\tau = 1/\sigma^2 > 0} the precision. \eqn{\eta_j} is the coordinate of
#' parameter \eqn{j} on the unconstrained scale of its link, and \eqn{h_j' =
#' \partial \theta_j / \partial \eta_j} the chain-rule factor onto it.
#'
#' @seealso [distrib_hessian.Gaussian3Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Gaussian3Distrib()] for their expectation,
#'   [distrib_gradient.Gaussian1Distrib()] for the same score in the standard
#'   deviation, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, tau = 0.25)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' r <- y - 1
#' all.equal(g$mu, 0.25 * r)
#' all.equal(g$tau, 1 / (2 * 0.25) - r^2 / 2)
#'
#' # The summed score vanishes at the maximum likelihood estimate, where the
#' # precision is the reciprocal of the variance with divisor n.
#' set.seed(1)
#' z <- rnorm(200, mean = 3, sd = 2)
#' mle <- list(mu = mean(z), tau = 1 / mean((z - mean(z))^2))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#'
#' # log(tau) is minus twice log(sigma), so the link-scale score here is minus
#' # half gaussian1's.
#' distrib_gradient(d, y, th, scale = "link")$tau /
#'   distrib_gradient(gaussian1_distrib(), y, list(mu = 1, sigma = 2),
#'                    scale = "link")$sigma
S7::method(distrib_gradient, Gaussian3Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian3_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Observed Hessian in Mean and Precision
#' @name distrib_hessian.Gaussian3Distrib
#' @description
#' Computes the three distinct second derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{\tau}, one value per observation, in
#' closed form. With \eqn{r = y - \mu},
#' \deqn{\ell^{(\mu\mu)} = -\tau, \qquad
#'       \ell^{(\mu\tau)} = r, \qquad
#'       \ell^{(\tau\tau)} = -\dfrac{1}{2\tau^2}.}
#' Two of the three are free of the data, so the observed and the expected
#' Hessians differ in the mixed entry alone. That entry sums to
#' \eqn{\sum_i (y_i - \mu)}, which vanishes at \eqn{\hat\mu = \bar y}: the two
#' matrices agree exactly once the mean equation is solved and differ before
#' that, so Fisher scoring and Newton's method take the same step at the
#' maximum and different steps on the way there.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `tau` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_tau` and
#'   `tau_tau`, each of length `max(length(y), length(mu), length(tau))`. The
#'   three name the distinct entries of a symmetric \eqn{2 \times 2} matrix per
#'   observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density with respect
#' to parameters \eqn{i} and \eqn{j}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_gradient.Gaussian3Distrib()] for the score,
#'   [distrib_expected_hessian.Gaussian3Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Gaussian3Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, tau = 0.25)
#' h <- distrib_hessian(d, y, th)
#'
#' # Two of the three are constants; only the mixed entry carries the data.
#' h$mu_mu
#' h$tau_tau
#' all.equal(h$mu_tau, y - 1)
#'
#' # The mixed entry sums to zero at the estimated mean and not elsewhere.
#' set.seed(3)
#' z <- distrib_rng(d, 500, list(mu = 3, tau = 0.25))
#' c(at_mle = sum(distrib_hessian(d, z, list(mu = mean(z), tau = 0.25))$mu_tau),
#'   at_2   = sum(distrib_hessian(d, z, list(mu = 2, tau = 0.25))$mu_tau))
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 1, tau = 0.25 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 1, tau = 0.25 - eps))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_tau, tolerance = 1e-6)
S7::method(distrib_hessian, Gaussian3Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian3_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Expected Hessian in Mean and Precision
#' @name distrib_expected_hessian.Gaussian3Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\tau, \qquad
#'       \mathbb{E}\left[\ell^{(\mu\tau)}\right] = 0, \qquad
#'       \mathbb{E}\left[\ell^{(\tau\tau)}\right] = -\dfrac{1}{2\tau^2}.}
#' They follow from \eqn{\mathbb{E}[Y-\mu] = 0}, the other two entries being
#' free of the data already. The negative of this matrix is the Fisher
#' information for one observation; the zero off-diagonal says the mean and the
#' precision are orthogonal, so their estimates are asymptotically independent.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `tau` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available. Accepted so that the
#'   signature matches the generic's, where it selects between the Bartlett,
#'   quadrature, Monte Carlo and outer-product routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_tau` and
#'   `tau_tau`, each of length `max(length(y), length(mu), length(tau))` and
#'   constant within itself when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The Gaussian is
#' a regular family, so the second Bartlett identity holds and this equals the
#' variance of the score.
#'
#' @seealso [distrib_hessian.Gaussian3Distrib()] for the observed quantity this
#'   is the expectation of, [distrib_expected_hessian.Gaussian1Distrib()] for
#'   the same information in the standard deviation, [fisher_scoring()], which
#'   inverts it at each step, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' th <- list(mu = 1, tau = 0.25)
#'
#' # The three constants, one value per observation.
#' lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
#'
#' # The information transforms by the delta method: with dtau/dsigma =
#' # -2/sigma^3, the precision entry carries onto gaussian1's -2/sigma^2.
#' distrib_expected_hessian(d, 0, th)$tau_tau * (2 / 2^3)^2
#' distrib_expected_hessian(gaussian1_distrib(), 0,
#'                          list(mu = 1, sigma = 2))$sigma_sigma
S7::method(distrib_expected_hessian, Gaussian3Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ..., threads = 1L) {
  gaussian3_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Third-Order Derivatives in Mean and Precision
#' @name distrib_deriv3.Gaussian3Distrib
#' @description
#' Computes the four distinct third derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{\tau}, in closed form:
#' \deqn{\ell^{(\mu\mu\mu)} = 0, \qquad
#'       \ell^{(\mu\mu\tau)} = -1, \qquad
#'       \ell^{(\mu\tau\tau)} = 0, \qquad
#'       \ell^{(\tau\tau\tau)} = \dfrac{1}{\tau^3}.}
#' Every one of them is free of the response, so the observed and the expected
#' values coincide and `expected` selects nothing: the same kernel runs either
#' way and the two results are identical to the bit. `approx` and `nsim` are
#' ignored for the same reason.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected third derivatives being the same numbers. Defaults to
#'   `FALSE`.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `tau` must
#'   be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_tau`,
#'   `mu_tau_tau` and `tau_tau_tau`, each of length `length(y)`. The names
#'   enumerate the distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.Gaussian3Distrib()] for the order below and
#'   [distrib_deriv4.Gaussian3Distrib()] for the order above;
#'   [distrib_deriv3.Gaussian1Distrib()] for the same order in the standard
#'   deviation, where two components do carry the residual; and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, tau = 0.25)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # Two of the four are non-zero, and both are constants.
#' lapply(d3, unique)
#' 1 / 0.25^3
#'
#' # Free of the response, so asking for the expectation changes nothing.
#' identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 1, tau = 0.25 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 1, tau = 0.25 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_tau, tolerance = 1e-6)
S7::method(distrib_deriv3, Gaussian3Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  gaussian3_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Fourth-Order Derivatives in Mean and Precision
#' @name distrib_deriv4.Gaussian3Distrib
#' @description
#' Computes the five distinct fourth derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{\tau}, in closed form. Four of the five
#' are exactly zero and the fifth is a constant:
#' \deqn{\ell^{(\tau\tau\tau\tau)} = -\dfrac{3}{\tau^4}.}
#' This is the flattest of the three parametrizations of the Gaussian at fourth
#' order; the log-density is quadratic in \eqn{\mu} and its \eqn{\tau} part is
#' \eqn{\tfrac{1}{2}\log\tau} plus a term linear in \eqn{\tau}, so every mixed
#' component past \eqn{\ell^{(\mu\mu\tau)}} vanishes. Like the third order,
#' the values are free of the response, so `expected`, `approx` and `nsim` are
#' all without effect.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected fourth derivatives being the same numbers. Defaults to
#'   `FALSE`.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `tau` must
#'   be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_tau`, `mu_mu_tau_tau`, `mu_tau_tau_tau` and `tau_tau_tau_tau`,
#'   each of length `length(y)`. Only the last is non-zero.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.Gaussian3Distrib()] for the order below,
#'   [distrib_deriv4.Gaussian1Distrib()] for the same order in the standard
#'   deviation, where three of the five are non-zero, and [distrib_deriv4()]
#'   for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, tau = 0.25)
#' d4 <- distrib_deriv4(d, y, th)
#'
#' # One non-zero component out of five, at -3/tau^4.
#' lapply(d4, unique)
#' -3 / 0.25^4
#'
#' # Free of the response, so asking for the expectation changes nothing.
#' identical(d4, distrib_deriv4(d, y, th, expected = TRUE))
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 1, tau = 0.25 + eps))$tau_tau_tau
#' dn <- distrib_deriv3(d, y, list(mu = 1, tau = 0.25 - eps))$tau_tau_tau
#' all.equal((up - dn) / (2 * eps), d4$tau_tau_tau_tau, tolerance = 1e-4)
S7::method(distrib_deriv4, Gaussian3Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  gaussian3_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian First Derivative in the Response, Mean and Precision
#' @name distrib_grad_y.Gaussian3Distrib
#' @description
#' Computes the first derivative of the Gaussian log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\tau(y - \mu),}
#' in closed form. The Gaussian is a location family in \eqn{\mu}, so the
#' response enters the log-density only through \eqn{y - \mu} and this
#' derivative is the negative of the score in \eqn{\mu}. Quantile residuals and
#' the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `tau` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(tau))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Gaussian3Distrib()] for the second derivative in
#'   the response, [distrib_gradient.Gaussian3Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic and the
#'   finite-difference fallback a family without a closed form takes.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, tau = 0.25)
#'
#' all.equal(distrib_grad_y(d, y, th), -0.25 * (y - 1))
#'
#' # A location family: the derivative in the response is minus the score in
#' # the location.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-5
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, Gaussian3Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * (y - theta[[1]])
}

#' @title Gaussian Second Derivative in the Response, Mean and Precision
#' @name distrib_hess_y.Gaussian3Distrib
#' @description
#' Computes the second derivative of the Gaussian log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\tau,}
#' in closed form. It does not depend on \eqn{y} or on \eqn{\mu}, so the value
#' is constant within a parameter setting and is recycled to the length of `y`.
#' Being a location family, the Gaussian has the same curvature in the response
#' as in its location, and this equals the `mu_mu` component of
#' [distrib_hessian.Gaussian3Distrib()]. In this parametrization that curvature
#' is the parameter itself, which is the sense in which the precision measures
#' how sharply the log-density peaks.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `tau`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `tau` must
#'   be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `length(y)`, every entry \eqn{-\tau}.
#'
#' @seealso [distrib_grad_y.Gaussian3Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Gaussian3Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- gaussian3_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, tau = 0.25)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same curvature as in the location, and here that
#' # curvature is minus the parameter.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Negative everywhere, so the log-density is concave in the response.
#' all(distrib_hess_y(d, y, th) < 0)
S7::method(distrib_hess_y, Gaussian3Distrib) <- function(distrib, y, theta, ...) {
  rep(-theta[[2]], length.out = length(y))
}


#' Gaussian Distribution, Mean and Precision
#'
#' @description
#' Builds the distribution object for the Gaussian (normal) family parametrized
#' by its mean \eqn{\mu} and its precision \eqn{\tau = 1/\sigma^2 > 0}. The
#' returned object carries closed-form derivatives of the log-density to fourth
#' order, in the parameters and in the response, and closed-form moments, so
#' every generic of the toolkit answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. The defaults are the identity for
#' the mean, which is already free, and the logarithm for the precision, which
#' keeps it positive at every predictor.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the mean ranging
#'   over the whole line already.
#' @param link_tau A `link` object from `linkfunctions7` for the precision
#'   \eqn{\tau}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \tau) = \sqrt{\dfrac{\tau}{2\pi}}
#'       \exp\left\{-\dfrac{\tau(y-\mu)^{2}}{2}\right\},}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\tau \in (0, \infty)}. The
#' mean is \eqn{\mu}, the variance \eqn{1/\tau}, and both the skewness and the
#' excess kurtosis are 0.
#'
#' This is the same law as [gaussian1_distrib()] in different coordinates, and
#' a separate family for the reason [gaussian2_distrib()] is: the parameter
#' here *is* the precision, and that is what the estimate, the standard error
#' and the interval describe. The precision is the parametrization a Bayesian
#' conjugate analysis uses, the gamma being conjugate for \eqn{\tau} at known
#' \eqn{\mu}.
#'
#' # Derivatives
#'
#' Writing \eqn{r = y - \mu}, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \tau r, \qquad
#'       \dfrac{\partial \ell}{\partial \tau} = \dfrac{1}{2\tau} - \dfrac{r^2}{2},}
#' the observed Hessian
#' \deqn{\ell^{(\mu\mu)} = -\tau, \quad
#'       \ell^{(\mu\tau)} = r, \quad
#'       \ell^{(\tau\tau)} = -\dfrac{1}{2\tau^2},}
#' and its expectation
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\tau, \quad
#'       \mathbb{E}\left[\ell^{(\mu\tau)}\right] = 0, \quad
#'       \mathbb{E}\left[\ell^{(\tau\tau)}\right] = -\dfrac{1}{2\tau^2}.}
#' The zero off-diagonal makes the mean and the precision orthogonal, so their
#' maximum likelihood estimates are asymptotically independent.
#'
#' It is the flattest of the three parametrizations. The log-density is
#' quadratic in \eqn{\mu} and, in \eqn{\tau}, is \eqn{\tfrac{1}{2}\log\tau}
#' plus a term linear in \eqn{\tau}, so every third and fourth derivative is
#' free of the response: at third order only \eqn{\ell^{(\mu\mu\tau)} = -1} and
#' \eqn{\ell^{(\tau\tau\tau)} = 1/\tau^3} survive, and at fourth order only
#' \eqn{\ell^{(\tau\tau\tau\tau)} = -3/\tau^4}. Asking those methods for the
#' expectation returns the same numbers.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale and reaches
#' the closed-form estimates \eqn{\hat\mu = \bar y} and
#' \eqn{\hat\tau = 1 / \{n^{-1}\sum (y_i - \bar y)^2\}}, the reciprocal of the
#' maximum likelihood variance with divisor \eqn{n}. The example below checks
#' both against the sample.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\tau = 1/\sigma^2 > 0} the precision. \eqn{\ell^{(ij)}} is a second
#' derivative of \eqn{\ell} in parameters \eqn{i} and \eqn{j}. \eqn{\eta} is a
#' parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Gaussian3Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"gaussian3"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "tau")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [gaussian1_distrib()] and [gaussian2_distrib()] for the same law in the
#' standard deviation and in the variance; [gamma1_distrib()] for the conjugate
#' prior of \eqn{\tau} at known \eqn{\mu}; [student_t1_distrib()] and
#' [laplace_distrib()] for heavier tails; [fit_distrib()] to estimate the
#' parameters; [check_distrib()] to validate a family of your own against the
#' same battery this one passes; [Gaussian3Distrib] for the class.
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. and
#' Rubin, D. B. (2013). *Bayesian Data Analysis*, 3rd edition, Chapter 2.
#' Chapman and Hall/CRC, Boca Raton.
#'
#' @examples
#' d <- gaussian3_distrib()
#' d
#'
#' # The same law as gaussian1 at sigma = 1/sqrt(tau).
#' y <- c(-1.2, 0.3, 2.5)
#' all.equal(distrib_pdf(d, y, list(mu = 1, tau = 0.25)),
#'           distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#'
#' # Moments in closed form: the variance is the reciprocal of the parameter.
#' th <- list(mu = 1, tau = 0.25)
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # Fitting recovers the closed-form maximum likelihood estimates.
#' set.seed(7)
#' z <- distrib_rng(d, 400, list(mu = 3, tau = 0.25))
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = mean(z), tau = 1 / mean((z - mean(z))^2)))
#'
#' # The estimate is a precision, so the interval it reports is an interval for
#' # the precision and stays positive.
#' confint(fit)
#'
#' @export
gaussian3_distrib <- function(link_mu = identity_link(),
                              link_tau = log_link()) {
  Gaussian3Distrib(
    distrib_name = "gaussian3",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "tau"),
    params_interpretation = c(mu = "mean", tau = "precision"),
    n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), tau = c(0, Inf)),
    link_params = list(mu = link_mu, tau = link_tau)
  )
}
