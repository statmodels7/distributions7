#' @include distrib.R generics.R
NULL

#' @title Gaussian Distribution Class, Mean and Variance
#' @name Gaussian2Distrib
#'
#' @description
#' The S7 class of the Gaussian (normal) family parametrized by its mean
#' \eqn{\mu} and its variance \eqn{\sigma^2 > 0}, with density
#' \eqn{f(y) = (2\pi\sigma^2)^{-1/2}\exp\{-(y-\mu)^2/(2\sigma^2)\}} on the
#' whole real line. It inherits from `continuous_distrib`, so it answers every
#' generic of the `distrib` contract; the eleven methods listed below are
#' registered on it directly and everything else comes from the parent.
#'
#' Build one with [gaussian2_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Gaussian2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [gaussian2_distrib()] they hold `"gaussian2"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma2")`, the interpretations
#'   `c(mu = "mean", sigma2 = "variance")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [gaussian2_distrib()] to build one;
#'   [gaussian1_distrib()] for the same law in mean and standard deviation and
#'   [gaussian3_distrib()] for mean and precision;
#'   [distrib_pdf.Gaussian2Distrib()] and [distrib_gradient.Gaussian2Distrib()]
#'   for the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Gaussian2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gaussian2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gaussian2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gaussian2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gaussian2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gaussian2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gaussian2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gaussian2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gaussian2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gaussian2Distrib],
#'   [`distrib_rng()`][distrib_rng.Gaussian2Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- gaussian2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@params_bounds
#'
#' # The second parameter is the variance, so the second interpretation and
#' # the fitted standard error both describe the variance.
#' d@params_interpretation[["sigma2"]]
Gaussian2Distrib <- S7::new_class("Gaussian2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Probability Density Function in Mean and Variance
#' @name distrib_pdf.Gaussian2Distrib
#' @description
#' Computes the Gaussian density
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{\sqrt{2\pi\sigma^2}}
#'       \exp\left\{-\dfrac{(y-\mu)^2}{2\sigma^2}\right\}}
#' by calling [stats::dnorm()] at `mean = mu` and `sd = sqrt(sigma2)`, so the
#' accuracy and the underflow behavior are R's own. With `log = TRUE` the
#' logarithm is formed inside `dnorm()` and stays finite far into the tails,
#' where the density itself underflows to zero.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive; a zero or negative value
#'   makes `sqrt(sigma2)` zero or `NaN` and [stats::dnorm()] warns.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_cdf.Gaussian2Distrib()] for the distribution function,
#'   [distrib_gradient.Gaussian2Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf.Gaussian1Distrib()] for the same density in the
#'   standard deviation, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#'
#' # The method is stats::dnorm at the square root of the variance.
#' all.equal(distrib_pdf(d, y, list(mu = 1, sigma2 = 4)),
#'           dnorm(y, mean = 1, sd = 2))
#'
#' # Same law as gaussian1 at sigma = sqrt(sigma2), to the last bit.
#' all.equal(distrib_pdf(d, y, list(mu = 1, sigma2 = 4)),
#'           distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0, 1, 2), sigma2 = c(1, 2.25, 4)))
#'
#' # In the far tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 40, list(mu = 0, sigma2 = 1))
#' distrib_pdf(d, 40, list(mu = 0, sigma2 = 1), log = TRUE)
S7::method(distrib_pdf, Gaussian2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(y, mean = theta[[1]], sd = sqrt(theta[[2]]), log = log)
}

#' @title Gaussian Cumulative Distribution Function in Mean and Variance
#' @name distrib_cdf.Gaussian2Distrib
#' @description
#' Computes the Gaussian distribution function
#' \deqn{F(q; \mu, \sigma^2) = \Phi\left(\dfrac{q-\mu}{\sigma}\right),
#'       \qquad \sigma = \sqrt{\sigma^2},}
#' with \eqn{\Phi} the standard normal distribution function, by calling
#' [stats::pnorm()]. Both tails are available exactly: `lower.tail = FALSE`
#' evaluates \eqn{1 - F} without forming the difference, and `log.p = TRUE`
#' returns a logarithm that stays finite where the probability itself
#' underflows to zero.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma2))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Gaussian2Distrib()] for the inverse,
#'   [distrib_pdf.Gaussian2Distrib()] for the density,
#'   [distrib_grad_cdf()] for the derivatives of this function in the
#'   parameters, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # The method is stats::pnorm at sd = sqrt(sigma2).
#' all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
#'           pnorm(c(-1.2, 0.3, 2.5), mean = 1, sd = 2))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
#'
#' # Forty standard deviations out the upper tail underflows; its log does not.
#' distrib_cdf(d, 40, list(mu = 0, sigma2 = 1), lower.tail = FALSE)
#' distrib_cdf(d, 40, list(mu = 0, sigma2 = 1), lower.tail = FALSE,
#'             log.p = TRUE)
S7::method(distrib_cdf, Gaussian2Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  stats::pnorm(q, mean = theta[[1]], sd = sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Quantile Function in Mean and Variance
#' @name distrib_quantile.Gaussian2Distrib
#' @description
#' Computes the Gaussian quantile function
#' \deqn{Q(p; \mu, \sigma^2) = \mu + \sqrt{\sigma^2}\,\Phi^{-1}(p)}
#' by calling [stats::qnorm()]. The Gaussian distribution function is strictly
#' increasing on the whole line, so \eqn{Q} is its exact inverse and the round
#' trip through [distrib_cdf.Gaussian2Distrib()] returns `p`.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; the endpoints give `-Inf` and `Inf`.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities, which is how a quantile deep in a tail is
#'   requested without the probability underflowing. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of length
#'   `max(length(p), length(mu), length(sigma2))`.
#'
#' @seealso [distrib_cdf.Gaussian2Distrib()], which this inverts;
#'   [distrib_rng.Gaussian2Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # The median is mu and the quartiles are symmetric about it.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # A quantile 40 standard deviations out, asked for on the log scale
#' # because the probability itself is zero in double precision.
#' distrib_quantile(d, pnorm(-40, log.p = TRUE), list(mu = 0, sigma2 = 1),
#'                  log.p = TRUE)
S7::method(distrib_quantile, Gaussian2Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  stats::qnorm(p, mean = theta[[1]], sd = sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Random Number Generator in Mean and Variance
#' @name distrib_rng.Gaussian2Distrib
#' @description
#' Draws `n` independent Gaussian variates by calling [stats::rnorm()] at
#' `sd = sqrt(sigma2)`, so the draws come from R's own normal generator and
#' depend on `.Random.seed` in the usual way. The inverse-transform fallback
#' the base class supplies is bypassed.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting.
#'   `sigma2` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.Gaussian2Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#'
#' # Same generator as stats::rnorm at sd = sqrt(sigma2).
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = c(0, 1, 2), sigma2 = c(1, 2.25, 4)))
#' set.seed(2)
#' identical(a, rnorm(3, mean = c(0, 1, 2), sd = c(1, 1.5, 2)))
#'
#' # The sample moments recover the parameters, the variance directly.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, sigma2 = 4))
#' c(mean = mean(z), var = var(z))
S7::method(distrib_rng, Gaussian2Distrib) <- function(distrib, n, theta, ...) {
  stats::rnorm(n, mean = theta[[1]], sd = sqrt(theta[[2]]))
}

#' @title Gaussian Score in Mean and Variance
#' @name distrib_gradient.Gaussian2Distrib
#' @description
#' Computes the first derivatives of the Gaussian log-density with respect to
#' \eqn{\mu} and \eqn{v = \sigma^2}, one value per observation, in closed form.
#' With \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{v},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial v} = \dfrac{r^2 - v}{2v^2}.}
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning, giving \eqn{\partial \ell / \partial \eta_j
#' = h_j'(\eta_j)\, \partial \ell / \partial \theta_j}. This method always
#' returns the parameter scale; the transformation happens in the generic.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma2`, each of
#'   length `max(length(y), length(mu), length(sigma2))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{v = \sigma^2 > 0} the variance. \eqn{\eta_j} is the coordinate of
#' parameter \eqn{j} on the unconstrained scale of its link, and \eqn{h_j' =
#' \partial \theta_j / \partial \eta_j} the chain-rule factor onto it.
#'
#' @seealso [distrib_hessian.Gaussian2Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Gaussian2Distrib()] for their expectation,
#'   [distrib_gradient.Gaussian1Distrib()] for the same score in the standard
#'   deviation, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, sigma2 = 4)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' r <- y - 1
#' all.equal(g$mu, r / 4)
#' all.equal(g$sigma2, (r^2 - 4) / (2 * 4^2))
#'
#' # The summed score vanishes at the maximum likelihood estimate, where the
#' # variance carries the divisor n.
#' set.seed(1)
#' z <- rnorm(200, mean = 3, sd = 2)
#' mle <- list(mu = mean(z), sigma2 = mean((z - mean(z))^2))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#'
#' # Both parametrizations ride a log link on the spread, and log(sigma^2) is
#' # twice log(sigma), so the link-scale score here is half gaussian1's.
#' distrib_gradient(d, y, th, scale = "link")$sigma2 /
#'   distrib_gradient(gaussian1_distrib(), y, list(mu = 1, sigma = 2),
#'                    scale = "link")$sigma
S7::method(distrib_gradient, Gaussian2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian2_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Observed Hessian in Mean and Variance
#' @name distrib_hessian.Gaussian2Distrib
#' @description
#' Computes the three distinct second derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{v = \sigma^2}, one value per observation,
#' in closed form. With \eqn{r = y - \mu},
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{1}{v}, \qquad
#'       \ell^{(\mu v)} = -\dfrac{r}{v^2}, \qquad
#'       \ell^{(vv)} = \dfrac{1}{2v^2} - \dfrac{r^2}{v^3}.}
#' Only the curvature in \eqn{\mu} is free of the data; the other two vary with
#' the residual, and their expectations are
#' [distrib_expected_hessian.Gaussian2Distrib()].
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_sigma2` and
#'   `sigma2_sigma2`, each of length
#'   `max(length(y), length(mu), length(sigma2))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density with respect
#' to parameters \eqn{i} and \eqn{j}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_gradient.Gaussian2Distrib()] for the score,
#'   [distrib_expected_hessian.Gaussian2Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Gaussian2Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, sigma2 = 4)
#' h <- distrib_hessian(d, y, th)
#'
#' # The curvature in mu is constant at -1/v; the other two are not.
#' h$mu_mu
#' r <- y - 1
#' all.equal(h$mu_sigma2, -r / 4^2)
#' all.equal(h$sigma2_sigma2, 1 / (2 * 4^2) - r^2 / 4^3)
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-5
#' up <- distrib_gradient(d, y, list(mu = 1, sigma2 = 4 + eps))$sigma2
#' dn <- distrib_gradient(d, y, list(mu = 1, sigma2 = 4 - eps))$sigma2
#' all.equal((up - dn) / (2 * eps), h$sigma2_sigma2, tolerance = 1e-6)
S7::method(distrib_hessian, Gaussian2Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian2_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Expected Hessian in Mean and Variance
#' @name distrib_expected_hessian.Gaussian2Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. With \eqn{v = \sigma^2},
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{v}, \qquad
#'       \mathbb{E}\left[\ell^{(\mu v)}\right] = 0, \qquad
#'       \mathbb{E}\left[\ell^{(vv)}\right] = -\dfrac{1}{2v^2}.}
#' They follow from \eqn{\mathbb{E}[(Y-\mu)^2] = v} and
#' \eqn{\mathbb{E}[Y-\mu] = 0}. The negative of this matrix is the Fisher
#' information for one observation; the zero off-diagonal says the mean and the
#' variance are orthogonal, so their estimates are asymptotically independent.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
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
#' @return A named list of three numeric vectors, `mu_mu`, `mu_sigma2` and
#'   `sigma2_sigma2`, each of length
#'   `max(length(y), length(mu), length(sigma2))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The Gaussian is
#' a regular family, so the second Bartlett identity holds and this equals the
#' variance of the score.
#'
#' @seealso [distrib_hessian.Gaussian2Distrib()] for the observed quantity this
#'   is the expectation of, [distrib_expected_hessian.Gaussian1Distrib()] for
#'   the same information in the standard deviation, [fisher_scoring()], which
#'   inverts it at each step, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' th <- list(mu = 1, sigma2 = 4)
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
#' # The information transforms by the delta method: with dv/dsigma = 2 sigma,
#' # the variance entry carries onto gaussian1's -2/sigma^2.
#' distrib_expected_hessian(d, 0, th)$sigma2_sigma2 * (2 * 2)^2
#' distrib_expected_hessian(gaussian1_distrib(), 0,
#'                          list(mu = 1, sigma = 2))$sigma_sigma
S7::method(distrib_expected_hessian, Gaussian2Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ..., threads = 1L) {
  gaussian2_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Third-Order Derivatives in Mean and Variance
#' @name distrib_deriv3.Gaussian2Distrib
#' @description
#' Computes the four distinct third derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{v = \sigma^2}, in closed form. Writing
#' \eqn{r = y - \mu},
#' \deqn{\ell^{(\mu\mu\mu)} = 0, \qquad
#'       \ell^{(\mu\mu v)} = \dfrac{1}{v^2}, \qquad
#'       \ell^{(\mu v v)} = \dfrac{2r}{v^3}, \qquad
#'       \ell^{(vvv)} = \dfrac{3r^2}{v^4} - \dfrac{1}{v^3}.}
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{r} with 0 and \eqn{r^2} with \eqn{v}: the component odd in \eqn{r}
#' vanishes and \eqn{\ell^{(vvv)}} becomes \eqn{2/v^3}. Both routes are closed
#' form, so no quadrature is run and `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
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
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma2`,
#'   `mu_sigma2_sigma2` and `sigma2_sigma2_sigma2`, each of length
#'   `max(length(y), length(mu), length(sigma2))`. The names enumerate the
#'   distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.Gaussian2Distrib()] for the order below and
#'   [distrib_deriv4.Gaussian2Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic and for the numerical route a family
#'   without a closed form takes.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, sigma2 = 4)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # The log-density is quadratic in mu, so the third derivative there is 0.
#' d3$mu_mu_mu
#'
#' # The other three, written out.
#' r <- y - 1
#' all.equal(d3$mu_mu_sigma2, rep(1 / 4^2, 3))
#' all.equal(d3$mu_sigma2_sigma2, 2 * r / 4^3)
#' all.equal(d3$sigma2_sigma2_sigma2, 3 * r^2 / 4^4 - 1 / 4^3)
#'
#' # Expected values: the component odd in the residual vanishes and the pure
#' # variance component becomes 2/v^3.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#' 2 / 4^3
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 1, sigma2 = 4 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 1, sigma2 = 4 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_sigma2, tolerance = 1e-6)
S7::method(distrib_deriv3, Gaussian2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian2_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian2_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian Fourth-Order Derivatives in Mean and Variance
#' @name distrib_deriv4.Gaussian2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{v = \sigma^2}, in closed form. Writing
#' \eqn{r = y - \mu},
#' \deqn{\ell^{(\mu\mu\mu\mu)} = \ell^{(\mu\mu\mu v)} = 0, \qquad
#'       \ell^{(\mu\mu vv)} = -\dfrac{2}{v^3}, \qquad
#'       \ell^{(\mu vvv)} = -\dfrac{6r}{v^4}, \qquad
#'       \ell^{(vvvv)} = \dfrac{3}{v^4} - \dfrac{12 r^2}{v^5}.}
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{r} with 0 and \eqn{r^2} with \eqn{v}, which leaves \eqn{-9/v^4} in the
#' last component and zero in the two odd ones. Both routes are closed form, so
#' `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
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
#'   `mu_mu_mu_sigma2`, `mu_mu_sigma2_sigma2`, `mu_sigma2_sigma2_sigma2` and
#'   `sigma2_sigma2_sigma2_sigma2`, each of length
#'   `max(length(y), length(mu), length(sigma2))`.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.Gaussian2Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.Gaussian2Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, sigma2 = 4)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Quadratic in mu, so the two components with three or more mu are zero.
#' c(d4$mu_mu_mu_mu[1], d4$mu_mu_mu_sigma2[1])
#'
#' # The mixed second-second component is constant at -2/v^3.
#' all.equal(d4$mu_mu_sigma2_sigma2, rep(-2 / 4^3, 3))
#'
#' # Expected values: -9/v^4 in the pure variance component.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#' -9 / 4^4
S7::method(distrib_deriv4, Gaussian2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian2_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian2_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian First Derivative in the Response, Mean and Variance
#' @name distrib_grad_y.Gaussian2Distrib
#' @description
#' Computes the first derivative of the Gaussian log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{y - \mu}{v},
#'       \qquad v = \sigma^2,}
#' in closed form. The Gaussian is a location family in \eqn{\mu}, so the
#' response enters the log-density only through \eqn{y - \mu} and this
#' derivative is the negative of the score in \eqn{\mu}. Quantile residuals and
#' the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Gaussian2Distrib()] for the second derivative in
#'   the response, [distrib_gradient.Gaussian2Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic and the
#'   finite-difference fallback a family without a closed form takes.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, sigma2 = 4)
#'
#' all.equal(distrib_grad_y(d, y, th), -(y - 1) / 4)
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
S7::method(distrib_grad_y, Gaussian2Distrib) <- function(distrib, y, theta, ...) {
  -(y - theta[[1]]) / theta[[2]]
}

#' @title Gaussian Second Derivative in the Response, Mean and Variance
#' @name distrib_hess_y.Gaussian2Distrib
#' @description
#' Computes the second derivative of the Gaussian log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{1}{v},
#'       \qquad v = \sigma^2,}
#' in closed form. It does not depend on \eqn{y} or on \eqn{\mu}, so the value
#' is constant within a parameter setting and is recycled to the length of `y`.
#' Being a location family, the Gaussian has the same curvature in the response
#' as in its location, and this equals the `mu_mu` component of
#' [distrib_hessian.Gaussian2Distrib()].
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `sigma2`
#'   must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `length(y)`, every entry \eqn{-1/v}.
#'
#' @seealso [distrib_grad_y.Gaussian2Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Gaussian2Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- gaussian2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 1, sigma2 = 4)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same curvature as in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Negative everywhere, so the log-density is concave in the response.
#' all(distrib_hess_y(d, y, th) < 0)
S7::method(distrib_hess_y, Gaussian2Distrib) <- function(distrib, y, theta, ...) {
  rep(-1 / theta[[2]], length.out = length(y))
}


#' Gaussian Distribution, Mean and Variance
#'
#' @description
#' Builds the distribution object for the Gaussian (normal) family parametrized
#' by its mean \eqn{\mu} and its variance \eqn{\sigma^2 > 0}. The returned
#' object carries closed-form derivatives of the log-density to fourth order,
#' in the parameters and in the response, and closed-form moments, so every
#' generic of the toolkit answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. The defaults are the identity for
#' the mean, which is already free, and the logarithm for the variance, which
#' keeps it positive at every predictor.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the mean ranging
#'   over the whole line already.
#' @param link_sigma2 A `link` object from `linkfunctions7` for the variance
#'   \eqn{\sigma^2}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{\sqrt{2\pi\sigma^{2}}}
#'       \exp\left\{-\dfrac{(y-\mu)^{2}}{2\sigma^{2}}\right\},}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\sigma^2 \in (0, \infty)}. The
#' mean is \eqn{\mu}, the variance \eqn{\sigma^2}, and both the skewness and
#' the excess kurtosis are 0. The numbering follows the literature where it has
#' one: this parametrization is `NO2` in the gamlss family catalog.
#'
#' This is the same law as [gaussian1_distrib()] in different coordinates,
#' \eqn{\sigma^2} here being the square of the \eqn{\sigma} there. The two are
#' separate families and not one family under a link, because a link changes
#' the scale a parameter is *modeled* on and leaves the parameter what it was,
#' while here the parameter, its interpretation, its standard error and its
#' confidence interval are all about the variance.
#'
#' # Derivatives
#'
#' Writing \eqn{r = y - \mu} and \eqn{v = \sigma^2}, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{v}, \qquad
#'       \dfrac{\partial \ell}{\partial v} = \dfrac{r^2 - v}{2v^2},}
#' the observed Hessian
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{1}{v}, \quad
#'       \ell^{(\mu v)} = -\dfrac{r}{v^2}, \quad
#'       \ell^{(vv)} = \dfrac{1}{2v^2} - \dfrac{r^2}{v^3},}
#' and its expectation
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{v}, \quad
#'       \mathbb{E}\left[\ell^{(\mu v)}\right] = 0, \quad
#'       \mathbb{E}\left[\ell^{(vv)}\right] = -\dfrac{1}{2v^2}.}
#' The zero off-diagonal makes the mean and the variance orthogonal, so their
#' maximum likelihood estimates are asymptotically independent. Orthogonality
#' holds in all three parametrizations of this family, the mean being
#' orthogonal to any smooth function of the spread alone.
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.Gaussian2Distrib()] and [distrib_deriv4.Gaussian2Distrib()],
#' as are the derivatives in the response,
#' [distrib_grad_y.Gaussian2Distrib()] and [distrib_hess_y.Gaussian2Distrib()].
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale and reaches
#' the closed-form estimates \eqn{\hat\mu = \bar y} and
#' \eqn{\hat\sigma^2 = n^{-1}\sum (y_i - \bar y)^2}, the maximum likelihood
#' estimate with divisor \eqn{n}. The example below checks both against the
#' sample.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\sigma^2 > 0} the variance. \eqn{\ell^{(ij)}} is a second derivative of
#' \eqn{\ell} in parameters \eqn{i} and \eqn{j}. \eqn{\eta} is a parameter on
#' the unconstrained scale of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Gaussian2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"gaussian2"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma2")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [gaussian1_distrib()] and [gaussian3_distrib()] for the same law in the
#' standard deviation and in the precision; [lognormal1_distrib()] for the
#' Gaussian on a log scale; [student_t1_distrib()] and [laplace_distrib()] for
#' heavier tails; [fit_distrib()] to estimate the parameters;
#' [check_distrib()] to validate a family of your own against the same battery
#' this one passes; [Gaussian2Distrib] for the class.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. *Journal of the Royal Statistical Society:
#' Series C* **54**, 507-554.
#'
#' @examples
#' d <- gaussian2_distrib()
#' d
#'
#' # The same law as gaussian1 at sigma = sqrt(sigma2).
#' y <- c(-1.2, 0.3, 2.5)
#' all.equal(distrib_pdf(d, y, list(mu = 1, sigma2 = 4)),
#'           distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#'
#' # Moments in closed form: the variance is the parameter itself.
#' th <- list(mu = 1, sigma2 = 4)
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # Fitting recovers the closed-form maximum likelihood estimates, with the
#' # variance divided by n.
#' set.seed(7)
#' z <- distrib_rng(d, 400, list(mu = 3, sigma2 = 4))
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = mean(z), sigma2 = mean((z - mean(z))^2)))
#'
#' # The estimate is a variance, so the interval it reports is an interval for
#' # the variance and stays positive.
#' confint(fit)
#'
#' @export
gaussian2_distrib <- function(link_mu = identity_link(),
                              link_sigma2 = log_link()) {
  Gaussian2Distrib(
    distrib_name = "gaussian2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma2"),
    params_interpretation = c(mu = "mean", sigma2 = "variance"),
    n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
    link_params = list(mu = link_mu, sigma2 = link_sigma2)
  )
}
