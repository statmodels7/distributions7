#' @include distrib.R generics.R
NULL

#' @title Gaussian Distribution Class, Mean and Standard Deviation
#' @name Gaussian1Distrib
#'
#' @description
#' The S7 class of the Gaussian (normal) family parametrized by its mean
#' \eqn{\mu} and its standard deviation \eqn{\sigma > 0}, with density
#' \eqn{f(y) = (2\pi)^{-1/2}\sigma^{-1}\exp\{-(y-\mu)^2/(2\sigma^2)\}} on the
#' whole real line. It inherits from `continuous_distrib`, so it answers every
#' generic of the `distrib` contract; the eleven methods listed below are
#' registered on it directly and everything else comes from the parent.
#'
#' Build one with [gaussian1_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Gaussian1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [gaussian1_distrib()] they hold `"gaussian1"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma")`, the interpretations
#'   `c(mu = "mean", sigma = "standard deviation")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [gaussian1_distrib()] to build one;
#'   [gaussian2_distrib()] for the same law in mean and variance and
#'   [gaussian3_distrib()] for mean and precision;
#'   [distrib_pdf.Gaussian1Distrib()] and [distrib_gradient.Gaussian1Distrib()]
#'   for the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Gaussian1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gaussian1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gaussian1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gaussian1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gaussian1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gaussian1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gaussian1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gaussian1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gaussian1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gaussian1Distrib],
#'   [`distrib_rng()`][distrib_rng.Gaussian1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- gaussian1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_bounds
#' d@params_interpretation
#'
#' # The links carry each parameter to the unconstrained scale an optimizer
#' # moves on: the mean is already free, the standard deviation rides a log.
#' vapply(d@link_params, function(l) l@link_name, character(1))
Gaussian1Distrib <- S7::new_class("Gaussian1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Probability Density Function
#' @name distrib_pdf.Gaussian1Distrib
#' @description
#' Computes the Gaussian density
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\,\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\}}
#' by calling [stats::dnorm()] at `mean = mu` and `sd = sigma`, so the accuracy
#' and the underflow behavior are R's own. With `log = TRUE` the logarithm is
#' formed inside `dnorm()` and stays finite far into the tails, where the
#' density itself underflows to zero.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive; a zero or negative value
#'   gives `NaN` with a warning from [stats::dnorm()].
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.Gaussian1Distrib()] for the distribution function,
#'   [distrib_gradient.Gaussian1Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf()] for the generic and
#'   [gaussian1_distrib()] for the family.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#'
#' # The method is stats::dnorm at this parametrization.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)),
#'           dnorm(y, mean = 0.4, sd = 1.5))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0, 1, 2), sigma = c(1, 1.5, 2)))
#'
#' # In the far tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 40, list(mu = 0, sigma = 1))
#' distrib_pdf(d, 40, list(mu = 0, sigma = 1), log = TRUE)
S7::method(distrib_pdf, Gaussian1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(
    x = y,
    mean = theta[[1]],
    sd = theta[[2]],
    log = log
  )
}

#' @title Gaussian Cumulative Distribution Function
#' @name distrib_cdf.Gaussian1Distrib
#' @description
#' Computes the Gaussian distribution function
#' \deqn{F(q; \mu, \sigma) = \Phi\left(\dfrac{q-\mu}{\sigma}\right)}
#' with \eqn{\Phi} the standard normal distribution function, by calling
#' [stats::pnorm()]. Both tails are available exactly: `lower.tail = FALSE`
#' evaluates \eqn{1 - F} without forming the difference, and `log.p = TRUE`
#' returns a logarithm that stays finite where the probability itself
#' underflows to zero.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
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
#' @seealso [distrib_quantile.Gaussian1Distrib()] for the inverse,
#'   [distrib_pdf.Gaussian1Distrib()] for the density,
#'   [distrib_grad_cdf()] for the derivatives of this function in the
#'   parameters, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' # The method is stats::pnorm at this parametrization.
#' all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
#'           pnorm(c(-1.2, 0.3, 2.5), mean = 0.4, sd = 1.5))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
#'
#' # Forty standard deviations out the upper tail underflows; its log does not.
#' distrib_cdf(d, 40, list(mu = 0, sigma = 1), lower.tail = FALSE)
#' distrib_cdf(d, 40, list(mu = 0, sigma = 1), lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, Gaussian1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pnorm(
    q = q,
    mean = theta[[1]],
    sd = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gaussian Quantile Function
#' @name distrib_quantile.Gaussian1Distrib
#' @description
#' Computes the Gaussian quantile function
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma\,\Phi^{-1}(p)}
#' by calling [stats::qnorm()]. The Gaussian distribution function is strictly
#' increasing on the whole line, so \eqn{Q} is its exact inverse and the round
#' trip through [distrib_cdf.Gaussian1Distrib()] returns `p`.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; the endpoints give `-Inf` and `Inf`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities, which is how a quantile deep in a tail is
#'   requested without the probability underflowing. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles, of length
#'   `max(length(p), length(mu), length(sigma))`.
#'
#' @seealso [distrib_cdf.Gaussian1Distrib()], which this inverts;
#'   [distrib_rng.Gaussian1Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
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
#' distrib_quantile(d, pnorm(-40, log.p = TRUE), list(mu = 0, sigma = 1),
#'                  log.p = TRUE)
S7::method(distrib_quantile, Gaussian1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qnorm(
    p = p,
    mean = theta[[1]],
    sd = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gaussian Random Number Generator
#' @name distrib_rng.Gaussian1Distrib
#' @description
#' Draws `n` independent Gaussian variates by calling [stats::rnorm()], so the
#' draws come from R's own normal generator and depend on `.Random.seed` in the
#' usual way. The inverse-transform fallback the base class supplies is
#' bypassed.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting.
#'   `sigma` must be strictly positive.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.Gaussian1Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#'
#' # Same generator as stats::rnorm, so the same seed gives the same draws.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = c(0, 1, 2), sigma = c(1, 1.5, 2)))
#' set.seed(2)
#' identical(a, rnorm(3, mean = c(0, 1, 2), sd = c(1, 1.5, 2)))
#'
#' # The sample moments recover the parameters.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
#' c(mean = mean(z), sd = sd(z))
S7::method(distrib_rng, Gaussian1Distrib) <- function(distrib, n, theta) {
  stats::rnorm(
    n = n,
    mean = theta[[1]],
    sd = theta[[2]]
  )
}

#' @title Gaussian Score
#' @name distrib_gradient.Gaussian1Distrib
#' @description
#' Computes the first derivatives of the Gaussian log-density with respect to
#' \eqn{\mu} and \eqn{\sigma}, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}.}
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning, giving \eqn{\partial \ell / \partial \eta_j
#' = h_j'(\eta_j)\, \partial \ell / \partial \theta_j}. This method always
#' returns the parameter scale; the transformation happens in the generic.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
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
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\sigma > 0} the standard deviation. \eqn{\eta_j} is the coordinate of
#' parameter \eqn{j} on the unconstrained scale of its link, and \eqn{h_j' =
#' \partial \theta_j / \partial \eta_j} the chain-rule factor onto it.
#'
#' @seealso [distrib_hessian.Gaussian1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Gaussian1Distrib()] for their expectation,
#'   [distrib_grad_y.Gaussian1Distrib()] for the derivative in the response,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' all.equal(g$mu, (y - 0.4) / 1.5^2)
#' all.equal(g$sigma, ((y - 0.4)^2 - 1.5^2) / 1.5^3)
#'
#' # The summed score vanishes at the maximum likelihood estimate.
#' set.seed(1)
#' z <- rnorm(200, mean = 3, sd = 2)
#' mle <- list(mu = mean(z), sigma = sqrt(mean((z - mean(z))^2)))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#'
#' # On the link scale the sigma component is multiplied by h' = sigma,
#' # the derivative of the inverse log link; mu rides the identity and is
#' # unchanged.
#' distrib_gradient(d, y, th, scale = "link")$sigma / g$sigma
S7::method(distrib_gradient, Gaussian1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Observed Hessian
#' @name distrib_hessian.Gaussian1Distrib
#' @description
#' Computes the three distinct second derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, one value per observation, in
#' closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y - \mu)^2}{\sigma^4},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma} = -\dfrac{2(y - \mu)}{\sigma^3}.}
#' Only the curvature in \eqn{\mu} is free of the data; the other two components
#' vary with the residual, and their expectations are
#' [distrib_expected_hessian.Gaussian1Distrib()].
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
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
#'   `mu_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @seealso [distrib_gradient.Gaussian1Distrib()] for the score,
#'   [distrib_expected_hessian.Gaussian1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Gaussian1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' h <- distrib_hessian(d, y, th)
#'
#' # The curvature in mu is constant at -1/sigma^2; the other two are not.
#' h$mu_mu
#' all.equal(h$sigma_sigma, (1.5^2 - 3 * (y - 0.4)^2) / 1.5^4)
#' all.equal(h$mu_sigma, -2 * (y - 0.4) / 1.5^3)
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-5
#' up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
S7::method(distrib_hessian, Gaussian1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Expected Hessian
#' @name distrib_expected_hessian.Gaussian1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}\right] = 0.}
#' They follow from \eqn{\mathbb{E}[(Y-\mu)^2] = \sigma^2} and
#' \eqn{\mathbb{E}[Y-\mu] = 0}. The negative of this matrix is the Fisher
#' information for one observation; the zero off-diagonal says the mean and the
#' standard deviation are orthogonal, so their estimates are asymptotically
#' independent.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
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
#' @return A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
#'   `mu_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The **expected information** is \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]},
#' the expectation of the **observed information** under the model. The Gaussian
#' is a regular family, so the second Bartlett identity holds and this equals
#' the variance of the score.
#'
#' @seealso [distrib_hessian.Gaussian1Distrib()] for the observed quantity this
#'   is the expectation of, [fit_distrib()] and [fisher_scoring()], which invert
#'   it at each step, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#' th <- list(mu = 0.4, sigma = 1.5)
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
#' # The mean and the standard deviation are orthogonal: the mixed entry is 0.
#' distrib_expected_hessian(d, 0, th)$mu_sigma
S7::method(distrib_expected_hessian, Gaussian1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  gaussian_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Third-Order Derivatives
#' @name distrib_deriv3.Gaussian1Distrib
#' @description
#' Computes the four distinct third derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, in closed form. Writing
#' \eqn{r = y - \mu},
#' \deqn{\ell^{(\mu\mu\mu)} = 0, \qquad
#'       \ell^{(\mu\mu\sigma)} = \dfrac{2}{\sigma^3}, \qquad
#'       \ell^{(\mu\sigma\sigma)} = \dfrac{6r}{\sigma^4}, \qquad
#'       \ell^{(\sigma\sigma\sigma)} = -\dfrac{2}{\sigma^3} + \dfrac{12 r^2}{\sigma^5}.}
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{r} with 0 and \eqn{r^2} with \eqn{\sigma^2}: the two components carrying
#' an odd power of \eqn{r} vanish and \eqn{\ell^{(\sigma\sigma\sigma)}} becomes
#' \eqn{10/\sigma^3}. Both routes are closed form, so no quadrature is run and
#' `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
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
#' \eqn{\ell^{(i j k)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives; a subscript on \eqn{\ell} never does.
#'
#' @seealso [distrib_hessian.Gaussian1Distrib()] for the order below and
#'   [distrib_deriv4.Gaussian1Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic and for the numerical route a family
#'   without a closed form takes.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # The log-density is quadratic in mu, so the third derivative there is 0.
#' d3$mu_mu_mu
#'
#' # The other three, written out.
#' all.equal(d3$mu_mu_sigma, rep(2 / 1.5^3, 3))
#' all.equal(d3$mu_sigma_sigma, 6 * (y - 0.4) / 1.5^4)
#'
#' # Expected values: the components odd in the residual vanish.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_sigma, tolerance = 1e-6)
S7::method(distrib_deriv3, Gaussian1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian Fourth-Order Derivatives
#' @name distrib_deriv4.Gaussian1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the Gaussian log-density
#' with respect to \eqn{\mu} and \eqn{\sigma}, in closed form. Writing
#' \eqn{r = y - \mu},
#' \deqn{\ell^{(\mu\mu\mu\mu)} = \ell^{(\mu\mu\mu\sigma)} = 0, \qquad
#'       \ell^{(\mu\mu\sigma\sigma)} = -\dfrac{6}{\sigma^4}, \qquad
#'       \ell^{(\mu\sigma\sigma\sigma)} = -\dfrac{24 r}{\sigma^5}, \qquad
#'       \ell^{(\sigma\sigma\sigma\sigma)} = \dfrac{6}{\sigma^4} - \dfrac{60 r^2}{\sigma^6}.}
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{r} with 0 and \eqn{r^2} with \eqn{\sigma^2}, which leaves
#' \eqn{-54/\sigma^4} in the last component and zero in the two odd ones. Both
#' routes are closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
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
#' @seealso [distrib_deriv3.Gaussian1Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.Gaussian1Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Quadratic in mu, so the two components with three or more mu are zero.
#' c(d4$mu_mu_mu_mu[1], d4$mu_mu_mu_sigma[1])
#'
#' # The mixed second-second component is constant at -6/sigma^4.
#' all.equal(d4$mu_mu_sigma_sigma, rep(-6 / 1.5^4, 3))
#'
#' # Expected values: -54/sigma^4 in the pure sigma component.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#' -54 / 1.5^4
S7::method(distrib_deriv4, Gaussian1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian First Derivative in the Response
#' @name distrib_grad_y.Gaussian1Distrib
#' @description
#' Computes the first derivative of the Gaussian log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{y - \mu}{\sigma^2},}
#' in closed form. The Gaussian is a location family in \eqn{\mu}, so the
#' response enters the log-density only through \eqn{y - \mu} and this
#' derivative is the negative of the score in \eqn{\mu}. Quantile residuals and
#' the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Gaussian1Distrib()] for the second derivative in
#'   the response, [distrib_gradient.Gaussian1Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic and the finite-difference
#'   fallback a family without a closed form takes.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' all.equal(distrib_grad_y(d, y, th), -(y - 0.4) / 1.5^2)
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
S7::method(distrib_grad_y, Gaussian1Distrib) <- function(distrib, y, theta) {
  -(y - theta[[1]]) / theta[[2]]^2
}

#' @title Gaussian Second Derivative in the Response
#' @name distrib_hess_y.Gaussian1Distrib
#' @description
#' Computes the second derivative of the Gaussian log-density with respect to
#' the response,
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{1}{\sigma^2},}
#' in closed form. It does not depend on \eqn{y} or on \eqn{\mu}, so the value
#' is constant within a parameter setting and is recycled to the length of `y`.
#' Being a location family, the Gaussian has the same curvature in the response
#' as in its location, and this equals the `mu_mu` component of
#' [distrib_hessian.Gaussian1Distrib()].
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `sigma` must
#'   be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `length(y)`, every entry
#'   \eqn{-1/\sigma^2}.
#'
#' @seealso [distrib_grad_y.Gaussian1Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Gaussian1Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same curvature as in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Negative everywhere, so the log-density is concave in the response.
#' all(distrib_hess_y(d, y, th) < 0)
S7::method(distrib_hess_y, Gaussian1Distrib) <- function(distrib, y, theta) {
  rep(-1 / theta[[2]]^2, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Gaussian Distribution, Mean and Standard Deviation
#'
#' @description
#' Builds the distribution object for the Gaussian (normal) family parametrized
#' by its mean \eqn{\mu} and its standard deviation \eqn{\sigma > 0}. The
#' returned object carries closed-form derivatives of the log-density to fourth
#' order, in the parameters and in the response, and closed-form moments, so
#' every generic of the toolkit answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. The defaults are the identity for
#' the mean, which is already free, and the logarithm for the standard
#' deviation, which keeps it positive at every predictor.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the mean ranging
#'   over the whole line already.
#' @param link_sigma A `link` object from `linkfunctions7` for the standard
#'   deviation \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which
#'   maps \eqn{(0, \infty)} onto the line and so keeps every fitted value
#'   positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\,\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\},}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\sigma \in (0, \infty)}. The
#' distribution function is \eqn{F(q) = \Phi((q-\mu)/\sigma)} and the quantile
#' function \eqn{Q(p) = \mu + \sigma \Phi^{-1}(p)}, with \eqn{\Phi} the standard
#' normal distribution function.
#'
#' The mean is \eqn{\mu}, the variance \eqn{\sigma^2}, and both the skewness and
#' the excess kurtosis are 0. Two sibling parametrizations of the same law are
#' [gaussian2_distrib()], in the mean and the variance, and
#' [gaussian3_distrib()], in the mean and the precision.
#'
#' # Derivatives
#'
#' The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3},}
#' the observed Hessian
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y-\mu)^2}{\sigma^4}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma} = -\dfrac{2(y-\mu)}{\sigma^3},}
#' and its expectation
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0.}
#' The zero off-diagonal makes the mean and the standard deviation orthogonal,
#' so their maximum likelihood estimates are asymptotically independent.
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.Gaussian1Distrib()] and [distrib_deriv4.Gaussian1Distrib()],
#' as are the derivatives in the response,
#' [distrib_grad_y.Gaussian1Distrib()] and [distrib_hess_y.Gaussian1Distrib()].
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale and reaches
#' the closed-form estimates \eqn{\hat\mu = \bar y} and
#' \eqn{\hat\sigma^2 = n^{-1}\sum (y_i - \bar y)^2}, the maximum likelihood
#' estimate of the variance with divisor \eqn{n}. The example below checks
#' both against the sample.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\sigma > 0} the standard deviation. \eqn{\eta} is a parameter on the
#' unconstrained scale of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Gaussian1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"gaussian1"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [gaussian2_distrib()] and [gaussian3_distrib()] for the same law in the
#' variance and in the precision; [lognormal1_distrib()] for the Gaussian on a
#' log scale; [student_t1_distrib()] and [laplace_distrib()] for heavier tails;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [Gaussian1Distrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 13.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dnorm pnorm qnorm rnorm
#'
#' @examples
#' d <- gaussian1_distrib()
#' d
#'
#' # The density and the distribution function are R's own.
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, sigma = 1.5)
#' all.equal(distrib_pdf(d, y, th), dnorm(y, 0.4, 1.5))
#'
#' # Moments in closed form: skewness and excess kurtosis are 0.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # Fitting recovers the closed-form maximum likelihood estimates, with the
#' # variance divided by n.
#' set.seed(7)
#' z <- distrib_rng(d, 400, list(mu = 3, sigma = 2))
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = mean(z), sigma = sqrt(mean((z - mean(z))^2))))
#'
#' # A different link changes the scale the optimizer moves on, not the answer.
#' d2 <- gaussian1_distrib(link_sigma = linkfunctions7::sqrt_link())
#' all.equal(coef(fit_distrib(d2, z)), coef(fit), tolerance = 1e-6)
#'
#' @export
gaussian1_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  
  Gaussian1Distrib(
    distrib_name = "gaussian1",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean", sigma = "standard deviation"),
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
