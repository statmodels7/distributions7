#' @include distrib.R generics.R
NULL

#' @title Lognormal Distribution Class
#' @name Lognormal1Distrib
#'
#' @description
#' The S7 class of the lognormal family on \eqn{(0, \infty)}: the law of a
#' variable whose logarithm is Gaussian with mean \eqn{\mu} and variance
#' \eqn{\sigma^2 > 0}. Both parameters live **on the log scale**, so neither is
#' the mean or the variance of \eqn{Y} itself. It inherits from
#' `continuous_distrib`, so it answers every generic of the `distrib`
#' contract; the eleven methods listed below are registered on it directly and
#' everything else comes from the parent.
#'
#' Build one with [lognormal1_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Lognormal1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [lognormal1_distrib()] they hold `"lognormal1"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "sigma2")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [lognormal1_distrib()] to build one;
#'   [lognormal2_distrib()] for the same law parametrized by the mean and the
#'   variance of \eqn{Y}; [gaussian2_distrib()], which this becomes at
#'   \eqn{\log y}; [distrib_pdf.Lognormal1Distrib()] and
#'   [distrib_gradient.Lognormal1Distrib()] for the closed forms this class
#'   supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Lognormal1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Lognormal1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Lognormal1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Lognormal1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Lognormal1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Lognormal1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Lognormal1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Lognormal1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Lognormal1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Lognormal1Distrib],
#'   [`distrib_rng()`][distrib_rng.Lognormal1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_bounds
#' d@bounds
#'
#' # mu is the mean of log(Y), so exp(mu) is the MEDIAN of Y and the mean sits
#' # above it.
#' th <- list(mu = 0.5, sigma2 = 0.36)
#' c(exp_mu = exp(0.5), median = distrib_quantile(d, 0.5, th),
#'   mean = mean(d, th))
Lognormal1Distrib <- S7::new_class("Lognormal1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Lognormal Probability Density Function
#' @name distrib_pdf.Lognormal1Distrib
#' @description
#' Computes the lognormal density
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}}
#'       \exp\left\{-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\},
#'       \qquad y > 0,}
#' by calling [stats::dlnorm()] at `meanlog = mu` and `sdlog = sqrt(sigma2)`.
#' The factor \eqn{1/y} is the Jacobian of the log transformation; it carries
#' no parameter, which is why every derivative in \eqn{\mu} and \eqn{\sigma^2}
#' is the Gaussian's read at \eqn{\log y}.
#'
#' With `log = TRUE` the logarithm is formed inside `dlnorm()` and stays finite
#' where the density itself underflows.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of observations. The support is
#'   \eqn{(0, \infty)}; a value at or below zero gives 0.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_cdf.Lognormal1Distrib()] for the distribution function,
#'   [distrib_gradient.Lognormal1Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf.Gaussian2Distrib()] for the law this is the
#'   exponential of, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#'
#' # The method is stats::dlnorm at meanlog = mu, sdlog = sqrt(sigma2).
#' all.equal(distrib_pdf(d, y, th), dlnorm(y, 0.5, sqrt(0.36)))
#'
#' # Times the Jacobian y, it is the Gaussian density at log y.
#' all.equal(distrib_pdf(d, y, th) * y,
#'           distrib_pdf(gaussian2_distrib(), log(y), th))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0, 0.5, 1), sigma2 = 0.36))
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e12, th)
#' distrib_pdf(d, 1e12, th, log = TRUE)
S7::method(distrib_pdf, Lognormal1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dlnorm(
    x = y,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    log = log
  )
}

#' @title Lognormal Cumulative Distribution Function
#' @name distrib_cdf.Lognormal1Distrib
#' @description
#' Computes the lognormal distribution function
#' \deqn{F(q; \mu, \sigma^2) = \Phi\left(\dfrac{\log q - \mu}{\sigma}\right),
#'       \qquad \sigma = \sqrt{\sigma^2},}
#' with \eqn{\Phi} the standard normal distribution function, by calling
#' [stats::plnorm()]. The log transformation is monotone, so this is the
#' Gaussian distribution function evaluated at \eqn{\log q} and nothing else.
#' Both tails are available exactly, and `log.p = TRUE` returns a logarithm
#' that stays finite where the probability itself underflows.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param q A numeric vector of quantiles. A value at or below zero gives a
#'   lower-tail probability of 0.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma2))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Lognormal1Distrib()] for the inverse,
#'   [distrib_pdf.Lognormal1Distrib()] for the density, [distrib_grad_cdf()]
#'   for the derivatives of this function in the parameters, which are the
#'   Gaussian's at \eqn{\log q} and so closed form, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' th <- list(mu = 0.5, sigma2 = 0.36)
#'
#' # The method is stats::plnorm at meanlog = mu, sdlog = sqrt(sigma2).
#' all.equal(distrib_cdf(d, c(0.5, 1.6, 4), th),
#'           plnorm(c(0.5, 1.6, 4), 0.5, sqrt(0.36)))
#'
#' # A monotone transformation, so it is the normal cdf at log q.
#' all.equal(distrib_cdf(d, c(0.5, 1.6, 4), th),
#'           pnorm(log(c(0.5, 1.6, 4)), 0.5, sqrt(0.36)))
#'
#' # Half the mass lies below exp(mu), which is the median.
#' distrib_cdf(d, exp(0.5), th)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 1e12, th, lower.tail = FALSE)
#' distrib_cdf(d, 1e12, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, Lognormal1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::plnorm(
    q = q,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Lognormal Quantile Function
#' @name distrib_quantile.Lognormal1Distrib
#' @description
#' Computes the lognormal quantile function
#' \deqn{Q(p; \mu, \sigma^2) = \exp\left\{\mu + \sigma\,\Phi^{-1}(p)\right\},
#'       \qquad \sigma = \sqrt{\sigma^2},}
#' by calling [stats::qlnorm()]. It is closed form, the log transformation
#' being monotone, so the round trip through
#' [distrib_cdf.Lognormal1Distrib()] returns `p` exactly. The median is
#' \eqn{e^{\mu}}, which sits **below** the mean \eqn{e^{\mu + \sigma^2/2}}.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities, which is how a quantile deep in a tail is
#'   requested without the probability underflowing. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty]}, of length
#'   `max(length(p), length(mu), length(sigma2))`.
#'
#' @seealso [distrib_cdf.Lognormal1Distrib()], which this inverts;
#'   [distrib_rng.Lognormal1Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' th <- list(mu = 0.5, sigma2 = 0.36)
#'
#' # A central 95 percent interval, multiplicatively symmetric about exp(mu).
#' q <- distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#' q
#' c(q[2] / q[1], q[3] / q[2])
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The median is exp(mu) and the mean is larger.
#' c(median = distrib_quantile(d, 0.5, th), exp_mu = exp(0.5),
#'   mean = mean(d, th))
S7::method(distrib_quantile, Lognormal1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qlnorm(
    p = p,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Lognormal Random Number Generator
#' @name distrib_rng.Lognormal1Distrib
#' @description
#' Draws `n` independent lognormal variates by calling [stats::rlnorm()] at
#' `meanlog = mu` and `sdlog = sqrt(sigma2)`, so the draws come from R's own
#' generator and depend on `.Random.seed` in the usual way. The
#' ratio-of-uniforms fallback the base class supplies is bypassed.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting.
#'   `sigma2` must be strictly positive.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_quantile.Lognormal1Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#'
#' # Same generator as stats::rlnorm at this parametrization.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0.5, sigma2 = 0.36))
#' set.seed(2)
#' identical(a, rlnorm(3, 0.5, sqrt(0.36)))
#'
#' # The parameters are moments of the LOGARITHM, so that is where the sample
#' # recovers them.
#' set.seed(6)
#' z <- distrib_rng(d, 2e4, list(mu = 0.5, sigma2 = 0.36))
#' c(mu = mean(log(z)), sigma2 = var(log(z)))
#'
#' # On the original scale the sample mean is exp(mu + sigma2/2).
#' c(sample = mean(z), theory = exp(0.5 + 0.36 / 2))
S7::method(distrib_rng, Lognormal1Distrib) <- function(distrib, n, theta) {
  stats::rlnorm(
    n = n,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]])
  )
}

#' @title Lognormal Score
#' @name distrib_gradient.Lognormal1Distrib
#' @description
#' Computes the first derivatives of the lognormal log-density with respect to
#' \eqn{\mu} and \eqn{\sigma^2}, one value per observation, in closed form.
#' With \eqn{r = \log y - \mu},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma^2} =
#'         \dfrac{r^2 - \sigma^2}{2\sigma^4}.}
#'
#' **These are exactly the Gaussian's, read at \eqn{\log y}.** The lognormal
#' log-density is the Gaussian's at \eqn{\log y} minus \eqn{\log y}, and that
#' last term carries no parameter, so it disappears from every derivative in
#' \eqn{\mu} and \eqn{\sigma^2}. The same holds at all four orders and for the
#' expected values; what differs from [gaussian2_distrib()] is only the
#' derivatives in the response, where the Jacobian does enter.
#'
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of strictly positive observations.
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
#' \eqn{\sigma^2 > 0} the variance **of \eqn{\log Y}**, not of \eqn{Y}.
#' \eqn{r = \log y - \mu} is the residual on the log scale.
#'
#' @seealso [distrib_hessian.Lognormal1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Lognormal1Distrib()] for their expectation,
#'   [distrib_gradient.Gaussian2Distrib()], which returns the same numbers at
#'   \eqn{\log y}, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out on the log scale.
#' r <- log(y) - 0.5
#' all.equal(g$mu, r / 0.36)
#' all.equal(g$sigma2, (r^2 - 0.36) / (2 * 0.36^2))
#'
#' # Identical to the Gaussian's score at log y, component for component.
#' all.equal(g, distrib_gradient(gaussian2_distrib(), log(y), th))
#'
#' # The summed score vanishes at the closed-form estimates, which are the
#' # sample moments of the logarithm.
#' set.seed(6)
#' z <- distrib_rng(d, 2000, th)
#' mle <- list(mu = mean(log(z)),
#'             sigma2 = mean((log(z) - mean(log(z)))^2))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
S7::method(distrib_gradient, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  lognormal_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Lognormal Observed Hessian
#' @name distrib_hessian.Lognormal1Distrib
#' @description
#' Computes the three distinct second derivatives of the lognormal log-density
#' with respect to \eqn{\mu} and \eqn{\sigma^2}, one value per observation, in
#' closed form. With \eqn{r = \log y - \mu},
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{1}{\sigma^2}, \qquad
#'       \ell^{(\mu\sigma^2)} = -\dfrac{r}{\sigma^4}, \qquad
#'       \ell^{(\sigma^2\sigma^2)} = \dfrac{1}{2\sigma^4}
#'         - \dfrac{r^2}{\sigma^6}.}
#' These are the Gaussian's at \eqn{\log y}, the Jacobian of the log
#' transformation carrying no parameter. Only the curvature in \eqn{\mu} is
#' free of the data; the other two vary with the residual on the log scale, and
#' their expectations are
#' [distrib_expected_hessian.Lognormal1Distrib()].
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of strictly positive observations.
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
#' @return A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
#'   `mu_sigma2`, in that order, each of length
#'   `max(length(y), length(mu), length(sigma2))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#' \eqn{\sigma^2} is the variance of \eqn{\log Y}.
#'
#' @seealso [distrib_gradient.Lognormal1Distrib()] for the score,
#'   [distrib_expected_hessian.Lognormal1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Lognormal1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#' h <- distrib_hessian(d, y, th)
#'
#' # The curvature in mu is constant at -1/sigma2; the other two are not.
#' h$mu_mu
#' r <- log(y) - 0.5
#' all.equal(h$mu_sigma2, -r / 0.36^2)
#' all.equal(h$sigma2_sigma2, 1 / (2 * 0.36^2) - r^2 / 0.36^3)
#'
#' # Identical to the Gaussian's at log y.
#' all.equal(h, distrib_hessian(gaussian2_distrib(), log(y), th))
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 0.5, sigma2 = 0.36 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.5, sigma2 = 0.36 - eps))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_sigma2, tolerance = 1e-5)
S7::method(distrib_hessian, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  lognormal_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Lognormal Expected Hessian
#' @name distrib_expected_hessian.Lognormal1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\sigma^2}, \qquad
#'       \mathbb{E}\left[\ell^{(\sigma^2\sigma^2)}\right] =
#'         -\dfrac{1}{2\sigma^4}, \qquad
#'       \mathbb{E}\left[\ell^{(\mu\sigma^2)}\right] = 0.}
#' They follow from \eqn{\mathbb{E}[\log Y] = \mu} and
#' \eqn{\mathbb{E}[(\log Y - \mu)^2] = \sigma^2}, \eqn{\log Y} being Gaussian
#' by construction. The zero off-diagonal says the two parameters are
#' orthogonal, so their estimates are asymptotically independent.
#'
#' The information does not depend on the mean of \eqn{Y} at all, only on the
#' variance of its logarithm. Because the values do not depend on the data,
#' `approx` and `nsim` are ignored and `y` is read only for its length.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
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
#' @return A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
#'   `mu_sigma2`, in that order, each of length
#'   `max(length(y), length(mu), length(sigma2))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The lognormal
#' is a regular family, so the second Bartlett identity holds and this equals
#' the variance of the score.
#'
#' @seealso [distrib_hessian.Lognormal1Distrib()] for the observed quantity
#'   this is the expectation of,
#'   [distrib_expected_hessian.Gaussian2Distrib()], which returns the same
#'   numbers, [fisher_scoring()], which inverts it at each step, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' th <- list(mu = 0.5, sigma2 = 0.36)
#'
#' # The three constants, one value per observation.
#' lapply(distrib_expected_hessian(d, c(0.5, 1.6, 4), th), unique)
#' c(-1 / 0.36, -1 / (2 * 0.36^2), 0)
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
#'
#' # The information does not move with the mean of Y, only with sigma2.
#' vapply(c(-2, 0, 5),
#'        function(m) distrib_expected_hessian(d, 0,
#'                      list(mu = m, sigma2 = 0.36))$sigma2_sigma2,
#'        numeric(1))
S7::method(distrib_expected_hessian, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  lognormal_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Lognormal Third-Order Derivatives
#' @name distrib_deriv3.Lognormal1Distrib
#' @description
#' Computes the four distinct third derivatives of the lognormal log-density
#' with respect to \eqn{\mu} and \eqn{\sigma^2}, in closed form. With
#' \eqn{r = \log y - \mu},
#' \deqn{\ell^{(\mu\mu\mu)} = 0, \qquad
#'       \ell^{(\mu\mu\sigma^2)} = \dfrac{1}{\sigma^4}, \qquad
#'       \ell^{(\mu\sigma^2\sigma^2)} = \dfrac{2r}{\sigma^6}, \qquad
#'       \ell^{(\sigma^2\sigma^2\sigma^2)} = \dfrac{3r^2}{\sigma^8}
#'         - \dfrac{1}{\sigma^6}.}
#' The log-density is quadratic in \eqn{\mu}, so the third derivative there is
#' zero. As at the lower orders these are the Gaussian's read at \eqn{\log y}.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{r} with 0 and \eqn{r^2} with \eqn{\sigma^2}: the component odd in
#' \eqn{r} vanishes and the pure variance one becomes \eqn{2/\sigma^6}. Both
#' routes are closed form, so no quadrature is run and `approx` and `nsim` are
#' ignored.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
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
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_hessian.Lognormal1Distrib()] for the order below and
#'   [distrib_deriv4.Lognormal1Distrib()] for the order above;
#'   [distrib_deriv3.Gaussian2Distrib()], which returns the same numbers at
#'   \eqn{\log y}; [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # Quadratic in mu, so the third derivative there is 0.
#' d3$mu_mu_mu
#'
#' # The other three, written out on the log scale.
#' r <- log(y) - 0.5
#' all.equal(d3$mu_mu_sigma2, rep(1 / 0.36^2, 3))
#' all.equal(d3$mu_sigma2_sigma2, 2 * r / 0.36^3)
#'
#' # Expected values: the component odd in r vanishes.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#' 2 / 0.36^3
#'
#' # Identical to the Gaussian's at log y.
#' all.equal(d3, distrib_deriv3(gaussian2_distrib(), log(y), th))
S7::method(distrib_deriv3, Lognormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) lognormal_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else lognormal_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Lognormal Fourth-Order Derivatives
#' @name distrib_deriv4.Lognormal1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the lognormal log-density
#' with respect to \eqn{\mu} and \eqn{\sigma^2}, in closed form. With
#' \eqn{r = \log y - \mu},
#' \deqn{\ell^{(\mu\mu\mu\mu)} = \ell^{(\mu\mu\mu\sigma^2)} = 0, \qquad
#'       \ell^{(\mu\mu\sigma^2\sigma^2)} = -\dfrac{2}{\sigma^6}, \qquad
#'       \ell^{(\mu\sigma^2\sigma^2\sigma^2)} = -\dfrac{6r}{\sigma^8}, \qquad
#'       \ell^{(\sigma^{2\cdot 4})} = \dfrac{3}{\sigma^8}
#'         - \dfrac{12 r^2}{\sigma^{10}}.}
#' As at the lower orders these are the Gaussian's read at \eqn{\log y}.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{r} with 0 and \eqn{r^2} with \eqn{\sigma^2}, which leaves
#' \eqn{-9/\sigma^8} in the last component and zero in the two odd ones. Both
#' routes are closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
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
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.Lognormal1Distrib()] for the order below,
#'   [distrib_deriv4.Gaussian2Distrib()], which returns the same numbers at
#'   \eqn{\log y}, and [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Quadratic in mu, so the two components with three or more mu are zero.
#' c(d4$mu_mu_mu_mu[1], d4$mu_mu_mu_sigma2[1])
#'
#' # The mixed second-second component is constant at -2/sigma2^3.
#' all.equal(d4$mu_mu_sigma2_sigma2, rep(-2 / 0.36^3, 3))
#'
#' # Expected values: -9/sigma2^4 in the pure variance component.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#' -9 / 0.36^4
S7::method(distrib_deriv4, Lognormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) lognormal_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else lognormal_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Lognormal First Derivative in the Response
#' @name distrib_grad_y.Lognormal1Distrib
#' @description
#' Computes the first derivative of the lognormal log-density with respect to
#' the response, in closed form. With \eqn{r = \log y - \mu},
#' \deqn{\dfrac{\partial \ell}{\partial y} =
#'       -\dfrac{1}{y}\left(1 + \dfrac{r}{\sigma^2}\right).}
#' This is where the family parts company with the Gaussian: the Jacobian
#' \eqn{1/y} of the log transformation carries no parameter and so leaves every
#' derivative in \eqn{\mu} and \eqn{\sigma^2} alone, but it is a function of
#' the response and does enter here.
#'
#' The derivative vanishes at \eqn{r = -\sigma^2}, that is at
#' \eqn{y = e^{\mu - \sigma^2}}, which is the mode of the density and lies
#' below both the median \eqn{e^{\mu}} and the mean
#' \eqn{e^{\mu + \sigma^2/2}}.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of strictly positive observations. The value
#'   diverges as `y` approaches zero.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Lognormal1Distrib()] for the second derivative in
#'   the response, [distrib_gradient.Lognormal1Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#'
#' # Written out.
#' r <- log(y) - 0.5
#' all.equal(distrib_grad_y(d, y, th), -(1 + r / 0.36) / y)
#'
#' # Zero at the mode exp(mu - sigma2), which lies below the median exp(mu).
#' mode <- exp(0.5 - 0.36)
#' c(mode = mode, median = exp(0.5), at_mode = distrib_grad_y(d, mode, th))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-7
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, Lognormal1Distrib) <- function(distrib, y, theta) {
  r <- log(y) - theta[[1]]
  s2 <- theta[[2]]
  -(1 + r / s2) / y
}

#' @title Lognormal Second Derivative in the Response
#' @name distrib_hess_y.Lognormal1Distrib
#' @description
#' Computes the second derivative of the lognormal log-density with respect to
#' the response, in closed form. With \eqn{r = \log y - \mu},
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} =
#'       \dfrac{1}{y^2}\left(1 + \dfrac{r - 1}{\sigma^2}\right).}
#' The sign changes at \eqn{r = 1 - \sigma^2}, that is at
#' \eqn{y = e^{\mu + 1 - \sigma^2}}: the log-density is **concave** in the
#' response below that point and **convex** above it, so a lognormal
#' log-likelihood is not a concave function of an observation over the whole
#' support.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param y A numeric vector of strictly positive observations. The value
#'   diverges as `y` approaches zero.
#' @param theta A named list with components `mu` and `sigma2`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `sigma2` must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma2))`, one value per observation.
#'
#' @seealso [distrib_grad_y.Lognormal1Distrib()] for the first derivative in
#'   the response, [distrib_hessian.Lognormal1Distrib()] for the curvature in
#'   the parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- lognormal1_distrib()
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#'
#' r <- log(y) - 0.5
#' all.equal(distrib_hess_y(d, y, th), (1 + (r - 1) / 0.36) / y^2)
#'
#' # Concave below exp(mu + 1 - sigma2) and convex above it.
#' cut <- exp(0.5 + 1 - 0.36)
#' c(cut = cut,
#'   below = distrib_hess_y(d, cut / 2, th),
#'   at = distrib_hess_y(d, cut, th),
#'   above = distrib_hess_y(d, 2 * cut, th))
S7::method(distrib_hess_y, Lognormal1Distrib) <- function(distrib, y, theta) {
  r <- log(y) - theta[[1]]
  s2 <- theta[[2]]
  (1 + (r - 1) / s2) / y^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Lognormal Distribution, Log-Scale Parametrization
#'
#' @description
#' Builds the distribution object for the lognormal family on
#' \eqn{(0, \infty)}: the law of a variable whose logarithm is Gaussian with
#' mean \eqn{\mu} and variance \eqn{\sigma^2 > 0}. Both parameters are moments
#' **of the logarithm**, so neither is the mean or the variance of \eqn{Y}. The
#' returned object carries closed-form derivatives of the log-density to fourth
#' order, in the parameters and in the response, and closed-form moments, so
#' every generic of the toolkit answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. The defaults are the identity for
#' \eqn{\mu}, which ranges over the whole line, and the logarithm for
#' \eqn{\sigma^2}, which keeps it positive at every predictor.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the location
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], \eqn{\mu} being
#'   a mean on the log scale and so already free.
#' @param link_sigma2 A `link` object from `linkfunctions7` for the variance
#'   \eqn{\sigma^2} of the logarithm. Defaults to
#'   [linkfunctions7::log_link()], which maps \eqn{(0, \infty)} onto the line.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, \infty)} is
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}}
#'       \exp\left\{-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\},}
#' the distribution function \eqn{F(q) = \Phi\{(\log q - \mu)/\sigma\}} and the
#' quantile function \eqn{Q(p) = \exp\{\mu + \sigma\Phi^{-1}(p)\}}, both closed
#' form because the log transformation is monotone.
#'
#' On the original scale the mean is \eqn{e^{\mu + \sigma^2/2}}, the variance
#' \eqn{(e^{\sigma^2}-1)e^{2\mu+\sigma^2}}, the skewness
#' \eqn{(e^{\sigma^2}+2)\sqrt{e^{\sigma^2}-1}} and the excess kurtosis
#' \eqn{e^{4\sigma^2}+2e^{3\sigma^2}+3e^{2\sigma^2}-6}. Three quantities are
#' worth keeping apart: the mode is \eqn{e^{\mu-\sigma^2}}, the median
#' \eqn{e^{\mu}} and the mean \eqn{e^{\mu+\sigma^2/2}}, in that order.
#'
#' [lognormal2_distrib()] carries the same law parametrized by the mean and the
#' variance of \eqn{Y} itself.
#'
#' # Every derivative in the parameters is the Gaussian's
#'
#' The log-density is the Gaussian's evaluated at \eqn{\log y}, minus
#' \eqn{\log y}. That last term is the Jacobian of the transformation and
#' carries no parameter, so it disappears from every derivative in \eqn{\mu}
#' and \eqn{\sigma^2}. Writing \eqn{r = \log y - \mu}, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma^2} =
#'         \dfrac{r^2 - \sigma^2}{2\sigma^4},}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\ell^{(\sigma^2\sigma^2)}\right] =
#'         -\dfrac{1}{2\sigma^4}, \quad
#'       \mathbb{E}\left[\ell^{(\mu\sigma^2)}\right] = 0.}
#' Every one of these agrees with [gaussian2_distrib()] at \eqn{\log y},
#' component for component and at all four orders, observed and expected. The
#' two parameters are orthogonal and the information does not move with the
#' mean of \eqn{Y}.
#'
#' The derivatives **in the response** are where the two families differ, the
#' Jacobian being a function of \eqn{y}. They are still closed form; see
#' [distrib_grad_y.Lognormal1Distrib()] and
#' [distrib_hess_y.Lognormal1Distrib()], and note that the log-density is
#' convex in the response above \eqn{e^{\mu+1-\sigma^2}}.
#'
#' # Estimation
#'
#' Both maximum likelihood estimates are closed form and are the sample moments
#' of the logarithm:
#' \deqn{\hat\mu = \dfrac{1}{n}\sum_i \log y_i, \qquad
#'       \hat\sigma^2 = \dfrac{1}{n}\sum_i (\log y_i - \hat\mu)^2,}
#' the second with divisor \eqn{n}. [fit_distrib()] reaches them on the link
#' scale, and the example below checks both against the sample.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean and
#' \eqn{\sigma^2 > 0} the variance **of \eqn{\log Y}**. \eqn{\Phi} is the
#' standard normal distribution function and \eqn{r = \log y - \mu} the
#' residual on the log scale. \eqn{\eta} is a parameter on the unconstrained
#' scale of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Lognormal1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"lognormal1"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "sigma2")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @seealso
#' [lognormal2_distrib()] for the same law in the mean and the variance of
#' \eqn{Y}; [gaussian2_distrib()] for the law of its logarithm;
#' [gamma1_distrib()] and [invgauss1_distrib()] for other positive families
#' with a multiplicative variance function; [transformation()] for the general
#' wrapper this family is a named instance of; [fit_distrib()] to estimate the
#' parameters; [check_distrib()] to validate a family of your own against the
#' same battery this one passes; [Lognormal1Distrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' *Continuous Univariate Distributions*, Volume 1, 2nd edition, Chapter 14.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dlnorm plnorm qlnorm rlnorm
#'
#' @examples
#' d <- lognormal1_distrib()
#' d
#'
#' # The density is stats::dlnorm at meanlog = mu, sdlog = sqrt(sigma2).
#' y <- c(0.5, 1.6, 4)
#' th <- list(mu = 0.5, sigma2 = 0.36)
#' all.equal(distrib_pdf(d, y, th), dlnorm(y, 0.5, sqrt(0.36)))
#'
#' # Mode, median and mean, in that order.
#' c(mode = exp(0.5 - 0.36), median = exp(0.5), mean = mean(d, th))
#'
#' # Moments on the original scale, all exponential in sigma2.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # Every derivative in the parameters is the Gaussian's at log y.
#' all.equal(distrib_gradient(d, y, th),
#'           distrib_gradient(gaussian2_distrib(), log(y), th))
#'
#' # Fitting recovers the closed-form estimates, the sample moments of the
#' # logarithm with divisor n.
#' set.seed(6)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = mean(log(z)),
#'                  sigma2 = mean((log(z) - mean(log(z)))^2)))
#'
#' @export
lognormal1_distrib <- function(link_mu = identity_link(), link_sigma2 = log_link()) {
  
  Lognormal1Distrib(
    distrib_name = "lognormal1",
    dimension = "univariate",
    bounds = c(0, Inf),
    
    params = c("mu", "sigma2"),
    params_interpretation = c(mu = "mean (log scale)", sigma2 = "variance (log scale)"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma2 = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      sigma2 = link_sigma2
    )
  )
  
}
