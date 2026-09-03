#' @include distrib.R generics.R
NULL

#' @title Inverse Gaussian Distribution Class, Mean and Dispersion
#' @name InvGauss1Distrib
#'
#' @description
#' The S7 class of the inverse Gaussian family on \eqn{(0, \infty)}
#' parametrized by its mean \eqn{\mu > 0} and a dispersion \eqn{\phi > 0}, so
#' that \eqn{\operatorname{Var}(Y) = \phi\mu^3}. It inherits from
#' `continuous_distrib`, so it answers every generic of the `distrib` contract;
#' the eleven methods listed below are registered on it directly and everything
#' else comes from the parent.
#'
#' Build one with [invgauss1_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `InvGauss1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [invgauss1_distrib()] they hold `"invgauss1"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "phi")`, the interpretations
#'   `c(mu = "mean", phi = "dispersion")`, `2`, the domain \eqn{(0, \infty)}
#'   for both parameters, and the two links.
#'
#' @seealso [invgauss1_distrib()] to build one;
#'   [invgauss2_distrib()] for the same law in mean and variance;
#'   [distrib_pdf.InvGauss1Distrib()] and
#'   [distrib_gradient.InvGauss1Distrib()] for the closed forms this class
#'   supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.InvGauss1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.InvGauss1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.InvGauss1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.InvGauss1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.InvGauss1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.InvGauss1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.InvGauss1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.InvGauss1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.InvGauss1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.InvGauss1Distrib],
#'   [`distrib_rng()`][distrib_rng.InvGauss1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- invgauss1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # The variance function is mu^3, so the spread grows faster with the mean
#' # than a gamma's, whose variance function is mu^2.
#' vapply(c(1, 2, 4), function(m) variance(d, list(mu = m, phi = 2)),
#'        numeric(1))
InvGauss1Distrib <- S7::new_class("InvGauss1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Inverse Gaussian Probability Density Function in Mean and Dispersion
#' @name distrib_pdf.InvGauss1Distrib
#' @description
#' Computes the inverse Gaussian density
#' \deqn{f(y; \mu, \phi) = \sqrt{\dfrac{1}{2\pi\phi y^3}}
#'       \exp\left\{-\dfrac{(y-\mu)^2}{2\phi\mu^2 y}\right\}, \qquad y > 0,}
#' by calling [statmod::dinvgauss()] at `mean = mu` and `dispersion = phi`.
#' With `log = TRUE` the logarithm is formed inside that function and stays
#' finite where the density itself underflows, which happens quickly: the
#' density falls off like \eqn{\exp\{-y/(2\phi\mu^2)\}} in the right tail and
#' like \eqn{\exp\{-1/(2\phi y)\}} at the origin.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of observations. The support is
#'   \eqn{(0, \infty)}; a value at or below zero gives 0.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(phi))`, one value per observation.
#'
#' @seealso [distrib_cdf.InvGauss1Distrib()] for the distribution function,
#'   [distrib_gradient.InvGauss1Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf.InvGauss2Distrib()] for the same density in the
#'   variance, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#'
#' # The method is statmod::dinvgauss at this parametrization.
#' all.equal(distrib_pdf(d, y, th),
#'           statmod::dinvgauss(y, mean = 1, dispersion = 2))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0.5, 1, 2), phi = 2))
#'
#' # The density vanishes at the origin faster than any power.
#' distrib_pdf(d, c(1e-2, 1e-3, 1e-4), th)
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e4, th)
#' distrib_pdf(d, 1e4, th, log = TRUE)
S7::method(distrib_pdf, InvGauss1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  statmod::dinvgauss(
    x = y,
    mean = theta[[1]],
    dispersion = theta[[2]],
    log = log
  )
}

#' @title Inverse Gaussian Cumulative Distribution Function in Mean and Dispersion
#' @name distrib_cdf.InvGauss1Distrib
#' @description
#' Computes the inverse Gaussian distribution function, which is elementary in
#' the standard normal distribution function \eqn{\Phi}:
#' \deqn{F(q; \mu, \phi) = \Phi\left\{\sqrt{\dfrac{1}{\phi q}}
#'         \left(\dfrac{q}{\mu} - 1\right)\right\}
#'       + e^{2/(\phi\mu)}\,
#'       \Phi\left\{-\sqrt{\dfrac{1}{\phi q}}
#'         \left(\dfrac{q}{\mu} + 1\right)\right\},}
#' by calling [statmod::pinvgauss()] at `mean = mu` and `dispersion = phi`.
#' Both tails are available exactly, and `log.p = TRUE` returns a logarithm
#' that stays finite where the probability itself underflows.
#'
#' The exponential factor is the reason this expression is evaluated on the
#' log scale rather than as written: \eqn{e^{2/(\phi\mu)}} overflows at
#' ordinary settings, at \eqn{\mu = 0.01} and \eqn{\phi = 0.1} the exponent
#' already being 2000, and it multiplies a normal tail that underflows by the
#' same amount.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param q A numeric vector of quantiles. A value at or below zero gives a
#'   lower-tail probability of 0.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(phi))`. With `log.p = TRUE` the values
#'   are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.InvGauss1Distrib()] for the inverse,
#'   [distrib_pdf.InvGauss1Distrib()] for the density, [distrib_grad_cdf()] for
#'   the derivatives of this function in the parameters, which are closed form
#'   here because the expression above is elementary, and [distrib_cdf()] for
#'   the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' th <- list(mu = 1, phi = 2)
#' q <- c(0.5, 1, 2)
#'
#' # The method is statmod::pinvgauss at this parametrization.
#' all.equal(distrib_cdf(d, q, th),
#'           statmod::pinvgauss(q, mean = 1, dispersion = 2))
#'
#' # The closed form above, evaluated directly at these safe values.
#' a <- sqrt(1 / (2 * q)) * (q / 1 - 1)
#' b <- -sqrt(1 / (2 * q)) * (q / 1 + 1)
#' pnorm(a) + exp(2 / (2 * 1)) * pnorm(b)
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 2, th) + distrib_cdf(d, 2, th, lower.tail = FALSE)
#'
#' # The law is heavily right skewed, so most of the mass sits below the mean.
#' distrib_cdf(d, 1, th)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 1e4, th, lower.tail = FALSE)
#' distrib_cdf(d, 1e4, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, InvGauss1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  statmod::pinvgauss(
    q = q,
    mean = theta[[1]],
    dispersion = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Inverse Gaussian Quantile Function in Mean and Dispersion
#' @name distrib_quantile.InvGauss1Distrib
#' @description
#' Computes the inverse Gaussian quantile function by calling
#' [statmod::qinvgauss()] at `mean = mu` and `dispersion = phi`. There is no
#' closed form: the distribution function is elementary but not invertible in
#' elementary terms, so `qinvgauss()` inverts it numerically. The distribution
#' function is strictly increasing on \eqn{(0, \infty)}, so the round trip
#' through [distrib_cdf.InvGauss1Distrib()] returns `p`.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities, which is how a quantile deep in a tail is
#'   requested without the probability underflowing. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty]}, of length
#'   `max(length(p), length(mu), length(phi))`.
#'
#' @seealso [distrib_cdf.InvGauss1Distrib()], which this inverts;
#'   [distrib_rng.InvGauss1Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' th <- list(mu = 1, phi = 2)
#'
#' # A central 95 percent interval, extremely asymmetric about the mean of 1.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The median falls well below the mean at this dispersion.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, InvGauss1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  statmod::qinvgauss(
    p = p,
    mean = theta[[1]],
    dispersion = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Inverse Gaussian Random Number Generator in Mean and Dispersion
#' @name distrib_rng.InvGauss1Distrib
#' @description
#' Draws `n` independent inverse Gaussian variates by calling
#' [statmod::rinvgauss()] at `mean = mu` and `dispersion = phi`, so the draws
#' come from that package's generator and depend on `.Random.seed` in the usual
#' way. The ratio-of-uniforms fallback the base class supplies is bypassed.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. Both
#'   must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_quantile.InvGauss1Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#'
#' # Same generator as statmod::rinvgauss at this parametrization.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 1, phi = 2))
#' set.seed(2)
#' identical(a, statmod::rinvgauss(3, mean = 1, dispersion = 2))
#'
#' # Both maximum likelihood estimates are closed form here: the mean is the
#' # sample mean and the dispersion the mean of 1/y minus 1/ybar.
#' set.seed(3)
#' z <- distrib_rng(d, 2e4, list(mu = 1, phi = 2))
#' c(mu = mean(z), phi = mean(1 / z) - 1 / mean(z))
S7::method(distrib_rng, InvGauss1Distrib) <- function(distrib, n, theta, ...) {
  statmod::rinvgauss(
    n = n,
    mean = theta[[1]],
    dispersion = theta[[2]]
  )
}

#' @title Inverse Gaussian Score in Mean and Dispersion
#' @name distrib_gradient.InvGauss1Distrib
#' @description
#' Computes the first derivatives of the inverse Gaussian log-density with
#' respect to \eqn{\mu} and \eqn{\phi}, one value per observation, in closed
#' form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^3},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \phi} =
#'         \dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2}.}
#' The mean component is the score of an inverse Gaussian generalized linear
#' model, the residual divided by the variance function \eqn{\phi\mu^3}. The
#' arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `phi`, each of length
#'   `max(length(y), length(mu), length(phi))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\phi > 0} the dispersion, with
#' \eqn{\operatorname{Var}(Y) = \phi\mu^3}.
#'
#' @seealso [distrib_hessian.InvGauss1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.InvGauss1Distrib()] for their expectation,
#'   [distrib_grad_y.InvGauss1Distrib()] for the derivative in the response,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#' g <- distrib_gradient(d, y, th)
#'
#' # The generalized linear model score: residual over the variance function.
#' all.equal(g$mu, (y - 1) / (2 * 1^3))
#' all.equal(g$phi, ((y - 1)^2 - y * 1^2 * 2) / (2 * y * 2^2 * 1^2))
#'
#' # The mean component vanishes at y = mu, whatever the dispersion.
#' distrib_gradient(d, 1, list(mu = 1, phi = c(0.5, 2, 8)))$mu
#'
#' # Both estimating equations solve in closed form, so the summed score
#' # vanishes at the sample mean and mean(1/y) - 1/mean(y).
#' set.seed(3)
#' z <- distrib_rng(d, 2000, th)
#' mle <- list(mu = mean(z), phi = mean(1 / z) - 1 / mean(z))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
S7::method(distrib_gradient, InvGauss1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  invgauss_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Observed Hessian in Mean and Dispersion
#' @name distrib_hessian.InvGauss1Distrib
#' @description
#' Computes the three distinct second derivatives of the inverse Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\phi}, one value per
#' observation, in closed form:
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{3y - 2\mu}{\phi\mu^4}, \qquad
#'       \ell^{(\phi\phi)} =
#'         \dfrac{\phi - 2(y-\mu)^2/(\mu^2 y)}{2\phi^3}, \qquad
#'       \ell^{(\mu\phi)} = -\dfrac{y - \mu}{\phi^2\mu^3}.}
#'
#' **Neither diagonal entry is negative at every observation.** The curvature
#' in \eqn{\mu} turns positive wherever \eqn{y < 2\mu/3}, and the curvature in
#' \eqn{\phi} is negative only where \eqn{\phi < 2(y-\mu)^2/(\mu^2 y)}, so at
#' \eqn{y = \mu} it is \eqn{+1/(2\phi^2)}. An observation near the mean
#' therefore contributes positive curvature in the dispersion, and a Newton
#' step taken on the observed Hessian can move the wrong way. The expected
#' Hessian is negative definite everywhere, which is why
#' [fisher_scoring()] is the more stable route here; see
#' [distrib_expected_hessian.InvGauss1Distrib()].
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `phi_phi` and
#'   `mu_phi`, in that order, each of length
#'   `max(length(y), length(mu), length(phi))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_gradient.InvGauss1Distrib()] for the score,
#'   [distrib_expected_hessian.InvGauss1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.InvGauss1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # The three closed forms, written out.
#' all.equal(h$mu_mu, -(3 * y - 2 * 1) / (2 * 1^4))
#' all.equal(h$phi_phi, (2 - 2 * (y - 1)^2 / (1^2 * y)) / (2 * 2^3))
#' all.equal(h$mu_phi, -(y - 1) / (2^2 * 1^3))
#'
#' # Positive curvature in mu below 2 mu/3, and in phi at y = mu.
#' c(mu_mu_at_half = distrib_hessian(d, 0.5, th)$mu_mu,
#'   phi_phi_at_mu = distrib_hessian(d, 1, th)$phi_phi,
#'   one_over_2phi2 = 1 / (2 * 2^2))
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 1, phi = 2 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 1, phi = 2 - eps))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_phi, tolerance = 1e-5)
S7::method(distrib_hessian, InvGauss1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  invgauss_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Expected Hessian in Mean and Dispersion
#' @name distrib_expected_hessian.InvGauss1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\phi\mu^3},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\phi\phi)}\right] = -\dfrac{1}{2\phi^2},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\phi)}\right] = 0.}
#' They follow from \eqn{\mathbb{E}[Y] = \mu} and
#' \eqn{\mathbb{E}[(Y-\mu)^2/Y] = \phi\mu^2}. The second identity is the one
#' that gives the score in \eqn{\phi} mean zero.
#'
#' Both diagonal entries are negative at every parameter setting, so the
#' information is positive definite everywhere, where the observed Hessian of
#' [distrib_hessian.InvGauss1Distrib()] is not. The zero off-diagonal says the
#' mean and the dispersion are orthogonal, so the mean equation can be fitted
#' with the dispersion held at any value without biasing it: this is the
#' generalized linear model parametrization, with variance function
#' \eqn{V(\mu) = \mu^3}.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
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
#' @return A named list of three numeric vectors, `mu_mu`, `phi_phi` and
#'   `mu_phi`, in that order, each of length
#'   `max(length(y), length(mu), length(phi))` and constant within itself when
#'   the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The inverse
#' Gaussian is a regular family, so the second Bartlett identity holds and this
#' equals the variance of the score.
#'
#' @seealso [distrib_hessian.InvGauss1Distrib()] for the observed quantity this
#'   is the expectation of, [fisher_scoring()], which inverts it at each step,
#'   and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' th <- list(mu = 1, phi = 2)
#'
#' # The three constants, one value per observation.
#' lapply(distrib_expected_hessian(d, c(0.5, 1, 2), th), unique)
#' c(-1 / (2 * 1^3), -1 / (2 * 2^2), 0)
#'
#' # Negative definite at every setting, where the observed Hessian is
#' # positive in phi at y = mu.
#' c(expected = distrib_expected_hessian(d, 1, th)$phi_phi,
#'   observed = distrib_hessian(d, 1, th)$phi_phi)
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_expected_hessian, InvGauss1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
                                       threads = 1L) {
  invgauss_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Third-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv3.InvGauss1Distrib
#' @description
#' Computes the four distinct third derivatives of the inverse Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\phi}, in closed form. The
#' log-density is a linear combination of \eqn{y}, \eqn{1/y} and \eqn{\log y}
#' with coefficients rational in the parameters, so every derivative is a
#' rational function of \eqn{(\mu, \phi)} times one of those three statistics.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{Y} with \eqn{\mu} and \eqn{1/Y} with \eqn{1/\mu + \phi}. Both routes
#' are closed form, so no quadrature is run and `approx` and `nsim` are
#' ignored.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
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
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_phi`,
#'   `mu_phi_phi` and `phi_phi_phi`, each of length
#'   `max(length(y), length(mu), length(phi))`. The names enumerate the
#'   distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name derivatives.
#'
#' @seealso [distrib_hessian.InvGauss1Distrib()] for the order below and
#'   [distrib_deriv4.InvGauss1Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic and for the numerical route a family
#'   without a closed form takes.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # The expected values are constants at a fixed parameter setting, and the
#' # mixed mu-phi-phi component is exactly zero, the score in mu being linear
#' # in y and the two parameters entering it as a product.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 1 + eps, phi = 2))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 1 - eps, phi = 2))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
S7::method(distrib_deriv3, InvGauss1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) invgauss_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else invgauss_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Fourth-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv4.InvGauss1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the inverse Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\phi}, in closed form, by the
#' same route the lower orders take: the log-density is a linear combination of
#' \eqn{y}, \eqn{1/y} and \eqn{\log y} with coefficients rational in the
#' parameters.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{Y} with \eqn{\mu} and \eqn{1/Y} with \eqn{1/\mu + \phi}. Both routes
#' are closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
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
#'   `mu_mu_mu_phi`, `mu_mu_phi_phi`, `mu_phi_phi_phi` and `phi_phi_phi_phi`,
#'   each of length `max(length(y), length(mu), length(phi))`.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.InvGauss1Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.InvGauss1Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 1 + eps, phi = 2))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 1 - eps, phi = 2))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
S7::method(distrib_deriv4, InvGauss1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) invgauss_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else invgauss_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian First Derivative in the Response
#' @name distrib_grad_y.InvGauss1Distrib
#' @description
#' Computes the first derivative of the inverse Gaussian log-density with
#' respect to the response, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{3}{2y}
#'       - \dfrac{y^2 - \mu^2}{2\phi\mu^2 y^2}.}
#' The first term comes from the \eqn{y^{-3/2}} factor of the density and the
#' second from the exponent. It changes sign at the mode, which for this family
#' is \eqn{\mu\{(1 + 9\phi^2\mu^2/4)^{1/2} - 3\phi\mu/2\}}, always strictly
#' below the mean.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of strictly positive observations. The value
#'   diverges as `y` approaches zero.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(phi))`, one value per observation.
#'
#' @seealso [distrib_hess_y.InvGauss1Distrib()] for the second derivative in
#'   the response, [distrib_gradient.InvGauss1Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#'
#' # Written out.
#' all.equal(distrib_grad_y(d, y, th),
#'           -1.5 / y - (y^2 - 1^2) / (2 * 2 * 1^2 * y^2))
#'
#' # Zero at the mode, which lies below the mean.
#' mode <- 1 * (sqrt(1 + 9 * 2^2 * 1^2 / 4) - 3 * 2 * 1 / 2)
#' c(mode = mode, mean = 1, at_mode = distrib_grad_y(d, mode, th))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, InvGauss1Distrib) <- function(distrib, y, theta, ...) {
  mu <- theta[[1]]; phi <- theta[[2]]
  -1.5 / y - (y^2 - mu^2) / (2 * phi * mu^2 * y^2)
}

#' @title Inverse Gaussian Second Derivative in the Response
#' @name distrib_hess_y.InvGauss1Distrib
#' @description
#' Computes the second derivative of the inverse Gaussian log-density with
#' respect to the response, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = \dfrac{3}{2y^2}
#'       - \dfrac{1}{\phi y^3}.}
#' The mean drops out, the log-density depending on \eqn{\mu} only through
#' terms linear in \eqn{y} and constant in it. The sign changes at
#' \eqn{y = 2/(3\phi)}: the log-density is concave in the response below that
#' point and **convex above it**, so an inverse Gaussian log-likelihood is not
#' a concave function of an observation over the whole support.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param y A numeric vector of strictly positive observations. The value
#'   diverges as `y` approaches zero.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `phi` must
#'   be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `max(length(y), length(phi))`, one value
#'   per observation.
#'
#' @seealso [distrib_grad_y.InvGauss1Distrib()] for the first derivative in the
#'   response, [distrib_hessian.InvGauss1Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- invgauss1_distrib()
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#'
#' all.equal(distrib_hess_y(d, y, th), 1.5 / y^2 - 1 / (2 * y^3))
#'
#' # The mean does not enter.
#' identical(distrib_hess_y(d, y, th),
#'           distrib_hess_y(d, y, list(mu = 300, phi = 2)))
#'
#' # Concave below 2/(3 phi) and convex above it.
#' cut <- 2 / (3 * 2)
#' c(below = distrib_hess_y(d, cut / 2, th),
#'   at = distrib_hess_y(d, cut, th),
#'   above = distrib_hess_y(d, 2 * cut, th))
S7::method(distrib_hess_y, InvGauss1Distrib) <- function(distrib, y, theta, ...) {
  phi <- theta[[2]]
  1.5 / y^2 - 1 / (phi * y^3)
}

# --- CONSTRUCTOR WRAPPER ---

#' Inverse Gaussian Distribution, Mean and Dispersion
#'
#' @description
#' Builds the distribution object for the inverse Gaussian family on
#' \eqn{(0, \infty)} parametrized by its mean \eqn{\mu > 0} and a dispersion
#' \eqn{\phi > 0}, so that \eqn{\operatorname{Var}(Y) = \phi\mu^3}. This is the
#' generalized linear model parametrization, with variance function
#' \eqn{V(\mu) = \mu^3}. The returned object carries closed-form derivatives of
#' the log-density to fourth order, in the parameters and in the response, and
#' closed-form moments, so every generic of the toolkit answers without a
#' numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. Both default to the logarithm,
#' both parameters being positive.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted mean positive.
#' @param link_phi A `link` object from `linkfunctions7` for the dispersion
#'   \eqn{\phi}. Defaults to [linkfunctions7::log_link()], for the same reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, \infty)} is
#' \deqn{f(y; \mu, \phi) = \sqrt{\dfrac{1}{2\pi\phi y^3}}
#'       \exp\left\{-\dfrac{(y-\mu)^2}{2\phi\mu^2 y}\right\},}
#' the mean is \eqn{\mu}, the variance \eqn{\phi\mu^3}, the skewness
#' \eqn{3\sqrt{\phi\mu}} and the excess kurtosis \eqn{15\phi\mu}. The
#' distribution function is elementary in the standard normal one,
#' \deqn{F(q) = \Phi\left\{\sqrt{\dfrac{1}{\phi q}}
#'         \left(\dfrac{q}{\mu} - 1\right)\right\}
#'       + e^{2/(\phi\mu)}\,\Phi\left\{-\sqrt{\dfrac{1}{\phi q}}
#'         \left(\dfrac{q}{\mu} + 1\right)\right\},}
#' and the quantile function is its numerical inverse.
#'
#' The law is the first-passage time of a Brownian motion with drift, which is
#' why it is skewed however small the dispersion and why its variance grows
#' with the cube of the mean. Beside the gamma, whose variance function is
#' \eqn{\mu^2}, it puts more mass both near zero and far out.
#'
#' This is the same law as [invgauss2_distrib()], which carries the mean and
#' the *variance*, the two being related by \eqn{\sigma^2 = \phi\mu^3}. They
#' are separate families because the second parameter is a different quantity
#' in each, with its own interpretation, standard error and interval.
#'
#' # Derivatives
#'
#' The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^3},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \phi} =
#'         \dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2},}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\phi\mu^3},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\phi\phi)}\right] = -\dfrac{1}{2\phi^2},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\phi)}\right] = 0.}
#' The zero off-diagonal makes the mean and the dispersion orthogonal, so the
#' mean equation can be fitted with the dispersion held at any value without
#' biasing it. Note that the expected curvature in \eqn{\phi} is
#' \eqn{-1/(2\phi^2)} whatever the mean, so the dispersion is estimated with
#' the same precision at every scale of the response.
#'
#' The observed Hessian is a different matter, and the method page says so at
#' length: neither diagonal entry is negative at every observation. That is the
#' practical reason to fit this family by [fisher_scoring()].
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.InvGauss1Distrib()] and [distrib_deriv4.InvGauss1Distrib()],
#' as are the derivatives in the response,
#' [distrib_grad_y.InvGauss1Distrib()] and [distrib_hess_y.InvGauss1Distrib()].
#' The derivatives of the *distribution* function in the parameters are closed
#' form too, the expression above being elementary.
#'
#' # Estimation
#'
#' Both maximum likelihood estimates are available in closed form, which few
#' two-parameter families offer:
#' \deqn{\hat\mu = \bar y, \qquad
#'       \hat\phi = \dfrac{1}{n}\sum_i \dfrac{1}{y_i} - \dfrac{1}{\bar y}.}
#' [fit_distrib()] reaches them numerically on the link scale, and the example
#' below checks both against the sample.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\phi > 0} the dispersion, with \eqn{\operatorname{Var}(Y) = \phi\mu^3}.
#' \eqn{\Phi} is the standard normal distribution function. \eqn{\eta} is a
#' parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `InvGauss1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"invgauss1"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "phi")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [invgauss2_distrib()] for the same law in the mean and the variance;
#' [gamma1_distrib()] for the other positive generalized linear model family,
#' with variance function \eqn{\mu^2}; [lognormal1_distrib()] for a third;
#' [pig1_distrib()] for the counts this law mixes a Poisson into;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [InvGauss1Distrib] for the class.
#'
#' @references
#' Chhikara, R. S. and Folks, J. L. (1989). *The Inverse Gaussian
#' Distribution: Theory, Methodology, and Applications*. Marcel Dekker,
#' New York.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom statmod dinvgauss pinvgauss qinvgauss rinvgauss
#'
#' @examples
#' d <- invgauss1_distrib()
#' d
#'
#' # The density is statmod::dinvgauss at this parametrization.
#' y <- c(0.5, 1, 2)
#' th <- list(mu = 1, phi = 2)
#' all.equal(distrib_pdf(d, y, th),
#'           statmod::dinvgauss(y, mean = 1, dispersion = 2))
#'
#' # Moments: variance phi mu^3, skewness 3 sqrt(phi mu), kurtosis 15 phi mu.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(2 * 1^3, 3 * sqrt(2 * 1), 15 * 2 * 1)
#'
#' # Fitting recovers the closed-form maximum likelihood estimates.
#' set.seed(3)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = mean(z), phi = mean(1 / z) - 1 / mean(z)))
#'
#' # The mean and the dispersion are orthogonal: the mixed entry is 0.
#' distrib_expected_hessian(d, 1, th)$mu_phi
#'
#' @export
invgauss1_distrib <- function(link_mu = log_link(), link_phi = log_link()) {
  InvGauss1Distrib(
    distrib_name = "invgauss1", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu", "phi"), params_interpretation = c(mu = "mean", phi = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), phi = c(0, Inf)),
    link_params = list(mu = link_mu, phi = link_phi)
  )
}
