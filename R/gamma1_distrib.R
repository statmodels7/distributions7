#' @include distrib.R generics.R
NULL

#' @title Gamma Distribution Class, Mean and Dispersion
#' @name Gamma1Distrib
#'
#' @description
#' The S7 class of the gamma family parametrized by its mean \eqn{\mu > 0} and
#' a dispersion \eqn{\phi > 0}, so that \eqn{\operatorname{Var}(Y) = \phi\mu^2}.
#' The shape is \eqn{1/\phi} and the rate \eqn{1/(\phi\mu)}. It inherits from
#' `continuous_distrib`, so it answers every generic of the `distrib`
#' contract; the eleven methods listed below are registered on it directly and
#' everything else comes from the parent.
#'
#' Build one with [gamma1_distrib()], which supplies the two link functions and
#' fills the properties in. This page documents the raw S7 constructor, which
#' takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Gamma1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [gamma1_distrib()] they hold `"gamma1"`, `"univariate"`,
#'   `c(0, Inf)`, `c("mu", "phi")`, the interpretations
#'   `c(mu = "mean", phi = "dispersion")`, `2`, the domain \eqn{(0, \infty)}
#'   for both parameters, and the two links.
#'
#' @seealso [gamma1_distrib()] to build one;
#'   [gamma2_distrib()] for the same law in mean and variance;
#'   [distrib_pdf.Gamma1Distrib()] and [distrib_gradient.Gamma1Distrib()] for
#'   the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Gamma1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gamma1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gamma1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gamma1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gamma1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gamma1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gamma1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gamma1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gamma1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gamma1Distrib],
#'   [`distrib_rng()`][distrib_rng.Gamma1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- gamma1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # The dispersion multiplies the variance function V(mu) = mu^2, so the
#' # coefficient of variation is sqrt(phi) whatever the mean.
#' sqrt(variance(d, list(mu = 3, phi = 0.5))) / 3
#' sqrt(variance(d, list(mu = 300, phi = 0.5))) / 300
Gamma1Distrib <- S7::new_class("Gamma1Distrib", parent = continuous_distrib)

#' The Shape and Rate a Mean and Dispersion Imply
#'
#' @description
#' Converts the `gamma1` parameters to the shape and rate [stats::dgamma()]
#' takes: \eqn{a = 1/\phi} and \eqn{b = 1/(\phi\mu)}. Under that pairing the
#' mean \eqn{a/b} is \eqn{\mu} and the variance \eqn{a/b^2} is \eqn{\phi\mu^2}.
#' Nothing is validated; a non-positive `phi` or `mu` propagates as an infinite
#' or negative value to the caller.
#'
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector. The two are used elementwise, so components of different lengths
#'   recycle in the usual way.
#'
#' @return A list of two numeric vectors, `shape` and `rate`, of the lengths
#'   the arithmetic produces.
#'
#' @seealso [gamma1_distrib()] for the family and
#'   [distrib_pdf.Gamma1Distrib()] for the density this feeds.
#'
#' @keywords internal
gamma1_shape_rate <- function(theta) {
  list(shape = 1 / theta[[2]], rate = 1 / (theta[[2]] * theta[[1]]))
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gamma Probability Density Function in Mean and Dispersion
#' @name distrib_pdf.Gamma1Distrib
#' @description
#' Computes the gamma density
#' \deqn{f(y; \mu, \phi) = \dfrac{y^{1/\phi - 1}\,e^{-y/(\phi\mu)}}
#'       {(\phi\mu)^{1/\phi}\,\Gamma(1/\phi)}, \qquad y > 0,}
#' by calling [stats::dgamma()] at shape \eqn{a = 1/\phi} and rate
#' \eqn{b = 1/(\phi\mu)}. With `log = TRUE` the logarithm is formed inside
#' `dgamma()` and stays finite where the density itself underflows.
#'
#' The density is unbounded at the origin when \eqn{\phi > 1}, where the shape
#' falls below one; it is flat there at \eqn{\phi = 1}, the exponential case,
#' and vanishes at the origin for \eqn{\phi < 1}.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
#' @param y A numeric vector of observations. The support is \eqn{(0, \infty)};
#'   a negative value gives 0 and `y = 0` gives 0, `Inf` or the rate according
#'   to whether \eqn{\phi} is below, above or equal to 1.
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
#' @seealso [distrib_cdf.Gamma1Distrib()] for the distribution function,
#'   [distrib_gradient.Gamma1Distrib()] for the derivatives of the log-density,
#'   [gamma1_shape_rate()] for the conversion this uses, and [distrib_pdf()]
#'   for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#'
#' # The method is stats::dgamma at shape 1/phi and rate 1/(phi mu).
#' all.equal(distrib_pdf(d, y, th),
#'           dgamma(y, shape = 1 / 0.5, rate = 1 / (0.5 * 3)))
#'
#' # At phi = 1 the shape is 1 and the gamma is the exponential.
#' all.equal(distrib_pdf(d, y, list(mu = 3, phi = 1)),
#'           distrib_pdf(exponential_distrib(), y, list(mu = 3)))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(1, 3, 9), phi = 0.5))
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e4, th)
#' distrib_pdf(d, 1e4, th, log = TRUE)
S7::method(distrib_pdf, Gamma1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  sr <- gamma1_shape_rate(theta)
  stats::dgamma(y, shape = sr$shape, rate = sr$rate, log = log)
}

#' @title Gamma Cumulative Distribution Function in Mean and Dispersion
#' @name distrib_cdf.Gamma1Distrib
#' @description
#' Computes the gamma distribution function, the regularized incomplete gamma
#' function
#' \deqn{F(q; \mu, \phi) = \dfrac{1}{\Gamma(a)}\int_0^{bq} t^{a-1}e^{-t}\,dt,
#'       \qquad a = \dfrac{1}{\phi}, \quad b = \dfrac{1}{\phi\mu},}
#' by calling [stats::pgamma()] at that shape and rate. Both tails are
#' available exactly: `lower.tail = FALSE` evaluates \eqn{1 - F} without
#' forming the difference, and `log.p = TRUE` returns a logarithm that stays
#' finite where the probability itself underflows to zero.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' @seealso [distrib_quantile.Gamma1Distrib()] for the inverse,
#'   [distrib_pdf.Gamma1Distrib()] for the density,
#'   [distrib_grad_cdf()] for the derivatives of this function in the
#'   parameters, which the gamma takes by finite difference because the
#'   derivative of an incomplete gamma in its shape is hypergeometric, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' th <- list(mu = 3, phi = 0.5)
#'
#' # The method is stats::pgamma at the implied shape and rate.
#' all.equal(distrib_cdf(d, c(1, 3, 5), th),
#'           pgamma(c(1, 3, 5), shape = 2, rate = 1 / 1.5))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 5, th) + distrib_cdf(d, 5, th, lower.tail = FALSE)
#'
#' # The gamma is right skewed, so less than half the mass lies below the mean.
#' distrib_cdf(d, 3, th)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 3000, th, lower.tail = FALSE)
#' distrib_cdf(d, 3000, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, Gamma1Distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log.p = FALSE, ...) {
  sr <- gamma1_shape_rate(theta)
  stats::pgamma(q, shape = sr$shape, rate = sr$rate,
                lower.tail = lower.tail, log.p = log.p)
}

#' @title Gamma Quantile Function in Mean and Dispersion
#' @name distrib_quantile.Gamma1Distrib
#' @description
#' Computes the gamma quantile function, the inverse of the regularized
#' incomplete gamma function in its argument, by calling [stats::qgamma()] at
#' shape \eqn{a = 1/\phi} and rate \eqn{b = 1/(\phi\mu)}. The gamma
#' distribution function is strictly increasing on \eqn{(0, \infty)}, so the
#' round trip through [distrib_cdf.Gamma1Distrib()] returns `p`.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' @seealso [distrib_cdf.Gamma1Distrib()], which this inverts;
#'   [distrib_rng.Gamma1Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' th <- list(mu = 3, phi = 0.5)
#'
#' # A central 95 percent interval, visibly asymmetric about the mean of 3.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The median falls below the mean, the gamma being right skewed.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, Gamma1Distrib) <- function(distrib, p, theta,
                                                         lower.tail = TRUE,
                                                         log.p = FALSE, ...) {
  sr <- gamma1_shape_rate(theta)
  stats::qgamma(p, shape = sr$shape, rate = sr$rate,
                lower.tail = lower.tail, log.p = log.p)
}

#' @title Gamma Random Number Generator in Mean and Dispersion
#' @name distrib_rng.Gamma1Distrib
#' @description
#' Draws `n` independent gamma variates by calling [stats::rgamma()] at shape
#' \eqn{a = 1/\phi} and rate \eqn{b = 1/(\phi\mu)}, so the draws come from R's
#' own gamma generator and depend on `.Random.seed` in the usual way. The
#' ratio-of-uniforms fallback the base class supplies is bypassed.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. Both
#'   must be strictly positive.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_quantile.Gamma1Distrib()] for the inverse-transform route,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#'
#' # Same generator as stats::rgamma at the implied shape and rate.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 3, phi = 0.5))
#' set.seed(2)
#' identical(a, rgamma(3, shape = 2, rate = 1 / 1.5))
#'
#' # The sample moments recover the parameters: the mean directly, and the
#' # dispersion as the squared coefficient of variation.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, phi = 0.5))
#' c(mu = mean(z), phi = var(z) / mean(z)^2)
S7::method(distrib_rng, Gamma1Distrib) <- function(distrib, n, theta) {
  sr <- gamma1_shape_rate(theta)
  stats::rgamma(n, shape = sr$shape, rate = sr$rate)
}

#' @title Gamma Score in Mean and Dispersion
#' @name distrib_gradient.Gamma1Distrib
#' @description
#' Computes the first derivatives of the gamma log-density with respect to
#' \eqn{\mu} and \eqn{\phi}, one value per observation, in closed form. With
#' \eqn{s = 1/\phi} the shape and \eqn{z = y/\mu},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^2},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \phi} =
#'       -s^2\left\{\log s + 1 - \psi(s) + \log z - z\right\},}
#' with \eqn{\psi} the digamma function. The first is the score of a gamma
#' generalized linear model, the residual divided by the variance function
#' \eqn{\phi\mu^2}.
#'
#' The dispersion component is a difference of two quantities that agree to
#' leading order as \eqn{s} grows: \eqn{\log s - \psi(s)} and
#' \eqn{\log z - (z - 1)} each go to zero at the boundary this family tends
#' towards. The kernel computes each of them as a polygamma minus its own
#' asymptote, so the digits survive at large shape where the direct difference
#' loses them.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' \eqn{\phi > 0} the dispersion, with \eqn{\operatorname{Var}(Y) = \phi\mu^2}.
#' \eqn{\psi} is the digamma function, \eqn{\psi = (\log\Gamma)'}.
#'
#' @seealso [distrib_hessian.Gamma1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Gamma1Distrib()] for their expectation,
#'   [distrib_grad_y.Gamma1Distrib()] for the derivative in the response, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The generalized linear model score: residual over the variance function.
#' all.equal(g$mu, (y - 3) / (0.5 * 3^2))
#'
#' # The dispersion component, written out with the digamma function.
#' s <- 1 / 0.5
#' z <- y / 3
#' all.equal(g$phi, -s^2 * (log(s) + 1 - digamma(s) + log(z) - z))
#'
#' # The mean component vanishes at y = mu, whatever the dispersion.
#' distrib_gradient(d, 3, list(mu = 3, phi = c(0.1, 0.5, 2)))$mu
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' set.seed(4)
#' zz <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, zz)
#' vapply(distrib_gradient(d, zz, as.list(coef(fit))), sum, numeric(1))
S7::method(distrib_gradient, Gamma1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"), ...,
                                                         threads = 1L) {
  gamma1_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Observed Hessian in Mean and Dispersion
#' @name distrib_hessian.Gamma1Distrib
#' @description
#' Computes the three distinct second derivatives of the gamma log-density
#' with respect to \eqn{\mu} and \eqn{\phi}, one value per observation, in
#' closed form. With \eqn{s = 1/\phi} and \eqn{z = y/\mu},
#' \deqn{\ell^{(\mu\mu)} = \dfrac{s(1 - 2z)}{\mu^2}, \qquad
#'       \ell^{(\mu\phi)} = \dfrac{s^2(1 - z)}{\mu}, \qquad
#'       \ell^{(\phi\phi)} = s^4\left\{\dfrac{1}{s} - \psi'(s)\right\}
#'         + 2s^3\left\{\log s + 1 - \psi(s) + \log z - z\right\}.}
#' Every derivative in \eqn{\phi} is the corresponding derivative in \eqn{s}
#' carried across by the one-variable chain rule, with
#' \eqn{s' = -s^2} and \eqn{s'' = 2s^3}, so each polygamma function is
#' evaluated once.
#'
#' The curvature in \eqn{\mu} turns **positive** wherever \eqn{y < \mu/2}, so
#' the observed information of a gamma is not positive definite at every data
#' point. Its expectation is, which is one reason a gamma model is fitted by
#' Fisher scoring; see [distrib_expected_hessian.Gamma1Distrib()].
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' @return A named list of three numeric vectors, `mu_mu`, `mu_phi` and
#'   `phi_phi`, each of length `max(length(y), length(mu), length(phi))`. The
#'   three name the distinct entries of a symmetric \eqn{2 \times 2} matrix per
#'   observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives. \eqn{\psi}
#' and \eqn{\psi'} are the digamma and trigamma functions.
#'
#' @seealso [distrib_gradient.Gamma1Distrib()] for the score,
#'   [distrib_expected_hessian.Gamma1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Gamma1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#' h <- distrib_hessian(d, y, th)
#'
#' # The curvature in mu, written out, and positive at y = 1 < mu/2 = 1.5.
#' s <- 1 / 0.5
#' z <- y / 3
#' all.equal(h$mu_mu, s * (1 - 2 * z) / 3^2)
#' h$mu_mu
#'
#' # The mixed entry vanishes at y = mu.
#' distrib_hessian(d, 3, th)$mu_phi
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 3, phi = 0.5 + eps))$phi
#' dn <- distrib_gradient(d, y, list(mu = 3, phi = 0.5 - eps))$phi
#' all.equal((up - dn) / (2 * eps), h$phi_phi, tolerance = 1e-5)
S7::method(distrib_hessian, Gamma1Distrib) <- function(distrib, y, theta,
                                                        scale = c("parameter", "link"), ...,
                                                        threads = 1L) {
  gamma1_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Expected Hessian in Mean and Dispersion
#' @name distrib_expected_hessian.Gamma1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. With \eqn{s = 1/\phi},
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\phi\mu^2},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\phi)}\right] = 0,
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\phi\phi)}\right] =
#'         s^4\left\{\dfrac{1}{s} - \psi'(s)\right\}.}
#' They follow from \eqn{\mathbb{E}[Y] = \mu} and
#' \eqn{\mathbb{E}[\log(Y/\mu)] = \psi(s) - \log s}, the second of which is
#' exactly what makes the score in \eqn{\phi} have mean zero. The quantity
#' \eqn{1/s - \psi'(s)} is negative for every \eqn{s > 0}, so the pure
#' dispersion entry is negative and the information is positive definite
#' everywhere, where the observed Hessian is not.
#'
#' The zero off-diagonal says the mean and the dispersion are orthogonal. That
#' orthogonality is the reason this parametrization is the one a generalized
#' linear model uses: the mean equation can be fitted with the dispersion held
#' at any value without biasing it.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' @return A named list of three numeric vectors, `mu_mu`, `mu_phi` and
#'   `phi_phi`, each of length `max(length(y), length(mu), length(phi))` and
#'   constant within itself when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The gamma is a
#' regular family, so the second Bartlett identity holds and this equals the
#' variance of the score. \eqn{\psi'} is the trigamma function.
#'
#' @seealso [distrib_hessian.Gamma1Distrib()] for the observed quantity this is
#'   the expectation of, [fisher_scoring()], which inverts it at each step, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' th <- list(mu = 3, phi = 0.5)
#'
#' # The three constants, one value per observation.
#' lapply(distrib_expected_hessian(d, c(1, 3, 5), th), unique)
#'
#' # The dispersion entry, written out with the trigamma function.
#' s <- 1 / 0.5
#' s^4 * (1 / s - trigamma(s))
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
#'
#' # Negative definite at every parameter setting, where the observed Hessian
#' # is positive in mu wherever y < mu/2.
#' phi <- c(0.05, 0.5, 5)
#' vapply(phi, function(p) {
#'   s <- 1 / p
#'   s^4 * (1 / s - trigamma(s))
#' }, numeric(1))
S7::method(distrib_expected_hessian, Gamma1Distrib) <- function(distrib, y, theta,
                                                                 scale = c("parameter", "link"),
                                                                 approx = c("bartlett", "integrate", "mc", "opg"),
                                                                 nsim = 10000, ...,
                                                                 threads = 1L) {
  gamma1_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Third-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv3.Gamma1Distrib
#' @description
#' Computes the four distinct third derivatives of the gamma log-density with
#' respect to \eqn{\mu} and \eqn{\phi}, in closed form. With \eqn{s = 1/\phi}
#' and \eqn{z = y/\mu},
#' \deqn{\ell^{(\mu\mu\mu)} = \dfrac{s(6z - 2)}{\mu^3}, \qquad
#'       \ell^{(\mu\mu\phi)} = -\dfrac{s^2(1 - 2z)}{\mu^2}, \qquad
#'       \ell^{(\mu\phi\phi)} = \dfrac{2s^3(z - 1)}{\mu},}
#' and the pure dispersion component follows from the derivatives in \eqn{s} by
#' Faa di Bruno on \eqn{s(\phi) = 1/\phi},
#' \deqn{\ell^{(\phi\phi\phi)} = f_3 (s')^3 + 3 f_2 s' s'' + f_1 s''',}
#' with \eqn{s' = -s^2}, \eqn{s'' = 2s^3}, \eqn{s''' = -6s^4} and \eqn{f_k} the
#' \eqn{k}th derivative of \eqn{\ell} in \eqn{s}.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{z} with 1 and \eqn{\log z} with \eqn{\psi(s) - \log s}: \eqn{f_1} and
#' the two components odd in \eqn{z - 1} vanish. Both routes are closed form,
#' so no quadrature is run and `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\psi} is the digamma function.
#'
#' @seealso [distrib_hessian.Gamma1Distrib()] for the order below and
#'   [distrib_deriv4.Gamma1Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic and for the numerical route a family
#'   without a closed form takes.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # The three components involving mu, written out.
#' s <- 1 / 0.5
#' z <- y / 3
#' all.equal(d3$mu_mu_mu, s * (6 * z - 2) / 3^3)
#' all.equal(d3$mu_mu_phi, -s^2 * (1 - 2 * z) / 3^2)
#' all.equal(d3$mu_phi_phi, 2 * s^3 * (z - 1) / 3)
#'
#' # Expected values: the components odd in z - 1 vanish.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 3, phi = 0.5 + eps))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 3, phi = 0.5 - eps))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_phi, tolerance = 1e-5)
S7::method(distrib_deriv3, Gamma1Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                       scale = c("parameter", "link"),
                                                       approx = c("integrate", "bartlett", "mc", "opg"),
                                                       nsim = 10000, ...,
                                                       threads = 1L) {
  if (expected) {
    gamma1_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gamma1_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gamma Fourth-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv4.Gamma1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the gamma log-density with
#' respect to \eqn{\mu} and \eqn{\phi}, in closed form. With \eqn{s = 1/\phi}
#' and \eqn{z = y/\mu},
#' \deqn{\ell^{(\mu\mu\mu\mu)} = \dfrac{s(6 - 24z)}{\mu^4}, \qquad
#'       \ell^{(\mu\mu\mu\phi)} = -\dfrac{s^2(6z - 2)}{\mu^3}, \qquad
#'       \ell^{(\mu\mu\phi\phi)} = \dfrac{2s^3(1 - 2z)}{\mu^2}, \qquad
#'       \ell^{(\mu\phi\phi\phi)} = -\dfrac{6s^4(z - 1)}{\mu},}
#' and the pure dispersion component follows from the derivatives in \eqn{s} by
#' Faa di Bruno on \eqn{s(\phi) = 1/\phi},
#' \deqn{\ell^{(\phi^4)} = f_4 (s')^4 + 6 f_3 (s')^2 s'' + 3 f_2 (s'')^2
#'       + 4 f_2 s' s''' + f_1 s'''',}
#' with \eqn{s' = -s^2}, \eqn{s'' = 2s^3}, \eqn{s''' = -6s^4},
#' \eqn{s'''' = 24s^5} and \eqn{f_k} the \eqn{k}th derivative of \eqn{\ell} in
#' \eqn{s}. The four \eqn{f_k} are polygamma functions of \eqn{s}, each
#' computed as a polygamma minus its own leading asymptote so that the digits
#' survive at large shape.
#'
#' With `expected = TRUE` the expectations are returned, obtained by replacing
#' \eqn{z} with 1 and \eqn{\log z} with \eqn{\psi(s) - \log s}. Both routes are
#' closed form, so `approx` and `nsim` are ignored.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
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
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\psi} is the digamma
#' function.
#'
#' @seealso [distrib_deriv3.Gamma1Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.Gamma1Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The four components involving mu, written out.
#' s <- 1 / 0.5
#' z <- y / 3
#' all.equal(d4$mu_mu_mu_mu, s * (6 - 24 * z) / 3^4)
#' all.equal(d4$mu_phi_phi_phi, -6 * s^4 * (z - 1) / 3)
#'
#' # Expected values: the components odd in z - 1 vanish.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 3, phi = 0.5 + eps))$mu_mu_phi
#' dn <- distrib_deriv3(d, y, list(mu = 3, phi = 0.5 - eps))$mu_mu_phi
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_phi_phi, tolerance = 1e-4)
S7::method(distrib_deriv4, Gamma1Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                       scale = c("parameter", "link"),
                                                       approx = c("integrate", "bartlett", "mc", "opg"),
                                                       nsim = 10000, ...,
                                                       threads = 1L) {
  if (expected) {
    gamma1_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gamma1_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gamma First Derivative in the Response, Mean and Dispersion
#' @name distrib_grad_y.Gamma1Distrib
#' @description
#' Computes the first derivative of the gamma log-density with respect to the
#' response, in closed form at the implied shape \eqn{a = 1/\phi} and rate
#' \eqn{b = 1/(\phi\mu)}:
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{a - 1}{y} - b.}
#' It changes sign at the mode \eqn{y = (a-1)/b}, so it is positive below the
#' mode and negative above it. At \eqn{\phi = 1} the shape is 1, the first term
#' drops out and the derivative is the constant \eqn{-1/\mu} of an exponential.
#' At \eqn{\phi > 1} the shape falls below 1, the density has no interior mode
#' and the derivative is negative throughout.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
#' @param y A numeric vector of strictly positive observations. At `y = 0` the
#'   value is infinite unless the shape is exactly 1.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(phi))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Gamma1Distrib()] for the second derivative in the
#'   response, [distrib_gradient.Gamma1Distrib()] for the score in the
#'   parameters, [gamma1_shape_rate()] for the conversion this uses, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#'
#' # Written out at the implied shape and rate.
#' all.equal(distrib_grad_y(d, y, th), (2 - 1) / y - 1 / 1.5)
#'
#' # Zero at the mode, positive below it and negative above.
#' mode <- (1 / 0.5 - 1) / (1 / (0.5 * 3))
#' c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#' distrib_grad_y(d, c(0.5, 5), th)
#'
#' # At phi = 1 the family is exponential and the derivative is constant.
#' distrib_grad_y(d, y, list(mu = 3, phi = 1))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, Gamma1Distrib) <- function(distrib, y, theta, ...) {
  sr <- gamma1_shape_rate(theta)
  (sr$shape - 1) / y - sr$rate
}

#' @title Gamma Second Derivative in the Response, Mean and Dispersion
#' @name distrib_hess_y.Gamma1Distrib
#' @description
#' Computes the second derivative of the gamma log-density with respect to the
#' response, in closed form at the implied shape \eqn{a = 1/\phi}:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{a - 1}{y^2}.}
#' The rate drops out, the log-density being linear in \eqn{y} apart from the
#' \eqn{(a-1)\log y} term. The sign follows the shape: the log-density is
#' concave in the response for \eqn{\phi < 1}, exactly flat at \eqn{\phi = 1},
#' where the family is exponential, and convex for \eqn{\phi > 1}.
#'
#' @param distrib A `Gamma1Distrib` object, from [gamma1_distrib()].
#' @param y A numeric vector of strictly positive observations. At `y = 0` the
#'   value is infinite unless the shape is exactly 1.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read, the rate
#'   having cancelled. `phi` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `max(length(y), length(phi))`, one value
#'   per observation.
#'
#' @seealso [distrib_grad_y.Gamma1Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Gamma1Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- gamma1_distrib()
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#'
#' all.equal(distrib_hess_y(d, y, th), -(1 / 0.5 - 1) / y^2)
#'
#' # The mean does not enter: only the shape survives the second derivative.
#' identical(distrib_hess_y(d, y, th),
#'           distrib_hess_y(d, y, list(mu = 300, phi = 0.5)))
#'
#' # Concave below phi = 1, flat at it, convex above.
#' vapply(c(0.5, 1, 2), function(p) distrib_hess_y(d, 3, list(mu = 3, phi = p)),
#'        numeric(1))
S7::method(distrib_hess_y, Gamma1Distrib) <- function(distrib, y, theta, ...) {
  sr <- gamma1_shape_rate(theta)
  -(sr$shape - 1) / y^2
}


#' Gamma Distribution, Mean and Dispersion
#'
#' @description
#' Builds the distribution object for the gamma family parametrized by its mean
#' \eqn{\mu > 0} and a dispersion \eqn{\phi > 0}, so that
#' \eqn{\operatorname{Var}(Y) = \phi\mu^2}. This is the parametrization a
#' generalized linear model uses. The returned object carries closed-form
#' derivatives of the log-density to fourth order, in the parameters and in the
#' response, and closed-form moments, so every generic of the toolkit answers
#' without a numerical fallback.
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
#' \deqn{f(y; \mu, \phi) = \dfrac{y^{1/\phi - 1}\,e^{-y/(\phi\mu)}}
#'       {(\phi\mu)^{1/\phi}\,\Gamma(1/\phi)},}
#' the gamma at shape \eqn{a = 1/\phi} and rate \eqn{b = 1/(\phi\mu)}. The mean
#' is \eqn{\mu}, the variance \eqn{\phi\mu^2}, the skewness \eqn{2\sqrt{\phi}}
#' and the excess kurtosis \eqn{6\phi}. The coefficient of variation is
#' \eqn{\sqrt{\phi}} at every mean, so \eqn{\phi} measures relative rather than
#' absolute spread.
#'
#' The variance function is \eqn{V(\mu) = \mu^2} and \eqn{\phi} is the
#' dispersion multiplying it, which makes the score in \eqn{\mu} the
#' generalized linear model score \eqn{(y-\mu)/(\phi\mu^2)}. At \eqn{\phi = 1}
#' the shape is 1 and the family is the exponential, [exponential_distrib()];
#' at \eqn{\phi \to 0} it tends to a Gaussian with variance \eqn{\phi\mu^2}.
#'
#' This is the same law as [gamma2_distrib()], which carries the mean and the
#' *variance*, the two being related by \eqn{\sigma^2 = \phi\mu^2}. They are
#' separate families because the second parameter is a different quantity in
#' each, with its own interpretation, standard error and interval.
#'
#' # Derivatives
#'
#' Writing \eqn{s = 1/\phi} and \eqn{z = y/\mu}, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\phi\mu^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \phi} =
#'       -s^2\left\{\log s + 1 - \psi(s) + \log z - z\right\},}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{1}{\phi\mu^2}, \quad
#'       \mathbb{E}\left[\ell^{(\mu\phi)}\right] = 0, \quad
#'       \mathbb{E}\left[\ell^{(\phi\phi)}\right] =
#'         s^4\left\{\dfrac{1}{s} - \psi'(s)\right\}.}
#' The zero off-diagonal makes the mean and the dispersion orthogonal, so the
#' mean equation can be fitted with the dispersion held at any value without
#' biasing it.
#'
#' Every derivative in \eqn{\phi} is the corresponding derivative in \eqn{s}
#' carried across by the one-variable chain rule, with \eqn{s' = -s^2},
#' \eqn{s'' = 2s^3}, \eqn{s''' = -6s^4} and \eqn{s'''' = 24s^5}, so each
#' polygamma function is evaluated once. Two quantities in those derivatives
#' cancel as \eqn{s} grows, \eqn{\log s - \psi(s)} and its polygamma analogues;
#' each is computed as a polygamma minus its own leading asymptote, so the
#' digits survive at the large shape a nearly Gaussian gamma reaches.
#'
#' Third and fourth orders are closed form as well, observed and expected, in
#' [distrib_deriv3.Gamma1Distrib()] and [distrib_deriv4.Gamma1Distrib()], as
#' are the derivatives in the response, [distrib_grad_y.Gamma1Distrib()] and
#' [distrib_hess_y.Gamma1Distrib()]. The derivatives of the *distribution*
#' function in the parameters have no elementary form here, the derivative of
#' an incomplete gamma in its shape being hypergeometric, and are taken by
#' finite difference on the analytic cdf.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. The mean has
#' the closed-form estimate \eqn{\hat\mu = \bar y}. The shape \eqn{\hat s} then
#' solves
#' \deqn{\log \hat s - \psi(\hat s) = -\dfrac{1}{n}\sum_i \log(y_i/\bar y),}
#' which has no closed form and is reached numerically, and
#' \eqn{\hat\phi = 1/\hat s}. The method of moments starting value is the
#' squared coefficient of variation, and the example below shows it landing
#' beside the estimate.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\phi > 0} the dispersion. \eqn{\psi} and \eqn{\psi'} are the digamma
#' and trigamma functions. \eqn{\eta} is a parameter on the unconstrained scale
#' of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Gamma1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"gamma1"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "phi")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [gamma2_distrib()] for the same law in the mean and the variance;
#' [exponential_distrib()] for the case \eqn{\phi = 1};
#' [chisq_distrib()] and [gengamma1_distrib()] for relatives;
#' [invgauss1_distrib()] and [lognormal1_distrib()] for other positive
#' families with a multiplicative variance function; [fit_distrib()] to
#' estimate the parameters; [check_distrib()] to validate a family of your own
#' against the same battery this one passes; [Gamma1Distrib] for the class.
#'
#' @references
#' McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
#' edition, Chapter 8. Chapman and Hall, London.
#'
#' @examples
#' d <- gamma1_distrib()
#' d
#'
#' # The density is stats::dgamma at shape 1/phi and rate 1/(phi mu).
#' y <- c(1, 3, 5)
#' th <- list(mu = 3, phi = 0.5)
#' all.equal(distrib_pdf(d, y, th), dgamma(y, shape = 2, rate = 1 / 1.5))
#'
#' # Moments in closed form: skewness 2 sqrt(phi), excess kurtosis 6 phi.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(2 * sqrt(0.5), 6 * 0.5)
#'
#' # Fitting recovers the parameters, the mean exactly at the sample mean.
#' set.seed(4)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted  = coef(fit),
#'       moments = c(mu = mean(z), phi = var(z) / mean(z)^2))
#'
#' # At phi = 1 the family is the exponential.
#' all.equal(distrib_pdf(d, y, list(mu = 3, phi = 1)),
#'           distrib_pdf(exponential_distrib(), y, list(mu = 3)))
#'
#' @export
gamma1_distrib <- function(link_mu = log_link(), link_phi = log_link()) {
  Gamma1Distrib(
    distrib_name = "gamma1",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "phi"),
    params_interpretation = c(mu = "mean", phi = "dispersion"),
    n_params = 2,
    params_bounds = list(mu = c(0, Inf), phi = c(0, Inf)),
    link_params = list(mu = link_mu, phi = link_phi)
  )
}
