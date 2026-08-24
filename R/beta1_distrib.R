#' @include distrib.R generics.R
NULL

#' @title Beta Distribution Class, Mean and Precision
#' @name Beta1Distrib
#'
#' @description
#' The S7 class of the beta family on \eqn{(0, 1)} parametrized by its mean
#' \eqn{\mu \in (0, 1)} and a precision \eqn{\phi > 0}, so that the shapes are
#' \eqn{\alpha = \mu\phi} and \eqn{\beta = (1-\mu)\phi} and the variance is
#' \eqn{\mu(1-\mu)/(\phi+1)}. It inherits from `continuous_distrib`, so it
#' answers every generic of the `distrib` contract; the eleven methods listed
#' below are registered on it directly and everything else comes from the
#' parent.
#'
#' Build one with [beta1_distrib()], which supplies the two link functions and
#' fills the properties in. This page documents the raw S7 constructor, which
#' takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Beta1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [beta1_distrib()] they hold `"beta1"`, `"univariate"`,
#'   `c(0, 1)`, `c("mu", "phi")`, the interpretations
#'   `c(mu = "mean", phi = "precision")`, `2`, the domains \eqn{(0, 1)} and
#'   \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [beta1_distrib()] to build one;
#'   [beta2_distrib()] for the same law in mean and variance;
#'   [distrib_pdf.Beta1Distrib()] and [distrib_gradient.Beta1Distrib()] for the
#'   closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Beta1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Beta1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Beta1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Beta1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Beta1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Beta1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Beta1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Beta1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Beta1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Beta1Distrib],
#'   [`distrib_rng()`][distrib_rng.Beta1Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- beta1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@params_bounds
#'
#' # phi is a precision: the variance falls as it grows, at a fixed mean.
#' vapply(c(1, 5, 50), function(p) variance(d, list(mu = 0.4, phi = p)),
#'        numeric(1))
Beta1Distrib <- S7::new_class("Beta1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta Probability Density Function in Mean and Precision
#' @name distrib_pdf.Beta1Distrib
#' @description
#' Computes the beta density
#' \deqn{f(y; \mu, \phi) = \dfrac{\Gamma(\phi)}{\Gamma(\mu\phi)\,
#'       \Gamma((1-\mu)\phi)}\, y^{\mu\phi - 1}(1-y)^{(1-\mu)\phi - 1},
#'       \qquad 0 < y < 1,}
#' by calling [stats::dbeta()] at `shape1 = mu * phi` and
#' `shape2 = (1 - mu) * phi`. With `log = TRUE` the logarithm is formed inside
#' `dbeta()` and stays finite where the density itself underflows.
#'
#' The shape of the density is governed by whether each shape parameter
#' exceeds one. It is unbounded at 0 when \eqn{\mu\phi < 1} and at 1 when
#' \eqn{(1-\mu)\phi < 1}; at \eqn{\mu = 1/2} and \eqn{\phi = 2} both shapes are
#' 1 and the density is the uniform.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations. The support is \eqn{(0, 1)}; a
#'   value outside \eqn{[0, 1]} gives 0, and the endpoints give 0, a finite
#'   value or `Inf` according to the shapes.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(phi))`, one value per observation.
#'
#' @seealso [distrib_cdf.Beta1Distrib()] for the distribution function,
#'   [distrib_gradient.Beta1Distrib()] for the derivatives of the log-density,
#'   [distrib_pdf.Beta2Distrib()] for the same density in the variance, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#'
#' # The method is stats::dbeta at shape1 = mu phi and shape2 = (1 - mu) phi.
#' all.equal(distrib_pdf(d, y, th), dbeta(y, 0.4 * 5, 0.6 * 5))
#'
#' # Both shapes are 1 at mu = 1/2, phi = 2, where the beta is the uniform.
#' distrib_pdf(d, y, list(mu = 0.5, phi = 2))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(0.2, 0.5, 0.8), phi = 5))
#'
#' # Near the boundary the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e-40, list(mu = 0.4, phi = 50))
#' distrib_pdf(d, 1e-40, list(mu = 0.4, phi = 50), log = TRUE)
S7::method(distrib_pdf, Beta1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbeta(
    x = y,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    log = log
  )
}

#' @title Beta Cumulative Distribution Function in Mean and Precision
#' @name distrib_cdf.Beta1Distrib
#' @description
#' Computes the beta distribution function, the regularized incomplete beta
#' function
#' \deqn{F(q; \mu, \phi) = I_q(\alpha, \beta), \qquad \alpha = \mu\phi, \quad
#'       \beta = (1-\mu)\phi,}
#' by calling [stats::pbeta()] at that pair of shapes. Both tails are available
#' exactly: `lower.tail = FALSE` evaluates \eqn{1 - F} without forming the
#' difference, and `log.p = TRUE` returns a logarithm that stays finite where
#' the probability itself underflows to zero.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param q A numeric vector of quantiles. A value at or below 0 gives a
#'   lower-tail probability of 0 and one at or above 1 gives 1.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(phi))`. With `log.p = TRUE` the values
#'   are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Beta1Distrib()] for the inverse,
#'   [distrib_pdf.Beta1Distrib()] for the density, [distrib_grad_cdf()] for the
#'   derivatives of this function in the parameters, which the beta takes by
#'   finite difference because the derivative of an incomplete beta in its
#'   shapes is hypergeometric, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' th <- list(mu = 0.4, phi = 5)
#'
#' # The method is stats::pbeta at the implied shapes.
#' all.equal(distrib_cdf(d, c(0.2, 0.5, 0.8), th),
#'           pbeta(c(0.2, 0.5, 0.8), 2, 3))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 0.5, th) + distrib_cdf(d, 0.5, th, lower.tail = FALSE)
#'
#' # At mu = 1/2 and phi = 2 the beta is the uniform, so F(q) = q.
#' distrib_cdf(d, c(0.2, 0.5, 0.8), list(mu = 0.5, phi = 2))
#'
#' # Near the lower boundary the probability underflows and its log does not.
#' distrib_cdf(d, 1e-30, list(mu = 0.4, phi = 50))
#' distrib_cdf(d, 1e-30, list(mu = 0.4, phi = 50), log.p = TRUE)
S7::method(distrib_cdf, Beta1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbeta(
    q = q,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Beta Quantile Function in Mean and Precision
#' @name distrib_quantile.Beta1Distrib
#' @description
#' Computes the beta quantile function, the inverse of the regularized
#' incomplete beta function in its argument, by calling [stats::qbeta()] at
#' shapes \eqn{\alpha = \mu\phi} and \eqn{\beta = (1-\mu)\phi}. There is no
#' elementary closed form; `qbeta()` inverts the distribution function
#' numerically. The beta distribution function is strictly increasing on
#' \eqn{(0, 1)}, so the round trip through [distrib_cdf.Beta1Distrib()]
#' returns `p`.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; the endpoints give 0 and 1.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles in \eqn{[0, 1]}, of length
#'   `max(length(p), length(mu), length(phi))`.
#'
#' @seealso [distrib_cdf.Beta1Distrib()], which this inverts;
#'   [distrib_rng.Beta1Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' th <- list(mu = 0.4, phi = 5)
#'
#' # A central 95 percent interval, asymmetric about the mean of 0.4.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # At mu = 1/2 and phi = 2 the beta is the uniform, so Q(p) = p.
#' distrib_quantile(d, p, list(mu = 0.5, phi = 2))
S7::method(distrib_quantile, Beta1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbeta(
    p = p,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Beta Random Number Generator in Mean and Precision
#' @name distrib_rng.Beta1Distrib
#' @description
#' Draws `n` independent beta variates by calling [stats::rbeta()] at shapes
#' \eqn{\alpha = \mu\phi} and \eqn{\beta = (1-\mu)\phi}, so the draws come from
#' R's own beta generator and depend on `.Random.seed` in the usual way. The
#' ratio-of-uniforms fallback the base class supplies is bypassed.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. `mu`
#'   must lie strictly in \eqn{(0, 1)} and `phi` must be strictly positive.
#'
#' @return A numeric vector of `n` draws in \eqn{(0, 1)}.
#'
#' @seealso [distrib_quantile.Beta1Distrib()] for the inverse-transform route,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#'
#' # Same generator as stats::rbeta at the implied shapes.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0.4, phi = 5))
#' set.seed(2)
#' identical(a, rbeta(3, 2, 3))
#'
#' # The sample moments recover the parameters: the mean directly, and the
#' # precision as mu(1 - mu)/var - 1.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 0.4, phi = 5))
#' c(mu = mean(z), phi = mean(z) * (1 - mean(z)) / var(z) - 1)
S7::method(distrib_rng, Beta1Distrib) <- function(distrib, n, theta) {
  stats::rbeta(
    n = n,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]]
  )
}

#' @title Beta Score in Mean and Precision
#' @name distrib_gradient.Beta1Distrib
#' @description
#' Computes the first derivatives of the beta log-density with respect to
#' \eqn{\mu} and \eqn{\phi}, one value per observation, in closed form. With
#' \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi} and \eqn{\psi} the digamma
#' function,
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \phi\left\{
#'         \log\dfrac{y}{1-y} - \psi(\alpha) + \psi(\beta)\right\},}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\alpha)
#'       - (1-\mu)\psi(\beta) + \mu\log y + (1-\mu)\log(1-y).}
#' The data enter only through \eqn{\log y} and \eqn{\log(1-y)}, the beta being
#' an exponential family in the shapes, and the mean component sees them only
#' through the log-odds \eqn{\log\{y/(1-y)\}}.
#'
#' The two expectations that make the score have mean zero are
#' \eqn{\mathbb{E}[\log Y] = \psi(\alpha) - \psi(\phi)} and
#' \eqn{\mathbb{E}[\log(1-Y)] = \psi(\beta) - \psi(\phi)}.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. An endpoint makes
#'   a logarithm infinite and the score non-finite.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
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
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu \in (0,1)} the
#' mean and \eqn{\phi > 0} the precision. \eqn{\psi} is the digamma function,
#' \eqn{\psi = (\log\Gamma)'}.
#'
#' @seealso [distrib_hessian.Beta1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.Beta1Distrib()] for their expectation,
#'   [distrib_grad_y.Beta1Distrib()] for the derivative in the response, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out with the digamma function.
#' a <- 0.4 * 5
#' b <- 0.6 * 5
#' all.equal(g$mu, 5 * (log(y / (1 - y)) - digamma(a) + digamma(b)))
#' all.equal(g$phi, digamma(5) - 0.4 * digamma(a) - 0.6 * digamma(b) +
#'                  0.4 * log(y) + 0.6 * log(1 - y))
#'
#' # The mean component vanishes where the log-odds equal psi(a) - psi(b).
#' distrib_gradient(d, plogis(digamma(a) - digamma(b)), th)$mu
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' set.seed(6)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
S7::method(distrib_gradient, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  beta_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Observed Hessian in Mean and Precision
#' @name distrib_hessian.Beta1Distrib
#' @description
#' Computes the three distinct second derivatives of the beta log-density with
#' respect to \eqn{\mu} and \eqn{\phi}, one value per observation, in closed
#' form. With \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi} and
#' \eqn{\psi_1} the trigamma function,
#' \deqn{\ell^{(\mu\mu)} = -\phi^2\left\{\psi_1(\alpha) + \psi_1(\beta)\right\},
#'       \qquad
#'       \ell^{(\phi\phi)} = \psi_1(\phi) - \mu^2\psi_1(\alpha)
#'         - (1-\mu)^2\psi_1(\beta),}
#' \deqn{\ell^{(\mu\phi)} = \log\dfrac{y}{1-y} - \psi(\alpha) + \psi(\beta)
#'       - \phi\left\{\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right\}.}
#'
#' Two of the three are free of the data. The mixed entry is not: it carries
#' the log-odds residual \eqn{\log\{y/(1-y)\} - \psi(\alpha) + \psi(\beta)},
#' which is the mean component of the score divided by \eqn{\phi} and has
#' expectation zero. So this Hessian and
#' [distrib_expected_hessian.Beta1Distrib()] agree in `mu_mu` and `phi_phi` and
#' differ in `mu_phi` by exactly that residual.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. Only the mixed
#'   entry reads it.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
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
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives. \eqn{\psi}
#' and \eqn{\psi_1} are the digamma and trigamma functions.
#'
#' @seealso [distrib_gradient.Beta1Distrib()] for the score,
#'   [distrib_expected_hessian.Beta1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.Beta1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#' h <- distrib_hessian(d, y, th)
#'
#' # Two entries are constant across the observations and one is not.
#' h
#'
#' # The mixed entry differs from its expectation by the log-odds residual.
#' a <- 0.4 * 5
#' b <- 0.6 * 5
#' all.equal(h$mu_phi - distrib_expected_hessian(d, y, th)$mu_phi,
#'           log(y / (1 - y)) - digamma(a) + digamma(b))
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 0.4, phi = 5 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.4, phi = 5 - eps))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_phi, tolerance = 1e-5)
S7::method(distrib_hessian, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  beta_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Expected Hessian in Mean and Precision
#' @name distrib_expected_hessian.Beta1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. With \eqn{\alpha = \mu\phi},
#' \eqn{\beta = (1-\mu)\phi} and \eqn{\psi_1} the trigamma function,
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] =
#'         -\phi^2\left\{\psi_1(\alpha) + \psi_1(\beta)\right\}, \qquad
#'       \mathbb{E}\left[\ell^{(\phi\phi)}\right] = \psi_1(\phi)
#'         - \mu^2\psi_1(\alpha) - (1-\mu)^2\psi_1(\beta),}
#' \deqn{\mathbb{E}\left[\ell^{(\mu\phi)}\right] =
#'         -\phi\left\{\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right\}.}
#'
#' The first two are the observed values themselves, being free of the data.
#' The third drops the log-odds residual the observed mixed entry carries,
#' whose expectation is zero because
#' \eqn{\mathbb{E}[\log\{Y/(1-Y)\}] = \psi(\alpha) - \psi(\beta)}.
#'
#' The mixed entry does not vanish, so the mean and the precision are not
#' orthogonal; their maximum likelihood estimates are asymptotically
#' correlated. It does vanish at \eqn{\mu = 1/2}, where the two shapes are
#' equal and the density is symmetric.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
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
#' expectation of the **observed information** under the model. The beta is a
#' regular family, so the second Bartlett identity holds and this equals the
#' variance of the score.
#'
#' @seealso [distrib_hessian.Beta1Distrib()] for the observed quantity this is
#'   the expectation of, [fisher_scoring()], which inverts it at each step, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' th <- list(mu = 0.4, phi = 5)
#' eh <- distrib_expected_hessian(d, c(0.2, 0.5, 0.8), th)
#' lapply(eh, unique)
#'
#' # Written out with the trigamma function.
#' a <- 0.4 * 5
#' b <- 0.6 * 5
#' c(-5^2 * (trigamma(a) + trigamma(b)),
#'   trigamma(5) - 0.4^2 * trigamma(a) - 0.6^2 * trigamma(b),
#'   -5 * (0.4 * trigamma(a) - 0.6 * trigamma(b)))
#'
#' # The observed mixed entry averages onto it over a large sample; the other
#' # two are equal to it observation by observation.
#' set.seed(9)
#' z <- distrib_rng(d, 5e5, th)
#' c(observed = mean(distrib_hessian(d, z, th)$mu_phi),
#'   expected = eh$mu_phi[1])
#'
#' # The mixed entry vanishes at mu = 1/2, where the density is symmetric.
#' distrib_expected_hessian(d, 0.5, list(mu = 0.5, phi = 5))$mu_phi
S7::method(distrib_expected_hessian, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  beta_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Third-Order Derivatives in Mean and Precision
#' @name distrib_deriv3.Beta1Distrib
#' @description
#' Computes the four distinct third derivatives of the beta log-density with
#' respect to \eqn{\mu} and \eqn{\phi}, in closed form as combinations of
#' \eqn{\psi_2}, the second derivative of the digamma function, at
#' \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi} and \eqn{\phi}.
#'
#' **Every component of this order is free of the response.** The data enter
#' the log-density only through \eqn{(\alpha-1)\log y + (\beta-1)\log(1-y)},
#' which is linear in the shapes and so at most quadratic in \eqn{(\mu, \phi)}
#' through the bilinear map \eqn{\alpha = \mu\phi}; three derivatives kill it.
#' The observed and expected values therefore coincide exactly, and `expected`
#' selects nothing: the same kernel runs either way and the two results are
#' identical to the bit. `approx` and `nsim` are ignored for the same reason.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected third derivatives being the same numbers. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_phi`,
#'   `mu_phi_phi` and `phi_phi_phi`, each of length `length(y)` and constant
#'   within itself when the parameters are. The names enumerate the distinct
#'   multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\psi} is the digamma function and
#' \eqn{\psi_m} its \eqn{m}th derivative.
#'
#' @seealso [distrib_hessian.Beta1Distrib()] for the order below, whose mixed
#'   entry is the last one to carry the response, and
#'   [distrib_deriv4.Beta1Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # Four constants: nothing at this order depends on the observation.
#' lapply(d3, unique)
#'
#' # So asking for the expectation changes nothing.
#' identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 0.4 + eps, phi = 5))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, phi = 5))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
S7::method(distrib_deriv3, Beta1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  beta_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Fourth-Order Derivatives in Mean and Precision
#' @name distrib_deriv4.Beta1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the beta log-density with
#' respect to \eqn{\mu} and \eqn{\phi}, in closed form as combinations of
#' \eqn{\psi_3}, the third derivative of the digamma function, at
#' \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi} and \eqn{\phi}.
#'
#' As at third order, every component is free of the response, so the observed
#' and expected values coincide exactly and `expected`, `approx` and `nsim` are
#' all without effect.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected fourth derivatives being the same numbers. Defaults to
#'   `FALSE`.
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
#'   `mu_mu_mu_phi`, `mu_mu_phi_phi`, `mu_phi_phi_phi` and `phi_phi_phi_phi`,
#'   each of length `length(y)` and constant within itself when the parameters
#'   are.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\psi} is the digamma
#' function and \eqn{\psi_m} its \eqn{m}th derivative.
#'
#' @seealso [distrib_deriv3.Beta1Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.Beta1Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#' d4 <- distrib_deriv4(d, y, th)
#'
#' # Five constants: nothing at this order depends on the observation.
#' lapply(d4, unique)
#'
#' # So asking for the expectation changes nothing.
#' identical(d4, distrib_deriv4(d, y, th, expected = TRUE))
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, phi = 5))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, phi = 5))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
S7::method(distrib_deriv4, Beta1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  beta_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta First Derivative in the Response, Mean and Precision
#' @name distrib_grad_y.Beta1Distrib
#' @description
#' Computes the first derivative of the beta log-density with respect to the
#' response, in closed form at the implied shapes \eqn{\alpha = \mu\phi} and
#' \eqn{\beta = (1-\mu)\phi}:
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\alpha - 1}{y}
#'       - \dfrac{\beta - 1}{1 - y}.}
#' Where both shapes exceed one it changes sign at the mode
#' \eqn{(\alpha-1)/(\alpha+\beta-2)}, so it is positive below the mode and
#' negative above it. Where a shape falls below one the density is unbounded at
#' the corresponding endpoint and the derivative does not change sign there.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. An endpoint makes
#'   the value infinite unless the corresponding shape is exactly 1.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(phi))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Beta1Distrib()] for the second derivative in the
#'   response, [distrib_gradient.Beta1Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#'
#' # Written out at the implied shapes.
#' a <- 0.4 * 5
#' b <- 0.6 * 5
#' all.equal(distrib_grad_y(d, y, th), (a - 1) / y - (b - 1) / (1 - y))
#'
#' # Zero at the mode (a - 1)/(a + b - 2), positive below and negative above.
#' mode <- (a - 1) / (a + b - 2)
#' c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#'
#' # At mu = 1/2 and phi = 2 both shapes are 1 and the density is flat.
#' distrib_grad_y(d, y, list(mu = 0.5, phi = 2))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, Beta1Distrib) <- function(distrib, y, theta) {
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  (a - 1) / y - (b - 1) / (1 - y)
}

#' @title Beta Second Derivative in the Response, Mean and Precision
#' @name distrib_hess_y.Beta1Distrib
#' @description
#' Computes the second derivative of the beta log-density with respect to the
#' response, in closed form at the implied shapes \eqn{\alpha = \mu\phi} and
#' \eqn{\beta = (1-\mu)\phi}:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{\alpha - 1}{y^2}
#'       - \dfrac{\beta - 1}{(1 - y)^2}.}
#' It is negative throughout, so the log-density is concave in the response,
#' whenever both shapes are at least one; a shape below one makes the
#' corresponding term positive and can turn the curvature positive near that
#' endpoint. At \eqn{\mu = 1/2} and \eqn{\phi = 2} both shapes are 1 and the
#' curvature is exactly zero, the density being the uniform.
#'
#' @param distrib A `Beta1Distrib` object, from [beta1_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. An endpoint makes
#'   the value infinite unless the corresponding shape is exactly 1.
#' @param theta A named list with components `mu` and `phi`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie strictly in \eqn{(0, 1)} and `phi` must be
#'   strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(phi))`, one value per observation.
#'
#' @seealso [distrib_grad_y.Beta1Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Beta1Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- beta1_distrib()
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#'
#' a <- 0.4 * 5
#' b <- 0.6 * 5
#' all.equal(distrib_hess_y(d, y, th),
#'           -(a - 1) / y^2 - (b - 1) / (1 - y)^2)
#'
#' # Concave everywhere while both shapes exceed one.
#' all(distrib_hess_y(d, y, th) < 0)
#'
#' # Exactly flat at mu = 1/2, phi = 2, where the beta is the uniform.
#' distrib_hess_y(d, y, list(mu = 0.5, phi = 2))
S7::method(distrib_hess_y, Beta1Distrib) <- function(distrib, y, theta) {
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  -(a - 1) / y^2 - (b - 1) / (1 - y)^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Beta Distribution, Mean and Precision
#'
#' @description
#' Builds the distribution object for the beta family on \eqn{(0, 1)}
#' parametrized by its mean \eqn{\mu \in (0, 1)} and a precision \eqn{\phi > 0},
#' so that the shapes are \eqn{\alpha = \mu\phi} and \eqn{\beta = (1-\mu)\phi}.
#' This is the parametrization beta regression uses. The returned object carries
#' closed-form derivatives of the log-density to fourth order, in the parameters
#' and in the response, and closed-form moments, so every generic of the toolkit
#' answers without a numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. The defaults are the logit for
#' the mean, which lives in \eqn{(0, 1)}, and the logarithm for the precision,
#' which is positive.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::logit_link()], which maps
#'   \eqn{(0, 1)} onto the line and so keeps every fitted mean inside the unit
#'   interval.
#' @param link_phi A `link` object from `linkfunctions7` for the precision
#'   \eqn{\phi}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, 1)} is
#' \deqn{f(y; \mu, \phi) = \dfrac{\Gamma(\phi)}{\Gamma(\alpha)\Gamma(\beta)}\,
#'       y^{\alpha-1}(1-y)^{\beta-1}, \qquad \alpha = \mu\phi, \quad
#'       \beta = (1-\mu)\phi,}
#' the distribution function \eqn{F(q) = I_q(\alpha, \beta)} with \eqn{I} the
#' regularized incomplete beta function, and the quantile function its
#' numerical inverse. The mean is \eqn{\mu}, the variance
#' \eqn{\mu(1-\mu)/(\phi+1)} and the skewness
#' \eqn{2(1-2\mu)\sqrt{\phi+1}/\{(\phi+2)\sqrt{\mu(1-\mu)}\}}.
#'
#' The name *precision* is earned: at a fixed mean the variance falls as
#' \eqn{1/(\phi+1)}, so \eqn{\phi} says how tightly the mass concentrates. Two
#' settings are worth recognizing. At \eqn{\mu = 1/2} and \eqn{\phi = 2} both
#' shapes are 1 and the density is the uniform. Where a shape falls below one,
#' which happens when \eqn{\phi} is small, the density is unbounded at the
#' corresponding endpoint.
#'
#' This is the same law as [beta2_distrib()], which carries the mean and the
#' *variance*, the two being related by \eqn{\sigma^2 = \mu(1-\mu)/(\phi+1)}.
#' They are separate families because the second parameter is a different
#' quantity in each, with its own interpretation, standard error and interval.
#'
#' # Derivatives
#'
#' With \eqn{\psi} the digamma function and \eqn{\psi_1} the trigamma, the
#' score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \phi\left\{
#'         \log\dfrac{y}{1-y} - \psi(\alpha) + \psi(\beta)\right\}, \qquad
#'       \dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\alpha)
#'         - (1-\mu)\psi(\beta) + \mu\log y + (1-\mu)\log(1-y),}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] =
#'         -\phi^2\left\{\psi_1(\alpha) + \psi_1(\beta)\right\}, \qquad
#'       \mathbb{E}\left[\ell^{(\phi\phi)}\right] = \psi_1(\phi)
#'         - \mu^2\psi_1(\alpha) - (1-\mu)^2\psi_1(\beta),}
#' \deqn{\mathbb{E}\left[\ell^{(\mu\phi)}\right] =
#'         -\phi\left\{\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right\}.}
#' The mixed entry is zero only at \eqn{\mu = 1/2}, so the mean and the
#' precision are in general not orthogonal.
#'
#' # Where the response stops entering
#'
#' The data reach the log-density only through
#' \eqn{(\alpha-1)\log y + (\beta-1)\log(1-y)}, which is linear in the shapes
#' and so at most quadratic in \eqn{(\mu, \phi)} through the bilinear map
#' \eqn{\alpha = \mu\phi}. Two consequences follow, and both are visible in the
#' methods:
#'
#' - the observed Hessian equals its expectation in `mu_mu` and `phi_phi` and
#'   differs in `mu_phi` by the log-odds residual
#'   \eqn{\log\{y/(1-y)\} - \psi(\alpha) + \psi(\beta)}, whose expectation is
#'   zero;
#' - **every** third and fourth derivative is free of the response, so
#'   [distrib_deriv3.Beta1Distrib()] and [distrib_deriv4.Beta1Distrib()] return
#'   the same numbers whether or not `expected = TRUE` is asked for.
#'
#' The derivatives of the *distribution* function in the parameters have no
#' elementary form, the derivative of an incomplete beta in its shapes being
#' hypergeometric, and are taken by finite difference on the analytic cdf.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' estimate is closed form: both estimating equations involve digamma
#' functions of the shapes. The method of moments supplies the starting values
#' \eqn{\hat\mu = \bar y} and \eqn{\hat\phi = \bar y(1-\bar y)/s^2 - 1}, and
#' the example below shows them landing beside the estimates.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu \in (0,1)} the
#' mean and \eqn{\phi > 0} the precision. \eqn{\alpha} and \eqn{\beta} are the
#' implied shapes. \eqn{\psi} is the digamma function and \eqn{\psi_m} its
#' \eqn{m}th derivative. \eqn{\eta} is a parameter on the unconstrained scale
#' of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Beta1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"beta1"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, 1)`, `params` `c("mu", "phi")`, `n_params`
#'   `2`, `params_bounds` the list of \eqn{(0, 1)} and \eqn{(0, \infty)}, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [beta2_distrib()] for the same law in the mean and the variance;
#' [betabinom1_distrib()] for the beta as a mixing law over a binomial
#' probability; [bernoulli_distrib()] for a response at the endpoints rather
#' than between them; [zero_adjusted()] and [truncated()] for the wrappers a
#' response with mass at 0 or 1 needs; [fit_distrib()] to estimate the
#' parameters; [check_distrib()] to validate a family of your own against the
#' same battery this one passes; [Beta1Distrib] for the class.
#'
#' @references
#' Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
#' rates and proportions. *Journal of Applied Statistics* **31**, 799-815.
#'
#' @importFrom linkfunctions7 logit_link log_link
#' @importFrom stats dbeta pbeta qbeta rbeta
#'
#' @examples
#' d <- beta1_distrib()
#' d
#'
#' # The density is stats::dbeta at shape1 = mu phi, shape2 = (1 - mu) phi.
#' y <- c(0.2, 0.5, 0.8)
#' th <- list(mu = 0.4, phi = 5)
#' all.equal(distrib_pdf(d, y, th), dbeta(y, 2, 3))
#'
#' # Moments: variance mu(1 - mu)/(phi + 1), so phi is a precision.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' 0.4 * 0.6 / (5 + 1)
#'
#' # Both shapes are 1 at mu = 1/2, phi = 2: the beta is the uniform.
#' distrib_pdf(d, y, list(mu = 0.5, phi = 2))
#'
#' # Fitting recovers the parameters; the moment estimates start it off.
#' set.seed(6)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted  = coef(fit),
#'       moments = c(mu = mean(z),
#'                   phi = mean(z) * (1 - mean(z)) / var(z) - 1))
#'
#' # Every third derivative is free of the response, so it is a constant.
#' lapply(distrib_deriv3(d, y, th), unique)
#'
#' @export
beta1_distrib <- function(link_mu = logit_link(), link_phi = log_link()) {
  
  Beta1Distrib(
    distrib_name = "beta1",
    dimension = "univariate",
    bounds = c(0, 1),
    
    params = c("mu", "phi"),
    params_interpretation = c(mu = "mean", phi = "precision"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(0, 1),
      phi = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      phi = link_phi
    )
  )
  
}
