#' @include distrib.R generics.R
NULL

#' @title Laplace Distribution Class, Location and Rate
#' @name Laplace2Distrib
#'
#' @description
#' The S7 class of the Laplace (double exponential) family written by its
#' **rate** \eqn{\lambda > 0} in place of its scale, with density
#' \eqn{f(y) = (\lambda/2)\exp(-\lambda|y-\mu|)} on the whole real line. It is
#' the same law as [laplace_distrib()] at \eqn{\lambda = 1/\sigma}; the
#' parametrization differs and the derivatives with it. It inherits from
#' `continuous_distrib`, so it answers every generic of the `distrib` contract.
#'
#' The rate form is the one a lasso penalty is written in: with \eqn{\mu = 0}
#' the negative log-density is \eqn{\lambda|y| - \log(\lambda/2)}, so
#' \eqn{\lambda} is the penalty's own tuning parameter and larger values shrink
#' harder. `penalties7::lasso_penalty()` is this family with the location held
#' at zero.
#'
#' Like [laplace_distrib()], this family has a **kink** at \eqn{y = \mu} and is
#' non-regular in its location; `params_smooth` is `c(mu = FALSE,
#' lambda = TRUE)`.
#'
#' Build one with [laplace2_distrib()]. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Laplace2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [laplace2_distrib()] they hold `"laplace2"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "lambda")`, the interpretations
#'   `c(mu = "location", lambda = "rate")`, `2`, the domains
#'   \eqn{(-\infty, \infty)} and \eqn{(0, \infty)}, the two links, and
#'   `c(mu = FALSE, lambda = TRUE)` for `params_smooth`.
#'
#' @seealso [laplace2_distrib()] to build one;
#'   [laplace_distrib()] for the same law in its scale;
#'   [enet_distrib()], which mixes this family with a Gaussian.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.Laplace2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Laplace2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Laplace2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Laplace2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Laplace2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Laplace2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Laplace2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Laplace2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Laplace2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Laplace2Distrib],
#'   [`distrib_rng()`][distrib_rng.Laplace2Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#' d@params
#' d@params_interpretation
#' d@params_smooth
#'
#' # The same law as laplace_distrib() at lambda = 1/sigma.
#' y <- c(-1.2, 0.3, 2.5)
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 1 / 1.5)),
#'           distrib_pdf(laplace_distrib(), y, list(mu = 0.4, sigma = 1.5)))
#'
#' # The variance is 2/lambda^2, so a larger rate is a tighter distribution.
#' vapply(c(0.5, 1, 2), function(l) variance(d, list(mu = 0, lambda = l)),
#'        numeric(1))
Laplace2Distrib <- S7::new_class("Laplace2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Laplace Probability Density Function, Rate Parametrization
#' @name distrib_pdf.Laplace2Distrib
#' @description
#' Computes the Laplace density in the rate parametrization,
#' \deqn{f(y; \mu, \lambda) = \dfrac{\lambda}{2} \exp\left(-\lambda|y - \mu|\right),}
#' from the log-density \eqn{\log(\lambda/2) - \lambda|y-\mu|}, which is formed
#' first and exponentiated only when `log = FALSE`. It equals the density of
#' [laplace_distrib()] at \eqn{\sigma = 1/\lambda}.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive; the arithmetic is performed
#'   as written, so a non-positive value gives `NaN` without a warning of its
#'   own.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(lambda))`, one value per observation.
#'
#' @seealso [distrib_pdf.LaplaceDistrib()] for the scale parametrization,
#'   [distrib_cdf.Laplace2Distrib()] for the distribution function, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#'
#' # The density, written out.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 2)),
#'           exp(-2 * abs(y - 0.4)) * 2 / 2)
#'
#' # The same law as the scale parametrization at lambda = 1/sigma.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 1 / 1.5)),
#'           distrib_pdf(laplace_distrib(), y, list(mu = 0.4, sigma = 1.5)))
#'
#' # A larger rate concentrates the mass: the peak is lambda/2.
#' vapply(c(0.5, 1, 2), function(l) distrib_pdf(d, 0.4,
#'                                              list(mu = 0.4, lambda = l)),
#'        numeric(1))
S7::method(distrib_pdf, Laplace2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  log_d <- base::log(lam / 2) - lam * abs(y - mu)
  if (log) log_d else exp(log_d)
}

#' @title Laplace Cumulative Distribution Function, Rate Parametrization
#' @name distrib_cdf.Laplace2Distrib
#' @description
#' Computes the Laplace distribution function in the rate parametrization,
#' one exponential on each side of the location:
#' \deqn{F(q; \mu, \lambda) = \begin{cases}
#'   \tfrac{1}{2}\exp(\lambda(q-\mu)), & q < \mu,\\[4pt]
#'   1 - \tfrac{1}{2}\exp(-\lambda(q-\mu)), & q \ge \mu.
#' \end{cases}}
#' Both arms are continuous at \eqn{q = \mu}, where the value is \eqn{1/2}.
#' With `lower.tail = FALSE` the complement is taken after the branch, and with
#' `log.p = TRUE` the logarithm is taken afterwards; neither avoids
#' cancellation, so a probability that has already underflowed is not
#' recovered.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   computed as \eqn{1 - F}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(lambda))`.
#'
#' @seealso [distrib_cdf.LaplaceDistrib()] for the scale parametrization,
#'   [distrib_quantile.Laplace2Distrib()] for the inverse, and [distrib_cdf()]
#'   for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' th <- list(mu = 0.4, lambda = 2)
#'
#' # One half at the location.
#' distrib_cdf(d, 0.4, th)
#'
#' # Each arm is an exponential, written out.
#' q <- c(-2, 3)
#' all.equal(distrib_cdf(d, q, th),
#'           ifelse(q < 0.4, 0.5 * exp(2 * (q - 0.4)),
#'                  1 - 0.5 * exp(-2 * (q - 0.4))))
#'
#' # The same probabilities as the scale parametrization at sigma = 1/lambda.
#' all.equal(distrib_cdf(d, q, list(mu = 0.4, lambda = 1 / 1.5)),
#'           distrib_cdf(laplace_distrib(), q, list(mu = 0.4, sigma = 1.5)))
S7::method(distrib_cdf, Laplace2Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  res <- ifelse(q < mu, 0.5 * exp(lam * (q - mu)), 1 - 0.5 * exp(-lam * (q - mu)))
  if (!lower.tail) res <- 1 - res
  if (log.p) base::log(res) else res
}

#' @title Laplace Quantile Function, Rate Parametrization
#' @name distrib_quantile.Laplace2Distrib
#' @description
#' Computes the Laplace quantile function in the rate parametrization,
#' \deqn{Q(p; \mu, \lambda) = \mu - \dfrac{\mathrm{sign}(p - \tfrac{1}{2})}{\lambda}\,\log\left(1 - 2\left|p - \tfrac{1}{2}\right|\right),}
#' one expression covering both arms. The median is \eqn{\mu} and the quartiles
#' are \eqn{\mu \pm \log(2)/\lambda}.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. The endpoints give `-Inf` and `Inf`.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)} and `p` is replaced by
#'   `1 - p` before the formula is applied.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are
#'   exponentiated first. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of length
#'   `max(length(p), length(mu), length(lambda))`.
#'
#' @seealso [distrib_cdf.Laplace2Distrib()], which this inverts;
#'   [distrib_rng.Laplace2Distrib()], which draws through it; and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' th <- list(mu = 0.4, lambda = 2)
#'
#' # The median is mu and the quartiles are mu -/+ log(2)/lambda.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#' 0.4 + c(-1, 0, 1) * log(2) / 2
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
S7::method(distrib_quantile, Laplace2Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  mu - sign(p - 0.5) * base::log(1 - 2 * abs(p - 0.5)) / lam
}

#' @title Laplace Random Number Generator, Rate Parametrization
#' @name distrib_rng.Laplace2Distrib
#' @description
#' Draws `n` independent Laplace variates by **inverse transform**: `n` uniform
#' variates from [stats::runif()] are passed through
#' [distrib_quantile.Laplace2Distrib()]. The quantile function is elementary
#' and exact, so this is cheaper and more accurate than the generalized
#' ratio-of-uniforms fallback the base class supplies.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param n A single positive integer, the number of draws. One uniform variate
#'   is consumed per draw.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled.
#'   `lambda` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.Laplace2Distrib()], through which the draws are
#'   made; [distrib_rng.LaplaceDistrib()] for the scale parametrization; and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#'
#' # Inverse transform, so the draws are the quantiles of the uniforms drawn.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 0.4, lambda = 2))
#' set.seed(2)
#' identical(a, distrib_quantile(d, runif(3), list(mu = 0.4, lambda = 2)))
#'
#' # The sample variance recovers 2/lambda^2.
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 3, lambda = 0.5))
#' c(mean = mean(z), var = var(z), two_over_lambda_sq = 2 / 0.5^2)
S7::method(distrib_rng, Laplace2Distrib) <- function(distrib, n, theta, ...) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Laplace Score, Rate Parametrization
#' @name distrib_gradient.Laplace2Distrib
#' @description
#' Computes the first derivatives of the Laplace log-density with respect to
#' \eqn{\mu} and \eqn{\lambda}, one value per observation, in closed form.
#' Writing \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \lambda\,\mathrm{sign}(r),
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{\lambda} - |r|.}
#'
#' Both are simpler than their scale-parametrization counterparts: the rate is
#' the natural parameter of the exponential family in \eqn{|r|}, so its score
#' is the difference between \eqn{1/\lambda} and the sufficient statistic. The
#' score in \eqn{\mu} again carries only the sign of the residual, so the
#' estimate of \eqn{\mu} is the sample median.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `lambda`, each of
#'   length `max(length(y), length(mu), length(lambda))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\lambda > 0} the rate, with variance \eqn{2/\lambda^2}. \eqn{r = y-\mu}
#' is the residual. Here \eqn{\lambda} is a rate; the same letter names a
#' penalty parameter above, and the two meet in the lasso, which is this
#' family with \eqn{\mu} held at zero.
#'
#' @seealso [distrib_gradient.LaplaceDistrib()] for the scale parametrization,
#'   [distrib_expected_hessian.Laplace2Distrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, lambda = 2)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' r <- y - 0.4
#' all.equal(g$mu, 2 * sign(r))
#' all.equal(g$lambda, 1 / 2 - abs(r))
#'
#' # The rate's score vanishes summed at 1/mean|r|, its estimate.
#' set.seed(12)
#' z <- distrib_rng(d, 2000, list(mu = 3, lambda = 0.5))
#' lam_hat <- 1 / mean(abs(z - median(z)))
#' sum(distrib_gradient(d, z, list(mu = median(z), lambda = lam_hat))$lambda)
S7::method(distrib_gradient, Laplace2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  r <- y - mu
  list(
    mu = lam * sign(r),
    lambda = 1 / lam - abs(r)
  )
}

#' @title Laplace Observed Hessian, Rate Parametrization
#' @name distrib_hessian.Laplace2Distrib
#' @description
#' Computes the three distinct second derivatives of the Laplace log-density
#' with respect to \eqn{\mu} and \eqn{\lambda}, one value per observation, in
#' closed form. Writing \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = 0,
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \lambda^2} = -\dfrac{1}{\lambda^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \lambda} = \mathrm{sign}(r).}
#'
#' The first entry is exactly zero for every observation, the log-density being
#' piecewise linear in \eqn{\mu}; see
#' [distrib_hessian.LaplaceDistrib()] for what that costs. The curvature in
#' \eqn{\lambda} carries no data at all, so in this parametrization the
#' observed and expected values agree in that entry.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `lambda_lambda` and
#'   `mu_lambda`, each of length `length(y)`. `mu_mu` is a vector of zeros and
#'   `lambda_lambda` is constant at \eqn{-1/\lambda^2}.
#'
#' @seealso [distrib_expected_hessian.Laplace2Distrib()] for the information,
#'   whose location entry differs from this one;
#'   [distrib_hessian.LaplaceDistrib()] for the scale parametrization; and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, lambda = 2)
#' h <- distrib_hessian(d, y, th)
#'
#' # Zero in the location, constant in the rate, a sign in the mixed entry.
#' h$mu_mu
#' h$lambda_lambda
#' all.equal(h$mu_lambda, sign(y - 0.4))
#'
#' # The rate entry already equals its own expectation; the location one
#' # does not, and the expected page says why.
#' rbind(observed = vapply(h, function(v) v[1], numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, y, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_hessian, Laplace2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  lam <- theta[[2]]
  n <- length(y)
  list(
    mu_mu = rep(0, n),
    lambda_lambda = rep(-1 / lam^2, length.out = n),
    mu_lambda = sign(y - mu) + rep(0, n)
  )
}

#' @title Laplace Expected Hessian, Rate Parametrization
#' @name distrib_expected_hessian.Laplace2Distrib
#' @description
#' Returns the negative of the Fisher information, in closed form and with no
#' quadrature or simulation:
#' \deqn{-I(\mu) = -\lambda^2, \qquad
#'       -I(\lambda) = -\dfrac{1}{\lambda^2}, \qquad
#'       -I(\mu, \lambda) = 0.}
#'
#' For \eqn{\mu} this is the **variance of the score**. It is not the
#' expectation of [distrib_hessian.Laplace2Distrib()], which is identically
#' zero: the
#' family has a kink at \eqn{y = \mu} and the second Bartlett identity fails
#' there, exactly as in [distrib_expected_hessian.LaplaceDistrib()]. With
#' \eqn{\partial\ell/\partial\mu = \lambda\,\mathrm{sign}(r)} the variance is
#' \eqn{\lambda^2}. The mixed entry vanishes because
#' \eqn{\mathbb{E}[\mathrm{sign}(r)] = 0} by symmetry, and the rate entry does
#' equal its observed counterpart, that one carrying no data.
#'
#' Because the values do not depend on the data, `approx` and `nsim` are
#' ignored. `y` is read only for its length.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available. Accepted so that the
#'   signature matches the generic's. As for [laplace_distrib()], a strategy
#'   averaging the observed second derivative would answer 0 in the location.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `lambda_lambda` and
#'   `mu_lambda`, each of length `length(y)`.
#'
#' @section Notation:
#' The **expected information** is the variance of the score. For a regular
#' family it also equals the expectation of the **observed information**, by
#' the second Bartlett identity; this family is not regular in \eqn{\mu} and
#' the two differ there.
#'
#' @seealso [distrib_hessian.Laplace2Distrib()], whose `mu_mu` is zero;
#'   [distrib_expected_hessian.LaplaceDistrib()] for the same argument in the
#'   scale parametrization; [fisher_scoring()], which inverts this matrix; and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' th <- list(mu = 0.4, lambda = 2)
#'
#' # -lambda^2 in the location, -1/lambda^2 in the rate, 0 mixed.
#' lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#' c(-2^2, -1 / 2^2)
#'
#' # The location entry is the variance of the score.
#' set.seed(12)
#' z <- distrib_rng(d, 1e5, th)
#' c(var_of_score = mean(distrib_gradient(d, z, th)$mu^2),
#'   information = -distrib_expected_hessian(d, 0, th)$mu_mu)
#'
#' # It agrees with the scale parametrization once the chain rule is applied:
#' # I(mu) is the same number, the location being shared.
#' -distrib_expected_hessian(laplace_distrib(), 0,
#'                           list(mu = 0.4, sigma = 1 / 2))$mu_mu
S7::method(distrib_expected_hessian, Laplace2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  lam <- theta[[2]]
  n <- length(y)
  list(
    mu_mu = rep(-lam^2, length.out = n),
    lambda_lambda = rep(-1 / lam^2, length.out = n),
    mu_lambda = rep(0, n)
  )
}

#' @title Laplace Third-Order Derivatives, Rate Parametrization
#' @name distrib_deriv3.Laplace2Distrib
#' @description
#' Computes the four distinct third derivatives of the Laplace log-density with
#' respect to \eqn{\mu} and \eqn{\lambda}. In this parametrization the
#' log-density is \eqn{\log(\lambda/2) - \lambda|y-\mu|}, which is **linear in
#' \eqn{\lambda}** apart from \eqn{\log\lambda}, and piecewise linear in
#' \eqn{\mu}. Every third derivative involving \eqn{\mu} at all is therefore
#' zero, and the only surviving component comes from \eqn{\log\lambda}:
#' \deqn{\dfrac{\partial^3 \ell}{\partial \lambda^3} = \dfrac{2}{\lambda^3},}
#' with `mu_mu_mu`, `mu_mu_lambda` and `mu_lambda_lambda` all zero.
#'
#' None of these depends on the data, so the observed and expected values
#' coincide: `expected`, `approx` and `nsim` are accepted and have no effect.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` is not read. `lambda`
#'   must be strictly positive.
#' @param expected Ignored: the observed values do not depend on the data, so
#'   they are already their own expectations. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, for the same reason as `expected`.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_lambda`,
#'   `mu_lambda_lambda` and `lambda_lambda_lambda`, each of length `length(y)`.
#'   The first three are zero and the last is constant at \eqn{2/\lambda^3}.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k)}} is the third derivative of the log-density with respect
#' to parameters \eqn{i}, \eqn{j} and \eqn{k}. Parenthesized superscripts name
#' derivatives.
#'
#' @seealso [distrib_deriv3.LaplaceDistrib()] for the scale parametrization,
#'   where the picture is less sparse; [distrib_deriv4.Laplace2Distrib()] for
#'   the order above; [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, lambda = 2)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # Only the pure-rate component survives.
#' vapply(d3, function(v) v[1], numeric(1))
#' 2 / 2^3
#'
#' # expected = TRUE changes nothing, the values carrying no data.
#' identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
S7::method(distrib_deriv3, Laplace2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  lam <- theta[[2]]
  n <- length(y)
  zero <- rep(0, n)
  list(
    mu_mu_mu = zero,
    mu_mu_lambda = zero,
    mu_lambda_lambda = zero,
    lambda_lambda_lambda = rep(2 / lam^3, length.out = n)
  )
}

#' @title Laplace Fourth-Order Derivatives, Rate Parametrization
#' @name distrib_deriv4.Laplace2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the Laplace log-density
#' with respect to \eqn{\mu} and \eqn{\lambda}. As at third order, every
#' component involving \eqn{\mu} is zero and the only survivor comes from
#' \eqn{\log\lambda}:
#' \deqn{\dfrac{\partial^4 \ell}{\partial \lambda^4} = -\dfrac{6}{\lambda^4}.}
#'
#' None of these depends on the data, so the observed and expected values
#' coincide: `expected`, `approx` and `nsim` are accepted and have no effect.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `lambda`. `mu` is not
#'   read. `lambda` must be strictly positive.
#' @param expected Ignored: the observed values do not depend on the data, so
#'   they are already their own expectations. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, for the same reason as `expected`.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_lambda`, `mu_mu_lambda_lambda`, `mu_lambda_lambda_lambda` and
#'   `lambda_lambda_lambda_lambda`, each of length `length(y)`. The first four
#'   are zero and the last is constant at \eqn{-6/\lambda^4}.
#'
#' @section Notation:
#' \eqn{\ell^{(i j k l)}} is the fourth derivative of the log-density with
#' respect to parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. Parenthesized
#' superscripts name derivatives.
#'
#' @seealso [distrib_deriv3.Laplace2Distrib()] for the order below;
#'   [distrib_deriv4.LaplaceDistrib()] for the scale parametrization;
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, lambda = 2)
#'
#' # Only the pure-rate component survives.
#' vapply(distrib_deriv4(d, y, th), function(v) v[1], numeric(1))
#' -6 / 2^4
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, y, list(mu = 0.4, lambda = 2 + eps))
#' dn <- distrib_deriv3(d, y, list(mu = 0.4, lambda = 2 - eps))
#' all.equal((up$lambda_lambda_lambda - dn$lambda_lambda_lambda) / (2 * eps),
#'           distrib_deriv4(d, y, th)$lambda_lambda_lambda_lambda,
#'           tolerance = 1e-6)
S7::method(distrib_deriv4, Laplace2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  lam <- theta[[2]]
  n <- length(y)
  zero <- rep(0, n)
  list(
    mu_mu_mu_mu = zero,
    mu_mu_mu_lambda = zero,
    mu_mu_lambda_lambda = zero,
    mu_lambda_lambda_lambda = zero,
    lambda_lambda_lambda_lambda = rep(-6 / lam^4, length.out = n)
  )
}

#' @title Laplace First Derivative in the Response, Rate Parametrization
#' @name distrib_grad_y.Laplace2Distrib
#' @description
#' Computes the first derivative of the Laplace log-density with respect to the
#' response,
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\lambda\,\mathrm{sign}(y - \mu),}
#' in closed form. The family is a location family in \eqn{\mu}, so this is the
#' negative of the score in \eqn{\mu}. It takes only the three values
#' \eqn{-\lambda}, 0 and \eqn{\lambda}; at \eqn{y = \mu} exactly, `sign(0)` is
#' 0 and the method returns 0, the midpoint of the subdifferential, the
#' derivative not existing there.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `lambda` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(lambda))`, taking the values
#'   \eqn{-\lambda}, 0 and \eqn{\lambda} only.
#'
#' @seealso [distrib_hess_y.Laplace2Distrib()] for the second derivative, which
#'   is zero; [distrib_grad_y.LaplaceDistrib()] for the scale parametrization;
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, lambda = 2)
#'
#' all.equal(distrib_grad_y(d, y, th), -2 * sign(y - 0.4))
#'
#' # A location family: minus the score in the location.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Three values only, and 0 at the kink itself.
#' distrib_grad_y(d, 0.4 + c(-100, -1e-9, 0, 1e-9, 100), th)
S7::method(distrib_grad_y, Laplace2Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * sign(y - theta[[1]])
}

#' @title Laplace Second Derivative in the Response, Rate Parametrization
#' @name distrib_hess_y.Laplace2Distrib
#' @description
#' Returns zero for every observation. The log-density
#' \eqn{\log(\lambda/2) - \lambda|y-\mu|} is linear in \eqn{y} on each side of
#' the location, so its second derivative in the response vanishes wherever it
#' exists. At \eqn{y = \mu} the first derivative drops from \eqn{\lambda} to
#' \eqn{-\lambda} and the second derivative does not exist; the returned zero
#' is the value away from that single point.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `lambda`. Neither is read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length `length(y)`.
#'
#' @seealso [distrib_grad_y.Laplace2Distrib()] for the first derivative, which
#'   drops at the location; [distrib_hess_y.LaplaceDistrib()] for the scale
#'   parametrization; [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- laplace2_distrib()
#' y <- c(-1.2, 0.3, 2.5)
#' th <- list(mu = 0.4, lambda = 2)
#'
#' distrib_hess_y(d, y, th)
#'
#' # A location family: the same as the curvature in the location.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # The first derivative drops by 2 lambda across the location.
#' diff(distrib_grad_y(d, 0.4 + c(-1e-9, 1e-9), th))
#' -2 * 2
S7::method(distrib_hess_y, Laplace2Distrib) <- function(distrib, y, theta, ...) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Laplace Distribution, Location and Rate
#'
#' @description
#' Builds the distribution object for the Laplace (double exponential) family
#' written by its **rate** \eqn{\lambda > 0}, with density
#' \eqn{(\lambda/2)\exp(-\lambda|y-\mu|)}. It is the same law as
#' [laplace_distrib()] at \eqn{\lambda = 1/\sigma}, and every quantity that
#' does not involve the second parameter is identical between them.
#'
#' The rate form is the one penalized regression uses: with \eqn{\mu = 0} the
#' negative log-density is \eqn{\lambda|y| - \log(\lambda/2)}, so \eqn{\lambda}
#' is a shrinkage parameter and larger values penalize harder.
#' `penalties7::lasso_penalty()` is this family with the location held at zero.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the location
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the location
#'   ranging over the whole line already.
#' @param link_lambda A `link` object from `linkfunctions7` for the rate
#'   \eqn{\lambda}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \lambda) = \dfrac{\lambda}{2}\exp\left(-\lambda|y-\mu|\right),}
#' with \eqn{\mu \in (-\infty, \infty)} and \eqn{\lambda \in (0, \infty)}. The
#' mean and median are \eqn{\mu}, the variance is \eqn{2/\lambda^2}, the
#' skewness is 0 and the excess kurtosis is 3.
#'
#' Against [laplace_distrib()] the map is \eqn{\lambda = 1/\sigma}. The
#' location is shared, so its score, its information and its derivatives in the
#' response are the same quantities read at the same point; the second
#' parameter's are not, and this parametrization's are the simpler of the two,
#' \eqn{\lambda} being the natural parameter of the exponential family in
#' \eqn{|y - \mu|}.
#'
#' # The kink
#'
#' \eqn{|y - \mu|} is not differentiable at \eqn{y = \mu}, so this family is
#' non-regular in its location exactly as [laplace_distrib()] is: the observed
#' second derivative in \eqn{\mu} is zero, the information is
#' \eqn{\lambda^2}, obtained as the variance of the score, and `params_smooth`
#' is `c(mu = FALSE, lambda = TRUE)`. The consequences are set out on
#' [distrib_expected_hessian.LaplaceDistrib()].
#'
#' A second consequence is visible at the higher orders. The log-density is
#' \eqn{\log(\lambda/2) - \lambda|r|}, linear in \eqn{\lambda} apart from
#' \eqn{\log\lambda} and piecewise linear in \eqn{\mu}, so **every** third and
#' fourth derivative except the pure-\eqn{\lambda} ones is exactly zero, and
#' none of them depends on the data. `expected = TRUE` therefore returns the
#' same values at those orders.
#'
#' # Estimation
#'
#' Both estimates are closed form: \eqn{\hat\mu} is the sample median and
#' \eqn{\hat\lambda} the reciprocal of the mean absolute deviation about it.
#' [fit_distrib()] reaches them by Fisher scoring, which inverts the
#' information above.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location and
#' \eqn{\lambda > 0} the rate, with variance \eqn{2/\lambda^2}. Here
#' \eqn{\lambda} is a distribution parameter; the same letter names a smoothing
#' parameter in `penalties7` and above, and the two meet in the lasso, which is
#' this family with \eqn{\mu} held at zero.
#'
#' @return An S7 object of class `Laplace2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"laplace2"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "lambda")`,
#'   `n_params` `2`, `params_bounds` the list of \eqn{(-\infty, \infty)} and
#'   \eqn{(0, \infty)}, `link_params` the two links given here, and
#'   `params_smooth` `c(mu = FALSE, lambda = TRUE)`.
#'
#' @seealso
#' [laplace_distrib()] for the scale parametrization;
#' [enet_distrib()], which mixes this family with a Gaussian and contains it as
#' the pure-\eqn{\ell_1} case; [fixed()], which holds the location at zero to
#' produce a prior; [fit_distrib()] to estimate the parameters;
#' [Laplace2Distrib] for the class.
#'
#' @references
#' Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001).
#' *The Laplace Distribution and Generalizations*. Birkhauser, Boston.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' *Journal of the Royal Statistical Society, Series B*, **58**(1), 267--288.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats runif
#'
#' @examples
#' d <- laplace2_distrib()
#' d
#'
#' # The same law as the scale parametrization at lambda = 1/sigma.
#' y <- c(-1.2, 0.3, 2.5)
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 1 / 1.5)),
#'           distrib_pdf(laplace_distrib(), y, list(mu = 0.4, sigma = 1.5)))
#'
#' # The variance is 2/lambda^2: a larger rate is a tighter distribution.
#' vapply(c(0.5, 1, 2), function(l) variance(d, list(mu = 0, lambda = l)),
#'        numeric(1))
#'
#' # Both estimates are closed form.
#' set.seed(12)
#' z <- distrib_rng(d, 2000, list(mu = 3, lambda = 0.5))
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = median(z),
#'                  lambda = 1 / mean(abs(z - median(z)))))
#'
#' # Holding the location at zero gives the prior a lasso penalty is written
#' # from: the negative log-density is lambda |y| up to a constant.
#' prior <- fixed(d, mu = 0)
#' prior@params
#' -distrib_pdf(prior, c(-2, 0, 2), list(lambda = 0.5), log = TRUE)
#' 0.5 * abs(c(-2, 0, 2)) - log(0.5 / 2)
#'
#' @export
laplace2_distrib <- function(link_mu = identity_link(), link_lambda = log_link()) {
  Laplace2Distrib(
    distrib_name = "laplace2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "lambda"),
    params_interpretation = c(mu = "location", lambda = "rate"),
    n_params = 2,

    params_bounds = list(
      mu = c(-Inf, Inf),
      lambda = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      lambda = link_lambda
    ),

    params_smooth = c(mu = FALSE, lambda = TRUE)
  )
}
