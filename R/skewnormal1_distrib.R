#' @include distrib.R generics.R
NULL

#' @title Skew Normal Distribution Class
#' @name SkewNormal1Distrib
#'
#' @description
#' The S7 class of Azzalini's skew normal family, a Gaussian carrying a shape
#' parameter that tilts it. With \eqn{z = (y-\mu)/\sigma} the density is
#' \eqn{f(y) = 2\phi(z)\Phi(\alpha z)/\sigma}, so the factor \eqn{2\Phi(\alpha z)}
#' is above one on the side where \eqn{\alpha z > 0} and below one on the other.
#' At \eqn{\alpha = 0} that factor is identically one and the family is the
#' Gaussian.
#'
#' The class inherits from `continuous_distrib`. Its parametrization is
#' \eqn{(\mu, \sigma, \alpha)} with \eqn{\mu} a location and not the mean, since
#' \eqn{E[Y] = \mu + \sigma\delta\sqrt{2/\pi}} with
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}}.
#'
#' Build one with [skewnormal1_distrib()], which supplies the three link
#' functions and fills the properties in. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `SkewNormal1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [skewnormal1_distrib()] they hold `"skew normal1"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma", "alpha")`, the
#'   interpretations `c(mu = "location", sigma = "scale", alpha = "shape")`,
#'   `3`, and the domains \eqn{(-\infty,\infty)}, \eqn{(0,\infty)} and
#'   \eqn{(-\infty,\infty)}.
#'
#' @section Methods:
#' Registered in this file:
#'   [`distrib_pdf()`][distrib_pdf.SkewNormal1Distrib],
#'   [`distrib_cdf()`][distrib_cdf.SkewNormal1Distrib],
#'   [`distrib_rng()`][distrib_rng.SkewNormal1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.SkewNormal1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.SkewNormal1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.SkewNormal1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.SkewNormal1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.SkewNormal1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewNormal1Distrib].
#'
#' Registered elsewhere in the package, all closed form: the four moments
#' [`mean()`][mean.SkewNormal1Distrib], [`variance()`][variance.SkewNormal1Distrib],
#' [`skewness()`][skewness.SkewNormal1Distrib] and
#' [`kurtosis()`][kurtosis.SkewNormal1Distrib] in `moments.R`; the mixed
#' response-parameter derivative [`distrib_cross_y()`][distrib_cross_y] in
#' `cross_derivatives_families.R`; the four orders of the distribution
#' function's parameter derivatives, [`distrib_grad_cdf()`][distrib_grad_cdf]
#' through [`distrib_deriv4_cdf()`][distrib_deriv4_cdf], in
#' `cdf_skewnormal_higher.R`; and the second-order response derivatives
#' [`distrib_cross2_y()`][distrib_cross2_y],
#' [`distrib_grad_y_hess()`][distrib_grad_y_hess] and
#' [`distrib_hess_y_hess()`][distrib_grad_y_hess] in `theta2_families.R`.
#'
#' The **quantile** comes from [continuous_distrib()], by root finding on the
#' distribution function. So does the **expected information**: this family has
#' none in elementary form, so [distrib_expected_hessian()] approximates it by
#' the strategy named in its `approx` argument.
#'
#' @seealso [skewnormal1_distrib()] to build one;
#'   [skewnormal2_distrib()] for the same law in Azzalini's centered
#'   parametrization, whose information is not singular at symmetry;
#'   [skewt_distrib()], which adds degrees of freedom and reaches a skewness
#'   this family cannot;
#'   [gaussian1_distrib()] for the case \eqn{\alpha = 0}.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_interpretation
#' d@params_bounds
#'
#' # The scale rides a log link; the location and the shape are unconstrained.
#' vapply(d@link_params, function(l) l@link_name, character(1))
#'
#' # mu is a location and not the mean.
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#' c(mu = th$mu, mean = mean(d, th))
SkewNormal1Distrib <- S7::new_class("SkewNormal1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Skew Normal Density
#' @name distrib_pdf.SkewNormal1Distrib
#'
#' @description
#' Computes the skew normal density, with \eqn{z = (y-\mu)/\sigma}:
#' \deqn{f(y; \mu, \sigma, \alpha) = \dfrac{2}{\sigma}\,\phi(z)\,\Phi(\alpha z).}
#' The Gaussian density is multiplied by \eqn{2\Phi(\alpha z)}, which exceeds
#' one where \eqn{\alpha z > 0} and falls below it where \eqn{\alpha z < 0}, so
#' a positive \eqn{\alpha} moves mass to the right. At \eqn{\alpha = 0} the
#' factor is one everywhere and the density is \eqn{\phi(z)/\sigma}.
#'
#' The two logarithms are taken separately and added, `dnorm(log = TRUE)`
#' beside `pnorm(log.p = TRUE)`, so the light tail of the skewed side returns a
#' large negative number where forming \eqn{\Phi(\alpha z)} first and then
#' taking its logarithm would return `-Inf`.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations, anywhere on the real line.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`; a component of length 1
#'   is recycled. `sigma` must be strictly positive, and `mu` and `alpha` may
#'   take any finite value.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma), length(alpha))`, one value per
#'   observation.
#'
#' @section Notation:
#' \eqn{\phi} and \eqn{\Phi} are the standard Gaussian density and distribution
#' function, \eqn{\mu} the location, \eqn{\sigma > 0} the scale and \eqn{\alpha}
#' the shape. \eqn{\mu} is not the mean unless \eqn{\alpha = 0}.
#'
#' @seealso [distrib_cdf.SkewNormal1Distrib()] for Owen's T identity,
#'   [distrib_gradient.SkewNormal1Distrib()] for the score,
#'   [skewnormal1_distrib()] for the family, and [distrib_pdf()] for the
#'   generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-2, -0.5, 0.5, 2)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#'
#' # The formula written out.
#' all.equal(distrib_pdf(d, y, th),
#'           2 * dnorm(y) * pnorm(3 * y))
#'
#' # It integrates to one.
#' integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#'
#' # At shape zero the tilting factor is one and the family is the Gaussian.
#' all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 0)), dnorm(y))
#'
#' # The light tail stays a number on the log scale where the density itself
#' # has underflowed to zero.
#' rbind(density = distrib_pdf(d, -40, th),
#'       log_density = distrib_pdf(d, -40, th, log = TRUE))
S7::method(distrib_pdf, SkewNormal1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - mu) / sigma
  # The two logarithms are taken separately, so the far tail of the skewed side
  # gives a large negative number rather than log(0).
  log_d <- log(2) - log(sigma) + stats::dnorm(z, log = TRUE) +
    stats::pnorm(alpha * z, log.p = TRUE)
  if (log) log_d else exp(log_d)
}

#' @title Skew Normal Distribution Function
#' @name distrib_cdf.SkewNormal1Distrib
#'
#' @description
#' Computes the skew normal distribution function through Owen's T function,
#' with \eqn{z = (q-\mu)/\sigma}:
#' \deqn{F(q; \mu, \sigma, \alpha) = \Phi(z) - 2\,T(z, \alpha).}
#' The identity is Azzalini's. Each evaluation costs one bounded
#' one-dimensional quadrature, through [numericals7::owen_t()], where the base
#' class would integrate the density over a semi-infinite range; the bounded
#' integrand is both cheaper and more accurate. A thousand quantiles cost about
#' 20 milliseconds.
#'
#' The result is clamped to \eqn{[0, 1]} before the tail and the logarithm are
#' applied, so rounding in the quadrature cannot return a probability outside
#' the unit interval.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param q A numeric vector of quantiles, anywhere on the real line.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `q`. `sigma` must be
#'   strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}, formed as the
#'   complement.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`. The logarithm is taken after
#'   the clamp, so a quantile far into the light tail returns `-Inf`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms with `log.p = TRUE`, of length
#'   `max(length(q), length(mu), length(sigma), length(alpha))`.
#'
#' @section Notation:
#' \eqn{\Phi} is the standard Gaussian distribution function and
#' \eqn{T(h, a) = (2\pi)^{-1}\int_0^a e^{-h^2(1+x^2)/2}(1+x^2)^{-1}dx} is
#' Owen's T. \eqn{\alpha} is the shape.
#'
#' @seealso [numericals7::owen_t()] for the quadrature,
#'   [distrib_pdf.SkewNormal1Distrib()] for the density it integrates,
#'   [distrib_grad_cdf()] for its parameter derivatives, and [distrib_cdf()]
#'   for the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' q <- c(-2, -0.5, 0.5, 2)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#'
#' # Owen's identity against a direct quadrature of the density.
#' rbind(owen = distrib_cdf(d, q, th),
#'       quadrature = vapply(q, function(u)
#'         integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#'
#' # At shape zero, T(z, 0) = 0 and the identity is the Gaussian's.
#' all.equal(distrib_cdf(d, q, list(mu = 0, sigma = 1, alpha = 0)), pnorm(q))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, q, th) + distrib_cdf(d, q, th, lower.tail = FALSE)
S7::method(distrib_cdf, SkewNormal1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  z <- (q - theta[[1]]) / theta[[2]]
  res <- stats::pnorm(z) - 2 * numericals7::owen_t(z, theta[[3]])
  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Skew Normal Random Generation
#' @name distrib_rng.SkewNormal1Distrib
#'
#' @description
#' Draws from the skew normal exactly, using Azzalini's stochastic
#' representation: with \eqn{U_0} and \eqn{U_1} independent standard Gaussians
#' and \eqn{\delta = \alpha/\sqrt{1+\alpha^2}},
#' \deqn{Z = \delta\,|U_0| + \sqrt{1-\delta^2}\;U_1}
#' is standard skew normal with shape \eqn{\alpha}, and \eqn{Y = \mu + \sigma Z}.
#' No inversion and no rejection is involved. The cost is two `rnorm` calls
#' whatever the shape is, and the draws follow the density exactly.
#'
#' The representation also reads off why the family cannot be very skewed: the
#' half-normal component enters with weight \eqn{\delta}, which saturates at one
#' as \eqn{\alpha \to \infty}, so the most skewed member is the half-normal
#' itself.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of length `n`; a component of length 1 is
#'   recycled, so a parameter varying by observation draws one value per
#'   observation from its own member of the family.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @section Notation:
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} is the correlation between \eqn{Z}
#' and the half-normal component; \eqn{\mu} is the location and \eqn{\sigma} the
#' scale.
#'
#' @seealso [distrib_pdf.SkewNormal1Distrib()] for the density the draws follow,
#'   [mean.SkewNormal1Distrib()] and [variance.SkewNormal1Distrib()] for the
#'   moments they should reproduce, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#'
#' set.seed(1)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # Three moments against their closed forms.
#' rbind(sample = c(mean(x), var(x), mean((x - mean(x))^3) / sd(x)^3),
#'       theory = c(mean(d, th), variance(d, th), skewness(d, th)))
#'
#' # The representation, written out at the same seed.
#' delta <- 3 / sqrt(1 + 9)
#' set.seed(1)
#' z <- delta * abs(rnorm(2e5)) + sqrt(1 - delta^2) * rnorm(2e5)
#' all.equal(x, z)
S7::method(distrib_rng, SkewNormal1Distrib) <- function(distrib, n, theta, ...) {
  alpha <- theta[[3]]
  delta <- alpha / sqrt(1 + alpha^2)
  z <- delta * abs(stats::rnorm(n)) + sqrt(1 - delta^2) * stats::rnorm(n)
  theta[[1]] + theta[[2]] * z
}

#' @title Skew Normal Score
#' @name distrib_gradient.SkewNormal1Distrib
#'
#' @description
#' Computes the three first derivatives of the log-density in closed form. With
#' \eqn{z = (y-\mu)/\sigma}, \eqn{t = \alpha z} and
#' \eqn{R(t) = \phi(t)/\Phi(t)} the inverse Mills ratio,
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{z - \alpha R}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{z^2 - 1 - \alpha z R}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \alpha} = z R.}
#' Every derivative of \eqn{\log\Phi(t)} is a polynomial in \eqn{t} and
#' \eqn{R}, because \eqn{R' = -R(t+R)} closes the recursion, so the whole
#' derivative surface of this family stays elementary.
#'
#' The ratio comes from [numericals7::mills_ratio()], which forms it on the log
#' scale. Below about \eqn{t = -38} both \eqn{\phi(t)} and \eqn{\Phi(t)}
#' underflow while their ratio is finite and close to \eqn{-t}: measured at
#' \eqn{t = -400} the ratio is 400.0025 where `dnorm(t)/pnorm(t)` is `NaN`.
#'
#' At \eqn{\alpha = 0} the shape score is \eqn{z\sqrt{2/\pi}} and the location
#' score is \eqn{z/\sigma}, so the two are exactly proportional. The expected
#' information is singular at symmetry for that reason.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`. `sigma` must be
#'   strictly positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. On the link
#'   scale the derivatives are taken with respect to the unconstrained
#'   coordinates \eqn{\eta} of the three link functions. The transformation is
#'   applied in the generic's body, so this method always returns the parameter
#'   scale and the argument is here to match the generic's signature.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu`, `sigma` and `alpha`,
#'   each of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\phi} and \eqn{\Phi}
#' the standard Gaussian density and distribution function, and
#' \eqn{R(t) = \phi(t)/\Phi(t)} the inverse Mills ratio.
#'
#' @seealso [distrib_hessian.SkewNormal1Distrib()] for the second derivatives,
#'   [distrib_grad_y.SkewNormal1Distrib()] for the derivative in the response,
#'   [numericals7::mills_ratio()] for the ratio, and [distrib_gradient()] for
#'   the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#' g <- distrib_gradient(d, y, th)
#'
#' # The shape component written out.
#' all.equal(g$alpha, y * numericals7::mills_ratio(3 * y)$r)
#'
#' # The score sums to nearly zero at the maximum likelihood estimate.
#' set.seed(3)
#' x <- distrib_rng(d, 4000, th)
#' fit <- fit_distrib(d, x)
#' vapply(distrib_gradient(d, x, as.list(coef(fit))), sum, 0) / 4000
#'
#' # At symmetry the shape score is a fixed multiple of the location score,
#' # which is where this parametrization loses rank.
#' g0 <- distrib_gradient(d, y, list(mu = 0, sigma = 1.4, alpha = 0))
#' c(ratio = unique(round(g0$alpha / g0$mu, 12)), sigma_root_2_pi = 1.4 * sqrt(2 / pi))
S7::method(distrib_gradient, SkewNormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - mu) / sigma
  m <- numericals7::mills_ratio(alpha * z)
  list(
    mu = (z - alpha * m$r) / sigma,
    sigma = (z^2 - 1 - alpha * z * m$r) / sigma,
    alpha = z * m$r
  )
}

#' @title Skew Normal Observed Hessian
#' @name distrib_hessian.SkewNormal1Distrib
#'
#' @description
#' Computes the six second derivatives of the log-density in closed form. In
#' the notation of [distrib_gradient.SkewNormal1Distrib()], with
#' \eqn{R' = -R(t+R)},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{\alpha^2 R' - 1}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{\alpha^2 z R' - 2z + \alpha R}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \alpha}
#'         = -\dfrac{R + \alpha z R'}{\sigma},}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = \dfrac{1 - 3z^2 + 2\alpha z R + \alpha^2 z^2 R'}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma \, \partial \alpha}
#'         = -\dfrac{z R + \alpha z^2 R'}{\sigma},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \alpha^2} = z^2 R'.}
#'
#' This is the **observed** curvature at the data. The family has no elementary
#' expected information, so [distrib_expected_hessian()] falls to the base
#' class and approximates it; see there for the strategies and their cost.
#'
#' @details
#' # Singularity at symmetry
#'
#' At \eqn{\alpha = 0} the expected information of this parametrization has rank
#' 2, not 3. The reason is in the score: the shape and location components are
#' exactly proportional there, so no data can separate them. Measured on the
#' approximated information at \eqn{\mu = 0}, \eqn{\sigma = 1}, its eigenvalues
#' are 2, 1.637 and \eqn{-5.6\times10^{-17}}, and the smallest one grows like
#' \eqn{\alpha^4} as the shape moves off zero: \eqn{4.4\times10^{-10}} at
#' \eqn{\alpha = 0.01} and \eqn{1.9\times10^{-3}} at \eqn{\alpha = 0.5}.
#'
#' The singularity belongs to the family as parametrized, not to the
#' implementation. Azzalini's centered parametrization removes it, and lives in
#' the separate object [skewnormal2_distrib()].
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`. `sigma` must be
#'   strictly positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors, in the order `mu_mu`,
#'   `sigma_sigma`, `alpha_alpha`, `mu_sigma`, `mu_alpha`, `sigma_alpha`, each
#'   of the length of the recycled inputs. The names are the diagonal first,
#'   then the off-diagonal pairs, which is [hess_names()]'s ordering.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{t = \alpha z}, \eqn{R(t) = \phi(t)/\Phi(t)}
#' and \eqn{R' = -R(t+R)}.
#'
#' @seealso [distrib_gradient.SkewNormal1Distrib()] for the score,
#'   [distrib_deriv3.SkewNormal1Distrib()] for the next order,
#'   [distrib_expected_hessian()] for the approximated expectation, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#' h <- distrib_hessian(d, y, th)
#'
#' # The shape component written out.
#' all.equal(h$alpha_alpha, y^2 * numericals7::mills_ratio(3 * y)$dr)
#'
#' # Against a central difference of the score.
#' eps <- 1e-5
#' rbind(analytic = h$mu_alpha,
#'       numeric = (distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 3 + eps))$mu -
#'                  distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 3 - eps))$mu) /
#'                 (2 * eps))
#'
#' # The information loses rank at symmetry, and recovers it as alpha^4.
#' rank_gap <- function(a) {
#'   e <- distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, alpha = a))
#'   M <- matrix(c(e$mu_mu, e$mu_sigma, e$mu_alpha,
#'                 e$mu_sigma, e$sigma_sigma, e$sigma_alpha,
#'                 e$mu_alpha, e$sigma_alpha, e$alpha_alpha), 3, 3)
#'   min(abs(eigen(-M, only.values = TRUE)$values))
#' }
#' vapply(c(0, 0.01, 0.5), rank_gap, 0)
S7::method(distrib_hessian, SkewNormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - mu) / sigma
  m <- numericals7::mills_ratio(alpha * z)
  r <- m$r
  dr <- m$dr
  s2 <- sigma^2
  list(
    mu_mu = (alpha^2 * dr - 1) / s2,
    sigma_sigma = (1 - 3 * z^2 + 2 * alpha * z * r + alpha^2 * z^2 * dr) / s2,
    alpha_alpha = z^2 * dr,
    mu_sigma = (alpha^2 * z * dr - 2 * z + alpha * r) / s2,
    mu_alpha = -(r + alpha * z * dr) / sigma,
    sigma_alpha = -(z * r + alpha * z^2 * dr) / sigma
  )
}

#' @title Skew Normal Third Derivatives
#' @name distrib_deriv3.SkewNormal1Distrib
#'
#' @description
#' Computes the ten third derivatives of the log-density in closed form, in a
#' compiled kernel. The family's whole derivative surface stays elementary for
#' one reason: with \eqn{t = \alpha z} and \eqn{R(t) = \phi(t)/\Phi(t)} the
#' inverse Mills ratio, \eqn{R' = -R(t+R)} closes the recursion, so every
#' derivative of \eqn{\log\Phi(t)} is a polynomial in \eqn{t} and \eqn{R} and
#' no new special function appears at any order.
#'
#' With `expected = TRUE` the value is an expectation instead, and there it is
#' **numerical**: the integrals are the ones that block the expected
#' information, so the call routes to `expected_derivative()` and `approx` and
#' `nsim` are read.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` its values
#'   are ignored and only its length matters, the result being one expectation
#'   repeated.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned from the compiled kernel. When `TRUE` the
#'   expectation is approximated by the strategy in `approx`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body, so this method always
#'   returns the parameter scale.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, the
#'   strategy for the expectation. Read only when `expected = TRUE`; the first
#'   is the default and quadratures the observed kernel against the density.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Ignored by every other strategy. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of ten numeric vectors, one per distinct third-order
#'   component, named as [deriv_names()] names them: `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_mu_alpha`, `mu_sigma_sigma`, `mu_sigma_alpha`, `mu_alpha_alpha`,
#'   `sigma_sigma_sigma`, `sigma_sigma_alpha`, `sigma_alpha_alpha` and
#'   `alpha_alpha_alpha`.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{t = \alpha z}, \eqn{R} the inverse Mills
#' ratio and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_hessian.SkewNormal1Distrib()] for the order below,
#'   [distrib_deriv4.SkewNormal1Distrib()] for the order above,
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Against a central difference of the analytic Hessian.
#' eps <- 1e-4
#' rbind(analytic = d3$mu_mu_alpha,
#'       numeric = (distrib_hessian(d, y, list(mu = 0, sigma = 1, alpha = 3 + eps))$mu_mu -
#'                  distrib_hessian(d, y, list(mu = 0, sigma = 1, alpha = 3 - eps))$mu_mu) /
#'                 (2 * eps))
#'
#' # At shape zero the third derivative in the location is the Gaussian's,
#' # which is exactly zero.
#' distrib_deriv3(d, y, list(mu = 0, sigma = 1, alpha = 0))$mu_mu_mu
#'
#' # The expectation is numerical here, and the strategy is read.
#' set.seed(4)
#' c(integrate = distrib_deriv3(d, 0, th, expected = TRUE)$alpha_alpha_alpha,
#'   mc = distrib_deriv3(d, 0, th, expected = TRUE, approx = "mc",
#'                       nsim = 2000)$alpha_alpha_alpha)
S7::method(distrib_deriv3, SkewNormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    skewnormal_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

#' @title Skew Normal Fourth Derivatives
#' @name distrib_deriv4.SkewNormal1Distrib
#'
#' @description
#' Computes the fifteen fourth derivatives of the log-density in closed form,
#' in a compiled kernel, in the notation of
#' [distrib_deriv3.SkewNormal1Distrib()]. The recursion
#' \eqn{R' = -R(t+R)} keeps every one of them a polynomial in \eqn{t} and the
#' inverse Mills ratio, so the fourth order needs nothing the third did not.
#'
#' With `expected = TRUE` the value is a numerical expectation, for the reason
#' it is at third order: the integrals are the ones that block the expected
#' information.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length matters.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned from the compiled kernel.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of fifteen numeric vectors, one per distinct
#'   fourth-order component, named as [deriv_names()] names them, from
#'   `mu_mu_mu_mu` to `alpha_alpha_alpha_alpha`.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{t = \alpha z} and \eqn{R} the inverse Mills
#' ratio.
#'
#' @seealso [distrib_deriv3.SkewNormal1Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and [skewnormal1_distrib()] for the
#'   family.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#' d4 <- distrib_deriv4(d, y, th)
#' length(d4)
#'
#' # Against a central difference of the third order.
#' eps <- 1e-4
#' rbind(analytic = d4$mu_mu_alpha_alpha,
#'       numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1, alpha = 3 + eps))$mu_mu_alpha -
#'                  distrib_deriv3(d, y, list(mu = 0, sigma = 1, alpha = 3 - eps))$mu_mu_alpha) /
#'                 (2 * eps))
#'
#' # The count is bit for bit independent of the thread count.
#' identical(distrib_deriv4(d, y, th, threads = 1L),
#'           distrib_deriv4(d, y, th, threads = 2L))
S7::method(distrib_deriv4, SkewNormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    skewnormal_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

#' @title Skew Normal Response Derivative
#' @name distrib_grad_y.SkewNormal1Distrib
#'
#' @description
#' Computes \eqn{\partial\ell/\partial y}, the derivative of the log-density in
#' the response. With \eqn{z = (y-\mu)/\sigma} and \eqn{R} the inverse Mills
#' ratio at \eqn{t = \alpha z},
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\alpha R - z}{\sigma}.}
#' This is minus the derivative in \eqn{\mu}, exactly, because the response and
#' the location enter the density only through their difference. The identity
#' holds for every location family, so this method is one sign away from
#' [distrib_gradient.SkewNormal1Distrib()].
#'
#' The quantity is what a censored likelihood and a quantile residual need, and
#' it is the first ingredient of [distrib_cross_y()].
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the length of the recycled inputs, one value per
#'   observation.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{z = (y-\mu)/\sigma}
#' and \eqn{R(t) = \phi(t)/\Phi(t)} the inverse Mills ratio.
#'
#' @seealso [distrib_hess_y.SkewNormal1Distrib()] for the second derivative,
#'   [distrib_cross_y()] for the mixed one,
#'   [distrib_gradient.SkewNormal1Distrib()] for the parameter derivatives, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#'
#' # It is minus the score in the location, exactly.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Against a central difference of the log-density.
#' eps <- 1e-6
#' rbind(analytic = distrib_grad_y(d, y, th),
#'       numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
#'                  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
S7::method(distrib_grad_y, SkewNormal1Distrib) <- function(distrib, y, theta, ...) {
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - theta[[1]]) / sigma
  (alpha * numericals7::mills_ratio(alpha * z)$r - z) / sigma
}

#' @title Skew Normal Second Response Derivative
#' @name distrib_hess_y.SkewNormal1Distrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y^2}, the curvature of the log-density
#' in the response. With \eqn{R' = -R(t+R)} at \eqn{t = \alpha z},
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2}
#'       = \dfrac{\alpha^2 R' - 1}{\sigma^2},}
#' which is the same expression as \eqn{\partial^2\ell/\partial\mu^2}: two
#' signs cancel where one did not at first order, so the response curvature
#' equals the location curvature instead of being its negative.
#'
#' Since \eqn{R'} is negative for every \eqn{t}, the value is strictly negative
#' at every observation and every shape. The skew normal log-density is
#' concave in the response, as the Gaussian's is.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the length of the recycled inputs, negative
#'   throughout.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{t = \alpha z}, \eqn{R} the inverse Mills
#' ratio and \eqn{R' = -R(t+R)} its derivative.
#'
#' @seealso [distrib_grad_y.SkewNormal1Distrib()] for the first derivative,
#'   [distrib_hessian.SkewNormal1Distrib()] for the parameter curvature it
#'   shares an expression with, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- skewnormal1_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#'
#' # It equals the curvature in the location, without a sign change.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Against a central difference of the response derivative.
#' eps <- 1e-5
#' rbind(analytic = distrib_hess_y(d, y, th),
#'       numeric = (distrib_grad_y(d, y + eps, th) -
#'                  distrib_grad_y(d, y - eps, th)) / (2 * eps))
#'
#' # Concave in the response at every shape, and bounded above by
#' # -1/sigma^2, which is the Gaussian's own curvature.
#' vapply(c(-20, -1, 0, 1, 20), function(a)
#'   max(distrib_hess_y(d, seq(-6, 6, by = 0.5), list(mu = 0, sigma = 1, alpha = a))), 0)
S7::method(distrib_hess_y, SkewNormal1Distrib) <- function(distrib, y, theta, ...) {
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  z <- (y - theta[[1]]) / sigma
  (alpha^2 * numericals7::mills_ratio(alpha * z)$dr - 1) / sigma^2
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Skew Normal Distribution Object
#'
#' @description
#' Builds a skew normal distribution object in Azzalini's direct
#' parametrization: a location \eqn{\mu}, a scale \eqn{\sigma > 0} and a shape
#' \eqn{\alpha} that tilts the Gaussian. With \eqn{z = (y-\mu)/\sigma} the
#' density is \eqn{2\phi(z)\Phi(\alpha z)/\sigma}, and \eqn{\alpha = 0} gives
#' the Gaussian exactly.
#'
#' The returned object carries the three link functions and every method of the
#' family. `mu` and `sigma` are a location and a scale, not the mean and the
#' standard deviation; [mean.SkewNormal1Distrib()] and
#' [variance.SkewNormal1Distrib()] give those.
#'
#' @param link_mu A `linkfunctions7` link object for the location \eqn{\mu},
#'   which is unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_sigma A link object for the scale \eqn{\sigma}, which must be
#'   strictly positive. Defaults to [linkfunctions7::log_link()], so that any
#'   real linear predictor maps to an admissible scale.
#' @param link_alpha A link object for the shape \eqn{\alpha}, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()].
#'
#' @details
#' # Density and distribution function
#'
#' \deqn{f(y; \mu, \sigma, \alpha) = \dfrac{2}{\sigma}\,\phi(z)\,\Phi(\alpha z),
#'       \qquad z = (y-\mu)/\sigma.}
#' The factor \eqn{2\Phi(\alpha z)} exceeds one where \eqn{\alpha z > 0} and
#' falls below it where \eqn{\alpha z < 0}, so a positive \eqn{\alpha} moves
#' mass to the right. The distribution function is Azzalini's identity
#' \eqn{F(q) = \Phi(z) - 2T(z, \alpha)} with \eqn{T} Owen's T, one bounded
#' quadrature per observation through [numericals7::owen_t()]. The quantile
#' function has no closed form and comes from the base class by root finding.
#'
#' # Derivatives
#'
#' The score and the observed Hessian are closed form, written in the inverse
#' Mills ratio \eqn{R(t) = \phi(t)/\Phi(t)} at \eqn{t = \alpha z}, and so are
#' the third and fourth orders, in compiled kernels. One identity does all of
#' it: \eqn{R' = -R(t+R)}, so every derivative of \eqn{\log\Phi(t)} is a
#' polynomial in \eqn{t} and \eqn{R}.
#'
#' The **expected** information has no elementary form, so none is registered
#' and [distrib_expected_hessian()] approximates it by the strategy named in
#' its `approx` argument. The expected third and fourth orders share that
#' obstruction.
#'
#' # Singularity at symmetry
#'
#' At \eqn{\alpha = 0} the expected information has rank 2 and not 3. The score
#' shows why: the shape component is \eqn{z\sqrt{2/\pi}} and the location
#' component \eqn{z/\sigma}, so the two are exactly proportional and no data
#' can separate them. Measured, the smallest eigenvalue is
#' \eqn{-5.6\times10^{-17}} against a largest of 2 at \eqn{\alpha = 0}, and it
#' grows like \eqn{\alpha^4} thereafter.
#'
#' The consequence for use is that a fit whose true shape is near zero
#' identifies \eqn{\alpha} weakly, and that a variance matrix computed at
#' exactly zero is not invertible. Azzalini and Capitanio's centered
#' parametrization removes the singularity; it is [skewnormal2_distrib()], a
#' separate object.
#'
#' # Moments, and the bound on the skewness
#'
#' With \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and \eqn{b = \sqrt{2/\pi}},
#' \deqn{E[Y] = \mu + \sigma b \delta, \qquad
#'       \mathrm{Var}(Y) = \sigma^2(1 - b^2\delta^2).}
#' All four moments are closed form. The skewness is bounded: \eqn{\delta}
#' saturates at one as \eqn{|\alpha|\to\infty}, so \eqn{\gamma_1} cannot leave
#' \eqn{(-0.9953, 0.9953)} whatever the shape is. Measured, it reaches 0.9556
#' at \eqn{\alpha = 10} and 0.99527 at \eqn{\alpha = 10^4}, and does not move
#' after that. A sample skewer than this needs [skewt_distrib()].
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (-\infty, \infty)}
#' - \eqn{\sigma \in (0, \infty)}
#' - \eqn{\alpha \in (-\infty, \infty)}
#'
#' @return An S7 object of class [SkewNormal1Distrib], inheriting from
#'   `continuous_distrib`. Its `params` are `c("mu", "sigma", "alpha")`, its
#'   `bounds` `c(-Inf, Inf)`, and its `link_params` the three links given here.
#'
#' @references
#' Azzalini, A. (1985). A class of distributions which includes the normal ones.
#' *Scandinavian Journal of Statistics* 12, 171-178.
#'
#' Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
#' Families*. Cambridge University Press.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dnorm pnorm rnorm integrate
#'
#' @examples
#' d <- skewnormal1_distrib()
#' d@params
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#'
#' # The location is not the mean, and the scale is not the standard deviation.
#' rbind(parameter = c(th$mu, th$sigma),
#'       moment = c(mean(d, th), sqrt(variance(d, th))))
#'
#' # Shape zero is the Gaussian, in the density and in the distribution function.
#' y <- c(-1, 0, 1)
#' th0 <- list(mu = 0, sigma = 1, alpha = 0)
#' c(density = max(abs(distrib_pdf(d, y, th0) - dnorm(y))),
#'   cdf = max(abs(distrib_cdf(d, y, th0) - pnorm(y))))
#'
#' # The skewness the family can reach is bounded, and saturates early.
#' vapply(c(1, 3, 10, 50, 1e4),
#'        function(a) skewness(d, list(mu = 0, sigma = 1, alpha = a)), 0)
#'
#' # A fit recovers all three parameters.
#' set.seed(5)
#' x <- distrib_rng(d, 3000, list(mu = 2, sigma = 1.5, alpha = 4))
#' coef(fit_distrib(d, x))
#'
#' @seealso [skewnormal2_distrib()] for the centered parametrization,
#'   [skewt_distrib()] for the four-parameter family that reaches a larger
#'   skewness, [gaussian1_distrib()] for the case \eqn{\alpha = 0}, and
#'   [SkewNormal1Distrib] for the class and its method list.
#' @export
skewnormal1_distrib <- function(link_mu = identity_link(),
                               link_sigma = log_link(),
                               link_alpha = identity_link()) {
  SkewNormal1Distrib(
    distrib_name = "skew normal1",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "sigma", "alpha"),
    params_interpretation = c(mu = "location", sigma = "scale", alpha = "shape"),
    n_params = 3,

    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf),
      alpha = c(-Inf, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma,
      alpha = link_alpha
    )
  )
}
