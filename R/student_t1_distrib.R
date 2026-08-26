#' @include distrib.R generics.R
NULL

#' @title Student t Distribution Class, Location, Scale and Degrees of Freedom
#' @name StudentT1Distrib
#'
#' @description
#' The S7 class of the location-scale Student t family, parametrized by a
#' location \eqn{\mu}, a scale \eqn{\sigma > 0} and degrees of freedom
#' \eqn{\nu > 0}, with density proportional to
#' \eqn{\{1 + (y-\mu)^2/(\nu\sigma^2)\}^{-(\nu+1)/2}} on the whole real line.
#' It inherits from `continuous_distrib`, so it answers every generic of the
#' `distrib` contract; the eleven methods listed below are registered on it
#' directly and everything else comes from the parent.
#'
#' Build one with [student_t1_distrib()], which supplies the three link
#' functions and fills the properties in. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `StudentT1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [student_t1_distrib()] they hold `"student t1"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma", "nu")`, the
#'   interpretations `c(mu = "location", sigma = "scale", nu = "shape")`, `3`,
#'   the domains \eqn{(-\infty, \infty)} and \eqn{(0, \infty)} twice, and the
#'   three links.
#'
#' @seealso [student_t1_distrib()] to build one;
#'   [student_t2_distrib()] for the same law parametrized by the standard
#'   deviation; [gaussian1_distrib()] for the limit at \eqn{\nu \to \infty};
#'   [cauchy_distrib()] for the case \eqn{\nu = 1};
#'   [distrib_gradient.StudentT1Distrib()] for the redescending score.
#'
#' @section Methods:
#' Registered on this class in this file:
#'   [`distrib_cdf()`][distrib_cdf.StudentT1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.StudentT1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.StudentT1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.StudentT1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.StudentT1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.StudentT1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.StudentT1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.StudentT1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.StudentT1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.StudentT1Distrib],
#'   [`distrib_rng()`][distrib_rng.StudentT1Distrib].
#'
#' Nine more are registered on the class from other files: the mixed
#' derivatives [`distrib_cross_y()`][distrib_cross_y.StudentT1Distrib] and
#'   [`distrib_cross2_y()`][distrib_cross2_y.StudentT1Distrib], the
#'   distribution-function derivatives
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.StudentT1Distrib] and
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.StudentT1Distrib],
#'   `distrib_grad_y_hess()` and `distrib_hess_y_hess()` in
#'   `theta2_families.R`, and the four moments
#'   [`mean()`][mean.StudentT1Distrib], [`variance()`][variance.StudentT1Distrib],
#'   [`skewness()`][skewness.StudentT1Distrib] and
#'   [`kurtosis()`][kurtosis.StudentT1Distrib] in `moments.R`.
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@params_bounds
#'
#' # The location is already free; the scale and the degrees of freedom ride
#' # a log, both being positive.
#' vapply(d@link_params, function(l) l@link_name, character(1))
StudentT1Distrib <- S7::new_class("StudentT1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Student t Probability Density Function
#' @name distrib_pdf.StudentT1Distrib
#' @description
#' Computes the location-scale Student t density
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\left(\dfrac{\nu}{2}\right)} \left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}}
#' by calling [stats::dt()] at the standardized value \eqn{(y-\mu)/\sigma} and
#' dividing by \eqn{\sigma}, the Jacobian of the standardization. With
#' `log = TRUE` the division becomes a subtraction of \eqn{\log\sigma}, so the
#' logarithm stays finite far into the tails.
#'
#' The tails are polynomial, of order \eqn{|y|^{-(\nu+1)}}, so the density
#' decays far more slowly than a Gaussian's and no value of \eqn{y} underflows
#' at an ordinary \eqn{\nu}.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive; a zero or
#'   negative value gives `NaN` with a warning from [stats::dt()].
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`, one value per
#'   observation.
#'
#' @seealso [distrib_cdf.StudentT1Distrib()] for the distribution function,
#'   [distrib_gradient.StudentT1Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf()] for the generic and
#'   [student_t1_distrib()] for the family.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#'
#' # The method is stats::dt at the standardized value, over sigma.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 5)),
#'           dt((y - 0.4) / 1.2, df = 5) / 1.2)
#'
#' # One degree of freedom is the Cauchy.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 1)),
#'           dcauchy(y, location = 0.4, scale = 1.2))
#'
#' # The tails are polynomial, so a value forty scales out still carries mass
#' # where a Gaussian's density has underflowed to zero.
#' c(t = distrib_pdf(d, 48, list(mu = 0, sigma = 1.2, nu = 5)),
#'   gaussian = dnorm(48, 0, 1.2))
S7::method(distrib_pdf, StudentT1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  val <- stats::dt(
    x = (y - theta[[1]]) / theta[[2]],
    df = theta[[3]],
    log = log
  )
  if (log) {
    val - log(theta[[2]])
  } else {
    val / theta[[2]]
  }
}

#' @title Student t Cumulative Distribution Function
#' @name distrib_cdf.StudentT1Distrib
#' @description
#' Computes the location-scale Student t distribution function
#' \deqn{F(q; \mu, \sigma, \nu) = T_\nu\!\left(\dfrac{q-\mu}{\sigma}\right)}
#' with \eqn{T_\nu} the standard Student t distribution function on \eqn{\nu}
#' degrees of freedom, by calling [stats::pt()] at the standardized value.
#' Both tails are available exactly: `lower.tail = FALSE` evaluates
#' \eqn{1 - F} without forming the difference, and `log.p = TRUE` returns a
#' logarithm that stays finite where the probability itself underflows.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `q`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma), length(nu))`. With
#'   `log.p = TRUE` the values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.StudentT1Distrib()] for the inverse,
#'   [distrib_pdf.StudentT1Distrib()] for the density,
#'   [distrib_grad_cdf.StudentT1Distrib()] for the derivatives of this function
#'   in the parameters, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#'
#' # The method is stats::pt at the standardized value.
#' all.equal(distrib_cdf(d, c(-2.5, 0.3, 1.8), th),
#'           pt((c(-2.5, 0.3, 1.8) - 0.4) / 1.2, df = 5))
#'
#' # The law is symmetric about mu, so the two tails at equal distance match.
#' c(distrib_cdf(d, 0.4 - 2, th),
#'   distrib_cdf(d, 0.4 + 2, th, lower.tail = FALSE))
#'
#' # The upper tail decays polynomially, so at forty scales out it is still
#' # representable where a Gaussian's has underflowed.
#' c(t = distrib_cdf(d, 48, list(mu = 0, sigma = 1.2, nu = 5),
#'                   lower.tail = FALSE),
#'   gaussian = pnorm(48, 0, 1.2, lower.tail = FALSE))
S7::method(distrib_cdf, StudentT1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  stats::pt(
    q = (q - theta[[1]]) / theta[[2]],
    df = theta[[3]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Student t Quantile Function
#' @name distrib_quantile.StudentT1Distrib
#' @description
#' Computes the location-scale Student t quantile function
#' \deqn{Q(p; \mu, \sigma, \nu) = \mu + \sigma\, T_\nu^{-1}(p)}
#' with \eqn{T_\nu^{-1}} the standard Student t quantile function on \eqn{\nu}
#' degrees of freedom, by calling [stats::qt()]. The distribution function is
#' strictly increasing on the whole line, so the inverse is exact and unique
#' and the root-finding fallback the base class supplies is bypassed. `Q(0)` is
#' `-Inf` and `Q(1)` is `Inf`.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN` with
#'   a warning from [stats::qt()].
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `p`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is read as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles on \eqn{[-\infty, \infty]}, of length
#'   `max(length(p), length(mu), length(sigma), length(nu))`.
#'
#' @seealso [distrib_cdf.StudentT1Distrib()] for the function inverted here,
#'   [distrib_rng.StudentT1Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#'
#' # The quartiles, and the round trip back through the distribution function.
#' q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#' q
#' all.equal(distrib_cdf(d, q, th), c(0.25, 0.5, 0.75))
#'
#' # The median is the location, the law being symmetric about it.
#' distrib_quantile(d, 0.5, th)
#'
#' # Heavy tails put the extreme quantiles far further out than a Gaussian's.
#' rbind(t = distrib_quantile(d, c(0.001, 0.999), th),
#'       gaussian = qnorm(c(0.001, 0.999), 0.4, 1.2))
S7::method(distrib_quantile, StudentT1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  theta[[1]] + theta[[2]] * stats::qt(
    p = p,
    df = theta[[3]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Student t Random Number Generator
#' @name distrib_rng.StudentT1Distrib
#' @description
#' Draws `n` independent location-scale Student t variates as
#' \eqn{\mu + \sigma T} with \eqn{T} from [stats::rt()], so the draws come from
#' R's own generator and depend on `.Random.seed` in the usual way. The
#' generalized ratio-of-uniforms fallback the base class supplies is bypassed.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of length `n`. A component of length 1 is
#'   recycled, so a vector of length `n` draws one variate per parameter
#'   setting. `sigma` and `nu` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.StudentT1Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#'
#' # The draws are stats::rt shifted and scaled, so the same seed reproduces
#' # them.
#' set.seed(12)
#' a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.2, nu = 5))
#' set.seed(12)
#' identical(a, 0.4 + 1.2 * rt(3, df = 5))
#'
#' # The sample variance recovers sigma^2 * nu / (nu - 2), not sigma^2.
#' set.seed(13)
#' z <- distrib_rng(d, 5e4, list(mu = 0.4, sigma = 1.2, nu = 5))
#' c(sample = var(z), theoretical = 1.2^2 * 5 / 3, sigma_sq = 1.2^2)
S7::method(distrib_rng, StudentT1Distrib) <- function(distrib, n, theta, ...) {
  theta[[1]] + theta[[2]] * stats::rt(
    n = n,
    df = theta[[3]]
  )
}

#' @title Student t Score
#' @name distrib_gradient.StudentT1Distrib
#' @description
#' Computes the first derivatives of the location-scale Student t log-density
#' with respect to \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}, one value per
#' observation, in closed form. Writing \eqn{r = y - \mu} and
#' \eqn{D = \nu\sigma^2 + r^2},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{(\nu+1)r}{D}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{\nu\left(r^2 - \sigma^2\right)}{\sigma D},}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left[
#'   -\dfrac{1}{\nu} - \psi\!\left(\dfrac{\nu}{2}\right)
#'   + \psi\!\left(\dfrac{\nu+1}{2}\right)
#'   + \dfrac{(\nu+1)r^2}{\nu D}
#'   - \log\!\left(1 + \dfrac{r^2}{\nu\sigma^2}\right)\right].}
#'
#' The location component **redescends**: it rises to \eqn{(\nu+1)/(2\sigma\sqrt{\nu})}
#' at \eqn{|r| = \sigma\sqrt{\nu}} and falls back towards zero beyond it, so a
#' gross outlier contributes almost nothing to the estimating equation. A
#' Gaussian's location score grows without bound instead, which is the whole
#' reason this family is used for robust location estimation.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale; the transformation happens in the generic. The arithmetic runs in a
#' compiled kernel decomposed over the elements of the output, so the result
#' does not depend on the thread count.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu`, `sigma` and `nu`, each
#'   of length `max(length(y), length(mu), length(sigma), length(nu))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale and \eqn{\nu > 0} the degrees of freedom.
#' \eqn{\psi} is the digamma function, `digamma()` in R.
#'
#' @section Large degrees of freedom:
#' Every component is written as a ratio in \eqn{z = r/\sigma} and
#' \eqn{u = z^2/\nu} instead of as a quotient by powers of \eqn{D}, because
#' \eqn{\nu\sigma^2} overflows well before the log link's own clamp is reached.
#' All three components stay finite at every \eqn{\nu} the chart can produce,
#' up to `.Machine$double.xmax`.
#'
#' @seealso [distrib_hessian.StudentT1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.StudentT1Distrib()] for their expectation,
#'   [distrib_grad_y.StudentT1Distrib()] for the derivative in the response,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The location and scale components, written out.
#' r <- y - 0.4; D <- 5 * 1.2^2 + r^2
#' all.equal(g$mu, 6 * r / D)
#' all.equal(g$sigma, 5 * (r^2 - 1.2^2) / (1.2 * D))
#'
#' # The location score redescends: it peaks at |r| = sigma * sqrt(nu) = 2.68
#' # and falls back, where a Gaussian's grows without bound.
#' rr <- c(0.5, 1, 2, 4, 8, 16)
#' rbind(residual = rr,
#'       t = 6 * rr / (5 * 1.2^2 + rr^2),
#'       gaussian = rr / 1.2^2)
#'
#' # The summed score vanishes at the maximum likelihood estimate.
#' set.seed(4)
#' z <- distrib_rng(d, 3000, list(mu = 1, sigma = 2, nu = 4))
#' mle <- as.list(coef(fit_distrib(d, z)))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
S7::method(distrib_gradient, StudentT1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  student_t_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Student t Observed Hessian
#' @name distrib_hessian.StudentT1Distrib
#' @description
#' Computes the six distinct second derivatives of the location-scale Student t
#' log-density with respect to \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}, one value
#' per observation, in closed form. With \eqn{r = y - \mu} and
#' \eqn{D = \nu\sigma^2 + r^2},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{(\nu+1)\left(r^2 - \nu\sigma^2\right)}{D^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\nu\left\{\nu\sigma^4 - (3\nu+1)\sigma^2 r^2 - r^4\right\}}{\sigma^2 D^2},}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4}\left[-\psi_1\!\left(\dfrac{\nu}{2}\right) + \psi_1\!\left(\dfrac{\nu+1}{2}\right) + \dfrac{2\left(\nu\sigma^4 + r^4\right)}{\nu D^2}\right],}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma} = -\dfrac{2\nu(\nu+1)\sigma r}{D^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \nu} = \dfrac{r\left(r^2 - \sigma^2\right)}{D^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma \, \partial \nu} = \dfrac{r^2\left(r^2 - \sigma^2\right)}{\sigma D^2}.}
#'
#' The curvature in \eqn{\mu} **turns positive** wherever \eqn{|r| >
#' \sigma\sqrt{\nu}}, the same point at which the score peaks, so the observed
#' information is indefinite in that direction at an outlying observation while
#' its expectation is negative definite everywhere.
#'
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of six numeric vectors, `mu_mu`, `sigma_sigma`,
#'   `nu_nu`, `mu_sigma`, `mu_nu` and `sigma_nu`, each of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`. The six name the
#'   distinct entries of a symmetric \eqn{3 \times 3} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale and \eqn{\nu > 0} the degrees of freedom.
#' \eqn{\psi_1} is the trigamma function, `trigamma()` in R.
#'
#' @seealso [distrib_gradient.StudentT1Distrib()] for the score,
#'   [distrib_expected_hessian.StudentT1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.StudentT1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # The location and mixed location-scale components, written out.
#' r <- y - 0.4; D <- 5 * 1.2^2 + r^2
#' all.equal(h$mu_mu, 6 * (r^2 - 5 * 1.2^2) / D^2)
#' all.equal(h$mu_sigma, -2 * 5 * 6 * 1.2 * r / D^2)
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-5
#' up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 5))$mu
#' dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 5))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#'
#' # The curvature in mu is positive beyond |r| = sigma * sqrt(nu) = 2.68.
#' rr <- c(0.5, 1, 2, 4, 8, 16)
#' rbind(residual = rr,
#'       mu_mu = 6 * (rr^2 - 5 * 1.2^2) / (5 * 1.2^2 + rr^2)^2)
S7::method(distrib_hessian, StudentT1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  student_t_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Student t Expected Hessian
#' @name distrib_expected_hessian.StudentT1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2\nu}{\sigma^2(\nu+3)},}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \nu^2}\right] = \dfrac{1}{4}\left[\psi_1\!\left(\dfrac{\nu+1}{2}\right) - \psi_1\!\left(\dfrac{\nu}{2}\right)\right] + \dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma \, \partial \nu}\right] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}.}
#'
#' The location is **orthogonal** to both other parameters: the two mixed
#' entries containing \eqn{\mu} are exactly zero, the law being symmetric about
#' \eqn{\mu} while the other two components are even in \eqn{r}. So
#' \eqn{\hat\mu} is asymptotically independent of \eqn{\hat\sigma} and
#' \eqn{\hat\nu}, which the scale and the degrees of freedom are not of each
#' other.
#'
#' Because a closed form exists, `approx` and `nsim` are ignored: every
#' strategy returns the same six numbers. The arithmetic runs in a compiled
#' kernel, so the result does not depend on the thread count.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations. Only its length is read, the
#'   expectation not depending on the data; the values are ignored.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being exact. Accepted so that
#'   the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of six numeric vectors, `mu_mu`, `sigma_sigma`,
#'   `nu_nu`, `mu_sigma`, `mu_nu` and `sigma_nu`, each of length `length(y)`
#'   and each constant along it. `mu_sigma` and `mu_nu` are exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the degrees of freedom and
#' \eqn{\psi_1} the trigamma function.
#'
#' @section Large degrees of freedom:
#' The \eqn{\nu\nu} entry is a difference of two trigammas half a unit apart
#' plus a term that cancels it to leading order, so writing it directly loses
#' every digit at a large \eqn{\nu}: before the asymptotic branches were added
#' it changed sign at \eqn{\nu \approx 3\times10^5} and read \eqn{-n/2} on the
#' link scale. The shipped kernel switches to a series at measured crossovers
#' and stays correctly signed to `.Machine$double.xmax`.
#'
#' @seealso [distrib_hessian.StudentT1Distrib()] for the quantity this is the
#'   expectation of, [distrib_gradient.StudentT1Distrib()] for the score, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#' eh <- distrib_expected_hessian(d, y, th)
#' vapply(eh, function(z) z[1], numeric(1))
#'
#' # The four non-zero closed forms, written out.
#' c(mu_mu = -6 / (1.2^2 * 8),
#'   sigma_sigma = -10 / (1.2^2 * 8),
#'   nu_nu = (trigamma(3) - trigamma(2.5)) / 4 + 10 / (2 * 5 * 6 * 8),
#'   sigma_nu = 2 / (1.2 * 6 * 8))
#'
#' # Averaging the observed Hessian over draws reaches the same six numbers,
#' # the two containing mu going to zero.
#' set.seed(1)
#' z <- distrib_rng(d, 4e5, th)
#' vapply(distrib_hessian(d, z, th), mean, numeric(1))
#'
#' # The strategy argument is inert, the expectation being exact.
#' identical(eh, distrib_expected_hessian(d, y, th, approx = "mc", nsim = 50))
S7::method(distrib_expected_hessian, StudentT1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  student_t_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Student t Third-Order Derivatives
#' @name distrib_deriv3.StudentT1Distrib
#' @description
#' Computes the ten distinct third derivatives of the location-scale Student t
#' log-density in \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}. The observed values
#' are closed form and run in a compiled kernel decomposed over the elements of
#' the output, so they do not depend on the thread count.
#'
#' **The expected values have no closed form.** With `expected = TRUE` the
#' method calls [expected_derivative()], which integrates the observed
#' derivatives against the density by the strategy `approx` names. That is the
#' one place on this page where `approx` and `nsim` are read; on the observed
#' branch both are ignored.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is read.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`, the strategy [expected_derivative()] uses. Read only when
#'   `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Read only on the observed branch. Defaults to `1L`.
#'
#' @return A named list of ten numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_mu_nu`, `mu_sigma_sigma`, `mu_sigma_nu`, `mu_nu_nu`,
#'   `sigma_sigma_sigma`, `sigma_sigma_nu`, `sigma_nu_nu` and `nu_nu_nu`, each
#'   of length `max(length(y), length(mu), length(sigma), length(nu))`.
#'
#' @section Large degrees of freedom:
#' Every component is divided by \eqn{D^3} with \eqn{D = \nu\sigma^2 + r^2},
#' and \eqn{D^3} overflows at \eqn{5.6\times10^{102}} where the log link
#' reaches \eqn{1.8\times10^{308}}. The shipped kernel is written in
#' \eqn{z = r/\sigma}, \eqn{u = z^2/\nu} and \eqn{t = 1/(1+u)} instead, so all
#' ten stay finite to `.Machine$double.xmax`. **On the link scale they do
#' not**: the chain rule forms \eqn{(h')^k} against a component of order
#' \eqn{\nu^{-k}}, and one of the ten ceases to be finite at
#' \eqn{\nu = 10^{150}}. That regime is where the family is a Gaussian in all
#' but name.
#'
#' @seealso [distrib_hessian.StudentT1Distrib()] for the order below,
#'   [distrib_deriv4.StudentT1Distrib()] for the order above,
#'   [expected_derivative()] for the numerical expectation, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-location
#' # component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 5))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 5))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#'
#' # The expected branch is a quadrature, and averaging the observed one over
#' # draws reaches it; the components odd in the residual go to zero.
#' set.seed(2)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(expected = vapply(distrib_deriv3(d, y, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       averaged = vapply(distrib_deriv3(d, z, th), mean, numeric(1)))
S7::method(distrib_deriv3, StudentT1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    student_t_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

#' @title Student t Fourth-Order Derivatives
#' @name distrib_deriv4.StudentT1Distrib
#' @description
#' Computes the fifteen distinct fourth derivatives of the location-scale
#' Student t log-density in \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}. The observed
#' values are closed form and run in a compiled kernel decomposed over the
#' elements of the output, so they do not depend on the thread count.
#'
#' **The expected values have no closed form.** With `expected = TRUE` the
#' method calls [expected_derivative()], which integrates the observed
#' derivatives against the density by the strategy `approx` names. That is the
#' one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is read.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`, the strategy [expected_derivative()] uses. Read only when
#'   `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Read only on the observed branch. Defaults to `1L`.
#'
#' @return A named list of fifteen numeric vectors named for the multi-index
#'   they carry, from `mu_mu_mu_mu` to `nu_nu_nu_nu`, each of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`.
#'
#' @section Large degrees of freedom:
#' Unlike the third order, the fourth is **not** rewritten in the ratio
#' variables and ceases to be finite at a large \eqn{\nu}: measured at
#' \eqn{\sigma = 1.2}, eight of the fifteen components are finite at
#' \eqn{\nu = 10^{150}}, five at \eqn{10^{300}} and two at
#' `.Machine$double.xmax`. A `NaN` there is a loud failure and is preferable to
#' a plausible wrong number, and the regime is one in which the family is a
#' Gaussian in all but name. An outer criterion that reads this order at a
#' \eqn{\nu} run to its clamp is reported as having no finite gradient rather
#' than being given one.
#'
#' @seealso [distrib_deriv3.StudentT1Distrib()] for the order below,
#'   [distrib_hessian.StudentT1Distrib()] for the second order,
#'   [expected_derivative()] for the numerical expectation, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#' d4 <- distrib_deriv4(d, y, th)
#' length(d4)
#' names(d4)[1:4]
#'
#' # A central difference of the third order reproduces the pure-location
#' # component.
#' eps <- 1e-4
#' up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 5))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 5))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#'
#' # At a degrees of freedom the log link can produce, part of the order is
#' # not representable and says so.
#' big <- distrib_deriv4(d, y, list(mu = 0.4, sigma = 1.2, nu = 1e300))
#' sum(vapply(big, function(v) is.finite(v[1]), logical(1)))
S7::method(distrib_deriv4, StudentT1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    student_t_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

#' @title Student t First Derivative in the Response
#' @name distrib_grad_y.StudentT1Distrib
#' @description
#' Computes \eqn{\partial \ell / \partial y}, the derivative of the Student t
#' log-density with respect to the response, in closed form. With
#' \eqn{r = y - \mu} and \eqn{D = \nu\sigma^2 + r^2},
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{(\nu+1)r}{D}.}
#' The family is a location family in \eqn{\mu}, so this is exactly the
#' negative of the location score
#' [distrib_gradient.StudentT1Distrib()]`$mu`, and it redescends with it: it is
#' largest in size at \eqn{|r| = \sigma\sqrt{\nu}} and falls back towards zero
#' beyond that. This quantity is what a quantile residual's delta-method
#' standard error and a change of variable in the response both need.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`, one value per
#'   observation.
#'
#' @seealso [distrib_hess_y.StudentT1Distrib()] for the second derivative in
#'   the response, [distrib_cross_y.StudentT1Distrib()] for the mixed
#'   derivative in the response and the parameters,
#'   [distrib_gradient.StudentT1Distrib()] for the derivatives in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#'
#' # The closed form, written out.
#' r <- y - 0.4
#' all.equal(distrib_grad_y(d, y, th), -6 * r / (5 * 1.2^2 + r^2))
#'
#' # A location family, so this is minus the location score.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # It is the derivative of the log-density, so a central difference of the
#' # log-density in y reproduces it.
#' eps <- 1e-6
#' all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
#'            distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
#'           distrib_grad_y(d, y, th), tolerance = 1e-6)
S7::method(distrib_grad_y, StudentT1Distrib) <- function(distrib, y, theta, ...) {
  r <- y - theta[[1]]
  nu <- theta[[3]]
  -(nu + 1) * r / (nu * theta[[2]]^2 + r^2)
}

#' @title Student t Second Derivative in the Response
#' @name distrib_hess_y.StudentT1Distrib
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2}, the second derivative of the
#' Student t log-density with respect to the response, in closed form. With
#' \eqn{r = y - \mu} and \eqn{D = \nu\sigma^2 + r^2},
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2}
#'   = \dfrac{(\nu+1)\left(r^2 - \nu\sigma^2\right)}{D^2}.}
#' The family is a location family, so this equals the pure-location entry of
#' [distrib_hessian.StudentT1Distrib()], with no sign change; two derivatives
#' in \eqn{\mu} carry two factors of \eqn{-1}. It is negative near the mode and
#' **positive** beyond \eqn{|r| = \sigma\sqrt{\nu}}, so the log-density is not
#' concave in the response, unlike a Gaussian's.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`, one value per
#'   observation.
#'
#' @seealso [distrib_grad_y.StudentT1Distrib()] for the first derivative in the
#'   response, [distrib_hessian.StudentT1Distrib()] for the second derivatives
#'   in the parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- student_t1_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#'
#' # The closed form, written out.
#' r <- y - 0.4; D <- 5 * 1.2^2 + r^2
#' all.equal(distrib_hess_y(d, y, th), 6 * (r^2 - 5 * 1.2^2) / D^2)
#'
#' # A location family, so this is the pure-location entry of the Hessian.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # A central difference of the first derivative reproduces it.
#' eps <- 1e-5
#' all.equal((distrib_grad_y(d, y + eps, th) -
#'            distrib_grad_y(d, y - eps, th)) / (2 * eps),
#'           distrib_hess_y(d, y, th), tolerance = 1e-6)
#'
#' # Positive beyond |r| = sigma * sqrt(nu) = 2.68, so the log-density is
#' # convex in the response out there.
#' distrib_hess_y(d, 0.4 + c(1, 2, 4, 8), th)
S7::method(distrib_hess_y, StudentT1Distrib) <- function(distrib, y, theta, ...) {
  r <- y - theta[[1]]
  nu <- theta[[3]]
  vs2 <- nu * theta[[2]]^2
  (nu + 1) * (r^2 - vs2) / (vs2 + r^2)^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Student t Distribution, Location, Scale and Degrees of Freedom
#'
#' @description
#' Builds the distribution object for the location-scale Student t family,
#' parametrized by a location \eqn{\mu}, a scale \eqn{\sigma > 0} and degrees
#' of freedom \eqn{\nu > 0}. The returned object carries closed-form
#' derivatives of the log-density to fourth order in the parameters, closed
#' first and second derivatives in the response, and a closed expected Hessian;
#' the expected third and fourth orders are the only quantities that go through
#' a numerical route.
#'
#' The family is the standard heavy-tailed alternative to a Gaussian. Its
#' location score redescends, so a gross outlier contributes almost nothing to
#' the estimating equation, and \eqn{\nu} is estimated from the data rather
#' than set: a large fitted \eqn{\nu} reports that the sample looks Gaussian.
#'
#' The three arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the location
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the location
#'   ranging over the whole line already.
#' @param link_sigma A `link` object from `linkfunctions7` for the scale
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#' @param link_nu A `link` object from `linkfunctions7` for the degrees of
#'   freedom \eqn{\nu}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\left(\dfrac{\nu}{2}\right)} \left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}},}
#' the distribution function \eqn{F(q) = T_\nu\{(q-\mu)/\sigma\}} and the
#' quantile function \eqn{Q(p) = \mu + \sigma T_\nu^{-1}(p)}, with \eqn{T_\nu}
#' the standard Student t distribution function.
#'
#' \eqn{\sigma} is the **scale** and not the standard deviation, which is
#' \eqn{\sigma\sqrt{\nu/(\nu-2)}} and exists only above \eqn{\nu = 2}. Keeping
#' the two apart lets the family be fitted where the second moment does not
#' exist; [student_t2_distrib()] is the parametrization by the standard
#' deviation, for a reader who wants one.
#'
#' Two shapes are named families: \eqn{\nu = 1} is the Cauchy, and
#' \eqn{\nu \to \infty} the Gaussian with standard deviation \eqn{\sigma}. The
#' approach to the limit is \eqn{O(1/\nu)}, so at \eqn{\nu = 10^5} the density
#' already agrees with a Gaussian's to seven figures.
#'
#' # Derivatives
#'
#' Write \eqn{r = y - \mu} and \eqn{D = \nu\sigma^2 + r^2}. The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{(\nu+1)r}{D}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{\nu\left(r^2 - \sigma^2\right)}{\sigma D},}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left[
#'   -\dfrac{1}{\nu} - \psi\!\left(\dfrac{\nu}{2}\right)
#'   + \psi\!\left(\dfrac{\nu+1}{2}\right) + \dfrac{(\nu+1)r^2}{\nu D}
#'   - \log\!\left(1 + \dfrac{r^2}{\nu\sigma^2}\right)\right],}
#' with \eqn{\psi} the digamma function. The observed Hessian is in
#' [distrib_hessian.StudentT1Distrib()], and its expectation is
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2\nu}{\sigma^2(\nu+3)},}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \nu^2}\right] = \dfrac{1}{4}\left[\psi_1\!\left(\dfrac{\nu+1}{2}\right) - \psi_1\!\left(\dfrac{\nu}{2}\right)\right] + \dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma \, \partial \nu}\right] = \dfrac{2}{\sigma(\nu+1)(\nu+3)},}
#' with the two entries containing \eqn{\mu} exactly zero. The location is
#' therefore orthogonal to the other two, whose estimates are asymptotically
#' correlated with each other.
#'
#' # Redescent, and where it costs
#'
#' The location score \eqn{(\nu+1)r/D} rises to
#' \eqn{(\nu+1)/(2\sigma\sqrt{\nu})} at \eqn{|r| = \sigma\sqrt{\nu}} and falls
#' back towards zero, so an observation far from the location is downweighted
#' automatically. The price is that the curvature in \eqn{\mu} turns
#' **positive** past the same point, so the observed information is indefinite
#' at an outlying observation. Fisher scoring is unaffected, the expected
#' information being negative definite everywhere, which is why it is the
#' default in [fit_distrib()].
#'
#' # Moments
#'
#' Each exists only above its own threshold, and the family reports `NaN` or
#' `Inf` where it does not: the mean is \eqn{\mu} for \eqn{\nu > 1}, the
#' variance \eqn{\sigma^2\nu/(\nu-2)} for \eqn{\nu > 2}, the skewness 0 for
#' \eqn{\nu > 3} and the excess kurtosis \eqn{6/(\nu-4)} for \eqn{\nu > 4}. At
#' \eqn{\nu = 1.5}, for instance, the mean is \eqn{\mu} and the variance is
#' `Inf`.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale by Fisher
#' scoring. No estimate is closed form, the three estimating equations being
#' coupled through \eqn{D}. On Gaussian data \eqn{\hat\nu} runs towards its
#' upper boundary, which is the correct answer and is reported as a value near
#' the clamp of its log link rather than as `Inf`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale and \eqn{\nu > 0} the degrees of freedom.
#' \eqn{\psi} and \eqn{\psi_1} are the digamma and trigamma functions,
#' `digamma()` and `trigamma()` in R. \eqn{\eta} is a parameter on the
#' unconstrained scale of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `StudentT1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"student t1"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params`
#'   `c("mu", "sigma", "nu")`, `n_params` `3`, `params_bounds` the domains
#'   \eqn{(-\infty, \infty)}, \eqn{(0, \infty)} and \eqn{(0, \infty)}, and
#'   `link_params` the three links given here.
#'
#' @references
#' Lange, K. L., Little, R. J. A. and Taylor, J. M. G. (1989).
#' Robust statistical modeling using the t distribution.
#' *Journal of the American Statistical Association*, **84**(408), 881-896.
#'
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995).
#' *Continuous Univariate Distributions*, Volume 2, 2nd edition, Chapter 28.
#' Wiley, New York.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dt pt qt rt
#'
#' @examples
#' d <- student_t1_distrib()
#' d
#'
#' # The density is stats::dt at the standardized value, over sigma.
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 5)
#' all.equal(distrib_pdf(d, y, th), dt((y - 0.4) / 1.2, df = 5) / 1.2)
#'
#' # The scale is not the standard deviation; the moments carry thresholds.
#' c(scale = th$sigma, mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # Below its own threshold each moment reports that it does not exist.
#' t(vapply(c(0.8, 1.5, 2.5, 3.5),
#'          function(v) {
#'            p <- list(mu = 0.4, sigma = 1.2, nu = v)
#'            c(nu = v, mean = mean(d, p), var = variance(d, p),
#'              kurt = kurtosis(d, p))
#'          }, numeric(4)))
#'
#' # One degree of freedom is the Cauchy; a large one is the Gaussian.
#' all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 1)),
#'           dcauchy(y, 0.4, 1.2))
#' max(abs(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 1e5)) -
#'         dnorm(y, 0.4, 1.2)))
#'
#' # Fitting recovers all three parameters.
#' set.seed(4)
#' z <- distrib_rng(d, 3000, list(mu = 1, sigma = 2, nu = 4))
#' coef(fit_distrib(d, z))
#'
#' @seealso
#' [student_t2_distrib()] for the same law in the standard deviation;
#' [cauchy_distrib()] for \eqn{\nu = 1} and [gaussian1_distrib()] for the
#' limit; [skewt_distrib()] to add a shape parameter;
#' [pseudohuber_distrib()] and [laplace_distrib()] for other robust families;
#' [mvstudent_t_distrib()] for the multivariate version;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [StudentT1Distrib] for the class.
#' @export
student_t1_distrib <- function(link_mu = identity_link(), link_sigma = log_link(), link_nu = log_link()) {
  
  StudentT1Distrib(
    distrib_name = "student t1", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma", "nu"),
    params_interpretation = c(mu = "location", sigma = "scale", nu = "shape"),
    n_params = 3,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf), nu = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma, nu = link_nu)
  )
  
}
