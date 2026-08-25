#' @include distrib.R generics.R
NULL

#' @title von Mises Distribution Class
#' @name VonMises1Distrib
#'
#' @description
#' The S7 class of the von Mises family, the natural distribution for an angle,
#' parametrized by a mean direction \eqn{\mu} and a concentration
#' \eqn{\kappa > 0}, with density
#' \eqn{f(y) = e^{\kappa\cos(y-\mu)}/\{2\pi I_0(\kappa)\}} on
#' \eqn{[-\pi, \pi)}. It inherits from `continuous_distrib`; the nine methods
#' listed below are registered on it directly.
#'
#' This is the first family in the package whose support has the topology of a
#' circle: the two ends of the interval are the same point, so the density need
#' not vanish at either.
#'
#' Build one with [vonmises1_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `VonMises1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [vonmises1_distrib()] they hold `"von mises1"`,
#'   `"univariate"`, `c(-pi, pi)`, `c("mu", "kappa")`, the interpretations
#'   `c(mu = "mean direction", kappa = "concentration")`, `2`, and the domains
#'   \eqn{(-\pi, \pi)} and \eqn{(0, \infty)}.
#'
#' @seealso [vonmises1_distrib()] to build one;
#'   [vonmises2_distrib()] for the same law parametrized by the mean resultant
#'   length; [gaussian1_distrib()] for the analogous family on the line;
#'   [distrib_cdf.VonMises1Distrib()] for the Bessel series.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.VonMises1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.VonMises1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.VonMises1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.VonMises1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.VonMises1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.VonMises1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.VonMises1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.VonMises1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.VonMises1Distrib],
#'   [`distrib_rng()`][distrib_rng.VonMises1Distrib].
#'
#' The **quantile** comes from [continuous_distrib()], by root finding on the
#' distribution function; the distribution function itself is this class's own
#' and is a Bessel series, not the parent's quadrature. The four moments come
#' from [mean.distrib()] and its siblings, numerically, and are the ordinary
#' moments of \eqn{Y} as a number rather than the circular ones.
#'
#' @examples
#' d <- vonmises1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The support is the circle, written as a half-open interval.
#' d@bounds
#' d@params
#' d@params_interpretation
#'
#' # The direction rides a bounded link and the concentration a log.
#' vapply(d@link_params, function(l) l@link_name, character(1))
VonMises1Distrib <- S7::new_class("VonMises1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title von Mises Density
#' @name distrib_pdf.VonMises1Distrib
#' @description
#' Computes the von Mises density
#' \deqn{f(y; \mu, \kappa) = \dfrac{e^{\kappa \cos(y - \mu)}}{2\pi I_0(\kappa)},
#'       \qquad y \in [-\pi, \pi),}
#' with \eqn{I_0} the modified Bessel function of the first kind. Outside
#' \eqn{[-\pi, \pi)} the density is 0: the support is the declared interval,
#' and an angle is not wrapped into it.
#'
#' The normalizing constant goes through [numericals7::log_bessel_i()] and not
#' through [base::besselI()]. R's exponentially scaled `besselI` **underflows
#' to an exact zero** between \eqn{\kappa = 10^5} and \eqn{10^6}, where the
#' logarithm then returns `-Inf`; `log_bessel_i` stays finite wherever the
#' logarithm itself is representable. Measured at \eqn{\kappa = 10^6} it
#' returns 999992.17 where the base route returns `-Inf`.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles. A value outside \eqn{[-\pi, \pi)} is
#'   off the support and gives a density of 0, or `-Inf` with `log = TRUE`.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(-\pi, \pi)} and `kappa` be strictly
#'   positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads
#'   [numericals7::log_bessel_i()] may use. It is carried down because this is
#'   where the family spends its time: profiled at 80.8 per cent of a fit whose
#'   concentration is modelled, so that `kappa` is a vector. Defaults to `1L`.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(kappa))`, one value per observation.
#'
#' @section Notation:
#' \eqn{\mu} is the mean direction, \eqn{\kappa > 0} the concentration and
#' \eqn{I_0} the modified Bessel function of the first kind of order zero,
#' `besselI(x, 0)` in R.
#'
#' @seealso [distrib_cdf.VonMises1Distrib()] for the distribution function,
#'   [distrib_gradient.VonMises1Distrib()] for the derivatives of the
#'   log-density, [numericals7::log_bessel_i()] for the constant, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' th <- list(mu = 0.5, kappa = 2)
#' y <- c(-1, 0, 0.5, 2)
#'
#' # The formula written out.
#' all.equal(distrib_pdf(d, y, th),
#'           exp(2 * cos(y - 0.5)) / (2 * pi * besselI(2, 0)))
#'
#' # It integrates to one over the circle.
#' integrate(function(v) distrib_pdf(d, v, th), -pi, pi)$value
#'
#' # Outside the declared interval the density is zero; an angle is not
#' # wrapped into it.
#' distrib_pdf(d, c(-4, 3.5), th)
#'
#' # The constant survives a concentration at which the base route does not.
#' c(ours = numericals7::log_bessel_i(1e6, 0),
#'   base = log(besselI(1e6, 0, expon.scaled = TRUE)) + 1e6)
S7::method(distrib_pdf, VonMises1Distrib) <- function(distrib, y, theta,
                                                     log = FALSE, ...,
                                                     threads = 1L) {
  k <- theta[[2]]
  # log I_0 from numericals7: the scaled besselI underflows to an exact zero
  # between kappa = 1e5 and 1e6, where log() returns -Inf; log_bessel_i is
  # finite wherever the logarithm itself is representable. It is also where
  # this family spends its time -- profiled at 80.8 per cent of a fit whose
  # concentration is modelled, so the concentration is a vector -- which is
  # why the count is carried down to it.
  log_i0 <- numericals7::log_bessel_i(k, 0, threads)
  out <- k * cos(y - theta[[1]]) - log(2 * pi) - log_i0
  out[y < -pi | y >= pi] <- -Inf
  if (log) out else exp(out)
}

#' @title von Mises Random Generation
#' @name distrib_rng.VonMises1Distrib
#' @description
#' Draws `n` independent angles by the rejection algorithm of Best and Fisher
#' (1979), which proposes from a wrapped Cauchy envelope and accepts with a
#' probability that **involves no Bessel function at all**, so the generator is
#' cheap at any concentration: the normalizing constant, the expensive part of
#' the density, never enters.
#'
#' The accepted angles are drawn about zero and then shifted by \eqn{\mu} and
#' wrapped back into \eqn{[-\pi, \pi)}, so every draw lies in the declared
#' support. The loop over-proposes and repeats until `n` draws have been
#' accepted, so it consumes an unpredictable number of R's uniform streams.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1. `mu` must lie in \eqn{(-\pi, \pi)} and `kappa` be
#'   strictly positive. The envelope's constants are built once per call, so a
#'   parameter varying by observation is not supported here.
#'
#' @return A numeric vector of `n` angles in \eqn{[-\pi, \pi)}.
#'
#' @references
#' Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von Mises
#' distribution. *Journal of the Royal Statistical Society, Series C*,
#' **28**(2), 152-157.
#'
#' @seealso [distrib_pdf.VonMises1Distrib()] for the density,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' th <- list(mu = 0.5, kappa = 2)
#' set.seed(1)
#' z <- distrib_rng(d, 3e5, th)
#'
#' # Every draw is in the declared support.
#' range(z)
#'
#' # The circular mean recovers mu and the mean resultant length recovers
#' # A(kappa) = I_1 / I_0; the ordinary mean recovers neither.
#' c(circular_mean = atan2(mean(sin(z)), mean(cos(z))), mu = 0.5)
#' c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2),
#'   A = numericals7::bessel_i_ratio(2))
S7::method(distrib_rng, VonMises1Distrib) <- function(distrib, n, theta) {
  mu <- theta[[1]]
  k <- theta[[2]]
  a <- 1 + sqrt(1 + 4 * k * k)
  b <- (a - sqrt(2 * a)) / (2 * k)
  r <- (1 + b * b) / (2 * b)

  out <- numeric(0)
  while (length(out) < n) {
    m <- 2 * (n - length(out)) + 16
    u1 <- stats::runif(m); u2 <- stats::runif(m); u3 <- stats::runif(m)
    z <- cos(pi * u1)
    f <- (1 + r * z) / (r + z)
    c0 <- k * (r - f)
    ok <- (c0 * (2 - c0) - u2 > 0) | (log(c0 / u2) + 1 - c0 >= 0)
    ang <- sign(u3 - 0.5) * acos(pmin(pmax(f, -1), 1))
    out <- c(out, ang[ok])
  }
  # wrapped back into the declared support
  ((out[seq_len(n)] + mu + pi) %% (2 * pi)) - pi
}

#' @title von Mises Score
#' @name distrib_gradient.VonMises1Distrib
#' @description
#' Computes the first derivatives of the von Mises log-density with respect to
#' the mean direction \eqn{\mu} and the concentration \eqn{\kappa}, one value
#' per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \kappa \sin(y - \mu), \qquad
#'       \dfrac{\partial \ell}{\partial \kappa} = \cos(y - \mu) - A(\kappa),}
#' with \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} the derivative of
#' \eqn{\log I_0} and also the **mean resultant length** of the family.
#'
#' The ratio comes from [numericals7::bessel_i_ratio()], which switches to an
#' asymptotic expansion past \eqn{\kappa = 10^4}. R's own scaled `besselI`
#' underflows to an exact zero between \eqn{10^5} and \eqn{10^6}, so forming
#' the ratio from two calls gives `NaN` over part of that band.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries. This method always returns the parameter scale.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(-\pi, \pi)} and `kappa` be strictly
#'   positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `kappa`, each of
#'   length `max(length(y), length(mu), length(kappa))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\kappa > 0} the concentration, \eqn{I_m} the modified Bessel
#' function of the first kind of order \eqn{m}, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_hessian.VonMises1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.VonMises1Distrib()] for their expectation,
#'   [numericals7::bessel_i_ratio()] for \eqn{A}, and [distrib_gradient()] for
#'   the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, kappa = 2)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' all.equal(g$mu, 2 * sin(y - 0.5))
#' all.equal(g$kappa, cos(y - 0.5) - besselI(2, 1) / besselI(2, 0))
#'
#' # numDeriv on the summed log-density reproduces the summed score.
#' fn <- function(p)
#'   sum(distrib_pdf(d, y, list(mu = p[1], kappa = p[2]), log = TRUE))
#' rbind(numeric = numDeriv::grad(fn, c(0.5, 2)),
#'       closed = vapply(g, sum, numeric(1)))
#'
#' # The summed score vanishes at the maximum likelihood estimate.
#' set.seed(7)
#' z <- distrib_rng(d, 2000, list(mu = 0.8, kappa = 3))
#' mle <- as.list(coef(fit_distrib(d, z)))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
S7::method(distrib_gradient, VonMises1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  d <- y - theta[[1]]
  list(mu = theta[[2]] * sin(d), kappa = cos(d) - numericals7::bessel_i_ratio(theta[[2]]))
}

#' @title von Mises Observed Hessian
#' @name distrib_hessian.VonMises1Distrib
#' @description
#' Computes the three distinct second derivatives of the von Mises log-density,
#' one value per observation, in closed form:
#' \deqn{\ell^{(\mu\mu)} = -\kappa\cos(y-\mu), \qquad
#'       \ell^{(\mu\kappa)} = \sin(y-\mu), \qquad
#'       \ell^{(\kappa\kappa)} = -A'(\kappa),}
#' with \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}. The
#' pure-concentration entry is **free of the data**: the log-density is
#' \eqn{\kappa\cos(y-\mu) - \log\{2\pi I_0(\kappa)\}}, linear in \eqn{\kappa}
#' apart from the normalizing constant, so \eqn{\kappa} appears twice only
#' inside \eqn{\log I_0}. It therefore equals its own expectation at every
#' observation.
#'
#' \eqn{A'} comes from the Riccati recurrence
#' \eqn{A' = 1 - A/\kappa - A^2}, which [numericals7::bessel_i_ratio_derivs()]
#' runs, so no second Bessel evaluation is needed. It is the variance of
#' \eqn{\cos(Y-\mu)} and is positive, so the information is positive
#' definite.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(-\pi, \pi)} and `kappa` be strictly
#'   positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_kappa` and
#'   `kappa_kappa`, each of length
#'   `max(length(y), length(mu), length(kappa))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\kappa > 0} the concentration, \eqn{I_m} the modified Bessel
#' function of the first kind of order \eqn{m}, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_gradient.VonMises1Distrib()] for the score,
#'   [distrib_expected_hessian.VonMises1Distrib()] for the expectation of this
#'   quantity, [numericals7::bessel_i_ratio_derivs()] for \eqn{A'}, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, kappa = 2)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # The pure-concentration entry is one number, repeated.
#' unique(h$kappa_kappa)
#'
#' # And it is minus A', which the Riccati recurrence gives.
#' A <- numericals7::bessel_i_ratio(2)
#' c(riccati = 1 - A / 2 - A^2,
#'   supplied = numericals7::bessel_i_ratio_derivs(2)$d1)
#'
#' # numDeriv on the summed log-density reproduces the summed matrix.
#' fn <- function(p)
#'   sum(distrib_pdf(d, y, list(mu = p[1], kappa = p[2]), log = TRUE))
#' H <- numDeriv::hessian(fn, c(0.5, 2))
#' rbind(numeric = c(H[1, 1], H[2, 2], H[1, 2]),
#'       closed = c(sum(h$mu_mu), sum(h$kappa_kappa), sum(h$mu_kappa)))
S7::method(distrib_hessian, VonMises1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...) {
  d <- y - theta[[1]]
  k <- theta[[2]]
  list(mu_mu = -k * cos(d), mu_kappa = sin(d),
       kappa_kappa = rep_len(-numericals7::bessel_i_ratio_derivs(k)$d1, length(d)))
}

#' @title von Mises Expected Hessian
#' @name distrib_expected_hessian.VonMises1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. Under the model
#' \eqn{\mathbb{E}[\cos(Y-\mu)] = A(\kappa)} and
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0} by symmetry, so
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\kappa A(\kappa), \qquad
#'       \mathbb{E}[\ell^{(\mu\kappa)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\kappa\kappa)}] = -A'(\kappa),}
#' with \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' The zero off-diagonal makes the direction and the concentration
#' **orthogonal**: the expected information is diagonal, their estimates are
#' asymptotically independent, and Fisher scoring updates the two
#' independently. `approx` and `nsim` are ignored, the expectation being exact.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles. Only its length is read, the
#'   expectation not depending on the data; the values are ignored.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1. `kappa` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being exact. Accepted so that
#'   the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_kappa` and
#'   `kappa_kappa`, each of length `length(y)` and constant along it.
#'   `mu_kappa` is exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\kappa > 0} the concentration, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} the mean resultant length.
#'
#' @seealso [distrib_hessian.VonMises1Distrib()] for the quantity this is the
#'   expectation of, [distrib_gradient.VonMises1Distrib()] for the score, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' th <- list(mu = 0.5, kappa = 2)
#' eh <- distrib_expected_hessian(d, c(-1, 0, 0.5, 2), th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # The closed forms, written out; the off-diagonal is exactly zero, so the
#' # two parameters are orthogonal.
#' A <- numericals7::bessel_i_ratio(2)
#' c(mu_mu = -2 * A, mu_kappa = 0,
#'   kappa_kappa = -numericals7::bessel_i_ratio_derivs(2)$d1)
#'
#' # Averaging the observed Hessian over draws reaches the same three numbers.
#' set.seed(1)
#' z <- distrib_rng(d, 3e5, th)
#' vapply(distrib_hessian(d, z, th), mean, numeric(1))
#'
#' # The strategy argument is inert, the expectation being exact.
#' identical(eh, distrib_expected_hessian(d, c(-1, 0, 0.5, 2), th,
#'                                        approx = "mc", nsim = 50))
S7::method(distrib_expected_hessian, VonMises1Distrib) <- function(distrib, y, theta,
                                                                   scale = c("parameter", "link"),
                                                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                                                   nsim = 10000, ...) {
  k <- theta[[2]]
  a <- numericals7::bessel_i_ratio_derivs(k)
  n <- length(y)
  list(mu_mu = rep_len(-k * a$A, n), mu_kappa = rep_len(0, n),
       kappa_kappa = rep_len(-a$d1, n))
}

#' @title von Mises First Derivative in the Response
#' @name distrib_grad_y.VonMises1Distrib
#' @description
#' Computes \eqn{\partial \ell / \partial y}, the derivative of the von Mises
#' log-density with respect to the angle, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\kappa \sin(y - \mu).}
#' The family is a location family on the circle, so this is exactly the
#' negative of the direction score
#' [distrib_gradient.VonMises1Distrib()]`$mu`. It is **bounded by** \eqn{\kappa}
#' and vanishes at both the mode \eqn{y = \mu} and the antimode
#' \eqn{y = \mu \pm \pi}, as a periodic density must.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles. The expression is evaluated wherever it
#'   is given, including outside \eqn{[-\pi, \pi)}, where the density itself is
#'   zero.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `kappa` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(kappa))`, one value per observation.
#'
#' @seealso [distrib_hess_y.VonMises1Distrib()] for the second derivative in
#'   the angle, [distrib_gradient.VonMises1Distrib()] for the derivatives in
#'   the parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, kappa = 2)
#'
#' # The closed form, written out.
#' all.equal(distrib_grad_y(d, y, th), -2 * sin(y - 0.5))
#'
#' # A location family on the circle, so this is minus the direction score.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # A central difference of the log-density reproduces it.
#' eps <- 1e-6
#' all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
#'            distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
#'           distrib_grad_y(d, y, th), tolerance = 1e-6)
#'
#' # It vanishes at the mode and at the antimode, and is bounded by kappa.
#' distrib_grad_y(d, c(0.5, 0.5 - pi), th)
#' max(abs(distrib_grad_y(d, seq(-pi, pi, length.out = 401), th)))
S7::method(distrib_grad_y, VonMises1Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * sin(y - theta[[1]])
}

#' @title von Mises Second Derivative in the Response
#' @name distrib_hess_y.VonMises1Distrib
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2} in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\kappa \cos(y - \mu).}
#' The family is a location family on the circle, so this equals the
#' pure-direction entry of [distrib_hessian.VonMises1Distrib()], with no sign
#' change. It is negative near the mode and **positive** on the half of the
#' circle further than a quarter turn from it, the density having a minimum at
#' the antimode; a periodic density cannot be concave everywhere.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles. The expression is evaluated wherever it
#'   is given, including outside \eqn{[-\pi, \pi)}, where the density itself is
#'   zero.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `kappa` must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(kappa))`, one value per observation.
#'
#' @seealso [distrib_grad_y.VonMises1Distrib()] for the first derivative in the
#'   angle, [distrib_hessian.VonMises1Distrib()] for the second derivatives in
#'   the parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, kappa = 2)
#'
#' # The closed form, written out.
#' all.equal(distrib_hess_y(d, y, th), -2 * cos(y - 0.5))
#'
#' # A location family, so this is the pure-direction entry of the Hessian.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # A central difference of the first derivative reproduces it.
#' eps <- 1e-5
#' all.equal((distrib_grad_y(d, y + eps, th) -
#'            distrib_grad_y(d, y - eps, th)) / (2 * eps),
#'           distrib_hess_y(d, y, th), tolerance = 1e-6)
#'
#' # Negative at the mode and positive at the antimode: a periodic density is
#' # not concave everywhere.
#' distrib_hess_y(d, c(0.5, 0.5 - pi), th)
S7::method(distrib_hess_y, VonMises1Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * cos(y - theta[[1]])
}

# --- CONSTRUCTOR WRAPPER ---

#' von Mises Distribution, Mean Direction and Concentration
#'
#' @description
#' Builds the distribution object for the von Mises family, the natural
#' distribution for an angle, parametrized by a mean direction \eqn{\mu} and a
#' concentration \eqn{\kappa > 0}. The returned object carries closed-form
#' derivatives of the log-density to fourth order in the parameters and in the
#' response, a closed-form expected information, and a **distribution function
#' from a Bessel series** in place of the quadrature the base class would use.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean direction
#'   \eqn{\mu}. Defaults to `linkfunctions7::bounded_link(lwr = -pi, upr = pi)`,
#'   which maps the free scale onto \eqn{(-\pi, \pi)}. See the note on the
#'   chart below.
#' @param link_kappa A `link` object from `linkfunctions7` for the
#'   concentration \eqn{\kappa}. Defaults to [linkfunctions7::log_link()],
#'   which maps \eqn{(0, \infty)} onto the line and so keeps every fitted value
#'   positive.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in [-\pi, \pi)} is
#' \deqn{f(y; \mu, \kappa) = \dfrac{e^{\kappa\cos(y-\mu)}}{2\pi I_0(\kappa)},}
#' with \eqn{I_0} the modified Bessel function of the first kind. This is the
#' first family in the package whose support has the topology of a circle: the
#' two ends of the interval are the same point, so the density need not vanish
#' at either. [vonmises2_distrib()] is the same law parametrized by the mean
#' resultant length in place of the concentration.
#'
#' The normalizing constant goes through [numericals7::log_bessel_i()]. R's
#' exponentially scaled `besselI` underflows to an exact zero between
#' \eqn{\kappa = 10^5} and \eqn{10^6}, where the logarithm returns `-Inf`;
#' measured at \eqn{\kappa = 10^6} the toolkit's routine returns 999992.17.
#'
#' # Score, information, and orthogonality
#'
#' Writing \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} for the derivative of
#' \eqn{\log I_0},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
#'       \dfrac{\partial\ell}{\partial\kappa} = \cos(y-\mu) - A(\kappa),}
#' and every second derivative is closed form, with
#' \eqn{A'(\kappa) = 1 - A(\kappa)/\kappa - A(\kappa)^2} from the Bessel
#' recurrences and no further evaluation. Since
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0} by symmetry, **the two parameters are
#' orthogonal**: the expected information is diagonal, and Fisher scoring
#' updates the direction and the concentration independently.
#'
#' \eqn{A(\kappa)} is the mean resultant length, so \eqn{A'(\kappa)} is the
#' variance of \eqn{\cos(Y-\mu)} and is positive, which keeps the information
#' positive definite.
#'
#' The log-density is linear in \eqn{\kappa} apart from the normalizing
#' constant, so at orders three and four every component naming one \eqn{\mu}
#' and two or more \eqn{\kappa} is exactly zero.
#'
#' # The moments a generic returns are not the circular ones
#'
#' [mean.distrib()], [variance()], [skewness()] and [kurtosis()] are the
#' ordinary moments of \eqn{Y} as a number on \eqn{[-\pi, \pi)}, obtained
#' numerically. They are **not** the circular quantities and not what \eqn{\mu}
#' and \eqn{\kappa} describe. \eqn{\mu} is the mean *direction*, and
#' \eqn{\mathbb{E}[Y] \ne \mu} whenever \eqn{\mu \ne 0}, because the interval
#' is cut at \eqn{\pm\pi} rather than at \eqn{\mu \pm \pi} and the density is
#' not symmetric about \eqn{\mu} on it: at \eqn{\mu = 1.2} and \eqn{\kappa = 2}
#' the ordinary mean is 1.079.
#'
#' The circular mean is \eqn{\mu} and the mean resultant length is
#' \eqn{\rho = I_1(\kappa)/I_0(\kappa)}, both closed form. Neither is returned
#' by a generic whose name means something else; compute them from a sample as
#' `atan2(mean(sin(z)), mean(cos(z)))` and
#' `sqrt(mean(cos(z))^2 + mean(sin(z))^2)`.
#'
#' # The direction is carried on a bounded chart
#'
#' The default link maps the free scale onto \eqn{(-\pi, \pi)}, which keeps
#' \eqn{\mu} identified. The cost is that a fit cannot walk across the
#' boundary, so data concentrated near \eqn{\pm\pi} are better rotated before
#' fitting than handed to the optimizer as they are. Leaving \eqn{\mu}
#' unbounded would make the likelihood periodic and every maximum one of
#' infinitely many.
#'
#' # The distribution function
#'
#' The density has no elementary antiderivative. The base class would
#' integrate it once per observation, which measured 4.5 seconds at ten
#' thousand points; this class supplies [distrib_cdf.VonMises1Distrib()]
#' instead, a Fourier series integrated term by term whose term count is
#' \eqn{8.5\sqrt{\kappa} + 10}, measured. The **quantile** does come from
#' [continuous_distrib()], by root finding on that series.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\kappa > 0} the concentration, \eqn{I_m} the modified Bessel
#' function of the first kind of order \eqn{m}, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} the mean resultant length.
#' \eqn{\eta} is a parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `VonMises1Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"von mises1"`, `dimension`
#'   `"univariate"`, `bounds` `c(-pi, pi)`, `params` `c("mu", "kappa")`,
#'   `n_params` `2`, `params_bounds` the domains \eqn{(-\pi, \pi)} and
#'   \eqn{(0, \infty)}, and `link_params` the two links given here.
#'
#' @references
#' Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von Mises
#' distribution. *Journal of the Royal Statistical Society, Series C*,
#' **28**(2), 152-157.
#'
#' Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*, Chapter 3.
#' Wiley, Chichester.
#'
#' @importFrom linkfunctions7 bounded_link log_link
#' @importFrom stats runif
#'
#' @examples
#' d <- vonmises1_distrib()
#' d
#'
#' # The density integrates to one over the circle.
#' th <- list(mu = 0.5, kappa = 2)
#' integrate(function(v) distrib_pdf(d, v, th), -pi, pi)$value
#'
#' # The expected information is diagonal: direction and concentration are
#' # orthogonal, the sine having mean zero.
#' vapply(distrib_expected_hessian(d, 0, th), function(v) v[1], numeric(1))
#'
#' # The ordinary mean is not the mean direction, the interval being cut at
#' # plus or minus pi.
#' c(ordinary = mean(d, list(mu = 1.2, kappa = 2)), direction = 1.2)
#'
#' # A sample recovers the circular quantities, which no generic returns.
#' set.seed(1)
#' z <- distrib_rng(d, 3e5, th)
#' c(circular_mean = atan2(mean(sin(z)), mean(cos(z))), mu = 0.5)
#' c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2),
#'   A = numericals7::bessel_i_ratio(2))
#'
#' # Fitting recovers both parameters.
#' set.seed(7)
#' coef(fit_distrib(d, distrib_rng(d, 2000, list(mu = 0.8, kappa = 3))))
#'
#' @seealso
#' [vonmises2_distrib()] for the same law in the mean resultant length;
#' [gaussian1_distrib()] for the analogous family on the line, which the von
#' Mises approaches at a large concentration;
#' [numericals7::log_bessel_i()] and [numericals7::bessel_i_ratio()] for the
#' Bessel machinery; [fit_distrib()] to estimate the parameters;
#' [check_distrib()] to validate a family of your own against the same battery
#' this one passes; [VonMises1Distrib] for the class.
#' @export
vonmises1_distrib <- function(link_mu = bounded_link(lwr = -pi, upr = pi),
                             link_kappa = log_link()) {
  VonMises1Distrib(
    distrib_name = "von mises1", dimension = "univariate",
    bounds = c(-pi, pi),
    params = c("mu", "kappa"),
    params_interpretation = c(mu = "mean direction", kappa = "concentration"),
    n_params = 2,
    params_bounds = list(mu = c(-pi, pi), kappa = c(0, Inf)),
    link_params = list(mu = link_mu, kappa = link_kappa)
  )
}


#' @title von Mises Third-Order Derivatives
#' @name distrib_deriv3.VonMises1Distrib
#' @description
#' Computes the four distinct third derivatives of the von Mises log-density in
#' \eqn{\mu} and \eqn{\kappa}, in closed form. The log-density is
#' \eqn{\kappa\cos(y-\mu) - \log\{2\pi I_0(\kappa)\}}, **linear in**
#' \eqn{\kappa} apart from the normalizing constant, so a component naming one
#' \eqn{\mu} and two or more \eqn{\kappa} is **exactly zero**. What remains
#' cycles: the pure-\eqn{\mu} components run through
#' \eqn{\kappa\{\sin, -\cos, -\sin, \cos\}(y-\mu)} with the order, the
#' \eqn{\mu\mu\kappa} component is \eqn{-\cos(y-\mu)}, and the
#' pure-\eqn{\kappa} one is \eqn{-A''(\kappa)}, which
#' [numericals7::bessel_i_ratio_derivs()] supplies from the Riccati recursion
#' \eqn{A' = 1 - A/\kappa - A^2} differentiated.
#'
#' With `expected = TRUE` the method calls [expected_derivative()], which is
#' the one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. `kappa` must be strictly
#'   positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method. **Note the
#'   argument order**: `scale` precedes `expected` here, unlike on most
#'   families, so both are best given by name.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_kappa`,
#'   `mu_kappa_kappa` and `kappa_kappa_kappa`, each of length `length(y)`.
#'   `mu_kappa_kappa` is exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\kappa > 0} the concentration, \eqn{I_m} the modified Bessel
#' function of the first kind of order \eqn{m}, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_hessian.VonMises1Distrib()] for the order below,
#'   [distrib_deriv4.VonMises1Distrib()] for the order above,
#'   [numericals7::bessel_i_ratio_derivs()] for the derivatives of \eqn{A}, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, kappa = 2)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # The component with one mu and two kappa is exactly zero, the log-density
#' # being linear in kappa apart from the constant.
#' d3$mu_kappa_kappa
#'
#' # A central difference of the Hessian reproduces the pure-direction one.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.5 + eps, kappa = 2))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.5 - eps, kappa = 2))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
S7::method(distrib_deriv3, VonMises1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         expected = FALSE,
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  d <- y - theta[[1]]
  k <- theta[[2]]
  ad <- numericals7::bessel_i_ratio_derivs(k)
  n <- length(d)
  list(mu_mu_mu = -k * sin(d),
       mu_mu_kappa = -cos(d),
       mu_kappa_kappa = rep_len(0, n),
       kappa_kappa_kappa = rep_len(-ad$d2, n))
}

#' @title von Mises Fourth-Order Derivatives
#' @name distrib_deriv4.VonMises1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the von Mises log-density
#' in \eqn{\mu} and \eqn{\kappa}, in closed form, by the construction
#' [distrib_deriv3.VonMises1Distrib()] describes carried one order further. The
#' log-density is linear in \eqn{\kappa} apart from the normalizing constant,
#' so **two of the five are exactly zero**: `mu_mu_kappa_kappa` and
#' `mu_kappa_kappa_kappa`. The pure-direction component is
#' \eqn{\kappa\cos(y-\mu)}, the \eqn{\mu\mu\mu\kappa} one \eqn{-\sin(y-\mu)},
#' and the pure-concentration one \eqn{-A^{(3)}(\kappa)}.
#'
#' With `expected = TRUE` the method calls [expected_derivative()], which is
#' the one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `y`. `kappa` must be strictly
#'   positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method. **Note the
#'   argument order**: `scale` precedes `expected` here, so both are best given
#'   by name.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_kappa`, `mu_mu_kappa_kappa`, `mu_kappa_kappa_kappa` and
#'   `kappa_kappa_kappa_kappa`, each of length `length(y)`. The middle two are
#'   exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\kappa > 0} the concentration, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_deriv3.VonMises1Distrib()] for the order below and the
#'   construction, [numericals7::bessel_i_ratio_derivs()] for the derivatives
#'   of \eqn{A}, and [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, kappa = 2)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Two of the five are exactly zero.
#' c(d4$mu_mu_kappa_kappa[1], d4$mu_kappa_kappa_kappa[1])
#'
#' # A central difference of the third order reproduces the pure-direction
#' # component.
#' eps <- 1e-4
#' up <- distrib_deriv3(d, y, list(mu = 0.5 + eps, kappa = 2))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.5 - eps, kappa = 2))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
S7::method(distrib_deriv4, VonMises1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         expected = FALSE,
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  d <- y - theta[[1]]
  k <- theta[[2]]
  ad <- numericals7::bessel_i_ratio_derivs(k)
  n <- length(d)
  z <- rep_len(0, n)
  list(mu_mu_mu_mu = k * cos(d),
       mu_mu_mu_kappa = -sin(d),
       mu_mu_kappa_kappa = z,
       mu_kappa_kappa_kappa = z,
       kappa_kappa_kappa_kappa = rep_len(-ad$d3, n))
}


#' The Distribution Function of a von Mises by Its Bessel Series
#'
#' @description
#' Returns \eqn{F(x)} on \eqn{[-\pi, \pi)} from the Fourier expansion of the
#' density, integrated term by term, in place of a quadrature.
#'
#' @details
#' # Why a series
#'
#' The density has no elementary antiderivative, and the base class integrates
#' it numerically: one quadrature per observation, which measured 4.5 seconds
#' at ten thousand points and made this family a thousand times dearer than any
#' other. The series is rapidly convergent. From
#' \eqn{e^{\kappa\cos\theta} = I_0(\kappa) + 2\sum_{j\ge1} I_j(\kappa)
#' \cos(j\theta)}, integrating from \eqn{-\pi} to \eqn{x},
#' \deqn{F(x) = \frac{x + \pi}{2\pi} + \frac{1}{\pi I_0(\kappa)}
#'   \sum_{j\ge1} \frac{I_j(\kappa)}{j}
#'   \big[\sin(j(x - \mu)) + \sin(j(\pi + \mu))\big].}
#' The second sine is the lower limit, and it makes this the distribution
#' function of the family **as written**: the support is
#' \eqn{[-\pi, \pi)} with the direction inside it, not a variable wrapped
#' around the circle.
#'
#' Only the **ratios** \eqn{I_j/I_0} are needed, and
#' [numericals7::bessel_i_ratios()] gives them by a backward recurrence whose
#' loop runs over the series index and not over the data. The series is cheap
#' for that reason: \eqn{n} quadratures become a few dozen vectorized steps.
#'
#' # How many terms, measured
#'
#' Compared against the same series at four times the length, machine precision
#' is reached at 10 terms at \eqn{\kappa = 0.5}, 26 at 10, 90 at 100, 242 at
#' 1000 and 404 at 3000, always under \eqn{8.5\sqrt{\kappa} + 10}, which is the
#' rule used, with a floor of 20.
#'
#' The sum is accumulated over blocks of observations, because the natural
#' expression forms an \eqn{n \times m} matrix, which at a hundred thousand
#' points and a concentration of a hundred is already hundreds of megabytes.
#' The result is clamped to \eqn{[0, 1]}: the series is exact and its rounding
#' is not, and a distribution function cannot leave its own interval.
#'
#' @param y A numeric vector of quantiles. Below \eqn{-\pi} the value is 0 and
#'   at or above \eqn{\pi} it is 1. `NA` and non-finite values give `NA`.
#' @param mu A numeric vector of mean directions, recycled to the length of
#'   `y`.
#' @param kappa A numeric vector of concentrations, recycled to the length of
#'   `y`. A value at or below zero gives `NA`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of the recycled
#'   length of the inputs.
#'
#' @seealso [distrib_cdf.VonMises1Distrib()], which calls this;
#'   [numericals7::bessel_i_ratios()] for the recurrence; and
#'   [vonmises1_distrib()] for the family.
#'
#' @keywords internal
vm_cdf <- function(y, mu, kappa) {
  n <- max(length(y), length(mu), length(kappa))
  y <- rep_len(as.numeric(y), n)
  mu <- rep_len(as.numeric(mu), n)
  kappa <- rep_len(as.numeric(kappa), n)
  out <- rep(NA_real_, n)
  ok <- is.finite(y) & is.finite(mu) & is.finite(kappa) & kappa > 0
  out[!is.na(y) & y < -pi] <- 0
  out[!is.na(y) & y >= pi] <- 1
  ok <- ok & y >= -pi & y < pi
  if (!any(ok)) return(out)
  km <- max(kappa[ok])
  m <- max(20L, as.integer(ceiling(8.5 * sqrt(km))) + 10L)
  idx <- which(ok)
  # in blocks, so that nothing of size n by m is ever formed
  step <- max(1L, as.integer(ceiling(1e6 / m)))
  j <- seq_len(m)
  for (from in seq.int(1L, length(idx), by = step)) {
    ii <- idx[seq.int(from, min(from + step - 1L, length(idx)))]
    R <- numericals7::bessel_i_ratios(kappa[ii], m)
    s <- rowSums(R * (sin(outer(y[ii] - mu[ii], j)) +
                        sin(outer(pi + mu[ii], j))) /
                   rep(j, each = length(ii)))
    out[ii] <- (y[ii] + pi) / (2 * pi) + s / pi
  }
  # the series is exact and its rounding is not: a distribution function
  # cannot leave its own interval
  pmin(pmax(out, 0), 1)
}

#' @title von Mises Cumulative Distribution Function
#' @name distrib_cdf.VonMises1Distrib
#' @description
#' Computes \eqn{F(q) = P(Y \le q)} on \eqn{[-\pi, \pi)} from the Fourier
#' series of the density integrated term by term, through [vm_cdf()]. The
#' density has no elementary antiderivative, and the quadrature the base class
#' would use costs one integration per observation; the series replaces
#' \eqn{n} quadratures with a few dozen vectorized steps, using only the Bessel
#' **ratios** \eqn{I_j/I_0} from a backward recurrence.
#'
#' Below \eqn{-\pi} the value is 0 and at or above \eqn{\pi} it is 1, the
#' support being the declared interval. The result is clamped to
#' \eqn{[0, 1]}: the series is exact and its rounding is not.
#'
#' @param distrib A `VonMises1Distrib` object, from [vonmises1_distrib()].
#' @param q A numeric vector of angles.
#' @param theta A named list with components `mu` and `kappa`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(-\pi, \pi)} and `kappa` be strictly
#'   positive; a non-positive `kappa` gives `NA`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'   This method takes **no** `lower.tail` or `log.p`: the upper tail is
#'   `1 - F(q)` and the logarithm is `log(F(q))`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(kappa))`.
#'
#' @seealso [vm_cdf()] for the series and its measured term count,
#'   [distrib_pdf.VonMises1Distrib()] for the density,
#'   [distrib_quantile()], which inverts this by root finding, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- vonmises1_distrib()
#' th <- list(mu = 0.5, kappa = 2)
#' y <- c(-1, 0, 0.5, 2)
#'
#' # The series agrees with a direct quadrature of the density.
#' rbind(series = distrib_cdf(d, y, th),
#'       quadrature = vapply(y, function(v)
#'         integrate(function(u) distrib_pdf(d, u, th), -pi, v)$value,
#'         numeric(1)))
#'
#' # It runs from 0 to 1 across the declared support.
#' c(distrib_cdf(d, -pi, th), distrib_cdf(d, pi - 1e-12, th))
#'
#' # And the quantile function, which the base class obtains by root finding
#' # on this, inverts it.
#' q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
#' all.equal(distrib_cdf(d, q, th), c(0.25, 0.5, 0.75), tolerance = 1e-6)
S7::method(distrib_cdf, VonMises1Distrib) <- function(distrib, q,
                                                      theta, ...) {
  vm_cdf(q, theta[[1]], theta[[2]])
}
