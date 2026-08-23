#' @include distrib.R generics.R
NULL

#' @title S7 Class for the von Mises Distribution
#' @name VonMises1Distrib
#'
#' @description A subclass of `continuous_distrib` representing the von
#'   Mises distribution on the circle.
#' @inheritParams distrib
#' @return An object of class `VonMises1Distrib`.
#' @seealso [vonmises1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.VonMises1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.VonMises1Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.VonMises1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.VonMises1Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.VonMises1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.VonMises1Distrib],
#'   [`distrib_rng()`][distrib_rng.VonMises1Distrib]
#'
#' The distribution function and the quantile come from
#' [continuous_distrib()], by quadrature and root finding over the
#' bounded support.
VonMises1Distrib <- S7::new_class("VonMises1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title von Mises Density
#' @name distrib_pdf.VonMises1Distrib
#' @description
#' \deqn{f(y) = \dfrac{e^{\kappa \cos(y - \mu)}}{2\pi I_0(\kappa)},
#'       \qquad y \in [-\pi, \pi)}
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [vonmises1_distrib()]
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
#' The rejection algorithm of Best and Fisher (1979), which draws from a
#' wrapped Cauchy envelope and accepts with a probability that does not
#' involve the Bessel function at all.
#' @param distrib A `VonMises1Distrib` object.
#' @param n The number of draws.
#' @param theta A list containing `mu` and `kappa`.
#' @return A numeric vector of length `n`, in \eqn{[-\pi, \pi)}.
#' @seealso [vonmises1_distrib()]
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

#' @title von Mises Analytical Gradient
#' @name distrib_gradient.VonMises1Distrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \kappa \sin(y - \mu), \qquad
#'       \dfrac{\partial \ell}{\partial \kappa} = \cos(y - \mu) - A(\kappa)}
#' with \eqn{A = I_1/I_0}, the derivative of \eqn{\log I_0}.
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list with the `mu` and `kappa` components.
#' @seealso [vonmises1_distrib()]
S7::method(distrib_gradient, VonMises1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  d <- y - theta[[1]]
  list(mu = theta[[2]] * sin(d), kappa = cos(d) - numericals7::bessel_i_ratio(theta[[2]]))
}

#' @title von Mises Analytical Observed Hessian
#' @name distrib_hessian.VonMises1Distrib
#' @description
#' \deqn{\ell^{(\mu\mu)} = -\kappa\cos(y-\mu), \qquad
#'       \ell^{(\mu\kappa)} = \sin(y-\mu), \qquad
#'       \ell^{(\kappa\kappa)} = -A'(\kappa)}
#' the last free of the data, \eqn{\log I_0} being the only place
#' \eqn{\kappa} appears alone.
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso [vonmises1_distrib()]
S7::method(distrib_hessian, VonMises1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...) {
  d <- y - theta[[1]]
  k <- theta[[2]]
  list(mu_mu = -k * cos(d), mu_kappa = sin(d),
       kappa_kappa = rep_len(-numericals7::bessel_i_ratio_derivs(k)$d1, length(d)))
}

#' @title von Mises Analytical Expected Hessian
#' @name distrib_expected_hessian.VonMises1Distrib
#' @description
#' \eqn{\mathbb{E}[\cos(Y-\mu)] = A(\kappa)} and
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0} by symmetry, so
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\kappa A(\kappa), \qquad
#'       \mathbb{E}[\ell^{(\mu\kappa)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\kappa\kappa)}] = -A'(\kappa)}
#' The location and the concentration are therefore orthogonal.
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso [vonmises1_distrib()]
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

#' @title von Mises Response Gradient
#' @name distrib_grad_y.VonMises1Distrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\kappa \sin(y - \mu)}
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [vonmises1_distrib()]
S7::method(distrib_grad_y, VonMises1Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * sin(y - theta[[1]])
}

#' @title von Mises Response Hessian
#' @name distrib_hess_y.VonMises1Distrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\kappa \cos(y - \mu)}
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [vonmises1_distrib()]
S7::method(distrib_hess_y, VonMises1Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * cos(y - theta[[1]])
}

# --- CONSTRUCTOR WRAPPER ---

#' von Mises Distribution Object
#'
#' @description
#' Creates a distribution object for the von Mises distribution, the natural
#' family for an angle, parametrized by a mean direction \eqn{\mu} and a
#' concentration \eqn{\kappa}.
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   [linkfunctions7::bounded_link()] on \eqn{(-\pi, \pi)}.
#' @param link_kappa A link function object for \eqn{\kappa}. Defaults to
#'   [linkfunctions7::log_link()] to ensure positivity.
#'
#' @details
#' The observation is an angle and the support is a circle, written here as
#' \eqn{[-\pi, \pi)}. This is the first family in the package whose support has
#' that topology: the two ends of the interval are the same point, so a density
#' need not vanish at either.
#'
#' **Density:**
#' \deqn{f(y) = \dfrac{e^{\kappa\cos(y-\mu)}}{2\pi I_0(\kappa)}}
#' The normalizing constant is a modified Bessel function, and it is evaluated
#' exponentially scaled with the exponent added back, so a concentration past
#' \eqn{\kappa = 700} does not overflow.
#'
#' **Moments.** [mean()], [variance()],
#' [skewness()] and [kurtosis()] are the ordinary moments
#' of \eqn{Y} as a number on \eqn{[-\pi, \pi)}, and they are obtained
#' numerically. They are not the circular quantities and they are not what
#' \eqn{\mu} and \eqn{\kappa} describe: \eqn{\mu} is the mean *direction*,
#' and \eqn{\mathbb{E}[Y] \ne \mu} whenever \eqn{\mu \ne 0}, because the
#' interval is cut at \eqn{\pm\pi} rather than at \eqn{\mu \pm \pi} and the
#' density is not symmetric about \eqn{\mu} on it. At \eqn{\mu = 1.2} and
#' \eqn{\kappa = 2} the ordinary mean is 1.079. The circular mean is
#' \eqn{\mu} and the mean resultant length is
#' \eqn{\rho = I_1(\kappa)/I_0(\kappa)}, both closed form, and neither is
#' returned by a generic whose name means something else.
#'
#' **Score, observed and expected Hessian.** Writing
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} for the derivative of
#' \eqn{\log I_0},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
#'       \dfrac{\partial\ell}{\partial\kappa} = \cos(y-\mu) - A(\kappa),}
#' and every second derivative is closed form, with
#' \eqn{A'(\kappa) = 1 - A(\kappa)/\kappa - A(\kappa)^2} obtained from the
#' Bessel recurrences rather than from a further evaluation. Since
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0} by symmetry, **the two parameters
#' are orthogonal**: the expected information is diagonal, and Fisher scoring
#' updates the direction and the concentration independently.
#'
#' \eqn{A(\kappa)} is the mean resultant length, so \eqn{A'(\kappa)} is the
#' variance of \eqn{\cos(Y-\mu)} and is positive, which is what makes the
#' information positive definite.
#'
#' **The mean direction is carried on a bounded chart**, the default link
#' mapping the free scale onto \eqn{(-\pi, \pi)}. That keeps the parameter
#' identified, at the cost that a fit cannot walk across the boundary: data
#' concentrated near \eqn{\pm\pi} are better rotated before fitting than
#' handed to the optimizer as they are. Leaving \eqn{\mu} unbounded instead
#' would make the likelihood periodic and every maximum one of infinitely
#' many.
#'
#' **Parameter domains:**
#'
#' - \eqn{\mu \in (-\pi, \pi)}
#' - \eqn{\kappa \in (0, +\infty)}
#'
#' The distribution function has no elementary form and comes from the base
#' class by quadrature over the bounded support, with the quantile by root
#' finding on it.
#'
#' @return An S7 object of class `VonMises1Distrib`.
#'
#' @references
#' Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von Mises
#' distribution. *Applied Statistics* 28, 152-157.
#'
#' Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*. Wiley.
#'
#' @seealso [gaussian1_distrib()] for the analogous family on the line
#'
#' @importFrom linkfunctions7 bounded_link log_link
#' @importFrom stats runif
#' @examples
#' d <- vonmises1_distrib()
#' d@params
#'
#' theta <- list(mu = 0.5, kappa = 2)
#' distrib_pdf(d, c(-1, 0, 0.5, 2), theta)
#'
#' # the expected information is diagonal: direction and concentration are
#' # orthogonal, the sine having mean zero
#' distrib_expected_hessian(d, 0, theta)
#'
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


#' @title von Mises Third and Fourth Derivatives
#' @name distrib_deriv3.VonMises1Distrib
#' @description
#' Closed form at both orders. The log-density is
#' \eqn{\kappa\cos(y-\mu) - \log(2\pi I_0(\kappa))}, linear in
#' \eqn{\kappa} apart from the normalizing constant, so every component with
#' one \eqn{\mu} and two or more \eqn{\kappa} vanishes exactly. The pure
#' \eqn{\mu} components cycle through
#' \eqn{\kappa\{\sin, -\cos, -\sin, \cos\}(y-\mu)}, and the pure
#' \eqn{\kappa} ones are minus the derivatives of \eqn{A(\kappa)}, which
#' \pkg{numericals7} supplies from the Riccati recursion
#' \eqn{A' = 1 - A/\kappa - A^2}.
#' @param distrib A `VonMises1Distrib` object.
#' @param y A numeric vector of angles.
#' @param theta A list containing `mu` and `kappa`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @param approx The approximation used when `expected` is `TRUE`.
#' @param nsim Monte Carlo draws when `approx = "mc"`.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [vonmises1_distrib()]
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

#' @rdname distrib_deriv3.VonMises1Distrib
#' @name distrib_deriv4.VonMises1Distrib
#' @return A named list of fourth-derivative components.
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
#' \eqn{F(x)} on \eqn{[-\pi, \pi)}, from the Fourier expansion of the
#' density integrated term by term.
#'
#' @details
#' The density has no elementary antiderivative, and the base class
#' integrates it numerically: one quadrature per observation, which measured
#' 4.5 seconds at ten thousand points and made this family a thousand times
#' dearer than any other. It has a rapidly convergent series instead. From
#' \eqn{e^{\kappa\cos\theta} = I_0(\kappa) + 2\sum_{j\ge1} I_j(\kappa)
#' \cos(j\theta)}, integrating from \eqn{-\pi} to \eqn{x},
#' \deqn{F(x) = \frac{x + \pi}{2\pi} + \frac{1}{\pi I_0(\kappa)}
#'   \sum_{j\ge1} \frac{I_j(\kappa)}{j}
#'   \big[\sin(j(x - \mu)) + \sin(j(\pi + \mu))\big].}
#' The second sine is the lower limit and is what makes this the
#' distribution function of the family as written -- the support is
#' \eqn{[-\pi, \pi)} with the location inside it, not a variable wrapped
#' around the circle.
#'
#' Only the RATIOS \eqn{I_j/I_0} are needed and
#' [numericals7::bessel_i_ratios()] gives them by a backward
#' recurrence whose loop runs over the series index rather than over the
#' data. That is what makes the series cheaper: \eqn{n} quadratures become a
#' few dozen vectorized steps.
#'
#' HOW MANY TERMS is measured rather than assumed. Comparing against the same
#' series at four times the length, machine precision is reached at 10 terms
#' at \eqn{\kappa = 0.5}, 26 at 10, 90 at 100, 242 at 1000 and 404 at 3000 --
#' always under \eqn{8.5\sqrt{\kappa} + 10}, which is the rule used. The sum
#' is accumulated over blocks of observations because the natural expression
#' forms an \eqn{n \times m} matrix, which at a hundred thousand points and a
#' concentration of a hundred is already hundreds of megabytes.
#'
#' @param y The quantiles.
#' @param mu The mean directions.
#' @param kappa The concentrations.
#'
#' @return The distribution function at `y`.
#'
#' @seealso [numericals7::bessel_i_ratios()]
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

S7::method(distrib_cdf, VonMises1Distrib) <- function(distrib, q,
                                                      theta, ...) {
  vm_cdf(q, theta[[1]], theta[[2]])
}
