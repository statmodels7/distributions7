#' @include distrib.R generics.R numerical_functions.R negbin2_distrib.R pseudohuber_distrib.R laplace_distrib.R weibull1_distrib.R gumbel_distrib.R skewnormal1_distrib.R skewt_distrib.R
NULL

#' @title Raw and Central Moments of a Distribution
#' @name moment
#'
#' @description
#' Computes raw moments \eqn{E[Y^p]} or central moments \eqn{E[(Y-\mu)^p]} of a
#' distribution numerically, via \code{\link{expectation}} (numerical integration
#' for continuous distributions, series summation for discrete ones).
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param theta A named list of parameters. Vectors are supported (vectorized computation).
#' @param p Numeric. The order of the moment. Can be a vector (recycled against \code{theta}).
#' @param central Logical. If \code{TRUE}, computes the central moment \eqn{E[(Y-\mu)^p]},
#'   where \eqn{\mu = E[Y]} (computed numerically unless \code{mu} is supplied).
#' @param mu Optional numeric. The centering value(s) used when \code{central = TRUE}.
#'   If \code{NULL}, the mean is computed numerically.
#' @param ... Additional arguments passed to \code{\link{expectation}}.
#'
#' @return A numeric vector of moments, with length equal to the maximum length
#'   among \code{theta} components and \code{p}.
#'
#' @examples
#' \dontrun{
#' d <- gaussian1_distrib()
#' moment(d, list(mu = 2, sigma = 3), p = 1)                 # 2
#' moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE) # 9
#' }
#'
#' @export
moment <- function(distrib, theta, p = 1, central = FALSE, mu = NULL, ...) {
  if (central) {
    if (is.null(mu)) {
      mu <- moment(distrib, theta, p = 1, central = FALSE, ...)
    }
  } else {
    mu <- 0
  }

  # Dot-prefixed argument names cannot collide with distribution parameters
  expectation(
    distrib,
    function(y, theta, .pow, .center) (y - .center)^.pow,
    theta,
    .pow = p,
    .center = mu,
    ...
  )
}

#' Mean of a Distribution Object
#'
#' @name mean.distrib
#' @description
#' Computes the expected value \eqn{E[Y]} of a distribution object numerically via
#' \code{\link{moment}}. Distribution classes with a closed-form mean may override
#' this method with an analytical version.
#'
#' @param x An object inheriting from class \code{"distrib"}.
#' @param theta A named list of parameters. Vectors are supported.
#' @param ... Additional arguments passed to \code{\link{moment}}.
#' @return A numeric vector of means.
#' @keywords internal
S7::method(mean, distrib) <- function(x, theta, ...) {
  moment(x, theta, p = 1, central = FALSE, ...)
}

#' Variance of a Distribution or Sample
#'
#' @description
#' Computes the variance. For \code{distrib} objects the second central moment is
#' evaluated numerically (analytical methods may override this for specific
#' distributions); for numeric vectors the sample variance \code{\link[stats]{var}} is returned.
#'
#' @param x An object inheriting from class \code{"distrib"}, or a numeric vector.
#' @param ... For \code{distrib} objects: \code{theta} (a named list of parameters) and
#'   further arguments passed to \code{\link{moment}}. For numeric vectors: \code{na.rm}.
#' @return A numeric vector.
#' @examples
#' variance(gaussian1_distrib(), list(mu = 0, sigma = 2))
#' variance(poisson_distrib(), list(mu = 3))
#'
#' @export
variance <- S7::new_generic("variance", "x")

#' @title Variance of a Distribution
#' @name variance.distrib
#' @description Computes the second central moment through \code{\link{moment}}, the mean being evaluated first and passed as the centre.
#' @param x A \code{distrib} object.
#' @param theta A named list of parameters.
#' @param ... Passed to \code{\link{moment}}.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, distrib) <- function(x, theta, ...) {
  m <- moment(x, theta, p = 1, central = FALSE, ...)
  moment(x, theta, p = 2, central = TRUE, mu = m, ...)
}

#' @title Sample Variance
#' @name variance.numeric
#' @description The sample variance of a numeric vector, delegated to \code{\link[stats]{var}}.
#' @param x A numeric vector.
#' @param na.rm Remove missing values first?
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(variance, S7::class_numeric) <- function(x, na.rm = FALSE, ...) {
  stats::var(x, na.rm = na.rm)
}

#' Standard Deviation of a Distribution or Sample
#'
#' @description
#' Computes the standard deviation as the square root of \code{\link{variance}}.
#' For numeric vectors the sample standard deviation \code{\link[stats]{sd}} is returned.
#'
#' @param x An object inheriting from class \code{"distrib"}, or a numeric vector.
#' @param ... For \code{distrib} objects: \code{theta} and further arguments passed to
#'   \code{\link{moment}}. For numeric vectors: \code{na.rm}.
#' @return A numeric vector.
#' @examples
#' std_dev(gaussian1_distrib(), list(mu = 0, sigma = 2))
#'
#' @export
std_dev <- S7::new_generic("std_dev", "x")

#' @title Standard Deviation of a Distribution
#' @name std_dev.distrib
#' @description The square root of \code{\link{variance}}.
#' @param x A \code{distrib} object.
#' @param theta A named list of parameters.
#' @param ... Passed to \code{\link{moment}}.
#' @return A numeric vector.
#' @keywords internal
S7::method(std_dev, distrib) <- function(x, theta, ...) {
  sqrt(variance(x, theta, ...))
}

#' @title Sample Standard Deviation
#' @name std_dev.numeric
#' @description The sample standard deviation of a numeric vector, delegated to \code{\link[stats]{sd}}.
#' @param x A numeric vector.
#' @param na.rm Remove missing values first?
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(std_dev, S7::class_numeric) <- function(x, na.rm = FALSE, ...) {
  stats::sd(x, na.rm = na.rm)
}

#' Skewness of a Distribution or Sample
#'
#' @description
#' Computes the skewness (third standardized moment). For \code{distrib} objects it is
#' evaluated numerically via \code{\link{moment}}; for numeric vectors the sample
#' skewness (population denominator) is returned.
#'
#' @param x An object inheriting from class \code{"distrib"}, or a numeric vector.
#' @param ... For \code{distrib} objects: \code{theta} and further arguments passed to
#'   \code{\link{moment}}. For numeric vectors: \code{na.rm}.
#' @return A numeric vector.
#' @examples
#' skewness(gaussian1_distrib(), list(mu = 0, sigma = 1))
#' skewness(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#'
#' @export
skewness <- S7::new_generic("skewness", "x")

# distrib method: theta is a named list of parameters
#' @title Skewness of a Distribution
#' @name skewness.distrib
#' @description Computes the standardised third central moment through \code{\link{moment}}.
#' @param x A \code{distrib} object.
#' @param theta A named list of parameters.
#' @param ... Passed to \code{\link{moment}}.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, distrib) <- function(x, theta, ...) {
  m <- moment(x, theta, p = 1, central = FALSE, ...)
  m2 <- moment(x, theta, p = 2, central = TRUE, mu = m, ...)
  m3 <- moment(x, theta, p = 3, central = TRUE, mu = m, ...)
  m3 / m2^1.5
}

#' @title Sample Skewness
#' @name skewness.numeric
#' @description The sample skewness of a numeric vector, computed from the sample's central moments.
#' @param x A numeric vector.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(skewness, S7::class_numeric) <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x <- x[!is.na(x)]
  m <- base::mean(x)
  s <- sqrt(base::mean((x - m)^2))
  base::mean((x - m)^3) / s^3
}

#' Excess Kurtosis of a Distribution or Sample
#'
#' @description
#' Computes the excess kurtosis (fourth standardized moment minus 3). For
#' \code{distrib} objects it is evaluated numerically via \code{\link{moment}}; for
#' numeric vectors the sample excess kurtosis (population denominator) is returned.
#'
#' @param x An object inheriting from class \code{"distrib"}, or a numeric vector.
#' @param ... For \code{distrib} objects: \code{theta} and further arguments passed to
#'   \code{\link{moment}}. For numeric vectors: \code{na.rm}.
#' @return A numeric vector.
#' @examples
#' kurtosis(gaussian1_distrib(), list(mu = 0, sigma = 1))
#' kurtosis(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#'
#' @export
kurtosis <- S7::new_generic("kurtosis", "x")

# distrib method: theta is a named list of parameters
#' @title Kurtosis of a Distribution
#' @name kurtosis.distrib
#' @description Computes the excess kurtosis, the standardised fourth central moment minus three through \code{\link{moment}}.
#' @param x A \code{distrib} object.
#' @param theta A named list of parameters.
#' @param ... Passed to \code{\link{moment}}.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, distrib) <- function(x, theta, ...) {
  m <- moment(x, theta, p = 1, central = FALSE, ...)
  m2 <- moment(x, theta, p = 2, central = TRUE, mu = m, ...)
  m4 <- moment(x, theta, p = 4, central = TRUE, mu = m, ...)
  m4 / m2^2 - 3
}

#' @title Sample Kurtosis
#' @name kurtosis.numeric
#' @description The sample excess kurtosis of a numeric vector, computed from the sample's central moments.
#' @param x A numeric vector.
#' @param ... Unused.
#' @return A single number.
#' @keywords internal
S7::method(kurtosis, S7::class_numeric) <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x <- x[!is.na(x)]
  m <- base::mean(x)
  s2 <- base::mean((x - m)^2)
  base::mean((x - m)^4) / s2^2 - 3
}

# --- ANALYTICAL OVERRIDES ---
# Closed-form moments for distributions whose numerical fallbacks would be
# needlessly slow (or are used internally, e.g. the Pseudo-Huber quantile
# bracket relies on its analytical variance).

#' @title Mean of the Negative Binomial Distribution
#' @name mean.NegBin2Distrib
#' @description Closed form, replacing the numerical default: \eqn{E[Y] = \mu}.
#' @param x A \code{NegBin2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Variance of the Negative Binomial Distribution
#' @name variance.NegBin2Distrib
#' @description Closed form, replacing the numerical default: \eqn{Var(Y) = \mu + \mu^2/\theta}.
#' @param x A \code{NegBin2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + theta[[1]]^2 / theta[[2]]
}

#' @title Skewness of the Negative Binomial Distribution
#' @name skewness.NegBin2Distrib
#' @description Closed form, replacing the numerical default: \eqn{(\theta + 2\mu)/\sqrt{\mu\theta(\theta+\mu)}}.
#' @param x A \code{NegBin2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  (th + 2 * mu) / sqrt(mu * th * (th + mu))
}

#' @title Kurtosis of the Negative Binomial Distribution
#' @name kurtosis.NegBin2Distrib
#' @description Closed form, replacing the numerical default: \eqn{6/\theta + \theta/(\mu(\theta+\mu))}.
#' @param x A \code{NegBin2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  6 / th + th / (mu * (th + mu))
}

#' @title Mean of the Pseudo-Huber Distribution
#' @name mean.PseudoHuberDistrib
#' @description Closed form, replacing the numerical default: \eqn{E[Y] = \mu}, by symmetry.
#' @param x A \code{PseudoHuberDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, PseudoHuberDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(3)])))
}

#' @title Variance of the Pseudo-Huber Distribution
#' @name variance.PseudoHuberDistrib
#' @description Closed form, replacing the numerical default: a closed form in exponentially scaled Bessel functions.
#' @param x A \code{PseudoHuberDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, PseudoHuberDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  sq_nu <- sqrt(theta[[3]])
  # Scaled Bessel functions: the exponential factors cancel in the ratio
  ratio <- besselK(sq_nu, 2, expon.scaled = TRUE) / besselK(sq_nu, 1, expon.scaled = TRUE)
  theta[[2]]^2 * sq_nu * ratio
}

#' @title Skewness of the Pseudo-Huber Distribution
#' @name skewness.PseudoHuberDistrib
#' @description Closed form, replacing the numerical default: zero, by symmetry.
#' @param x A \code{PseudoHuberDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, PseudoHuberDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(0, length.out = max(lengths(theta[seq_len(3)])))
}

#' @title Kurtosis of the Pseudo-Huber Distribution
#' @name kurtosis.PseudoHuberDistrib
#' @description Closed form, replacing the numerical default: a closed form in exponentially scaled Bessel functions.
#' @param x A \code{PseudoHuberDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, PseudoHuberDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  sq_nu <- sqrt(theta[[3]])
  k1 <- besselK(sq_nu, 1, expon.scaled = TRUE)
  k2 <- besselK(sq_nu, 2, expon.scaled = TRUE)
  k3 <- besselK(sq_nu, 3, expon.scaled = TRUE)
  3 * (k3 * k1) / (k2^2) - 3
}

#' @title Mean of the Laplace Distribution
#' @name mean.LaplaceDistrib
#' @description Closed form, replacing the numerical default: \eqn{E[Y] = \mu}.
#' @param x A \code{LaplaceDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Variance of the Laplace Distribution
#' @name variance.LaplaceDistrib
#' @description Closed form, replacing the numerical default: \eqn{Var(Y) = 2b^2}.
#' @param x A \code{LaplaceDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  2 * theta[[2]]^2
}

#' @title Skewness of the Laplace Distribution
#' @name skewness.LaplaceDistrib
#' @description Closed form, replacing the numerical default: zero, by symmetry.
#' @param x A \code{LaplaceDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(0, length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Kurtosis of the Laplace Distribution
#' @name kurtosis.LaplaceDistrib
#' @description Closed form, replacing the numerical default: \eqn{3} (excess).
#' @param x A \code{LaplaceDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(3, length.out = max(lengths(theta[seq_len(2)])))
}

# --- WEIBULL ----------------------------------------------------------------
#
# Every moment of a Weibull is the scale to a power times a gamma function of
# the shape: with g_k = Gamma(1 + k/sigma), E[Y^k] = mu^k g_k. The four
# quantities below are the standardised combinations of those, so only the
# first depends on mu at all.

#' Gamma Factors of a Weibull's Moments
#'
#' @description
#' Returns \eqn{g_k = \Gamma(1 + k/\sigma)} for \eqn{k = 1, \ldots, 4}, from
#' which every moment of a Weibull follows as \eqn{E[Y^k] = \mu^k g_k}.
#'
#' @param sigma The shape parameter.
#' @param k How many factors to return.
#'
#' @return A list of numeric vectors, \code{g1} to \code{gk}.
#'
#' @keywords internal
weibull_gamma_factors <- function(sigma, k = 4L) {
  out <- lapply(seq_len(k), function(j) gamma(1 + j / sigma))
  names(out) <- paste0("g", seq_len(k))
  out
}

#' @title Mean of the Weibull Distribution
#' @name mean.Weibull1Distrib
#' @description Closed form: \eqn{\mu\,\Gamma(1 + 1/\sigma)}. The scale
#'   \eqn{\mu} is not the mean, which is what this method reports.
#' @param x A \code{\link{Weibull1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] * gamma(1 + 1 / theta[[2]])
}

#' @title Variance of the Weibull Distribution
#' @name variance.Weibull1Distrib
#' @description Closed form: \eqn{\mu^2\left(g_2 - g_1^2\right)} with
#'   \eqn{g_k = \Gamma(1 + k/\sigma)}.
#' @param x A \code{\link{Weibull1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  g <- weibull_gamma_factors(theta[[2]], 2L)
  theta[[1]]^2 * (g$g2 - g$g1^2)
}

#' @title Skewness of the Weibull Distribution
#' @name skewness.Weibull1Distrib
#' @description Closed form:
#'   \eqn{(g_3 - 3g_1g_2 + 2g_1^3)/(g_2 - g_1^2)^{3/2}}, free of the scale.
#' @param x A \code{\link{Weibull1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  g <- weibull_gamma_factors(theta[[2]], 3L)
  v <- g$g2 - g$g1^2
  (g$g3 - 3 * g$g1 * g$g2 + 2 * g$g1^3) / v^1.5
}

#' @title Kurtosis of the Weibull Distribution
#' @name kurtosis.Weibull1Distrib
#' @description Closed form, excess:
#'   \eqn{(g_4 - 4g_1g_3 + 6g_1^2g_2 - 3g_1^4)/(g_2 - g_1^2)^2 - 3}, free of
#'   the scale.
#' @param x A \code{\link{Weibull1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  g <- weibull_gamma_factors(theta[[2]], 4L)
  v <- g$g2 - g$g1^2
  (g$g4 - 4 * g$g1 * g$g3 + 6 * g$g1^2 * g$g2 - 3 * g$g1^4) / v^2 - 3
}


# --- GUMBEL -----------------------------------------------------------------
#
# A location-scale family with a FIXED shape: the third and fourth standardised
# moments are constants, which is the substantive statement here rather than an
# arithmetic convenience.

#' @title Mean of the Gumbel Distribution
#' @name mean.GumbelDistrib
#' @description Closed form: \eqn{\mu + \gamma\sigma}, with \eqn{\gamma} the
#'   Euler-Mascheroni constant.
#' @param x A \code{\link{GumbelDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + (-digamma(1)) * theta[[2]]
}

#' @title Variance of the Gumbel Distribution
#' @name variance.GumbelDistrib
#' @description Closed form: \eqn{\pi^2\sigma^2/6}.
#' @param x A \code{\link{GumbelDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  pi^2 * theta[[2]]^2 / 6
}

#' @title Skewness of the Gumbel Distribution
#' @name skewness.GumbelDistrib
#' @description Closed form: \eqn{12\sqrt{6}\,\zeta(3)/\pi^3 \approx 1.1395},
#'   the same for every location and scale.
#' @param x A \code{\link{GumbelDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  # zeta(3), Apery's constant, written out rather than reached for: no base R
  # function returns it and one more dependency for one number is not worth it.
  zeta3 <- 1.2020569031595942854
  rep(12 * sqrt(6) * zeta3 / pi^3, length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Kurtosis of the Gumbel Distribution
#' @name kurtosis.GumbelDistrib
#' @description Closed form, excess: \eqn{12/5}, the same for every location
#'   and scale.
#' @param x A \code{\link{GumbelDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(12 / 5, length.out = max(lengths(theta[seq_len(2)])))
}


# --- SKEW NORMAL ------------------------------------------------------------
#
# Every standardised moment is a function of delta = alpha/sqrt(1 + alpha^2)
# alone, which is the quantity the shape enters through. As |alpha| grows delta
# tends to 1, so the skewness and the kurtosis the family can reach are bounded
# -- that bound is the reason the skew t exists.

#' The Shape a Skew Normal's Moments Depend On
#'
#' @description
#' Returns \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and the product
#' \eqn{b\delta} with \eqn{b = \sqrt{2/\pi}}, from which every moment of a skew
#' normal follows.
#'
#' @param alpha The shape parameter.
#'
#' @return A list with \code{delta} and \code{bd}.
#'
#' @keywords internal
skewnormal_delta <- function(alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  list(delta = delta, bd = sqrt(2 / pi) * delta)
}

#' @title Mean of the Skew Normal Distribution
#' @name mean.SkewNormal1Distrib
#' @description Closed form: \eqn{\mu + \sigma b \delta} with
#'   \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and \eqn{b = \sqrt{2/\pi}}.
#' @param x A \code{\link{SkewNormal1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + theta[[2]] * skewnormal_delta(theta[[3]])$bd
}

#' @title Variance of the Skew Normal Distribution
#' @name variance.SkewNormal1Distrib
#' @description Closed form: \eqn{\sigma^2\left(1 - b^2\delta^2\right)}.
#' @param x A \code{\link{SkewNormal1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 * (1 - skewnormal_delta(theta[[3]])$bd^2)
}

#' @title Skewness of the Skew Normal Distribution
#' @name skewness.SkewNormal1Distrib
#' @description Closed form:
#'   \eqn{\tfrac{4-\pi}{2}(b\delta)^3/(1 - b^2\delta^2)^{3/2}}. Its range is
#'   \eqn{(-0.9953, 0.9953)}, whatever the shape.
#' @param x A \code{\link{SkewNormal1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  bd <- skewnormal_delta(theta[[3]])$bd
  ((4 - pi) / 2) * bd^3 / (1 - bd^2)^1.5
}

#' @title Kurtosis of the Skew Normal Distribution
#' @name kurtosis.SkewNormal1Distrib
#' @description Closed form, excess:
#'   \eqn{2(\pi - 3)(b\delta)^4/(1 - b^2\delta^2)^2}. It is non-negative and
#'   bounded above by about \eqn{0.8692}.
#' @param x A \code{\link{SkewNormal1Distrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  bd <- skewnormal_delta(theta[[3]])$bd
  2 * (pi - 3) * bd^4 / (1 - bd^2)^2
}


# --- SKEW t -----------------------------------------------------------------
#
# The moments exist only up to order nu, and each is NaN below its threshold.
# The density is perfectly well defined there, so a family that reported a
# number would be reporting one that does not exist.

#' The Quantities a Skew t's Moments Are Built From
#'
#' @description
#' Returns \eqn{\delta = \alpha/\sqrt{1+\alpha^2}}, the constant
#' \eqn{b_\nu = \sqrt{\nu/\pi}\,\Gamma\{(\nu-1)/2\}/\Gamma(\nu/2)}, the mean
#' \eqn{\mu_z = \delta b_\nu} and the variance \eqn{\sigma_z^2} of the
#' standardised variable.
#'
#' @details
#' \eqn{b_\nu} is finite only for \eqn{\nu > 1} and \eqn{\sigma_z^2} only for
#' \eqn{\nu > 2}; each is \code{NaN} otherwise, and the moments that use it
#' inherit that.
#'
#' @param alpha The shape parameter.
#' @param nu The degrees of freedom.
#'
#' @return A list with \code{delta}, \code{bnu}, \code{mz} and \code{vz}.
#'
#' @keywords internal
skewt_moment_pieces <- function(alpha, nu) {
  delta <- alpha / sqrt(1 + alpha^2)
  bnu <- ifelse(nu > 1, sqrt(nu / pi) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2)), NaN)
  mz <- delta * bnu
  vz <- ifelse(nu > 2, nu / (nu - 2) - mz^2, NaN)
  list(delta = delta, bnu = bnu, mz = mz, vz = vz)
}

#' @title Mean of the Skew t Distribution
#' @name mean.SkewTDistrib
#' @description Closed form for \eqn{\nu > 1}: \eqn{\mu + \sigma\delta b_\nu};
#'   \code{NaN} otherwise, the mean not existing there.
#' @param x A \code{\link{SkewTDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  p <- skewt_moment_pieces(theta[[3]], theta[[4]])
  theta[[1]] + theta[[2]] * p$mz
}

#' @title Variance of the Skew t Distribution
#' @name variance.SkewTDistrib
#' @description Closed form for \eqn{\nu > 2}:
#'   \eqn{\sigma^2\{\nu/(\nu-2) - \mu_z^2\}}; \code{NaN} otherwise.
#' @param x A \code{\link{SkewTDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 * skewt_moment_pieces(theta[[3]], theta[[4]])$vz
}

#' @title Skewness of the Skew t Distribution
#' @name skewness.SkewTDistrib
#' @description Closed form for \eqn{\nu > 3}, from Azzalini and Capitanio
#'   (2003); \code{NaN} otherwise. Unlike the skew normal's, it is unbounded.
#' @param x A \code{\link{SkewTDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  nu <- theta[[4]]
  p <- skewt_moment_pieces(theta[[3]], nu)
  val <- (p$mz / p$vz^1.5) *
    (nu * (3 - p$delta^2) / (nu - 3) - 3 * nu / (nu - 2) + 2 * p$mz^2)
  ifelse(nu > 3, val, NaN)
}

#' @title Kurtosis of the Skew t Distribution
#' @name kurtosis.SkewTDistrib
#' @description Closed form, excess, for \eqn{\nu > 4}, from Azzalini and
#'   Capitanio (2003); \code{NaN} otherwise.
#' @param x A \code{\link{SkewTDistrib}}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  nu <- theta[[4]]
  p <- skewt_moment_pieces(theta[[3]], nu)
  val <- (3 * nu^2 / ((nu - 2) * (nu - 4)) -
    4 * p$mz^2 * nu * (3 - p$delta^2) / (nu - 3) +
    6 * p$mz^2 * nu / (nu - 2) - 3 * p$mz^4) / p$vz^2 - 3
  ifelse(nu > 4, val, NaN)
}
