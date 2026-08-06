#' @include distrib.R generics.R numerical_functions.R negbin2_distrib.R pseudohuber_distrib.R laplace_distrib.R weibull1_distrib.R gumbel_distrib.R skewnormal1_distrib.R skewt_distrib.R gaussian1_distrib.R cauchy_distrib.R logistic_distrib.R student_t1_distrib.R gamma2_distrib.R exponential_distrib.R chisq_distrib.R lognormal1_distrib.R invgauss1_distrib.R beta1_distrib.R gpd_distrib.R gengamma1_distrib.R poisson_distrib.R bernoulli_distrib.R binomial_distrib.R geometric_distrib.R negbin1_distrib.R betabinom1_distrib.R gaussian2_distrib.R gaussian3_distrib.R gamma1_distrib.R invgauss2_distrib.R beta2_distrib.R betabinom2_distrib.R pig1_distrib.R pig2_distrib.R
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
#' @description Computes the second central moment through \code{\link{moment}}, the mean being evaluated first and passed as the center.
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
#' @description Computes the standardized third central moment through \code{\link{moment}}.
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
#' @description Computes the excess kurtosis, the standardized fourth central moment minus three through \code{\link{moment}}.
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
# quantities below are the standardized combinations of those, so only the
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
# A location-scale family with a FIXED shape: the third and fourth standardized
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
# Every standardized moment is a function of delta = alpha/sqrt(1 + alpha^2)
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
#' standardized variable.
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

# ---------------------------------------------------------------------------
# Closed-form moments.
#
# Every method here replaces the numerical default of moment(), which is a
# quadrature and carries its error: the excess kurtosis of a gaussian comes
# back as -4.8e-07 through it and as exactly 0 from the formula. Kurtosis is
# always the EXCESS, so a gaussian is 0 and a Laplace is 3.
#
# A moment that does not exist is NaN rather than an error. The alternative,
# letting the quadrature diverge, reports "maximum number of subdivisions
# reached", which names the arithmetic instead of the mathematics.
# ---------------------------------------------------------------------------

#' Recycle a Constant Moment to the Length of the Parameters
#'
#' @description
#' Returns \code{value} repeated to the length the parameters imply, which is
#' what a moment that does not depend on them still has to be.
#'
#' @param theta An aligned named list of parameters.
#' @param k The number of parameters to read the length from.
#' @param value The constant.
#'
#' @return A numeric vector.
#'
#' @keywords internal
moment_const <- function(theta, k, value) {
  rep(value, length.out = max(lengths(theta[seq_len(k)])))
}

# --- gaussian1 -------------------------------------------------------------

#' @title Mean of the Gaussian Distribution
#' @name mean.Gaussian1Distrib
#' @description Closed form: \eqn{E[Y] = \mu}.
#' @param x A \code{Gaussian1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Gaussian1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gaussian Distribution
#' @name variance.Gaussian1Distrib
#' @description Closed form: \eqn{\operatorname{Var}(Y) = \sigma^2}.
#' @param x A \code{Gaussian1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Gaussian1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 2L, 0)
}

#' @title Skewness of the Gaussian Distribution
#' @name skewness.Gaussian1Distrib
#' @description Closed form: zero, the density being symmetric.
#' @param x A \code{Gaussian1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Gaussian1Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Kurtosis of the Gaussian Distribution
#' @name kurtosis.Gaussian1Distrib
#' @description Closed form: zero excess, which is what the scale is set by.
#' @param x A \code{Gaussian1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Gaussian1Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

# --- cauchy ----------------------------------------------------------------

#' @title The Cauchy Distribution Has No Moments
#' @name mean.CauchyDistrib
#' @description
#' \code{NaN}. No moment of the Cauchy exists, the tails being too heavy for
#' \eqn{\int |y| f(y) \, dy} to converge, and that is the honest answer rather
#' than the divergent quadrature the numerical default would attempt.
#' @param x A \code{CauchyDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of \code{NaN}.
#' @keywords internal
S7::method(mean, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Variance
#' @name variance.CauchyDistrib
#' @description \code{NaN}; see \code{\link[=mean.CauchyDistrib]{mean()}}.
#' @param x A \code{CauchyDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of \code{NaN}.
#' @keywords internal
S7::method(variance, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Skewness
#' @name skewness.CauchyDistrib
#' @description \code{NaN}; see \code{\link[=mean.CauchyDistrib]{mean()}}.
#' @param x A \code{CauchyDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of \code{NaN}.
#' @keywords internal
S7::method(skewness, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Kurtosis
#' @name kurtosis.CauchyDistrib
#' @description \code{NaN}; see \code{\link[=mean.CauchyDistrib]{mean()}}.
#' @param x A \code{CauchyDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of \code{NaN}.
#' @keywords internal
S7::method(kurtosis, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

# --- logistic --------------------------------------------------------------

#' @title Mean of the Logistic Distribution
#' @name mean.LogisticDistrib
#' @description Closed form: \eqn{E[Y] = \mu}.
#' @param x A \code{LogisticDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, LogisticDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Logistic Distribution
#' @name variance.LogisticDistrib
#' @description Closed form: \eqn{\operatorname{Var}(Y) = \pi^2\sigma^2/3}.
#' @param x A \code{LogisticDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, LogisticDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  pi^2 * theta[[2]]^2 / 3 + moment_const(theta, 2L, 0)
}

#' @title Skewness of the Logistic Distribution
#' @name skewness.LogisticDistrib
#' @description Closed form: zero, the density being symmetric.
#' @param x A \code{LogisticDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, LogisticDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Kurtosis of the Logistic Distribution
#' @name kurtosis.LogisticDistrib
#' @description Closed form: excess \eqn{6/5}, free of the parameters.
#' @param x A \code{LogisticDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, LogisticDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 6 / 5)
}

# --- student_t1 ------------------------------------------------------------
#
# Every moment of order k exists only for nu > k, so each returns NaN below
# its threshold. The variance is infinite on 1 < nu <= 2 and the excess
# kurtosis on 2 < nu <= 4, which is a different statement from not existing.

#' @title Mean of the Student t Distribution
#' @name mean.StudentT1Distrib
#' @description Closed form: \eqn{\mu} for \eqn{\nu > 1}, and \code{NaN} below,
#'   where the first moment does not exist.
#' @param x A \code{StudentT1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, StudentT1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(3)]))
  mu <- rep(theta[[1]], length.out = n)
  nu <- rep(theta[[3]], length.out = n)
  ifelse(nu > 1, mu, NaN)
}

#' @title Variance of the Student t Distribution
#' @name variance.StudentT1Distrib
#' @description
#' Closed form: \eqn{\sigma^2\nu/(\nu-2)} for \eqn{\nu > 2}, infinite for
#' \eqn{1 < \nu \le 2}, and \code{NaN} at or below one.
#' @param x A \code{StudentT1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, StudentT1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(3)]))
  s <- rep(theta[[2]], length.out = n)
  nu <- rep(theta[[3]], length.out = n)
  out <- ifelse(nu > 2, s^2 * nu / (nu - 2), Inf)
  out[nu <= 1] <- NaN
  out
}

#' @title Skewness of the Student t Distribution
#' @name skewness.StudentT1Distrib
#' @description Closed form: zero for \eqn{\nu > 3} by symmetry, \code{NaN}
#'   below, where the third moment does not exist.
#' @param x A \code{StudentT1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, StudentT1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(3)]))
  nu <- rep(theta[[3]], length.out = n)
  ifelse(nu > 3, 0, NaN)
}

#' @title Kurtosis of the Student t Distribution
#' @name kurtosis.StudentT1Distrib
#' @description
#' Closed form: excess \eqn{6/(\nu-4)} for \eqn{\nu > 4}, infinite for
#' \eqn{2 < \nu \le 4}, and \code{NaN} at or below two.
#' @param x A \code{StudentT1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, StudentT1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(3)]))
  nu <- rep(theta[[3]], length.out = n)
  out <- ifelse(nu > 4, 6 / (nu - 4), Inf)
  out[nu <= 2] <- NaN
  out
}

# --- gamma2 ----------------------------------------------------------------
#
# The parametrization is the mean and the variance, so the first two moments
# are the parameters and the shape a = mu^2/sigma2 carries the rest.

#' @title Mean of the Gamma Distribution
#' @name mean.Gamma2Distrib
#' @description Closed form: \eqn{E[Y] = \mu}, a parameter of this
#'   parametrization.
#' @param x A \code{Gamma2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gamma Distribution
#' @name variance.Gamma2Distrib
#' @description Closed form: \eqn{\sigma^2}, a parameter of this
#'   parametrization.
#' @param x A \code{Gamma2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[2]]
}

#' @title Skewness of the Gamma Distribution
#' @name skewness.Gamma2Distrib
#' @description Closed form: \eqn{2/\sqrt{a}} with shape
#'   \eqn{a = \mu^2/\sigma^2}, hence \eqn{2\sigma/\mu}.
#' @param x A \code{Gamma2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  2 * sqrt(theta[[2]]) / theta[[1]]
}

#' @title Kurtosis of the Gamma Distribution
#' @name kurtosis.Gamma2Distrib
#' @description Closed form: excess \eqn{6/a = 6\sigma^2/\mu^2}.
#' @param x A \code{Gamma2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  6 * theta[[2]] / theta[[1]]^2
}

# --- exponential -----------------------------------------------------------

#' @title Mean of the Exponential Distribution
#' @name mean.ExponentialDistrib
#' @description Closed form: \eqn{E[Y] = \mu}, the parameter itself.
#' @param x An \code{ExponentialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, ExponentialDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Exponential Distribution
#' @name variance.ExponentialDistrib
#' @description Closed form: \eqn{\mu^2}; the coefficient of variation is one.
#' @param x An \code{ExponentialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, ExponentialDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]^2
}

#' @title Skewness of the Exponential Distribution
#' @name skewness.ExponentialDistrib
#' @description Closed form: 2, free of the parameter.
#' @param x An \code{ExponentialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, ExponentialDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 1L, 2)
}

#' @title Kurtosis of the Exponential Distribution
#' @name kurtosis.ExponentialDistrib
#' @description Closed form: excess 6, free of the parameter.
#' @param x An \code{ExponentialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, ExponentialDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 1L, 6)
}

# --- chisq -----------------------------------------------------------------
#
# Parametrized by its mean, which for a chi-squared is its degrees of freedom.

#' @title Mean of the Chi-Squared Distribution
#' @name mean.ChisqDistrib
#' @description Closed form: \eqn{\mu}, which for this family is also the
#'   degrees of freedom.
#' @param x A \code{ChisqDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, ChisqDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Chi-Squared Distribution
#' @name variance.ChisqDistrib
#' @description Closed form: \eqn{2\mu}.
#' @param x A \code{ChisqDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, ChisqDistrib) <- function(x, theta, ...) {
  2 * align_theta(x, theta)[[1]]
}

#' @title Skewness of the Chi-Squared Distribution
#' @name skewness.ChisqDistrib
#' @description Closed form: \eqn{\sqrt{8/\mu}}.
#' @param x A \code{ChisqDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, ChisqDistrib) <- function(x, theta, ...) {
  sqrt(8 / align_theta(x, theta)[[1]])
}

#' @title Kurtosis of the Chi-Squared Distribution
#' @name kurtosis.ChisqDistrib
#' @description Closed form: excess \eqn{12/\mu}.
#' @param x A \code{ChisqDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, ChisqDistrib) <- function(x, theta, ...) {
  12 / align_theta(x, theta)[[1]]
}

# --- lognormal1 ------------------------------------------------------------
#
# The parameters live on the log scale, so none of them is a moment of Y.

#' @title Mean of the Lognormal Distribution
#' @name mean.Lognormal1Distrib
#' @description Closed form: \eqn{\exp(\mu + \sigma^2/2)}. The parameters
#'   describe \eqn{\log Y}, so \eqn{\mu} is not the mean.
#' @param x A \code{Lognormal1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  exp(theta[[1]] + theta[[2]] / 2)
}

#' @title Variance of the Lognormal Distribution
#' @name variance.Lognormal1Distrib
#' @description Closed form:
#'   \eqn{(e^{\sigma^2}-1)\,e^{2\mu+\sigma^2}}.
#' @param x A \code{Lognormal1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  expm1(theta[[2]]) * exp(2 * theta[[1]] + theta[[2]])
}

#' @title Skewness of the Lognormal Distribution
#' @name skewness.Lognormal1Distrib
#' @description Closed form:
#'   \eqn{(e^{\sigma^2}+2)\sqrt{e^{\sigma^2}-1}}, free of \eqn{\mu}.
#' @param x A \code{Lognormal1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  s2 <- theta[[2]]
  (exp(s2) + 2) * sqrt(expm1(s2)) + moment_const(theta, 2L, 0)
}

#' @title Kurtosis of the Lognormal Distribution
#' @name kurtosis.Lognormal1Distrib
#' @description Closed form: excess
#'   \eqn{e^{4\sigma^2} + 2e^{3\sigma^2} + 3e^{2\sigma^2} - 6}, free of
#'   \eqn{\mu}.
#' @param x A \code{Lognormal1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  s2 <- theta[[2]]
  exp(4 * s2) + 2 * exp(3 * s2) + 3 * exp(2 * s2) - 6 + moment_const(theta, 2L, 0)
}

# --- invgauss1 -------------------------------------------------------------
#
# statmod parametrizes by the mean and a dispersion phi, so Var = phi mu^3.

#' @title Mean of the Inverse Gaussian Distribution
#' @name mean.InvGauss1Distrib
#' @description Closed form: \eqn{\mu}, a parameter of this parametrization.
#' @param x An \code{InvGauss1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Inverse Gaussian Distribution
#' @name variance.InvGauss1Distrib
#' @description Closed form: \eqn{\phi\mu^3}, the dispersion multiplying the
#'   cube of the mean.
#' @param x An \code{InvGauss1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]] * theta[[1]]^3
}

#' @title Skewness of the Inverse Gaussian Distribution
#' @name skewness.InvGauss1Distrib
#' @description Closed form: \eqn{3\sqrt{\phi\mu}}.
#' @param x An \code{InvGauss1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  3 * sqrt(theta[[2]] * theta[[1]])
}

#' @title Kurtosis of the Inverse Gaussian Distribution
#' @name kurtosis.InvGauss1Distrib
#' @description Closed form: excess \eqn{15\phi\mu}.
#' @param x An \code{InvGauss1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  15 * theta[[2]] * theta[[1]]
}

# --- beta1 -----------------------------------------------------------------
#
# The shapes are a = mu phi and b = (1 - mu) phi, so a + b = phi.

#' @title Mean of the Beta Distribution
#' @name mean.Beta1Distrib
#' @description Closed form: \eqn{\mu}, a parameter of this parametrization.
#' @param x A \code{Beta1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Beta Distribution
#' @name variance.Beta1Distrib
#' @description Closed form: \eqn{\mu(1-\mu)/(\phi+1)}, which is why
#'   \eqn{\phi} is called a precision.
#' @param x A \code{Beta1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] * (1 - theta[[1]]) / (theta[[2]] + 1)
}

#' @title Skewness of the Beta Distribution
#' @name skewness.Beta1Distrib
#' @description Closed form in the shapes \eqn{a = \mu\phi} and
#'   \eqn{b = (1-\mu)\phi}:
#'   \eqn{2(b-a)\sqrt{a+b+1}\,/\,\{(a+b+2)\sqrt{ab}\}}.
#' @param x A \code{Beta1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b))
}

#' @title Kurtosis of the Beta Distribution
#' @name kurtosis.Beta1Distrib
#' @description Closed form in the shapes \eqn{a} and \eqn{b}: excess
#'   \eqn{6\{(a-b)^2(a+b+1) - ab(a+b+2)\}/\{ab(a+b+2)(a+b+3)\}}.
#' @param x A \code{Beta1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  6 * ((a - b)^2 * (a + b + 1) - a * b * (a + b + 2)) /
    (a * b * (a + b + 2) * (a + b + 3))
}

# --- gpd -------------------------------------------------------------------
#
# A moment of order k exists only for xi < 1/k, which is the condition for
# the integrand to be integrable in the upper tail.

#' @title Mean of the Generalized Pareto Distribution
#' @name mean.GPDDistrib
#' @description Closed form: \eqn{\sigma/(1-\xi)} for \eqn{\xi < 1}, and
#'   infinite at or above one.
#' @param x A \code{GPDDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, GPDDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(2)]))
  s <- rep(theta[[1]], length.out = n)
  xi <- rep(theta[[2]], length.out = n)
  ifelse(xi < 1, s / (1 - xi), Inf)
}

#' @title Variance of the Generalized Pareto Distribution
#' @name variance.GPDDistrib
#' @description Closed form:
#'   \eqn{\sigma^2/\{(1-\xi)^2(1-2\xi)\}} for \eqn{\xi < 1/2}, and infinite
#'   at or above one half.
#' @param x A \code{GPDDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, GPDDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(2)]))
  s <- rep(theta[[1]], length.out = n)
  xi <- rep(theta[[2]], length.out = n)
  ifelse(xi < 0.5, s^2 / ((1 - xi)^2 * (1 - 2 * xi)), Inf)
}

#' @title Skewness of the Generalized Pareto Distribution
#' @name skewness.GPDDistrib
#' @description Closed form:
#'   \eqn{2(1+\xi)\sqrt{1-2\xi}/(1-3\xi)} for \eqn{\xi < 1/3}, and infinite
#'   at or above one third.
#' @param x A \code{GPDDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, GPDDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(2)]))
  xi <- rep(theta[[2]], length.out = n)
  ifelse(xi < 1 / 3, 2 * (1 + xi) * sqrt(1 - 2 * xi) / (1 - 3 * xi), Inf)
}

#' @title Kurtosis of the Generalized Pareto Distribution
#' @name kurtosis.GPDDistrib
#' @description Closed form: excess
#'   \eqn{3(1-2\xi)(2\xi^2+\xi+3)/\{(1-3\xi)(1-4\xi)\} - 3} for
#'   \eqn{\xi < 1/4}, and infinite at or above one quarter.
#' @param x A \code{GPDDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, GPDDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(2)]))
  xi <- rep(theta[[2]], length.out = n)
  ifelse(xi < 0.25,
         3 * (1 - 2 * xi) * (2 * xi^2 + xi + 3) /
           ((1 - 3 * xi) * (1 - 4 * xi)) - 3,
         Inf)
}

# --- gengamma1 -------------------------------------------------------------
#
# Every raw moment is a ratio of gamma functions,
#   E[Y^k] = a^k Gamma((d+k)/p) / Gamma(d/p),
# so the four central moments follow from four of them. The ratio is formed
# through lgamma, since Gamma((d+k)/p) overflows long before the ratio does.

#' Raw Moments of the Generalized Gamma
#'
#' @description
#' \eqn{E[Y^k] = a^k\,\Gamma((d+k)/p)/\Gamma(d/p)} for \eqn{k = 1, \dots, 4},
#' formed on the log scale.
#'
#' @param a,d,p The Stacy parameters.
#'
#' @return A list of the four raw moments.
#'
#' @keywords internal
gengamma_raw_moments <- function(a, d, p) {
  k0 <- d / p
  lapply(1:4, function(k) exp(k * log(a) + lgamma(k0 + k / p) - lgamma(k0)))
}

#' @title Mean of the Generalized Gamma Distribution
#' @name mean.GenGamma1Distrib
#' @description Closed form:
#'   \eqn{a\,\Gamma((d+1)/p)/\Gamma(d/p)}.
#' @param x A \code{GenGamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])[[1]]
}

#' @title Variance of the Generalized Gamma Distribution
#' @name variance.GenGamma1Distrib
#' @description Closed form, from the first two raw moments.
#' @param x A \code{GenGamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])
  m[[2]] - m[[1]]^2
}

#' @title Skewness of the Generalized Gamma Distribution
#' @name skewness.GenGamma1Distrib
#' @description Closed form, from the first three raw moments.
#' @param x A \code{GenGamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])
  v <- m[[2]] - m[[1]]^2
  (m[[3]] - 3 * m[[1]] * m[[2]] + 2 * m[[1]]^3) / v^1.5
}

#' @title Kurtosis of the Generalized Gamma Distribution
#' @name kurtosis.GenGamma1Distrib
#' @description Closed form, from the first four raw moments; excess.
#' @param x A \code{GenGamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])
  v <- m[[2]] - m[[1]]^2
  (m[[4]] - 4 * m[[1]] * m[[3]] + 6 * m[[1]]^2 * m[[2]] - 3 * m[[1]]^4) / v^2 - 3
}

# --- poisson ---------------------------------------------------------------

#' @title Mean of the Poisson Distribution
#' @name mean.PoissonDistrib
#' @description Closed form: \eqn{\mu}.
#' @param x A \code{PoissonDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, PoissonDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Poisson Distribution
#' @name variance.PoissonDistrib
#' @description Closed form: \eqn{\mu}, equal to the mean.
#' @param x A \code{PoissonDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, PoissonDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Skewness of the Poisson Distribution
#' @name skewness.PoissonDistrib
#' @description Closed form: \eqn{1/\sqrt{\mu}}.
#' @param x A \code{PoissonDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, PoissonDistrib) <- function(x, theta, ...) {
  1 / sqrt(align_theta(x, theta)[[1]])
}

#' @title Kurtosis of the Poisson Distribution
#' @name kurtosis.PoissonDistrib
#' @description Closed form: excess \eqn{1/\mu}.
#' @param x A \code{PoissonDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, PoissonDistrib) <- function(x, theta, ...) {
  1 / align_theta(x, theta)[[1]]
}

# --- bernoulli -------------------------------------------------------------

#' @title Mean of the Bernoulli Distribution
#' @name mean.BernoulliDistrib
#' @description Closed form: \eqn{\mu}, the success probability.
#' @param x A \code{BernoulliDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, BernoulliDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Bernoulli Distribution
#' @name variance.BernoulliDistrib
#' @description Closed form: \eqn{\mu(1-\mu)}.
#' @param x A \code{BernoulliDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, BernoulliDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  p * (1 - p)
}

#' @title Skewness of the Bernoulli Distribution
#' @name skewness.BernoulliDistrib
#' @description Closed form: \eqn{(1-2\mu)/\sqrt{\mu(1-\mu)}}.
#' @param x A \code{BernoulliDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, BernoulliDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  (1 - 2 * p) / sqrt(p * (1 - p))
}

#' @title Kurtosis of the Bernoulli Distribution
#' @name kurtosis.BernoulliDistrib
#' @description Closed form: excess
#'   \eqn{\{1 - 6\mu(1-\mu)\}/\{\mu(1-\mu)\}}.
#' @param x A \code{BernoulliDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, BernoulliDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  v <- p * (1 - p)
  (1 - 6 * v) / v
}

# --- binomial --------------------------------------------------------------

#' @title Mean of the Binomial Distribution
#' @name mean.BinomialDistrib
#' @description Closed form: \eqn{n\mu}.
#' @param x A \code{BinomialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, BinomialDistrib) <- function(x, theta, ...) {
  x@size * align_theta(x, theta)[[1]]
}

#' @title Variance of the Binomial Distribution
#' @name variance.BinomialDistrib
#' @description Closed form: \eqn{n\mu(1-\mu)}.
#' @param x A \code{BinomialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, BinomialDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  x@size * p * (1 - p)
}

#' @title Skewness of the Binomial Distribution
#' @name skewness.BinomialDistrib
#' @description Closed form: \eqn{(1-2\mu)/\sqrt{n\mu(1-\mu)}}.
#' @param x A \code{BinomialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, BinomialDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  (1 - 2 * p) / sqrt(x@size * p * (1 - p))
}

#' @title Kurtosis of the Binomial Distribution
#' @name kurtosis.BinomialDistrib
#' @description Closed form: excess
#'   \eqn{\{1 - 6\mu(1-\mu)\}/\{n\mu(1-\mu)\}}.
#' @param x A \code{BinomialDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, BinomialDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  v <- p * (1 - p)
  (1 - 6 * v) / (x@size * v)
}

# --- geometric -------------------------------------------------------------
#
# Counting failures before the first success, with p = 1/(1 + mu).

#' @title Mean of the Geometric Distribution
#' @name mean.GeometricDistrib
#' @description Closed form: \eqn{\mu}, the parameter itself.
#' @param x A \code{GeometricDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, GeometricDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Geometric Distribution
#' @name variance.GeometricDistrib
#' @description Closed form: \eqn{\mu(1+\mu)}.
#' @param x A \code{GeometricDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, GeometricDistrib) <- function(x, theta, ...) {
  mu <- align_theta(x, theta)[[1]]
  mu * (1 + mu)
}

#' @title Skewness of the Geometric Distribution
#' @name skewness.GeometricDistrib
#' @description Closed form: \eqn{(1+2\mu)/\sqrt{\mu(1+\mu)}}.
#' @param x A \code{GeometricDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, GeometricDistrib) <- function(x, theta, ...) {
  mu <- align_theta(x, theta)[[1]]
  (1 + 2 * mu) / sqrt(mu * (1 + mu))
}

#' @title Kurtosis of the Geometric Distribution
#' @name kurtosis.GeometricDistrib
#' @description Closed form: excess \eqn{6 + 1/\{\mu(1+\mu)\}}.
#' @param x A \code{GeometricDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, GeometricDistrib) <- function(x, theta, ...) {
  mu <- align_theta(x, theta)[[1]]
  6 + 1 / (mu * (1 + mu))
}

# --- negbin1 ---------------------------------------------------------------
#
# The variance is linear in the mean, mu(1 + theta), which is what fixes the
# size at mu/theta and the success probability at 1/(1 + theta). Every moment
# below reduces to the Poisson one as theta goes to zero.

#' @title Mean of the NB1 Negative Binomial Distribution
#' @name mean.NegBin1Distrib
#' @description Closed form: \eqn{\mu}.
#' @param x A \code{NegBin1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the NB1 Negative Binomial Distribution
#' @name variance.NegBin1Distrib
#' @description Closed form: \eqn{\mu(1+\theta)}, linear in the mean.
#' @param x A \code{NegBin1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] * (1 + theta[[2]])
}

#' @title Skewness of the NB1 Negative Binomial Distribution
#' @name skewness.NegBin1Distrib
#' @description Closed form: \eqn{(1+2\theta)/\sqrt{\mu(1+\theta)}}, which
#'   tends to the Poisson \eqn{1/\sqrt{\mu}} as \eqn{\theta \to 0}.
#' @param x A \code{NegBin1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  (1 + 2 * th) / sqrt(mu * (1 + th))
}

#' @title Kurtosis of the NB1 Negative Binomial Distribution
#' @name kurtosis.NegBin1Distrib
#' @description Closed form: excess
#'   \eqn{6\theta/\mu + 1/\{\mu(1+\theta)\}}, which tends to the Poisson
#'   \eqn{1/\mu} as \eqn{\theta \to 0}.
#' @param x A \code{NegBin1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  6 * th / mu + 1 / (mu * (1 + th))
}

# --- betabinom1 ------------------------------------------------------------
#
# The shapes are a = mu/sigma and b = (1 - mu)/sigma. The central moments come
# from the falling factorial moments,
#   E[Y (Y-1) ... (Y-k+1)] = n^(k) prod_{j<k} (a+j)/(a+b+j),
# which are one product each, rather than from a transcribed quartic.

#' Falling Factorial Moments of the Beta-Binomial
#'
#' @description
#' \eqn{E[Y^{(k)}] = n^{(k)}\prod_{j=0}^{k-1}(a+j)/(a+b+j)} for
#' \eqn{k = 1, \dots, 4}, from which the raw and then the central moments
#' follow.
#'
#' @param a,b The beta shapes.
#' @param n The number of trials.
#'
#' @return A list of the four factorial moments.
#'
#' @keywords internal
betabinom_factorial_moments <- function(a, b, n) {
  lapply(1:4, function(k) {
    j <- seq_len(k) - 1L
    nf <- prod(n - j)
    r <- 1
    for (i in j) r <- r * (a + i) / (a + b + i)
    nf * r
  })
}

#' Central Moments From Falling Factorial Moments
#'
#' @description
#' Converts \eqn{E[Y^{(k)}]} to the first four central moments, through the
#' raw moments \eqn{m_1 = f_1}, \eqn{m_2 = f_2 + f_1},
#' \eqn{m_3 = f_3 + 3f_2 + f_1} and \eqn{m_4 = f_4 + 6f_3 + 7f_2 + f_1}.
#'
#' @param f A list of the four falling factorial moments.
#'
#' @return A list with the mean and the second, third and fourth central
#'   moments.
#'
#' @keywords internal
central_from_factorial <- function(f) {
  m1 <- f[[1]]
  m2 <- f[[2]] + f[[1]]
  m3 <- f[[3]] + 3 * f[[2]] + f[[1]]
  m4 <- f[[4]] + 6 * f[[3]] + 7 * f[[2]] + f[[1]]
  list(mean = m1,
       c2 = m2 - m1^2,
       c3 = m3 - 3 * m1 * m2 + 2 * m1^3,
       c4 = m4 - 4 * m1 * m3 + 6 * m1^2 * m2 - 3 * m1^4)
}

#' Central Moments of a Beta-Binomial
#'
#' @description
#' The mean and the second, third and fourth central moments, vectorized over
#' the parameters.
#'
#' @param mu,sigma The mean proportion and the dispersion.
#' @param n The number of trials.
#'
#' @return A list of four numeric vectors.
#'
#' @keywords internal
betabinom_central <- function(mu, sigma, n) {
  len <- max(length(mu), length(sigma))
  mu <- rep(mu, length.out = len)
  sigma <- rep(sigma, length.out = len)
  out <- lapply(seq_len(len), function(i) {
    central_from_factorial(
      betabinom_factorial_moments(mu[i] / sigma[i], (1 - mu[i]) / sigma[i], n)
    )
  })
  list(mean = vapply(out, function(z) z$mean, numeric(1)),
       c2 = vapply(out, function(z) z$c2, numeric(1)),
       c3 = vapply(out, function(z) z$c3, numeric(1)),
       c4 = vapply(out, function(z) z$c4, numeric(1)))
}

#' @title Mean of the Beta-Binomial Distribution
#' @name mean.BetaBinom1Distrib
#' @description Closed form: \eqn{n\mu}, the shapes canceling.
#' @param x A \code{BetaBinom1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  betabinom_central(theta[[1]], theta[[2]], x@size)$mean
}

#' @title Variance of the Beta-Binomial Distribution
#' @name variance.BetaBinom1Distrib
#' @description
#' Closed form: \eqn{n\mu(1-\mu)(1+n\sigma)/(1+\sigma)}, the binomial variance
#' inflated by the dispersion.
#' @param x A \code{BetaBinom1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  betabinom_central(theta[[1]], theta[[2]], x@size)$c2
}

#' @title Skewness of the Beta-Binomial Distribution
#' @name skewness.BetaBinom1Distrib
#' @description Closed form, from the falling factorial moments.
#' @param x A \code{BetaBinom1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- betabinom_central(theta[[1]], theta[[2]], x@size)
  m$c3 / m$c2^1.5
}

#' @title Kurtosis of the Beta-Binomial Distribution
#' @name kurtosis.BetaBinom1Distrib
#' @description Closed form, from the falling factorial moments; excess.
#' @param x A \code{BetaBinom1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- betabinom_central(theta[[1]], theta[[2]], x@size)
  m$c4 / m$c2^2 - 3
}

# --- gaussian2, gaussian3, gamma1 ------------------------------------------
#
# Three reparametrizations of two families already here, so each moment is the
# same number the twin reports, written in the coordinates of this one.

#' @title Mean of the Gaussian in Mean and Variance
#' @name mean.Gaussian2Distrib
#' @description Closed form: \eqn{\mu}.
#' @param x A \code{Gaussian2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Gaussian2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gaussian in Mean and Variance
#' @name variance.Gaussian2Distrib
#' @description Closed form: \eqn{\sigma^2}, the parameter itself.
#' @param x A \code{Gaussian2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Gaussian2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[2]]
}

#' @title Skewness of the Gaussian in Mean and Variance
#' @name skewness.Gaussian2Distrib
#' @description Closed form: zero.
#' @param x A \code{Gaussian2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Gaussian2Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Kurtosis of the Gaussian in Mean and Variance
#' @name kurtosis.Gaussian2Distrib
#' @description Closed form: zero excess.
#' @param x A \code{Gaussian2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Gaussian2Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Mean of the Gaussian in Mean and Precision
#' @name mean.Gaussian3Distrib
#' @description Closed form: \eqn{\mu}.
#' @param x A \code{Gaussian3Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Gaussian3Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gaussian in Mean and Precision
#' @name variance.Gaussian3Distrib
#' @description Closed form: \eqn{1/\tau}.
#' @param x A \code{Gaussian3Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Gaussian3Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  1 / theta[[2]] + moment_const(theta, 2L, 0)
}

#' @title Skewness of the Gaussian in Mean and Precision
#' @name skewness.Gaussian3Distrib
#' @description Closed form: zero.
#' @param x A \code{Gaussian3Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Gaussian3Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Kurtosis of the Gaussian in Mean and Precision
#' @name kurtosis.Gaussian3Distrib
#' @description Closed form: zero excess.
#' @param x A \code{Gaussian3Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Gaussian3Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Mean of the Gamma in Mean and Dispersion
#' @name mean.Gamma1Distrib
#' @description Closed form: \eqn{\mu}, the parameter itself.
#' @param x A \code{Gamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Gamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gamma in Mean and Dispersion
#' @name variance.Gamma1Distrib
#' @description Closed form: \eqn{\phi\mu^2}, which is what the
#'   parametrization is defined by.
#' @param x A \code{Gamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Gamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]] * theta[[1]]^2
}

#' @title Skewness of the Gamma in Mean and Dispersion
#' @name skewness.Gamma1Distrib
#' @description Closed form: \eqn{2\sqrt{\phi}}, free of the mean.
#' @param x A \code{Gamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Gamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  2 * sqrt(theta[[2]]) + moment_const(theta, 2L, 0)
}

#' @title Kurtosis of the Gamma in Mean and Dispersion
#' @name kurtosis.Gamma1Distrib
#' @description Closed form: excess \eqn{6\phi}, free of the mean.
#' @param x A \code{Gamma1Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Gamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  6 * theta[[2]] + moment_const(theta, 2L, 0)
}

# --- invgauss2, beta2, betabinom2 ------------------------------------------

#' @title Mean of the Inverse Gaussian in Mean and Shape
#' @name mean.InvGauss2Distrib
#' @description Closed form: \eqn{\mu}.
#' @param x An \code{InvGauss2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, InvGauss2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Inverse Gaussian in Mean and Shape
#' @name variance.InvGauss2Distrib
#' @description Closed form: \eqn{\mu^3/\lambda}.
#' @param x An \code{InvGauss2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, InvGauss2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]]^3 / theta[[2]]
}

#' @title Skewness of the Inverse Gaussian in Mean and Shape
#' @name skewness.InvGauss2Distrib
#' @description Closed form: \eqn{3\sqrt{\mu/\lambda}}.
#' @param x An \code{InvGauss2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, InvGauss2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  3 * sqrt(theta[[1]] / theta[[2]])
}

#' @title Kurtosis of the Inverse Gaussian in Mean and Shape
#' @name kurtosis.InvGauss2Distrib
#' @description Closed form: excess \eqn{15\mu/\lambda}.
#' @param x An \code{InvGauss2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, InvGauss2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  15 * theta[[1]] / theta[[2]]
}

#' @title Mean of the Beta in Its Shapes
#' @name mean.Beta2Distrib
#' @description Closed form: \eqn{\alpha/(\alpha+\beta)}.
#' @param x A \code{Beta2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Beta2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] / (theta[[1]] + theta[[2]])
}

#' @title Variance of the Beta in Its Shapes
#' @name variance.Beta2Distrib
#' @description Closed form:
#'   \eqn{\alpha\beta/\{(\alpha+\beta)^2(\alpha+\beta+1)\}}.
#' @param x A \code{Beta2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, Beta2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  a <- theta[[1]]
  b <- theta[[2]]
  a * b / ((a + b)^2 * (a + b + 1))
}

#' @title Skewness of the Beta in Its Shapes
#' @name skewness.Beta2Distrib
#' @description Closed form:
#'   \eqn{2(\beta-\alpha)\sqrt{\alpha+\beta+1}/
#'        \{(\alpha+\beta+2)\sqrt{\alpha\beta}\}}.
#' @param x A \code{Beta2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, Beta2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  a <- theta[[1]]
  b <- theta[[2]]
  2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b))
}

#' @title Kurtosis of the Beta in Its Shapes
#' @name kurtosis.Beta2Distrib
#' @description Closed form; excess.
#' @param x A \code{Beta2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, Beta2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  a <- theta[[1]]
  b <- theta[[2]]
  6 * ((a - b)^2 * (a + b + 1) - a * b * (a + b + 2)) /
    (a * b * (a + b + 2) * (a + b + 3))
}

#' @title Mean of the Beta-Binomial in Its Shapes
#' @name mean.BetaBinom2Distrib
#' @description Closed form: \eqn{n\alpha/(\alpha+\beta)}.
#' @param x A \code{BetaBinom2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, BetaBinom2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  betabinom_central(theta[[1]] / (theta[[1]] + theta[[2]]),
                    1 / (theta[[1]] + theta[[2]]), x@size)$mean
}

#' @title Variance of the Beta-Binomial in Its Shapes
#' @name variance.BetaBinom2Distrib
#' @description Closed form, from the falling factorial moments.
#' @param x A \code{BetaBinom2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, BetaBinom2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  betabinom_central(theta[[1]] / (theta[[1]] + theta[[2]]),
                    1 / (theta[[1]] + theta[[2]]), x@size)$c2
}

#' @title Skewness of the Beta-Binomial in Its Shapes
#' @name skewness.BetaBinom2Distrib
#' @description Closed form, from the falling factorial moments.
#' @param x A \code{BetaBinom2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, BetaBinom2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- betabinom_central(theta[[1]] / (theta[[1]] + theta[[2]]),
                         1 / (theta[[1]] + theta[[2]]), x@size)
  m$c3 / m$c2^1.5
}

#' @title Kurtosis of the Beta-Binomial in Its Shapes
#' @name kurtosis.BetaBinom2Distrib
#' @description Closed form, from the falling factorial moments; excess.
#' @param x A \code{BetaBinom2Distrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, BetaBinom2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- betabinom_central(theta[[1]] / (theta[[1]] + theta[[2]]),
                         1 / (theta[[1]] + theta[[2]]), x@size)
  m$c4 / m$c2^2 - 3
}


# --- poisson-inverse gaussian ------------------------------------------------

#' @title Closed Moments of the Poisson-Inverse Gaussian
#' @name moments.Pig1Distrib
#' @description The cumulants of the mixture are those of the Poisson plus
#' those the inverse Gaussian rate contributes:
#' \eqn{\kappa_1 = \mu}, \eqn{\kappa_2 = \mu + \sigma\mu^2},
#' \eqn{\kappa_3 = \mu + 3\sigma\mu^2 + 3\sigma^2\mu^3},
#' \eqn{\kappa_4 = \mu + 7\sigma\mu^2 + 18\sigma^2\mu^3 + 15\sigma^3\mu^4};
#' skewness and excess kurtosis are the standardized ratios.
#' @param x A \code{Pig1Distrib} object.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Pig1Distrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @rdname moments.Pig1Distrib
#' @name variance.Pig1Distrib
S7::method(variance, Pig1Distrib) <- function(x, theta, ...) {
  th <- align_theta(x, theta)
  th[[1]] + th[[2]] * th[[1]]^2
}

#' @rdname moments.Pig1Distrib
#' @name skewness.Pig1Distrib
S7::method(skewness, Pig1Distrib) <- function(x, theta, ...) {
  th <- align_theta(x, theta)
  mu <- th[[1]]; s <- th[[2]]
  k2 <- mu + s * mu^2
  k3 <- mu + 3 * s * mu^2 + 3 * s^2 * mu^3
  k3 / k2^1.5
}

#' @rdname moments.Pig1Distrib
#' @name kurtosis.Pig1Distrib
S7::method(kurtosis, Pig1Distrib) <- function(x, theta, ...) {
  th <- align_theta(x, theta)
  mu <- th[[1]]; s <- th[[2]]
  k2 <- mu + s * mu^2
  k4 <- mu + 7 * s * mu^2 + 18 * s^2 * mu^3 + 15 * s^3 * mu^4
  k4 / k2^2
}

#' @title Closed Moments of the Orthogonal Poisson-Inverse Gaussian
#' @name moments.Pig2Distrib
#' @description The cumulants of \code{\link[=moments.Pig1Distrib]{pig1}}
#' at the dispersion \code{\link{pig2_sigma}} implies.
#' @param x A \code{Pig2Distrib} object.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, Pig2Distrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @rdname moments.Pig2Distrib
#' @name variance.Pig2Distrib
S7::method(variance, Pig2Distrib) <- function(x, theta, ...) {
  th <- align_theta(x, theta)
  mu <- th[[1]]; s <- pig2_sigma(mu, th[[2]])
  mu + s * mu^2
}

#' @rdname moments.Pig2Distrib
#' @name skewness.Pig2Distrib
S7::method(skewness, Pig2Distrib) <- function(x, theta, ...) {
  th <- align_theta(x, theta)
  mu <- th[[1]]; s <- pig2_sigma(mu, th[[2]])
  k2 <- mu + s * mu^2
  (mu + 3 * s * mu^2 + 3 * s^2 * mu^3) / k2^1.5
}

#' @rdname moments.Pig2Distrib
#' @name kurtosis.Pig2Distrib
S7::method(kurtosis, Pig2Distrib) <- function(x, theta, ...) {
  th <- align_theta(x, theta)
  mu <- th[[1]]; s <- pig2_sigma(mu, th[[2]])
  k2 <- mu + s * mu^2
  (mu + 7 * s * mu^2 + 18 * s^2 * mu^3 + 15 * s^3 * mu^4) / k2^2
}
