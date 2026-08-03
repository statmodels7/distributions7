#' @include distrib.R generics.R numerical_functions.R negbin_distrib.R pseudohuber_distrib.R laplace_distrib.R

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
#' d <- gaussian_distrib()
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
#' @name mean.NegBinDistrib
#' @description Closed form, replacing the numerical default: \eqn{E[Y] = \mu}.
#' @param x A \code{NegBinDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, NegBinDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Variance of the Negative Binomial Distribution
#' @name variance.NegBinDistrib
#' @description Closed form, replacing the numerical default: \eqn{Var(Y) = \mu + \mu^2/\theta}.
#' @param x A \code{NegBinDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, NegBinDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + theta[[1]]^2 / theta[[2]]
}

#' @title Skewness of the Negative Binomial Distribution
#' @name skewness.NegBinDistrib
#' @description Closed form, replacing the numerical default: \eqn{(\theta + 2\mu)/\sqrt{\mu\theta(\theta+\mu)}}.
#' @param x A \code{NegBinDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, NegBinDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  (th + 2 * mu) / sqrt(mu * th * (th + mu))
}

#' @title Kurtosis of the Negative Binomial Distribution
#' @name kurtosis.NegBinDistrib
#' @description Closed form, replacing the numerical default: \eqn{6/\theta + \theta/(\mu(\theta+\mu))}.
#' @param x A \code{NegBinDistrib}.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, NegBinDistrib) <- function(x, theta, ...) {
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
