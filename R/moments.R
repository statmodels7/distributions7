#' @include distrib.R generics.R numerical_functions.R negbin2_distrib.R pseudohuber_distrib.R laplace_distrib.R laplace2_distrib.R weibull1_distrib.R gumbel_distrib.R skewnormal1_distrib.R skewt_distrib.R gaussian1_distrib.R cauchy_distrib.R logistic_distrib.R student_t1_distrib.R gamma2_distrib.R exponential_distrib.R chisq_distrib.R lognormal1_distrib.R invgauss1_distrib.R beta1_distrib.R gpd_distrib.R gengamma1_distrib.R poisson_distrib.R bernoulli_distrib.R binomial_distrib.R geometric_distrib.R negbin1_distrib.R betabinom1_distrib.R gaussian2_distrib.R gaussian3_distrib.R gamma1_distrib.R invgauss2_distrib.R beta2_distrib.R betabinom2_distrib.R pig1_distrib.R pig2_distrib.R
NULL

#' @title Raw and Central Moments of a Distribution
#' @name moment
#'
#' @description
#' Computes the raw moment \eqn{E[Y^p]} or the central moment
#' \eqn{E[(Y-\mu)^p]} of a distribution at one or more parameter settings. The
#' expectation is taken numerically by [expectation()], which integrates over
#' the support of a continuous family and sums the series of a discrete one.
#' The order \eqn{p} need not be a whole number, so fractional and absolute
#' moments are reachable from the same function. Every setting supplied in
#' `theta` is evaluated in one batched call, and the result carries one value
#' per setting.
#'
#' @details
#' # What is computed
#'
#' With `central = FALSE`,
#' \deqn{m_p = \int y^p f(y \mid \theta)\, \mathrm{d}y
#'       \quad\text{or}\quad \sum_y y^p\, f(y \mid \theta),}
#' and with `central = TRUE` the same integral or sum with \eqn{y} replaced by
#' \eqn{y - \mu}. The centering value \eqn{\mu} is the mean, obtained by a
#' first call at `p = 1`, unless it is supplied through the `mu` argument.
#' Supplying `mu = 0` with `central = TRUE` therefore returns the raw moment,
#' and supplying a fitted or a theoretical mean saves one pass.
#'
#' # Accuracy and cost
#'
#' A continuous family's moment is a quadrature and a discrete family's is a
#' series, so neither is exact. On a Gaussian the second and fourth central
#' moments come back at \eqn{9 - 2.8\times10^{-14}} and
#' \eqn{3 - 8.9\times10^{-15}} against the exact 9 and 3; on a Poisson the
#' series is exact to the last bit at ordinary means. The price is the
#' evaluation: one numerical variance of a Gaussian costs 2.3 ms against
#' 24 microseconds for the closed form the family registers, a factor of about
#' 94, and one numerical skewness costs 3.5 ms against the same 24
#' microseconds. That gap is why 43 of the 45 shipped families answer
#' [variance()] and [skewness()] with a formula of their own; the two von Mises
#' families are the ones that reach this function.
#'
#' # A moment that does not exist
#'
#' A divergent integral does not announce itself. The quadrature returns
#' whatever its truncation gives, a number that moves with the panel layout and
#' looks like an estimate. A family whose moments fail to exist therefore
#' registers a method that returns `NaN` directly, as the Cauchy does through
#' [mean.CauchyDistrib()], and a family whose moments exist only above a
#' threshold returns `Inf` where they do not, as the Student t does at
#' \eqn{\nu \le 4} for the kurtosis. Reading a number back from this function on
#' a family with no analytical method is a statement about the quadrature and
#' not about the law.
#'
#' @section Notation:
#' \eqn{Y} is the response, \eqn{f(y \mid \theta)} its density or mass function,
#' \eqn{\theta} the parameter on its own scale, and \eqn{\mu = E[Y]} the mean.
#'
#' @param distrib An object inheriting from `distrib`, from any of the family
#'   constructors such as [gaussian1_distrib()] or [poisson_distrib()].
#' @param theta A named list of parameters on the parameter scale, with one
#'   component per parameter of `distrib`. Components may be vectors and are
#'   recycled against one another and against `p`, so several settings are
#'   evaluated in one call. The list is aligned and validated by name, so a
#'   missing or out-of-bounds component throws.
#' @param p The order of the moment. A numeric vector of length 1 or of the
#'   number of settings, recycled against the components of `theta`. Defaults
#'   to 1, the mean. Non-integer orders are accepted; a negative order is
#'   accepted too and diverges for any family whose support reaches zero.
#' @param central Should the moment be taken about the mean? A single logical,
#'   `FALSE` by default, which gives the raw moment and needs no extra pass.
#'   `TRUE` costs one further evaluation unless `mu` is supplied.
#' @param mu The centering value used when `central = TRUE`. A numeric vector
#'   of length 1 or of the number of settings, or `NULL` (the default), which
#'   computes the mean numerically. Read only when `central = TRUE`.
#' @param ... Passed to [expectation()], and from there to the quadrature or
#'   the series. Names here must not collide with the names in `theta`.
#'
#' @return A numeric vector of moments, of length equal to the longest of the
#'   components of `theta` and of `p`. `NaN` where the integrand is not
#'   defined, `Inf` where the integral diverges to infinity.
#'
#' @seealso [expectation()] for the quadrature and the series this rests on;
#'   [mean.distrib()], [variance()], [std_dev()], [skewness()] and [kurtosis()]
#'   for the four standard moments, each of which prefers a family's closed
#'   form when one is registered.
#'
#' @examples
#' d <- gaussian1_distrib()
#'
#' # The first four moments of a Gaussian, raw and central.
#' moment(d, list(mu = 2, sigma = 3), p = 1)                  # 2
#' moment(d, list(mu = 2, sigma = 3), p = 2)                  # mu^2 + sigma^2
#' moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE)  # sigma^2
#'
#' # p is recycled against theta, so one call covers a grid of orders.
#' moment(d, list(mu = 0, sigma = 1), p = 1:4, central = TRUE)   # 0, 1, 0, 3
#'
#' # Centering at zero recovers the raw moment.
#' all.equal(moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE, mu = 0),
#'           moment(d, list(mu = 2, sigma = 3), p = 2))
#'
#' # One value per parameter setting.
#' moment(d, list(mu = c(0, 1, 2), sigma = 1), p = 2, central = TRUE)
#'
#' # On a discrete family the expectation is an exact sum.
#' all.equal(moment(poisson_distrib(), list(mu = 3), p = 2, central = TRUE), 3)
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
#'
#' @description
#' Evaluates \eqn{E[Y]} by quadrature over the support of a continuous family
#' or by summing the series of a discrete one, through one call of [moment()]
#' at `p = 1`. This is the fallback that runs when the family has registered no
#' closed form. Of the 45 shipped families only the two von Mises reach it; the
#' rest carry a formula, which is a hundredfold cheaper and exact.
#'
#' @param x An object inheriting from `distrib`.
#' @param theta A named list of parameters on the parameter scale, one
#'   component per parameter of `x`. Components may be vectors, in which case
#'   one mean is returned per setting.
#' @param ... Passed to [moment()] and from there to [expectation()].
#'
#' @return A numeric vector of means, of length equal to the longest component
#'   of `theta`. `NaN` or `Inf` where the integral does not converge, with no
#'   warning: a family whose mean does not exist registers its own method
#'   instead of relying on this one.
#'
#' @seealso [moment()] for the quadrature, [mean.CauchyDistrib()] for a family
#'   that overrides this method because its mean does not exist,
#'   [variance()], [skewness()] and [kurtosis()] for the higher moments.
#'
#' @examples
#' # The von Mises is one of the two families with no closed-form mean, so this
#' # method is what answers for it. On (-pi, pi] with mu = 0 the mean is 0.
#' round(mean(vonmises1_distrib(), list(mu = 0, kappa = 2)), 12)
#'
#' # A family with a closed form never reaches here; the numbers agree anyway.
#' all.equal(moment(gaussian1_distrib(), list(mu = 2, sigma = 3), p = 1), 2)
#'
#' @keywords internal
S7::method(mean, distrib) <- function(x, theta, ...) {
  moment(x, theta, p = 1, central = FALSE, ...)
}

#' Variance of a Distribution or Sample
#'
#' @description
#' Computes \eqn{\operatorname{Var}(Y) = E[(Y - E[Y])^2]} for a distribution
#' object, or the sample variance for a numeric vector. Dispatch is on the
#' first argument: a `distrib` uses the family's closed form where one is
#' registered and a quadrature or series otherwise, and a numeric vector is
#' passed to [stats::var()]. Where the variance does not exist the value is
#' `NaN` or `Inf`, never a truncated quadrature.
#'
#' @details
#' # The two routes
#'
#' \deqn{\operatorname{Var}(Y) = \mathbb{E}\left[(Y - \mathbb{E}[Y])^{2}\right].}
#'
#' 43 of the 45 shipped families register a closed form, so the second central
#' moment is a formula in the parameters and costs about 24 microseconds; the
#' two von Mises families fall through to [moment()], which is a quadrature and
#' costs about 12 ms. The two routes agree to the rounding of the quadrature,
#' 3e-15 relative on a Gaussian.
#'
#' # What the sample method returns
#'
#' On a numeric vector the value is [stats::var()], which divides by
#' \eqn{n - 1}. That is the unbiased estimator, and it differs from the
#' convention [skewness()] and [kurtosis()] use for their sample versions,
#' which divide by \eqn{n}. The difference is deliberate: both follow the
#' commonest convention for the quantity in question.
#'
#' @param x An object inheriting from `distrib`, or a numeric vector.
#' @param ... For a `distrib`: `theta`, a named list of parameters, followed by
#'   any further arguments for [moment()]. For a numeric vector: `na.rm`, a
#'   single logical, `FALSE` by default.
#'
#' @return A numeric vector for a `distrib`, one value per parameter setting; a
#'   single number for a numeric vector. `NaN` for a family with no moments,
#'   `Inf` for one whose variance diverges.
#'
#' @seealso [std_dev()] for its square root, [skewness()] and [kurtosis()] for
#'   the standardized third and fourth moments, [moment()] for the numerical
#'   route, [expectation()] for the quadrature and the series.
#'
#' @examples
#' # A closed form, one value per setting.
#' variance(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))
#'
#' # A Poisson is equidispersed: the variance is the mean.
#' variance(poisson_distrib(), list(mu = 3))
#'
#' # A Student t has a variance only above two degrees of freedom.
#' variance(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(3, 10)))
#'
#' # On a numeric vector this is stats::var, with the n - 1 denominator.
#' set.seed(1)
#' y <- rnorm(50)
#' all.equal(variance(y), var(y))
#'
#' @export
variance <- S7::new_generic("variance", "x")

#' @title Variance of a Distribution
#' @name variance.distrib
#'
#' @description
#' Evaluates the second central moment numerically: one call of [moment()] at
#' `p = 1` gives the mean, and a second at `p = 2` takes the central moment
#' about it. Passing the mean forward keeps the second pass to a single
#' quadrature. This is the fallback for a family with no closed form, which
#' among the shipped families means the two von Mises.
#'
#' @param x A `distrib` object.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, giving one variance per setting.
#' @param ... Passed to [moment()] and from there to [expectation()].
#'
#' @return A numeric vector of variances, of length equal to the longest
#'   component of `theta`.
#'
#' @seealso [variance()] for the generic and the sample version,
#'   [moment()] for the quadrature, [std_dev.distrib()] for its square root.
#'
#' @examples
#' # The von Mises reaches this method; on (-pi, pi] its variance is finite.
#' variance(vonmises1_distrib(), list(mu = 0, kappa = 2))
#'
#' # Two passes of moment() give the same answer as the closed form elsewhere.
#' d <- gaussian1_distrib()
#' m <- moment(d, list(mu = 2, sigma = 3), p = 1)
#' all.equal(moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE, mu = m), 9)
#'
#' @keywords internal
S7::method(variance, distrib) <- function(x, theta, ...) {
  m <- moment(x, theta, p = 1, central = FALSE, ...)
  moment(x, theta, p = 2, central = TRUE, mu = m, ...)
}

#' @title Sample Variance
#' @name variance.numeric
#'
#' @description
#' The sample variance of a numeric vector, delegated to [stats::var()] and so
#' divided by \eqn{n - 1}. The method exists so that [variance()] reads the same
#' on data as on a distribution object; nothing is computed here that
#' [stats::var()] does not compute.
#'
#' @param x A numeric vector. A vector of length 1 gives `NA`, as
#'   [stats::var()] does, there being no degree of freedom left.
#' @param na.rm Should missing values be dropped before the variance is taken?
#'   A single logical, `FALSE` by default, which propagates `NA`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number, `NA` if `x` has fewer than two non-missing values.
#'
#' @seealso [variance()] for the generic, [std_dev.numeric()] for the square
#'   root, [skewness.numeric()] and [kurtosis.numeric()], which use the
#'   \eqn{n} denominator instead.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(50)
#' all.equal(variance(y), var(y))
#'
#' # The n - 1 denominator, not n.
#' all.equal(variance(y), sum((y - mean(y))^2) / (length(y) - 1))
#'
#' # Missing values propagate unless dropped.
#' c(variance(c(y, NA)), variance(c(y, NA), na.rm = TRUE))
#'
#' @keywords internal
S7::method(variance, S7::class_numeric) <- function(x, na.rm = FALSE, ...) {
  stats::var(x, na.rm = na.rm)
}

#' Standard Deviation of a Distribution or Sample
#'
#' @description
#' Computes \eqn{\operatorname{sd}(Y) = \sqrt{\operatorname{Var}(Y)}} for a
#' distribution object, or the sample standard deviation for a numeric vector.
#' No family registers a closed form of its own, so a `distrib` always reaches
#' the square root of [variance()], which is where a family's formula is
#' consulted. A numeric vector is passed to [stats::sd()].
#'
#' @details
#' \deqn{\operatorname{sd}(Y) = \sqrt{\operatorname{Var}(Y)}.}
#'
#' The square root is taken after the variance, so accuracy and cost are the
#' variance's: a closed form for 43 of the 45 shipped families, a quadrature
#' for the two von Mises. A family whose variance is `NaN` or `Inf` gives the
#' same here, and a negative variance cannot arise, so the root is always real.
#'
#' On a numeric vector the value is [stats::sd()], the root of the \eqn{n - 1}
#' variance.
#'
#' @param x An object inheriting from `distrib`, or a numeric vector.
#' @param ... For a `distrib`: `theta`, a named list of parameters, followed by
#'   any further arguments for [moment()]. For a numeric vector: `na.rm`, a
#'   single logical, `FALSE` by default.
#'
#' @return A numeric vector for a `distrib`, one value per parameter setting; a
#'   single number for a numeric vector.
#'
#' @seealso [variance()], of which this is the square root; [skewness()] and
#'   [kurtosis()], which standardize by it; [moment()] for the numerical route.
#'
#' @examples
#' # The scale of a Gaussian is its standard deviation, exactly.
#' std_dev(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))
#'
#' # The square root of the variance, on any family.
#' d <- gamma2_distrib()
#' all.equal(std_dev(d, list(mu = 2, sigma2 = 1)),
#'           sqrt(variance(d, list(mu = 2, sigma2 = 1))))
#'
#' # On a numeric vector this is stats::sd.
#' set.seed(1)
#' y <- rnorm(50)
#' all.equal(std_dev(y), sd(y))
#'
#' @export
std_dev <- S7::new_generic("std_dev", "x")

#' @title Standard Deviation of a Distribution
#' @name std_dev.distrib
#'
#' @description
#' The square root of [variance()], evaluated at the same parameters. This is
#' the only route to a distribution's standard deviation: no shipped family
#' registers a method of its own, so whatever a family does for its variance,
#' closed form or quadrature, is what happens here and the root is taken after.
#'
#' @param x A `distrib` object.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, giving one standard deviation per setting.
#' @param ... Passed to [variance()] and from there to [moment()].
#'
#' @return A numeric vector of standard deviations, of length equal to the
#'   longest component of `theta`. `NaN` where the variance is `NaN`, `Inf`
#'   where it is `Inf`.
#'
#' @seealso [std_dev()] for the generic and the sample version,
#'   [variance.distrib()] for the quantity under the root.
#'
#' @examples
#' d <- weibull1_distrib()
#' th <- list(mu = 2, sigma = 3)
#' all.equal(std_dev(d, th), sqrt(variance(d, th)))
#'
#' # A Student t below two degrees of freedom has an infinite variance.
#' std_dev(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(1.5, 5)))
#'
#' @keywords internal
S7::method(std_dev, distrib) <- function(x, theta, ...) {
  sqrt(variance(x, theta, ...))
}

#' @title Sample Standard Deviation
#' @name std_dev.numeric
#'
#' @description
#' The sample standard deviation of a numeric vector, delegated to
#' [stats::sd()] and so the root of the \eqn{n - 1} variance. The method exists
#' so that [std_dev()] reads the same on data as on a distribution object.
#'
#' @param x A numeric vector. A vector of length 1 gives `NA`, as
#'   [stats::sd()] does.
#' @param na.rm Should missing values be dropped before the standard deviation
#'   is taken? A single logical, `FALSE` by default, which propagates `NA`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number, `NA` if `x` has fewer than two non-missing values.
#'
#' @seealso [std_dev()] for the generic, [variance.numeric()] for the square,
#'   [skewness.numeric()] and [kurtosis.numeric()].
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(50)
#' all.equal(std_dev(y), sd(y))
#' all.equal(std_dev(y), sqrt(variance(y)))
#'
#' # Missing values propagate unless dropped.
#' c(std_dev(c(y, NA)), std_dev(c(y, NA), na.rm = TRUE))
#'
#' @keywords internal
S7::method(std_dev, S7::class_numeric) <- function(x, na.rm = FALSE, ...) {
  stats::sd(x, na.rm = na.rm)
}

#' Skewness of a Distribution or Sample
#'
#' @description
#' Computes the skewness, the third standardized central moment
#'
#' \deqn{\gamma_1 = \mathbb{E}\!\left[\left(\frac{Y - \mathbb{E}[Y]}{\operatorname{sd}(Y)}\right)^{3}\right],}
#'
#' for a distribution object, or the sample skewness for a numeric vector. It
#' is zero for a symmetric law, positive for a right tail and negative for a
#' left one. 43 of the 45 shipped families register a closed form; the two von
#' Mises reach a quadrature through [moment()].
#'
#' @details
#' # What the sign and the size mean
#'
#' \eqn{\gamma_1} is invariant to location and to positive scaling, so it is a
#' property of the shape alone: every Gaussian has \eqn{\gamma_1 = 0} and every
#' exponential has \eqn{\gamma_1 = 2}, whatever their parameters. A family with
#' a shape parameter moves along the axis with it, and that is why the quantity
#' is worth reporting: a Poisson at mean \eqn{\mu} has
#' \eqn{\gamma_1 = \mu^{-1/2}}, so its asymmetry vanishes as the counts grow.
#'
#' # When it fails to exist
#'
#' The third moment is needed, so a family with heavy tails may have no
#' skewness even where it has a mean. A Student t has one only above three
#' degrees of freedom and returns `Inf` at or below them; a Cauchy has none at
#' any parameter value and returns `NaN`.
#'
#' # What the sample method returns
#'
#' On a numeric vector the third and second central moments both divide by
#' \eqn{n}, so the value is \eqn{m_3 / m_2^{3/2}} with
#' \eqn{m_k = n^{-1}\sum_i (y_i - \bar y)^k}. This is the population
#' denominator, and it differs from the \eqn{n - 1} convention [variance()]
#' uses on a vector.
#'
#' @param x An object inheriting from `distrib`, or a numeric vector.
#' @param ... For a `distrib`: `theta`, a named list of parameters, followed by
#'   any further arguments for [moment()]. For a numeric vector: `na.rm`, a
#'   single logical, `FALSE` by default.
#'
#' @return A numeric vector for a `distrib`, one value per parameter setting; a
#'   single number for a numeric vector. `NaN` for a family with no moments,
#'   `Inf` where the third moment diverges.
#'
#' @seealso [kurtosis()] for the fourth standardized moment, [variance()] and
#'   [std_dev()] for the second, [moment()] for the numerical route.
#'
#' @examples
#' # Zero for a symmetric family, at every parameter value.
#' skewness(gaussian1_distrib(), list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#'
#' # A shape parameter moves it: a Poisson is mu^(-1/2).
#' all.equal(skewness(poisson_distrib(), list(mu = 4)), 0.5)
#'
#' # A gamma is 2 / sqrt(shape), so it flattens as the shape grows.
#' skewness(gamma2_distrib(), list(mu = 2, sigma2 = c(4, 1, 0.25)))
#'
#' # The sample version uses the n denominator.
#' set.seed(1)
#' y <- rgamma(200, shape = 2)
#' c(sample = skewness(y), theory = 2 / sqrt(2))
#'
#' @export
skewness <- S7::new_generic("skewness", "x")

# distrib method: theta is a named list of parameters
#' @title Skewness of a Distribution
#' @name skewness.distrib
#'
#' @description
#' Evaluates the third standardized central moment numerically. Three calls of
#' [moment()] are made: the mean at `p = 1`, then the second and third central
#' moments about it, and the value is \eqn{m_3 / m_2^{3/2}}. Passing the mean
#' forward keeps the count at three quadratures. This is the fallback for a
#' family with no closed form, which among the shipped families means the two
#' von Mises.
#'
#' @param x A `distrib` object.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, giving one skewness per setting.
#' @param ... Passed to [moment()] and from there to [expectation()].
#'
#' @return A numeric vector, of length equal to the longest component of
#'   `theta`.
#'
#' @seealso [skewness()] for the generic and the sample version,
#'   [moment()] for the quadrature, [kurtosis.distrib()] for the fourth order.
#'
#' @examples
#' # The von Mises reaches this method, and is symmetric about mu.
#' round(skewness(vonmises1_distrib(), list(mu = 0, kappa = 2)), 10)
#'
#' # The same three moments, assembled by hand, on a family that has a formula.
#' d <- gamma2_distrib(); th <- list(mu = 2, sigma2 = 1)
#' m <- moment(d, th, p = 1)
#' m2 <- moment(d, th, p = 2, central = TRUE, mu = m)
#' m3 <- moment(d, th, p = 3, central = TRUE, mu = m)
#' all.equal(m3 / m2^1.5, skewness(d, th))
#'
#' @keywords internal
S7::method(skewness, distrib) <- function(x, theta, ...) {
  m <- moment(x, theta, p = 1, central = FALSE, ...)
  m2 <- moment(x, theta, p = 2, central = TRUE, mu = m, ...)
  m3 <- moment(x, theta, p = 3, central = TRUE, mu = m, ...)
  m3 / m2^1.5
}

#' @title Sample Skewness
#' @name skewness.numeric
#'
#' @description
#' The sample skewness of a numeric vector, \eqn{m_3 / m_2^{3/2}} with
#' \eqn{m_k = n^{-1}\sum_i (y_i - \bar y)^k}. Both central moments divide by
#' \eqn{n}, so this is the population denominator and the estimator is biased
#' towards zero in small samples. [variance.numeric()] uses \eqn{n - 1}
#' instead, each following the commonest convention for its own quantity.
#'
#' @param x A numeric vector. A constant vector gives `NaN`, the standardizing
#'   denominator being zero.
#' @param na.rm Should missing values be dropped before the moments are taken?
#'   A single logical, `FALSE` by default, which propagates `NA`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @seealso [skewness()] for the generic and the distribution version,
#'   [kurtosis.numeric()] for the fourth order, [variance.numeric()], whose
#'   denominator is \eqn{n - 1}.
#'
#' @examples
#' set.seed(1)
#' y <- rgamma(500, shape = 2)
#'
#' # The n denominator, written out.
#' m <- mean(y)
#' all.equal(skewness(y),
#'           mean((y - m)^3) / mean((y - m)^2)^1.5)
#'
#' # It estimates the population value 2 / sqrt(shape).
#' c(sample = skewness(y), theory = 2 / sqrt(2))
#'
#' # Missing values propagate unless dropped.
#' c(skewness(c(y, NA)), skewness(c(y, NA), na.rm = TRUE))
#'
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
#' Computes the excess kurtosis, the fourth standardized central moment less
#' the three a Gaussian has:
#'
#' \deqn{\gamma_2 = \mathbb{E}\!\left[\left(\frac{Y - \mathbb{E}[Y]}{\operatorname{sd}(Y)}\right)^{4}\right] - 3.}
#'
#' The subtraction puts the Gaussian at zero, so the sign reads as a comparison
#' with it: positive for heavier tails and a sharper peak, negative for lighter
#' ones. 42 of the 45 shipped families register a closed form; the elastic net
#' and the two von Mises reach a quadrature through [moment()].
#'
#' @details
#' # Why the three is subtracted
#'
#' The raw fourth standardized moment is 3 for every Gaussian, whatever its
#' mean and scale, so it carries no information about the Gaussian itself.
#' Subtracting it makes the quantity a signed distance from normality and puts
#' the two commonest reference laws at recognizable values: 0 for the Gaussian,
#' 3 for the Laplace, \eqn{6/(\nu-4)} for a Student t.
#'
#' # When it fails to exist
#'
#' The fourth moment is needed, so the threshold is higher than for
#' [skewness()]: a Student t has an excess kurtosis only above four degrees of
#' freedom and returns `Inf` at or below them, and a Cauchy returns `NaN` at
#' every parameter value.
#'
#' # What the sample method returns
#'
#' On a numeric vector the fourth and second central moments both divide by
#' \eqn{n}, giving \eqn{m_4 / m_2^{2} - 3}. This is the population denominator,
#' the same convention [skewness()] uses. [variance()] divides by \eqn{n - 1}
#' instead.
#'
#' @param x An object inheriting from `distrib`, or a numeric vector.
#' @param ... For a `distrib`: `theta`, a named list of parameters, followed by
#'   any further arguments for [moment()]. For a numeric vector: `na.rm`, a
#'   single logical, `FALSE` by default.
#'
#' @return A numeric vector for a `distrib`, one value per parameter setting; a
#'   single number for a numeric vector. `NaN` for a family with no moments,
#'   `Inf` where the fourth moment diverges.
#'
#' @seealso [skewness()] for the third standardized moment, [variance()] and
#'   [std_dev()] for the second, [moment()] for the numerical route.
#'
#' @examples
#' # Zero for every Gaussian, by construction.
#' kurtosis(gaussian1_distrib(), list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#'
#' # Three for the Laplace, which is the textbook heavy-tailed comparison.
#' all.equal(kurtosis(laplace_distrib(), list(mu = 0, sigma = 2)), 3)
#'
#' # A Student t is 6 / (nu - 4), and infinite at or below four.
#' kurtosis(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(3, 5, 10)))
#'
#' # The sample version uses the n denominator.
#' set.seed(1)
#' y <- rt(500, df = 10)
#' c(sample = kurtosis(y), theory = 6 / (10 - 4))
#'
#' @export
kurtosis <- S7::new_generic("kurtosis", "x")

# distrib method: theta is a named list of parameters
#' @title Excess Kurtosis of a Distribution
#' @name kurtosis.distrib
#'
#' @description
#' Evaluates the excess kurtosis numerically. Three calls of [moment()] are
#' made: the mean at `p = 1`, then the second and fourth central moments about
#' it, and the value is \eqn{m_4 / m_2^{2} - 3}. Passing the mean forward keeps
#' the count at three quadratures. This is the fallback for a family with no
#' closed form, which among the shipped families means the elastic net and the
#' two von Mises.
#'
#' @param x A `distrib` object.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, giving one excess kurtosis per setting.
#' @param ... Passed to [moment()] and from there to [expectation()].
#'
#' @return A numeric vector, of length equal to the longest component of
#'   `theta`.
#'
#' @seealso [kurtosis()] for the generic and the sample version,
#'   [moment()] for the quadrature, [skewness.distrib()] for the third order.
#'
#' @examples
#' # The elastic net reaches this method; its excess kurtosis lies between the
#' # Gaussian's 0 and the Laplace's 3.
#' round(kurtosis(enet_distrib(), list(mu = 0, lambda = 1, alpha = 0.5)), 4)
#'
#' # The same three moments, assembled by hand, on a family with a formula.
#' d <- gamma2_distrib(); th <- list(mu = 2, sigma2 = 1)
#' m <- moment(d, th, p = 1)
#' m2 <- moment(d, th, p = 2, central = TRUE, mu = m)
#' m4 <- moment(d, th, p = 4, central = TRUE, mu = m)
#' all.equal(m4 / m2^2 - 3, kurtosis(d, th))
#'
#' @keywords internal
S7::method(kurtosis, distrib) <- function(x, theta, ...) {
  m <- moment(x, theta, p = 1, central = FALSE, ...)
  m2 <- moment(x, theta, p = 2, central = TRUE, mu = m, ...)
  m4 <- moment(x, theta, p = 4, central = TRUE, mu = m, ...)
  m4 / m2^2 - 3
}

#' @title Sample Excess Kurtosis
#' @name kurtosis.numeric
#'
#' @description
#' The sample excess kurtosis of a numeric vector, \eqn{m_4 / m_2^{2} - 3} with
#' \eqn{m_k = n^{-1}\sum_i (y_i - \bar y)^k}. Both central moments divide by
#' \eqn{n}, so this is the population denominator, the convention
#' [skewness.numeric()] also uses. The estimator is biased downwards in small
#' samples, and heavy tails make it slow to settle: it depends on the fourth
#' moment, so a handful of extreme points dominate it.
#'
#' @param x A numeric vector. A constant vector gives `NaN`, the standardizing
#'   denominator being zero.
#' @param na.rm Should missing values be dropped before the moments are taken?
#'   A single logical, `FALSE` by default, which propagates `NA`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @seealso [kurtosis()] for the generic and the distribution version,
#'   [skewness.numeric()] for the third order, [variance.numeric()], whose
#'   denominator is \eqn{n - 1}.
#'
#' @examples
#' set.seed(1)
#' y <- rt(1000, df = 10)
#'
#' # The n denominator, written out.
#' m <- mean(y)
#' all.equal(kurtosis(y), mean((y - m)^4) / mean((y - m)^2)^2 - 3)
#'
#' # It estimates 6 / (nu - 4), and needs a lot of data to get there.
#' c(sample = kurtosis(y), theory = 6 / (10 - 4))
#'
#' # Missing values propagate unless dropped.
#' c(kurtosis(c(y, NA)), kurtosis(c(y, NA), na.rm = TRUE))
#'
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
#'
#' @description
#' Closed form, replacing the numerical default: \eqn{E[Y] = \mu}. This
#' parametrization carries the mean as its first parameter, so the method reads
#' it off and recycles it to the length the two parameters imply. Nothing is
#' integrated and nothing is summed.
#'
#' @param x A `NegBin2Distrib`, from [negbin2_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `theta` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`. Aligned and validated by name, so a missing or out-of-bounds
#'   component throws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$theta))`.
#'
#' @seealso [variance.NegBin2Distrib()], which exceeds this by
#'   \eqn{\mu^2/\theta}; [skewness.NegBin2Distrib()] and
#'   [kurtosis.NegBin2Distrib()]; [negbin2_distrib()] for the family.
#'
#' @examples
#' d <- negbin2_distrib()
#'
#' # The mean is the first parameter, and the dispersion does not move it.
#' mean(d, list(mu = c(1, 2, 3), theta = 2))
#'
#' @keywords internal
S7::method(mean, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Variance of the Negative Binomial Distribution
#' @name variance.NegBin2Distrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \eqn{\operatorname{Var}(Y) = \mu + \mu^2/\theta}. The quadratic term is the
#' overdispersion, so the variance exceeds the mean at every finite \eqn{\theta}
#' and falls back onto it as \eqn{\theta} grows, which is the Poisson limit of
#' the family.
#'
#' @param x A `NegBin2Distrib`, from [negbin2_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `theta` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`. Small `theta` means heavy overdispersion; the variance diverges as it
#'   approaches zero.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$theta))`.
#'
#' @seealso [mean.NegBin2Distrib()], [skewness.NegBin2Distrib()],
#'   [kurtosis.NegBin2Distrib()]; [variance.PoissonDistrib()] for the
#'   equidispersed limit.
#'
#' @examples
#' d <- negbin2_distrib()
#'
#' # Mean 4, dispersion 2: the variance is 4 + 16/2.
#' all.equal(variance(d, list(mu = 4, theta = 2)), 12)
#'
#' # The overdispersion vanishes as theta grows, leaving the Poisson variance.
#' variance(d, list(mu = 4, theta = c(1, 10, 1e4)))
#'
#' @keywords internal
S7::method(variance, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + theta[[1]]^2 / theta[[2]]
}

#' @title Skewness of the Negative Binomial Distribution
#' @name skewness.NegBin2Distrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \eqn{(\theta + 2\mu)/\sqrt{\mu\theta(\theta+\mu)}}. It is positive at every
#' parameter value, a count distribution having a right tail and a floor at
#' zero, and it decreases towards the Poisson value \eqn{\mu^{-1/2}} as
#' \eqn{\theta} grows.
#'
#' @param x A `NegBin2Distrib`, from [negbin2_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `theta` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$theta))`.
#'
#' @seealso [kurtosis.NegBin2Distrib()], [variance.NegBin2Distrib()],
#'   [skewness.PoissonDistrib()] for the limit.
#'
#' @examples
#' d <- negbin2_distrib()
#'
#' # The published form, written out.
#' mu <- 4; th <- 2
#' all.equal(skewness(d, list(mu = mu, theta = th)),
#'           (th + 2 * mu) / sqrt(mu * th * (th + mu)))
#'
#' # It falls onto the Poisson's mu^(-1/2) = 0.5 as the dispersion grows.
#' skewness(d, list(mu = 4, theta = c(1, 100, 1e8)))
#'
#' @keywords internal
S7::method(skewness, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  (th + 2 * mu) / sqrt(mu * th * (th + mu))
}

#' @title Excess Kurtosis of the Negative Binomial Distribution
#' @name kurtosis.NegBin2Distrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \eqn{6/\theta + \theta/\{\mu(\theta+\mu)\}}, the excess over the Gaussian's
#' three. Both terms are positive, so the family is always leptokurtic, and the
#' second is what survives in the Poisson limit, where the expression tends to
#' \eqn{1/\mu}.
#'
#' @param x A `NegBin2Distrib`, from [negbin2_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `theta` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$theta))`.
#'
#' @seealso [skewness.NegBin2Distrib()], [variance.NegBin2Distrib()],
#'   [kurtosis.PoissonDistrib()] for the limit.
#'
#' @examples
#' d <- negbin2_distrib()
#'
#' # The published form, written out.
#' mu <- 4; th <- 2
#' all.equal(kurtosis(d, list(mu = mu, theta = th)),
#'           6 / th + th / (mu * (th + mu)))
#'
#' # The Poisson limit is 1 / mu = 0.25.
#' kurtosis(d, list(mu = 4, theta = c(1, 100, 1e8)))
#'
#' @keywords internal
S7::method(kurtosis, NegBin2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  6 / th + th / (mu * (th + mu))
}

#' @title Mean of the Pseudo-Huber Distribution
#' @name mean.PseudoHuberDistrib
#'
#' @description
#' Closed form, replacing the numerical default: \eqn{E[Y] = \mu}. The density
#' is symmetric about \eqn{\mu} at every scale and shape, so the location
#' parameter is the mean, the median and the mode at once.
#'
#' @param x A `PseudoHuberDistrib`, from [pseudohuber_distrib()].
#' @param theta A named list with components `mu` (the location), `sigma` (the
#'   scale, positive) and `nu` (the shape, positive), each a numeric vector of
#'   length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length equal to the longest of the
#'   three components.
#'
#' @seealso [variance.PseudoHuberDistrib()], which does depend on the shape;
#'   [skewness.PseudoHuberDistrib()], zero by the same symmetry;
#'   [pseudohuber_distrib()] for the family.
#'
#' @examples
#' d <- pseudohuber_distrib()
#'
#' # The location is the mean, whatever the scale and the shape.
#' mean(d, list(mu = c(-1, 0, 3), sigma = 2, nu = 4))
#'
#' @keywords internal
S7::method(mean, PseudoHuberDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(3)])))
}

#' @title Variance of the Pseudo-Huber Distribution
#' @name variance.PseudoHuberDistrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \deqn{\operatorname{Var}(Y) = \sigma^2 \sqrt{\nu}\,
#'       \frac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})},}
#' with \eqn{K_r} the modified Bessel function of the second kind. The scale
#' parameter is therefore not the standard deviation: the shape multiplies it,
#' and the factor grows roughly like \eqn{\sqrt{\nu}} at large \eqn{\nu}.
#'
#' @details
#' The two Bessel functions are evaluated exponentially scaled, through
#' `besselK(sqrt(nu), r, expon.scaled = TRUE)`. Both carry the same factor
#' \eqn{e^{\sqrt{\nu}}}, which cancels in the ratio, so the quotient is exact
#' where the unscaled functions would have overflowed; the terms are degree
#' homogeneous, and the scaling has been checked to \eqn{\nu = 2000}.
#'
#' The quantity is also used inside the package: the bracket
#' [distrib_quantile()] searches over for this family is built from it, so an
#' error here would surface as a failed root search several frames away.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale, \eqn{\nu > 0} the shape, and \eqn{K_r} the
#' modified Bessel function of the second kind of order \eqn{r}.
#'
#' @param x A `PseudoHuberDistrib`, from [pseudohuber_distrib()].
#' @param theta A named list with components `mu` (the location), `sigma` (the
#'   scale, positive) and `nu` (the shape, positive), each a numeric vector of
#'   length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, as long as the longest of `sigma`
#'   and `nu`. The location does not enter the value and does not lengthen it,
#'   so a setting that varies `mu` alone comes back of length 1.
#'
#' @seealso [kurtosis.PseudoHuberDistrib()], which is a ratio of the same
#'   Bessel functions; [mean.PseudoHuberDistrib()]; [pseudohuber_distrib()].
#'
#' @examples
#' d <- pseudohuber_distrib()
#'
#' # The Bessel ratio, written out against the method.
#' nu <- 4
#' all.equal(variance(d, list(mu = 0, sigma = 1, nu = nu)),
#'           sqrt(nu) * besselK(sqrt(nu), 2, expon.scaled = TRUE) /
#'                      besselK(sqrt(nu), 1, expon.scaled = TRUE))
#'
#' # sigma is not the standard deviation: the shape scales it.
#' variance(d, list(mu = 0, sigma = 1, nu = c(1, 4, 100)))
#'
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
#'
#' @description
#' Exactly zero at every parameter value. The density is symmetric about
#' \eqn{\mu}, so every odd central moment vanishes and the third standardized
#' one with it. The constant is returned directly, recycled to the length the
#' parameters imply, so no quadrature is run.
#'
#' @param x A `PseudoHuberDistrib`, from [pseudohuber_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths, but the list is still aligned and validated.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length equal to the longest of the
#'   three components.
#'
#' @seealso [kurtosis.PseudoHuberDistrib()], which is not constant;
#'   [mean.PseudoHuberDistrib()]; [pseudohuber_distrib()].
#'
#' @examples
#' d <- pseudohuber_distrib()
#'
#' # Zero at every shape, by symmetry, with one value per setting.
#' skewness(d, list(mu = 0, sigma = 1, nu = c(1, 2, 3)))
#'
#' @keywords internal
S7::method(skewness, PseudoHuberDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(0, length.out = max(lengths(theta[seq_len(3)])))
}

#' @title Excess Kurtosis of the Pseudo-Huber Distribution
#' @name kurtosis.PseudoHuberDistrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \deqn{\gamma_2 = 3\,\frac{K_3(\sqrt{\nu})\,K_1(\sqrt{\nu})}
#'                          {K_2(\sqrt{\nu})^2} - 3,}
#' with \eqn{K_r} the modified Bessel function of the second kind. It depends
#' on the shape alone, the location and the scale cancelling in a standardized
#' moment, and it interpolates between the two limits the family is built to
#' span: about 3, the Laplace's, as \eqn{\nu \to 0}, and 0, the Gaussian's, as
#' \eqn{\nu} grows.
#'
#' @details
#' All three Bessel functions are evaluated exponentially scaled. Each carries
#' the factor \eqn{e^{\sqrt{\nu}}}, and the ratio is arranged so that the
#' factors cancel exactly, so the expression stays usable at large \eqn{\nu}
#' where the unscaled functions underflow.
#'
#' @section Notation:
#' \eqn{\nu > 0} is the shape and \eqn{K_r} the modified Bessel function of the
#' second kind of order \eqn{r}.
#'
#' @param x A `PseudoHuberDistrib`, from [pseudohuber_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or `n`. Only `nu` enters the value; the other
#'   two are read for their lengths and are still validated.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, as long as `nu`. The location
#'   and the scale do not enter the value and do not lengthen it, so a setting
#'   that varies either alone comes back of length 1.
#'
#' @seealso [variance.PseudoHuberDistrib()] for the other Bessel ratio,
#'   [kurtosis.LaplaceDistrib()] and [kurtosis.Gaussian1Distrib()] for the two
#'   limits, [pseudohuber_distrib()] for the family.
#'
#' @examples
#' d <- pseudohuber_distrib()
#'
#' # The shape carries the tail weight, from the Laplace's 3 to the Gaussian's 0.
#' round(kurtosis(d, list(mu = 0, sigma = 1, nu = c(1e-4, 1, 1e2, 1e4))), 4)
#'
#' # Neither the location nor the scale enters a standardized moment.
#' kurtosis(d, list(mu = c(0, 5), sigma = c(1, 9), nu = 4))
#'
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
#'
#' @description
#' Closed form, replacing the numerical default: \eqn{E[Y] = \mu}. The density
#' is symmetric about \eqn{\mu}, so the location is the mean and the median at
#' once. The value is recycled to the length the two parameters imply, and no
#' quadrature is run: the kink at \eqn{y = \mu} makes the numerical route
#' needlessly awkward for a quantity available in one read.
#'
#' @param x A `LaplaceDistrib`, from [laplace_distrib()].
#' @param theta A named list with components `mu` (the location) and `sigma`
#'   (the scale, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [variance.LaplaceDistrib()], [kurtosis.LaplaceDistrib()],
#'   [mean.Laplace2Distrib()] for the same family written by its rate,
#'   [laplace_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#'
#' # The location is the mean, and the scale does not move it.
#' mean(d, list(mu = c(-1, 0, 4), sigma = 3))
#'
#' @keywords internal
S7::method(mean, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Variance of the Laplace Distribution
#' @name variance.LaplaceDistrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \eqn{\operatorname{Var}(Y) = 2\sigma^2}. The scale is not the standard
#' deviation here: it is smaller by \eqn{\sqrt2}, which is worth knowing before
#' a Laplace scale is compared with a Gaussian one.
#'
#' @param x A `LaplaceDistrib`, from [laplace_distrib()].
#' @param theta A named list with components `mu` (the location) and `sigma`
#'   (the scale, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, as long as `sigma`. The location does
#'   not enter the value and does not lengthen it, so a setting that varies `mu`
#'   alone comes back of length 1.
#'
#' @seealso [mean.LaplaceDistrib()], [kurtosis.LaplaceDistrib()],
#'   [variance.Laplace2Distrib()] for the rate parametrization,
#'   [laplace_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#'
#' # Twice the square of the scale.
#' all.equal(variance(d, list(mu = 0, sigma = 2)), 8)
#'
#' # The standard deviation is sqrt(2) times the scale.
#' std_dev(d, list(mu = 0, sigma = 1)) / 1
#'
#' @keywords internal
S7::method(variance, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  2 * theta[[2]]^2
}

#' @title Skewness of the Laplace Distribution
#' @name skewness.LaplaceDistrib
#'
#' @description
#' Exactly zero at every parameter value. The density is symmetric about
#' \eqn{\mu}, so every odd central moment vanishes. The constant is returned
#' directly, recycled to the length the parameters imply.
#'
#' @param x A `LaplaceDistrib`, from [laplace_distrib()].
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or `n`. The values are not read, only their lengths,
#'   and the list is still aligned and validated.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [kurtosis.LaplaceDistrib()], which is 3;
#'   [mean.LaplaceDistrib()]; [laplace_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#'
#' # Zero at every location and scale, with one value per setting.
#' skewness(d, list(mu = c(0, 5), sigma = c(1, 2)))
#'
#' @keywords internal
S7::method(skewness, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(0, length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Excess Kurtosis of the Laplace Distribution
#' @name kurtosis.LaplaceDistrib
#'
#' @description
#' Exactly 3 at every parameter value, the excess over the Gaussian's own
#' fourth standardized moment. A standardized moment is free of location and
#' scale, and the Laplace has no shape parameter, so the whole family sits at
#' one point of the kurtosis axis. That point is the reference the toolkit's
#' heavier-tailed families are read against.
#'
#' @param x A `LaplaceDistrib`, from [laplace_distrib()].
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or `n`. The values are not read, only their lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of 3s, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [skewness.LaplaceDistrib()], [variance.LaplaceDistrib()],
#'   [kurtosis.PseudoHuberDistrib()], which approaches this value as its shape
#'   goes to zero; [laplace_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#'
#' # Three, whatever the location and the scale.
#' kurtosis(d, list(mu = c(0, 5), sigma = c(1, 5)))
#'
#' @keywords internal
S7::method(kurtosis, LaplaceDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(3, length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Mean of the Laplace Distribution in Location and Rate
#' @name mean.Laplace2Distrib
#'
#' @description
#' Closed form, replacing the numerical default: \eqn{E[Y] = \mu}. This is the
#' same family as [laplace_distrib()] written by its rate
#' \eqn{\lambda = 1/\sigma}, and the mean does not see the difference: it is
#' the location in either parametrization.
#'
#' @param x A `Laplace2Distrib`, from [laplace2_distrib()].
#' @param theta A named list with components `mu` (the location) and `lambda`
#'   (the rate, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$lambda))`.
#'
#' @seealso [variance.Laplace2Distrib()], where the parametrization does show;
#'   [mean.LaplaceDistrib()] for the scale form; [laplace2_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#'
#' # The location is the mean, and the rate does not move it.
#' mean(d, list(mu = c(-1, 0, 4), lambda = 2))
#'
#' @keywords internal
S7::method(mean, Laplace2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Variance of the Laplace Distribution in Location and Rate
#' @name variance.Laplace2Distrib
#'
#' @description
#' Closed form, replacing the numerical default:
#' \eqn{\operatorname{Var}(Y) = 2/\lambda^2}. With \eqn{\lambda = 1/\sigma}
#' this is the \eqn{2\sigma^2} of [variance.LaplaceDistrib()], so the two
#' parametrizations agree on the law and differ only in which number is
#' reported. A larger rate is a tighter distribution.
#'
#' @param x A `Laplace2Distrib`, from [laplace2_distrib()].
#' @param theta A named list with components `mu` (the location) and `lambda`
#'   (the rate, positive), each a numeric vector of length 1 or `n`. The
#'   variance diverges as the rate approaches zero.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, as long as `lambda`. The location does
#'   not enter the value and does not lengthen it, so a setting that varies `mu`
#'   alone comes back of length 1.
#'
#' @seealso [variance.LaplaceDistrib()], the same quantity in the scale
#'   parametrization; [mean.Laplace2Distrib()]; [laplace2_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#'
#' # Two over the square of the rate.
#' all.equal(variance(d, list(mu = 0, lambda = 0.5)), 8)
#'
#' # The two parametrizations agree when lambda = 1 / sigma.
#' all.equal(variance(d, list(mu = 0, lambda = 1 / 3)),
#'           variance(laplace_distrib(), list(mu = 0, sigma = 3)))
#'
#' @keywords internal
S7::method(variance, Laplace2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  2 / theta[[2]]^2
}

#' @title Skewness of the Laplace Distribution in Location and Rate
#' @name skewness.Laplace2Distrib
#'
#' @description
#' Exactly zero at every parameter value, the density being symmetric about
#' \eqn{\mu}. A standardized moment is free of the parametrization as well as
#' of the location and the scale, so this agrees with
#' [skewness.LaplaceDistrib()] identically.
#'
#' @param x A `Laplace2Distrib`, from [laplace2_distrib()].
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or `n`. The values are not read, only their lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length
#'   `max(length(theta$mu), length(theta$lambda))`.
#'
#' @seealso [kurtosis.Laplace2Distrib()], [skewness.LaplaceDistrib()],
#'   [laplace2_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#'
#' # Zero at every location and rate, with one value per setting.
#' skewness(d, list(mu = c(0, 5), lambda = c(1, 2)))
#'
#' @keywords internal
S7::method(skewness, Laplace2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(0, length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Excess Kurtosis of the Laplace Distribution in Location and Rate
#' @name kurtosis.Laplace2Distrib
#'
#' @description
#' Exactly 3 at every parameter value, the excess over the Gaussian. A
#' standardized moment does not see a change of parametrization, so this is the
#' same number [kurtosis.LaplaceDistrib()] returns and it moves with neither
#' the location nor the rate.
#'
#' @param x A `Laplace2Distrib`, from [laplace2_distrib()].
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or `n`. The values are not read, only their lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of 3s, of length
#'   `max(length(theta$mu), length(theta$lambda))`.
#'
#' @seealso [skewness.Laplace2Distrib()], [kurtosis.LaplaceDistrib()],
#'   [laplace2_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#'
#' # Three, in either parametrization.
#' c(kurtosis(d, list(mu = 0, lambda = 2)),
#'   kurtosis(laplace_distrib(), list(mu = 0, sigma = 0.5)))
#'
#' @keywords internal
S7::method(kurtosis, Laplace2Distrib) <- function(x, theta, ...) {
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
#' Returns \eqn{g_k = \Gamma(1 + k/\sigma)} for \eqn{k = 1, \ldots, K}. Every
#' raw moment of a Weibull is one of these times a power of the scale,
#' \eqn{E[Y^k] = \mu^k g_k}, so the four moment methods of the family share
#' this one helper and differ only in which combination of the factors they
#' assemble.
#'
#' @param sigma The shape parameter, a positive numeric vector. Each factor
#'   diverges as the shape approaches zero, the corresponding moment not
#'   existing there.
#' @param k How many factors to return. A single whole number, 4 by default,
#'   the number the excess kurtosis needs; the variance needs 2 and the skewness
#'   3, and each caller asks for only what it uses.
#'
#' @return A named list of `k` numeric vectors, `g1` to `gk`, each the length
#'   of `sigma`.
#'
#' @seealso [variance.Weibull1Distrib()], [skewness.Weibull1Distrib()] and
#'   [kurtosis.Weibull1Distrib()] for the three consumers.
#'
#' @examples
#' # At shape 1 the family is exponential and g_k is k factorial.
#' distributions7:::weibull_gamma_factors(1, 4)
#'
#' @keywords internal
weibull_gamma_factors <- function(sigma, k = 4L) {
  out <- lapply(seq_len(k), function(j) gamma(1 + j / sigma))
  names(out) <- paste0("g", seq_len(k))
  out
}

#' @title Mean of the Weibull Distribution
#' @name mean.Weibull1Distrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu\,\Gamma(1 + 1/\sigma)}. The first parameter of
#' this parametrization is the **scale**, following `gamlss`'s `WEI`, so it is
#' not the mean; the gamma factor is what separates the two, and it is 1 only
#' at shape 1, where the family is exponential.
#'
#' @details
#' A mean parametrization was available and was not taken: it would make every
#' derivative of the family a derivative of the gamma function and of its
#' inverse. The cost of the choice is paid here, in one gamma evaluation per
#' setting.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the scale and \eqn{\sigma > 0} the shape, in the
#' parametrization [weibull1_distrib()] uses.
#'
#' @param x A `Weibull1Distrib`, from [weibull1_distrib()].
#' @param theta A named list with components `mu` (the scale, positive) and
#'   `sigma` (the shape, positive), each a numeric vector of length 1 or `n`.
#'   The mean diverges as the shape approaches zero.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [variance.Weibull1Distrib()], [skewness.Weibull1Distrib()],
#'   [weibull_gamma_factors()] for the shared factors, [weibull1_distrib()].
#'
#' @examples
#' d <- weibull1_distrib()
#'
#' # The scale times Gamma(1 + 1 / shape).
#' all.equal(mean(d, list(mu = 2, sigma = 3)), 2 * gamma(1 + 1 / 3))
#'
#' # At shape 1 the family is exponential and the scale is the mean.
#' mean(d, list(mu = 2, sigma = 1))
#'
#' @keywords internal
S7::method(mean, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] * gamma(1 + 1 / theta[[2]])
}

#' @title Variance of the Weibull Distribution
#' @name variance.Weibull1Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu^2 (g_2 - g_1^2)} with
#' \eqn{g_k = \Gamma(1 + k/\sigma)}. The scale enters as a square and the shape
#' through the two gamma factors, so the coefficient of variation
#' \eqn{\sqrt{g_2 - g_1^2}/g_1} is a function of the shape alone.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the scale, \eqn{\sigma > 0} the shape and
#' \eqn{g_k = \Gamma(1 + k/\sigma)}.
#'
#' @param x A `Weibull1Distrib`, from [weibull1_distrib()].
#' @param theta A named list with components `mu` (the scale, positive) and
#'   `sigma` (the shape, positive), each a numeric vector of length 1 or `n`.
#'   The variance diverges as the shape falls below \eqn{1/2}, the second
#'   moment failing to exist there.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [mean.Weibull1Distrib()], [skewness.Weibull1Distrib()],
#'   [weibull_gamma_factors()], [weibull1_distrib()].
#'
#' @examples
#' d <- weibull1_distrib()
#'
#' # The two gamma factors, written out.
#' g1 <- gamma(1 + 1 / 3); g2 <- gamma(1 + 2 / 3)
#' all.equal(variance(d, list(mu = 2, sigma = 3)), 4 * (g2 - g1^2))
#'
#' # At shape 1 the variance is the square of the scale, as for an exponential.
#' variance(d, list(mu = 2, sigma = 1))
#'
#' @keywords internal
S7::method(variance, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  g <- weibull_gamma_factors(theta[[2]], 2L)
  theta[[1]]^2 * (g$g2 - g$g1^2)
}

#' @title Skewness of the Weibull Distribution
#' @name skewness.Weibull1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\gamma_1 = \frac{g_3 - 3 g_1 g_2 + 2 g_1^3}{(g_2 - g_1^2)^{3/2}},
#'       \qquad g_k = \Gamma(1 + k/\sigma).}
#' The scale cancels, so the value is a function of the shape alone. It falls
#' through zero at a shape near 3.60235, so a Weibull is right-skewed below
#' that and left-skewed above it, which is one of the reasons the family covers
#' more shapes than a gamma.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the shape and \eqn{g_k = \Gamma(1 + k/\sigma)}.
#'
#' @param x A `Weibull1Distrib`, from [weibull1_distrib()].
#' @param theta A named list with components `mu` (the scale, positive) and
#'   `sigma` (the shape, positive), each a numeric vector of length 1 or `n`.
#'   The third moment requires a shape above \eqn{1/3}.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, as long as `sigma`. The scale does not enter the
#'   value and does not lengthen it, so a setting that varies `mu` alone comes
#'   back of length 1.
#'
#' @seealso [kurtosis.Weibull1Distrib()], also free of the scale;
#'   [variance.Weibull1Distrib()], which is not;
#'   [weibull_gamma_factors()], [weibull1_distrib()].
#'
#' @examples
#' d <- weibull1_distrib()
#'
#' # At shape 1 the family is exponential, whose skewness is 2.
#' all.equal(skewness(d, list(mu = 1, sigma = 1)), 2)
#'
#' # The scale does not enter a standardized moment.
#' skewness(d, list(mu = c(0.1, 1, 100), sigma = 2))
#'
#' # The sign changes at a shape of about 3.60235.
#' round(skewness(d, list(mu = 1, sigma = c(2, 3.60235, 10))), 6)
#'
#' @keywords internal
S7::method(skewness, Weibull1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  g <- weibull_gamma_factors(theta[[2]], 3L)
  v <- g$g2 - g$g1^2
  (g$g3 - 3 * g$g1 * g$g2 + 2 * g$g1^3) / v^1.5
}

#' @title Excess Kurtosis of the Weibull Distribution
#' @name kurtosis.Weibull1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\gamma_2 = \frac{g_4 - 4 g_1 g_3 + 6 g_1^2 g_2 - 3 g_1^4}
#'                       {(g_2 - g_1^2)^2} - 3,
#'       \qquad g_k = \Gamma(1 + k/\sigma).}
#' The scale cancels, so the value depends on the shape alone. It is 6 at shape
#' 1, where the family is exponential, dips below zero over a band of moderate
#' shapes, and climbs again as the shape grows.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the shape and \eqn{g_k = \Gamma(1 + k/\sigma)}.
#'
#' @param x A `Weibull1Distrib`, from [weibull1_distrib()].
#' @param theta A named list with components `mu` (the scale, positive) and
#'   `sigma` (the shape, positive), each a numeric vector of length 1 or `n`.
#'   The fourth moment requires a shape above \eqn{1/4}.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, as long as `sigma`. The scale
#'   does not enter the value and does not lengthen it, so a setting that varies
#'   `mu` alone comes back of length 1.
#'
#' @seealso [skewness.Weibull1Distrib()], [variance.Weibull1Distrib()],
#'   [weibull_gamma_factors()], [weibull1_distrib()].
#'
#' @examples
#' d <- weibull1_distrib()
#'
#' # At shape 1 the family is exponential, whose excess kurtosis is 6.
#' all.equal(kurtosis(d, list(mu = 1, sigma = 1)), 6)
#'
#' # It goes negative over a band of moderate shapes.
#' round(kurtosis(d, list(mu = 1, sigma = c(1, 2, 3.6, 10))), 4)
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu + \gamma\sigma}, with
#' \eqn{\gamma \approx 0.5772157} the Euler-Mascheroni constant. The location
#' is therefore not the mean: the distribution of a maximum is shifted to the
#' right of its location by a fixed multiple of the scale. The constant is
#' obtained as `-digamma(1)`, which is exact to the last bit.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale and \eqn{\gamma} the
#' Euler-Mascheroni constant.
#'
#' @param x A `GumbelDistrib`, from [gumbel_distrib()].
#' @param theta A named list with components `mu` (the location) and `sigma`
#'   (the scale, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [variance.GumbelDistrib()]; [skewness.GumbelDistrib()] and
#'   [kurtosis.GumbelDistrib()], which are constants; [gumbel_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#'
#' # The location plus the Euler-Mascheroni constant times the scale.
#' all.equal(mean(d, list(mu = 0, sigma = 1)), -digamma(1))
#'
#' # The gap between location and mean grows with the scale.
#' mean(d, list(mu = 0, sigma = c(1, 2, 5)))
#'
#' @keywords internal
S7::method(mean, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + (-digamma(1)) * theta[[2]]
}

#' @title Variance of the Gumbel Distribution
#' @name variance.GumbelDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \pi^2\sigma^2/6}. The scale is not
#' the standard deviation: it is smaller by \eqn{\pi/\sqrt6 \approx 1.2825}.
#' The location does not enter, the family being location-scale.
#'
#' @param x A `GumbelDistrib`, from [gumbel_distrib()].
#' @param theta A named list with components `mu` (the location) and `sigma`
#'   (the scale, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, as long as `sigma`. The location does
#'   not enter the value and does not lengthen it, so a setting that varies `mu`
#'   alone comes back of length 1.
#'
#' @seealso [mean.GumbelDistrib()], [kurtosis.GumbelDistrib()],
#'   [gumbel_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#'
#' # pi^2 sigma^2 / 6.
#' all.equal(variance(d, list(mu = 0, sigma = 2)), pi^2 * 4 / 6)
#'
#' # The standard deviation is pi / sqrt(6) times the scale.
#' std_dev(d, list(mu = 0, sigma = 1))
#'
#' @keywords internal
S7::method(variance, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  pi^2 * theta[[2]]^2 / 6
}

#' @title Skewness of the Gumbel Distribution
#' @name skewness.GumbelDistrib
#'
#' @description
#' Constant: \eqn{\gamma_1 = 12\sqrt6\,\zeta(3)/\pi^3 \approx 1.1395}, with
#' \eqn{\zeta(3)} Apery's constant. The Gumbel is location-scale with no shape
#' parameter, so its standardized moments are numbers and not functions: every
#' Gumbel has this same right skew, whatever its location and scale.
#'
#' @details
#' \eqn{\zeta(3) = 1.2020569031595942854} is written out in the body. Base R
#' has no function that returns it, and one number does not justify a
#' dependency.
#'
#' @section Notation:
#' \eqn{\zeta} is the Riemann zeta function; \eqn{\mu} and \eqn{\sigma > 0} are
#' the location and the scale, neither of which enters the value.
#'
#' @param x A `GumbelDistrib`, from [gumbel_distrib()].
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or `n`. The values are not read, only their lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, every element \eqn{12\sqrt6\,\zeta(3)/\pi^3}, of
#'   length `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [kurtosis.GumbelDistrib()], the other constant;
#'   [mean.GumbelDistrib()] and [variance.GumbelDistrib()], which do move with
#'   the parameters; [gumbel_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#'
#' # One number for the whole family.
#' skewness(d, list(mu = c(0, 5, -3), sigma = c(1, 7, 0.2)))
#'
#' # Against the published constant.
#' all.equal(skewness(d, list(mu = 0, sigma = 1)),
#'           12 * sqrt(6) * 1.2020569031595942854 / pi^3)
#'
#' @keywords internal
S7::method(skewness, GumbelDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  # zeta(3), Apery's constant, written out rather than reached for: no base R
  # function returns it and one more dependency for one number is not worth it.
  zeta3 <- 1.2020569031595942854
  rep(12 * sqrt(6) * zeta3 / pi^3, length.out = max(lengths(theta[seq_len(2)])))
}

#' @title Excess Kurtosis of the Gumbel Distribution
#' @name kurtosis.GumbelDistrib
#'
#' @description
#' Constant: \eqn{\gamma_2 = 12/5 = 2.4}, the excess over the Gaussian. Like
#' the skewness it is fixed for the whole family, there being no shape parameter
#' for a standardized moment to depend on.
#'
#' @param x A `GumbelDistrib`, from [gumbel_distrib()].
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or `n`. The values are not read, only their lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of 2.4s, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [skewness.GumbelDistrib()], the other constant;
#'   [variance.GumbelDistrib()]; [gumbel_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#'
#' # 12 / 5 for every Gumbel.
#' kurtosis(d, list(mu = c(0, 5), sigma = c(1, 7)))
#'
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
#' \eqn{b\delta} with \eqn{b = \sqrt{2/\pi}}. Every moment of a skew normal is
#' a function of \eqn{b\delta}, so the four moment methods of the family share
#' this helper and assemble different combinations of one number.
#'
#' @details
#' \eqn{\delta} is bounded: it runs over \eqn{(-1, 1)} as \eqn{\alpha} runs
#' over the whole line, and it saturates quickly, reaching 0.995 at
#' \eqn{\alpha = 10}. That bound is what caps the skewness and the kurtosis the
#' family can reach, and is why the skew t exists.
#'
#' @param alpha The shape parameter, a numeric vector of any sign. Zero gives
#'   \eqn{\delta = 0} and the symmetric Gaussian case; the magnitude saturates
#'   at one as the shape grows.
#'
#' @return A named list with `delta` and `bd`, each a numeric vector the length
#'   of `alpha`.
#'
#' @seealso [mean.SkewNormal1Distrib()], [variance.SkewNormal1Distrib()],
#'   [skewness.SkewNormal1Distrib()] and [kurtosis.SkewNormal1Distrib()] for
#'   the four consumers.
#'
#' @examples
#' # delta saturates at one as the shape grows.
#' distributions7:::skewnormal_delta(c(0, 1, 10, 1e6))$delta
#'
#' @keywords internal
skewnormal_delta <- function(alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  list(delta = delta, bd = sqrt(2 / pi) * delta)
}

#' @title Mean of the Skew Normal Distribution
#' @name mean.SkewNormal1Distrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu + \sigma b \delta}, with
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and \eqn{b = \sqrt{2/\pi}}. The
#' location is not the mean unless the shape is zero: a positive shape pulls
#' the mass to the right and the mean with it, by at most
#' \eqn{\sigma\sqrt{2/\pi} \approx 0.7979\sigma}.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale and \eqn{\alpha} the
#' shape, in Azzalini's direct parametrization.
#'
#' @param x A `SkewNormal1Distrib`, from [skewnormal1_distrib()].
#' @param theta A named list with components `mu` (the location), `sigma` (the
#'   scale, positive) and `alpha` (the shape, any sign), each a numeric vector
#'   of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length equal to the longest of the
#'   three components.
#'
#' @seealso [variance.SkewNormal1Distrib()], [skewness.SkewNormal1Distrib()],
#'   [skewnormal_delta()] for the shared factor, [skewnormal1_distrib()].
#'
#' @examples
#' d <- skewnormal1_distrib()
#'
#' # The location plus sigma b delta.
#' delta <- 3 / sqrt(1 + 9)
#' all.equal(mean(d, list(mu = 0, sigma = 1, alpha = 3)),
#'           sqrt(2 / pi) * delta)
#'
#' # At shape zero the family is Gaussian and the location is the mean.
#' mean(d, list(mu = 2, sigma = 1, alpha = 0))
#'
#' @keywords internal
S7::method(mean, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + theta[[2]] * skewnormal_delta(theta[[3]])$bd
}

#' @title Variance of the Skew Normal Distribution
#' @name variance.SkewNormal1Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \sigma^2(1 - b^2\delta^2)}, with
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and \eqn{b = \sqrt{2/\pi}}. The
#' bracket is at most 1 and at least \eqn{1 - 2/\pi \approx 0.3634}, so
#' skewing the family narrows it: the scale is an upper bound on the standard
#' deviation and is attained only at shape zero.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale and \eqn{\alpha} the shape, in Azzalini's
#' direct parametrization.
#'
#' @param x A `SkewNormal1Distrib`, from [skewnormal1_distrib()].
#' @param theta A named list with components `mu` (the location), `sigma` (the
#'   scale, positive) and `alpha` (the shape, any sign), each a numeric vector
#'   of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, as long as the longer of `sigma` and
#'   `alpha`. The location does not enter the value and does not lengthen it, so
#'   a setting that varies `mu` alone comes back of length 1.
#'
#' @seealso [mean.SkewNormal1Distrib()], [kurtosis.SkewNormal1Distrib()],
#'   [skewnormal_delta()], [skewnormal1_distrib()].
#'
#' @examples
#' d <- skewnormal1_distrib()
#'
#' # The variance falls from sigma^2 towards sigma^2 (1 - 2/pi) as the shape grows.
#' round(variance(d, list(mu = 0, sigma = 1, alpha = c(0, 1, 3, 1e6))), 6)
#' 1 - 2 / pi
#'
#' @keywords internal
S7::method(variance, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 * (1 - skewnormal_delta(theta[[3]])$bd^2)
}

#' @title Skewness of the Skew Normal Distribution
#' @name skewness.SkewNormal1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\gamma_1 = \frac{4-\pi}{2}\,
#'       \frac{(b\delta)^3}{(1 - b^2\delta^2)^{3/2}},}
#' with \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and \eqn{b = \sqrt{2/\pi}}. It
#' takes the sign of the shape and vanishes at zero. Because \eqn{\delta} is
#' bounded by 1, the skewness the family can reach is bounded too, by about
#' 0.9953 in absolute value: a sample skewed more than that cannot be fitted by
#' a skew normal at any shape, and the skew t is the family to reach for.
#'
#' @section Notation:
#' \eqn{\alpha} is the shape, \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and
#' \eqn{b = \sqrt{2/\pi}}. Neither the location nor the scale enters a
#' standardized moment.
#'
#' @param x A `SkewNormal1Distrib`, from [skewnormal1_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive) and
#'   `alpha` (any sign), each a numeric vector of length 1 or `n`. Only `alpha`
#'   enters the value.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, as long as `alpha`, always inside
#'   \eqn{(-0.9953, 0.9953)}. Neither the location nor the scale enters the
#'   value or lengthens it, so a setting that varies either alone comes back of
#'   length 1.
#'
#' @seealso [kurtosis.SkewNormal1Distrib()], bounded for the same reason;
#'   [skewness.SkewTDistrib()], which is not bounded;
#'   [skewnormal_delta()], [skewnormal1_distrib()].
#'
#' @examples
#' d <- skewnormal1_distrib()
#'
#' # Zero at shape zero, and it takes the sign of the shape.
#' skewness(d, list(mu = 0, sigma = 1, alpha = c(-3, 0, 3)))
#'
#' # The reachable range stops short of one.
#' skewness(d, list(mu = 0, sigma = 1, alpha = 1e6))
#'
#' @keywords internal
S7::method(skewness, SkewNormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  bd <- skewnormal_delta(theta[[3]])$bd
  ((4 - pi) / 2) * bd^3 / (1 - bd^2)^1.5
}

#' @title Excess Kurtosis of the Skew Normal Distribution
#' @name kurtosis.SkewNormal1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\gamma_2 = 2(\pi - 3)\,
#'       \frac{(b\delta)^4}{(1 - b^2\delta^2)^{2}},}
#' with \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and \eqn{b = \sqrt{2/\pi}}. The
#' fourth power makes it even in the shape and non-negative, so a skew normal
#' is never lighter-tailed than a Gaussian, and \eqn{\delta}'s bound caps it at
#' about 0.8692. A sample needing more excess kurtosis than that is outside the
#' family whatever the shape.
#'
#' @section Notation:
#' \eqn{\alpha} is the shape, \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and
#' \eqn{b = \sqrt{2/\pi}}.
#'
#' @param x A `SkewNormal1Distrib`, from [skewnormal1_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive) and
#'   `alpha` (any sign), each a numeric vector of length 1 or `n`. Only `alpha`
#'   enters the value.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, as long as `alpha`, always in
#'   \eqn{[0, 0.8692)}. Neither the location nor the scale enters the value or
#'   lengthens it, so a setting that varies either alone comes back of length 1.
#'
#' @seealso [skewness.SkewNormal1Distrib()], bounded for the same reason;
#'   [kurtosis.SkewTDistrib()], which is not bounded;
#'   [skewnormal_delta()], [skewnormal1_distrib()].
#'
#' @examples
#' d <- skewnormal1_distrib()
#'
#' # Even in the shape, and zero only at zero.
#' kurtosis(d, list(mu = 0, sigma = 1, alpha = c(-3, 0, 3)))
#'
#' # The reachable range stops short of 0.8692.
#' kurtosis(d, list(mu = 0, sigma = 1, alpha = 1e6))
#'
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
#' Returns the four pieces the family's four moment methods share:
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}}, the constant
#' \eqn{b_\nu = \sqrt{\nu/\pi}\,\Gamma\{(\nu-1)/2\}/\Gamma(\nu/2)}, the mean
#' \eqn{\mu_z = \delta b_\nu} of the standardized variable and its variance
#' \eqn{\sigma_z^2 = \nu/(\nu-2) - \mu_z^2}.
#'
#' @details
#' The gamma ratio in \eqn{b_\nu} is formed as
#' `exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))`, so it stays finite at degrees
#' of freedom where the two gamma functions themselves would overflow.
#'
#' Existence is enforced here, once for all four callers: \eqn{b_\nu} is `NaN`
#' at \eqn{\nu \le 1} and \eqn{\sigma_z^2} is `NaN` at \eqn{\nu \le 2}, so a
#' moment built from either inherits the `NaN` and no method has to test the
#' threshold twice.
#'
#' @param alpha The shape parameter, a numeric vector of any sign.
#' @param nu The degrees of freedom, a positive numeric vector. Values at or
#'   below 1 give `NaN` in `bnu` and `mz`; values at or below 2 give `NaN` in
#'   `vz` as well.
#'
#' @return A named list with `delta`, `bnu`, `mz` and `vz`, each a numeric
#'   vector recycled to the longer of `alpha` and `nu`.
#'
#' @seealso [mean.SkewTDistrib()], [variance.SkewTDistrib()],
#'   [skewness.SkewTDistrib()] and [kurtosis.SkewTDistrib()] for the four
#'   consumers; [skewnormal_delta()] for the same shape factor in the
#'   Gaussian-tailed family.
#'
#' @examples
#' # Below two degrees of freedom the variance piece is NaN.
#' distributions7:::skewt_moment_pieces(alpha = 2, nu = c(0.5, 1.5, 5))$vz
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
#'
#' @description
#' Closed form for \eqn{\nu > 1}: \eqn{E[Y] = \mu + \sigma\delta b_\nu}, with
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and
#' \eqn{b_\nu = \sqrt{\nu/\pi}\,\Gamma\{(\nu-1)/2\}/\Gamma(\nu/2)}. At
#' \eqn{\nu \le 1} the value is `NaN`: the density is perfectly well defined
#' there and its first moment is not, so a number would be a claim the law does
#' not support.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\alpha} the
#' shape and \eqn{\nu > 0} the degrees of freedom, in Azzalini and Capitanio's
#' parametrization.
#'
#' @param x A `SkewTDistrib`, from [skewt_distrib()].
#' @param theta A named list with components `mu` (the location), `sigma` (the
#'   scale, positive), `alpha` (the shape, any sign) and `nu` (the degrees of
#'   freedom, positive), each a numeric vector of length 1 or `n`. Settings
#'   with \eqn{\nu \le 1} give `NaN` in their own positions and do not affect
#'   the others.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length equal to the longest of the
#'   four components, `NaN` wherever \eqn{\nu \le 1}.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2003). Distributions generated by
#' perturbation of symmetry with emphasis on a multivariate skew t
#' distribution. *Journal of the Royal Statistical Society Series B* **65**,
#' 367-389.
#'
#' @seealso [variance.SkewTDistrib()], whose threshold is 2;
#'   [skewt_moment_pieces()] for the shared quantities;
#'   [mean.SkewNormal1Distrib()], the Gaussian-tailed limit; [skewt_distrib()].
#'
#' @examples
#' d <- skewt_distrib()
#'
#' # Symmetric at shape zero, so the location is the mean.
#' mean(d, list(mu = 0, sigma = 1, alpha = 0, nu = 5))
#'
#' # The mean does not exist at or below one degree of freedom.
#' mean(d, list(mu = 0, sigma = 1, alpha = 2, nu = c(1, 5)))
#'
#' @keywords internal
S7::method(mean, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  p <- skewt_moment_pieces(theta[[3]], theta[[4]])
  theta[[1]] + theta[[2]] * p$mz
}

#' @title Variance of the Skew t Distribution
#' @name variance.SkewTDistrib
#'
#' @description
#' Closed form for \eqn{\nu > 2}:
#' \eqn{\operatorname{Var}(Y) = \sigma^2\{\nu/(\nu-2) - \mu_z^2\}}, with
#' \eqn{\mu_z = \delta b_\nu} the mean of the standardized variable. At
#' \eqn{\nu \le 2} the value is `NaN`, the second moment not existing there.
#' Two effects push in opposite directions: the \eqn{t} factor
#' \eqn{\nu/(\nu-2)} inflates the variance above the scale, and the skewing
#' term \eqn{\mu_z^2} pulls it back.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale, \eqn{\alpha} the shape, \eqn{\nu} the degrees
#' of freedom, \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} and
#' \eqn{b_\nu = \sqrt{\nu/\pi}\,\Gamma\{(\nu-1)/2\}/\Gamma(\nu/2)}.
#'
#' @param x A `SkewTDistrib`, from [skewt_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive), `alpha`
#'   (any sign) and `nu` (positive), each a numeric vector of length 1 or `n`.
#'   Settings with \eqn{\nu \le 2} give `NaN`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, as long as the longest of `sigma`,
#'   `alpha` and `nu`, `NaN` wherever \eqn{\nu \le 2}. The location does not
#'   enter the value and does not lengthen it.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2003). Distributions generated by
#' perturbation of symmetry with emphasis on a multivariate skew t
#' distribution. *Journal of the Royal Statistical Society Series B* **65**,
#' 367-389.
#'
#' @seealso [mean.SkewTDistrib()], whose threshold is 1;
#'   [skewness.SkewTDistrib()], whose threshold is 3;
#'   [skewt_moment_pieces()], [skewt_distrib()].
#'
#' @examples
#' d <- skewt_distrib()
#'
#' # At shape zero this is the Student t's nu / (nu - 2).
#' all.equal(variance(d, list(mu = 0, sigma = 1, alpha = 0, nu = 5)), 5 / 3)
#'
#' # It does not exist at or below two degrees of freedom.
#' variance(d, list(mu = 0, sigma = 1, alpha = 2, nu = c(2, 5)))
#'
#' @keywords internal
S7::method(variance, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 * skewt_moment_pieces(theta[[3]], theta[[4]])$vz
}

#' @title Skewness of the Skew t Distribution
#' @name skewness.SkewTDistrib
#'
#' @description
#' Closed form for \eqn{\nu > 3}, from Azzalini and Capitanio (2003):
#' \deqn{\gamma_1 = \frac{\mu_z}{\sigma_z^{3}}
#'       \left\{\frac{\nu(3-\delta^2)}{\nu-3}
#'              - \frac{3\nu}{\nu-2} + 2\mu_z^2\right\},}
#' and `NaN` at \eqn{\nu \le 3}. Unlike the skew normal's, this skewness is
#' unbounded: heavy tails supply asymmetry the Gaussian-tailed family cannot
#' reach at any shape.
#'
#' @section Notation:
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}}, \eqn{\mu_z = \delta b_\nu} and
#' \eqn{\sigma_z^2 = \nu/(\nu-2) - \mu_z^2} are the standardized mean and
#' variance from [skewt_moment_pieces()]. Neither the location nor the scale
#' enters a standardized moment.
#'
#' @param x A `SkewTDistrib`, from [skewt_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive), `alpha`
#'   (any sign) and `nu` (positive), each a numeric vector of length 1 or `n`.
#'   Only `alpha` and `nu` enter the value; settings with \eqn{\nu \le 3} give
#'   `NaN`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, as long as the longer of `alpha` and `nu`, `NaN`
#'   wherever \eqn{\nu \le 3}. Neither the location nor the scale enters the
#'   value or lengthens it.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2003). Distributions generated by
#' perturbation of symmetry with emphasis on a multivariate skew t
#' distribution. *Journal of the Royal Statistical Society Series B* **65**,
#' 367-389.
#'
#' @seealso [kurtosis.SkewTDistrib()], whose threshold is 4;
#'   [skewness.SkewNormal1Distrib()], the bounded limit;
#'   [skewt_moment_pieces()], [skewt_distrib()].
#'
#' @examples
#' d <- skewt_distrib()
#'
#' # Undefined at or below three degrees of freedom, then large and falling.
#' round(skewness(d, list(mu = 0, sigma = 1, alpha = 2, nu = c(3, 4, 10))), 4)
#'
#' # As nu grows it approaches the skew normal's value at the same shape.
#' c(skewt = skewness(d, list(mu = 0, sigma = 1, alpha = 3, nu = 1e6)),
#'   skewnormal = skewness(skewnormal1_distrib(),
#'                         list(mu = 0, sigma = 1, alpha = 3)))
#'
#' @keywords internal
S7::method(skewness, SkewTDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  nu <- theta[[4]]
  p <- skewt_moment_pieces(theta[[3]], nu)
  val <- (p$mz / p$vz^1.5) *
    (nu * (3 - p$delta^2) / (nu - 3) - 3 * nu / (nu - 2) + 2 * p$mz^2)
  ifelse(nu > 3, val, NaN)
}

#' @title Excess Kurtosis of the Skew t Distribution
#' @name kurtosis.SkewTDistrib
#'
#' @description
#' Closed form for \eqn{\nu > 4}, from Azzalini and Capitanio (2003):
#' \deqn{\gamma_2 = \frac{1}{\sigma_z^{4}}
#'       \left\{\frac{3\nu^2}{(\nu-2)(\nu-4)}
#'              - \frac{4\mu_z^2\nu(3-\delta^2)}{\nu-3}
#'              + \frac{6\mu_z^2\nu}{\nu-2} - 3\mu_z^4\right\} - 3,}
#' and `NaN` at \eqn{\nu \le 4}. The threshold is the highest of the four
#' moments, the fourth being the one that fails first as the tails grow heavy.
#'
#' @section Notation:
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}}, \eqn{\mu_z = \delta b_\nu} and
#' \eqn{\sigma_z^2 = \nu/(\nu-2) - \mu_z^2} are the standardized mean and
#' variance from [skewt_moment_pieces()].
#'
#' @param x A `SkewTDistrib`, from [skewt_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive), `alpha`
#'   (any sign) and `nu` (positive), each a numeric vector of length 1 or `n`.
#'   Only `alpha` and `nu` enter the value; settings with \eqn{\nu \le 4} give
#'   `NaN`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, as long as the longer of
#'   `alpha` and `nu`, `NaN` wherever \eqn{\nu \le 4}. Neither the location nor
#'   the scale enters the value or lengthens it.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2003). Distributions generated by
#' perturbation of symmetry with emphasis on a multivariate skew t
#' distribution. *Journal of the Royal Statistical Society Series B* **65**,
#' 367-389.
#'
#' @seealso [skewness.SkewTDistrib()], whose threshold is 3;
#'   [kurtosis.SkewNormal1Distrib()], the bounded limit;
#'   [skewt_moment_pieces()], [skewt_distrib()].
#'
#' @examples
#' d <- skewt_distrib()
#'
#' # At shape zero this is the Student t's 6 / (nu - 4).
#' kurtosis(d, list(mu = 0, sigma = 1, alpha = 0, nu = c(4, 5, 10)))
#'
#' # As nu grows it approaches the skew normal's value at the same shape.
#' c(skewt = kurtosis(d, list(mu = 0, sigma = 1, alpha = 3, nu = 1e6)),
#'   skewnormal = kurtosis(skewnormal1_distrib(),
#'                         list(mu = 0, sigma = 1, alpha = 3)))
#'
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
#' Returns `value` repeated to the length the parameters imply. A moment that
#' does not depend on the parameters still has to come back one value per
#' setting, so that every moment method returns the same shape whatever the
#' family; this helper is what enforces that for the constants.
#'
#' @param theta An aligned named list of parameters, as [align_theta()] returns.
#' @param k How many of its components to read the length from. A single whole
#'   number, usually the number of parameters the family carries.
#' @param value The constant to recycle. A numeric vector of length 1, which is
#'   `NaN` for a moment that does not exist and a number for one that is fixed
#'   by the family.
#'
#' @return A numeric vector of length `max(lengths(theta[seq_len(k)]))`, every
#'   element `value`.
#'
#' @seealso [skewness.CauchyDistrib()], which uses it for a moment that does
#'   not exist, and [kurtosis.GumbelDistrib()], which uses it for one that is
#'   constant.
#'
#' @examples
#' # One value per parameter setting, even though the value is fixed.
#' skewness(gumbel_distrib(), list(mu = c(0, 1, 2), sigma = 1))
#'
#' @keywords internal
moment_const <- function(theta, k, value) {
  rep(value, length.out = max(lengths(theta[seq_len(k)])))
}

# --- gaussian1 -------------------------------------------------------------

#' @title Mean of the Gaussian Distribution
#' @name mean.Gaussian1Distrib
#' @description Closed form: \eqn{E[Y] = \mu}.
#' @param x A `Gaussian1Distrib`.
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
#' @param x A `Gaussian1Distrib`.
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
#' @param x A `Gaussian1Distrib`.
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
#' @param x A `Gaussian1Distrib`.
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
#' `NaN`. No moment of the Cauchy exists, the tails being too heavy for
#' \eqn{\int |y| f(y) \, dy} to converge, and `NaN` is returned rather
#' than the artifact of the divergent quadrature the numerical default would
#' attempt.
#' @param x A `CauchyDistrib`.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of `NaN`.
#' @keywords internal
S7::method(mean, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Variance
#' @name variance.CauchyDistrib
#' @description `NaN`; see [`mean()`][mean.CauchyDistrib].
#' @param x A `CauchyDistrib`.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of `NaN`.
#' @keywords internal
S7::method(variance, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Skewness
#' @name skewness.CauchyDistrib
#'
#' @description
#' Returns `NaN`. Skewness is the third standardized central moment, so it needs
#' the first three moments of \eqn{Y}, and the Cauchy has none of them. Its
#' density decays like \eqn{y^{-2}}, so \eqn{\int |y|^p f(y)\,\mathrm{d}y}
#' diverges for every \eqn{p \ge 1}; already the mean fails to exist, and the
#' second and third moments fail worse. `NaN` is the value of a quantity that is
#' not defined, and it is returned directly so that no quadrature is attempted.
#'
#' @details
#' The generic [skewness()] evaluates
#' \eqn{\gamma_1 = \mathbb{E}[((Y - \mathbb{E}Y)/\mathrm{sd}(Y))^3]}
#' for a `distrib` object by calling [moment()] three times. On a Cauchy each of
#' those integrals diverges, and a numerical quadrature over a divergent
#' integral does not fail loudly: it returns whatever its truncation happens to
#' give, a number that changes with the panel layout and looks like an estimate.
#' This method short-circuits that.
#'
#' The same holds for [mean.CauchyDistrib()], [variance.CauchyDistrib()] and
#' [kurtosis.CauchyDistrib()], and for the sample versions: a Cauchy sample mean
#' does not converge as \eqn{n} grows, it is itself Cauchy with the same scale,
#' so no amount of data settles it. The example below shows that.
#'
#' What the Cauchy does have is a **median**, equal to \eqn{\mu}, and a
#' half-interquartile range, equal to \eqn{\sigma}. Both are exact and finite at
#' every parameter value, and both are available from [distrib_quantile()].
#' Reach for those when a location and a spread are wanted.
#'
#' @section Notation:
#' \eqn{\mu} is the location and \eqn{\sigma > 0} the scale of the Cauchy, in
#' the parametrization `cauchy_distrib()` uses. Neither is a moment: \eqn{\mu}
#' is the median and \eqn{\sigma} the half-interquartile range.
#'
#' @param x A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param theta A named list of parameters with components `mu` and `sigma`,
#'   each a numeric vector of length 1 or `n`. The values are not read, only
#'   their lengths, but they are still aligned and validated, so a missing or
#'   out-of-bounds component throws exactly as it would for a family whose
#'   moments exist.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'   Arguments meaningful to [moment()] have no effect here, no quadrature
#'   being run.
#'
#' @return A numeric vector of `NaN`, of length
#'   `max(length(theta$mu), length(theta$sigma))`, matching the one-value-per-
#'   observation shape every moment method returns.
#'
#' @seealso [mean.CauchyDistrib()], [variance.CauchyDistrib()] and
#'   [kurtosis.CauchyDistrib()], which return `NaN` for the same reason;
#'   [distrib_quantile()] for the median and quartiles, which do exist;
#'   [skewness()] for the generic and the families that answer it with a number.
#'
#' @examples
#' d <- cauchy_distrib()
#'
#' # No moment of the Cauchy exists, so all four are NaN.
#' skewness(d, list(mu = 0, sigma = 1))
#' c(mean(d, list(mu = 0, sigma = 1)),
#'   variance(d, list(mu = 0, sigma = 1)),
#'   kurtosis(d, list(mu = 0, sigma = 1)))
#'
#' # One value per observation, as for a family whose moments do exist.
#' skewness(d, list(mu = c(0, 1, 2), sigma = 1))
#'
#' # Why: a Cauchy sample mean is itself Cauchy, so it never settles.
#' set.seed(1)
#' y <- distrib_rng(d, 1e6, list(mu = 0, sigma = 1))
#' vapply(c(1e3, 1e4, 1e5, 1e6), function(n) mean(y[1:n]), numeric(1))
#'
#' # The median and the half-IQR are exact, and are mu and sigma.
#' distrib_quantile(d, 0.5, list(mu = 3, sigma = 2))
#' q <- distrib_quantile(d, c(0.25, 0.75), list(mu = 3, sigma = 2))
#' diff(q) / 2
#'
#' @keywords internal
S7::method(skewness, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Kurtosis
#' @name kurtosis.CauchyDistrib
#' @description `NaN`; see [`mean()`][mean.CauchyDistrib].
#' @param x A `CauchyDistrib`.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of `NaN`.
#' @keywords internal
S7::method(kurtosis, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

# --- logistic --------------------------------------------------------------

#' @title Mean of the Logistic Distribution
#' @name mean.LogisticDistrib
#' @description Closed form: \eqn{E[Y] = \mu}.
#' @param x A `LogisticDistrib`.
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
#' @param x A `LogisticDistrib`.
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
#' @param x A `LogisticDistrib`.
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
#' @param x A `LogisticDistrib`.
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
#' @description Closed form: \eqn{\mu} for \eqn{\nu > 1}, and `NaN` below,
#'   where the first moment does not exist.
#' @param x A `StudentT1Distrib`.
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
#' \eqn{1 < \nu \le 2}, and `NaN` at or below one.
#' @param x A `StudentT1Distrib`.
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
#' @description Closed form: zero for \eqn{\nu > 3} by symmetry, `NaN`
#'   below, where the third moment does not exist.
#' @param x A `StudentT1Distrib`.
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
#' \eqn{2 < \nu \le 4}, and `NaN` at or below two.
#' @param x A `StudentT1Distrib`.
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
#' @param x A `Gamma2Distrib`.
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
#' @param x A `Gamma2Distrib`.
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
#' @param x A `Gamma2Distrib`.
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
#' @param x A `Gamma2Distrib`.
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
#' @param x An `ExponentialDistrib`.
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
#' @param x An `ExponentialDistrib`.
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
#' @param x An `ExponentialDistrib`.
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
#' @param x An `ExponentialDistrib`.
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
#' @param x A `ChisqDistrib`.
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
#' @param x A `ChisqDistrib`.
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
#' @param x A `ChisqDistrib`.
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
#' @param x A `ChisqDistrib`.
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
#' @param x A `Lognormal1Distrib`.
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
#' @param x A `Lognormal1Distrib`.
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
#' @param x A `Lognormal1Distrib`.
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
#' @param x A `Lognormal1Distrib`.
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
#' @param x An `InvGauss1Distrib`.
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
#' @param x An `InvGauss1Distrib`.
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
#' @param x An `InvGauss1Distrib`.
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
#' @param x An `InvGauss1Distrib`.
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
#' @param x A `Beta1Distrib`.
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
#' @param x A `Beta1Distrib`.
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
#' @param x A `Beta1Distrib`.
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
#' @param x A `Beta1Distrib`.
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
#' @param x A `GPDDistrib`.
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
#' @param x A `GPDDistrib`.
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
#' @param x A `GPDDistrib`.
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
#' @param x A `GPDDistrib`.
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
#' @param x A `GenGamma1Distrib`.
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
#' @param x A `GenGamma1Distrib`.
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
#' @param x A `GenGamma1Distrib`.
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
#' @param x A `GenGamma1Distrib`.
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
#' @param x A `PoissonDistrib`.
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
#' @param x A `PoissonDistrib`.
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
#' @param x A `PoissonDistrib`.
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
#' @param x A `PoissonDistrib`.
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
#' @param x A `BernoulliDistrib`.
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
#' @param x A `BernoulliDistrib`.
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
#' @param x A `BernoulliDistrib`.
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
#' @param x A `BernoulliDistrib`.
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
#' @param x A `BinomialDistrib`.
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
#' @param x A `BinomialDistrib`.
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
#' @param x A `BinomialDistrib`.
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
#' @param x A `BinomialDistrib`.
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
#' @param x A `GeometricDistrib`.
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
#' @param x A `GeometricDistrib`.
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
#' @param x A `GeometricDistrib`.
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
#' @param x A `GeometricDistrib`.
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
#' @param x A `NegBin1Distrib`.
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
#' @param x A `NegBin1Distrib`.
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
#' @param x A `NegBin1Distrib`.
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
#' @param x A `NegBin1Distrib`.
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
#' @param x A `BetaBinom1Distrib`.
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
#' @param x A `BetaBinom1Distrib`.
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
#' @param x A `BetaBinom1Distrib`.
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
#' @param x A `BetaBinom1Distrib`.
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
#' @param x A `Gaussian2Distrib`.
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
#' @param x A `Gaussian2Distrib`.
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
#' @param x A `Gaussian2Distrib`.
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
#' @param x A `Gaussian2Distrib`.
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
#' @param x A `Gaussian3Distrib`.
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
#' @param x A `Gaussian3Distrib`.
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
#' @param x A `Gaussian3Distrib`.
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
#' @param x A `Gaussian3Distrib`.
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
#' @param x A `Gamma1Distrib`.
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
#' @param x A `Gamma1Distrib`.
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
#' @param x A `Gamma1Distrib`.
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
#' @param x A `Gamma1Distrib`.
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
#' @param x An `InvGauss2Distrib`.
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
#' @param x An `InvGauss2Distrib`.
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
#' @param x An `InvGauss2Distrib`.
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
#' @param x An `InvGauss2Distrib`.
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
#' @param x A `Beta2Distrib`.
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
#' @param x A `Beta2Distrib`.
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
#' @param x A `Beta2Distrib`.
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
#' @param x A `Beta2Distrib`.
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
#' @param x A `BetaBinom2Distrib`.
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
#' @param x A `BetaBinom2Distrib`.
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
#' @param x A `BetaBinom2Distrib`.
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
#' @param x A `BetaBinom2Distrib`.
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
#' @param x A `Pig1Distrib` object.
#' @param theta A list containing `mu` and `sigma`.
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
#' @description The cumulants of [`pig1()`][moments.Pig1Distrib]
#' at the dispersion [pig2_sigma()] implies.
#' @param x A `Pig2Distrib` object.
#' @param theta A list containing `mu` and `alpha`.
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
