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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. The first parameter of this parametrization
#' is the mean itself, so the method reads it off and recycles it to the length
#' the two parameters imply.
#'
#' @param x A `Gaussian1Distrib`, from [gaussian1_distrib()].
#' @param theta A named list with components `mu` (the mean, any real value)
#'   and `sigma` (the standard deviation, positive), each a numeric vector of
#'   length 1 or `n`. Aligned and validated by name, so a missing or
#'   out-of-bounds component throws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [variance.Gaussian1Distrib()], [skewness.Gaussian1Distrib()] and
#'   [kurtosis.Gaussian1Distrib()], which are the family's other three moments;
#'   [gaussian1_distrib()] for the family.
#'
#' @examples
#' d <- gaussian1_distrib()
#'
#' # The first parameter is the mean, and the scale does not move it.
#' mean(d, list(mu = c(-1, 0, 4), sigma = 3))
#'
#' @keywords internal
S7::method(mean, Gaussian1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gaussian Distribution
#' @name variance.Gaussian1Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \sigma^2}. This parametrization
#' carries the standard deviation, so the variance is its square and
#' [std_dev.distrib()] returns the parameter back unchanged.
#'
#' @param x A `Gaussian1Distrib`, from [gaussian1_distrib()].
#' @param theta A named list with components `mu` (the mean) and `sigma` (the
#'   standard deviation, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [mean.Gaussian1Distrib()]; [variance.Gaussian2Distrib()] and
#'   [variance.Gaussian3Distrib()], the same law written by its variance and by
#'   its precision; [gaussian1_distrib()].
#'
#' @examples
#' d <- gaussian1_distrib()
#'
#' # The square of the second parameter.
#' variance(d, list(mu = 0, sigma = c(1, 2, 4)))
#'
#' # The standard deviation is the parameter itself.
#' all.equal(std_dev(d, list(mu = 0, sigma = 3)), 3)
#'
#' @keywords internal
S7::method(variance, Gaussian1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 2L, 0)
}

#' @title Skewness of the Gaussian Distribution
#' @name skewness.Gaussian1Distrib
#'
#' @description
#' Exactly zero at every parameter value. The density is symmetric about
#' \eqn{\mu}, so every odd central moment vanishes. The constant is recycled to
#' the length the parameters imply, and no quadrature is run: the numerical
#' default returns a small non-zero number here, and an exact zero is available.
#'
#' @param x A `Gaussian1Distrib`, from [gaussian1_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths, and the list is still aligned and validated.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [kurtosis.Gaussian1Distrib()], zero by the excess convention;
#'   [skewness()] for the generic; [gaussian1_distrib()].
#'
#' @examples
#' d <- gaussian1_distrib()
#'
#' # Zero at every location and scale, with one value per setting.
#' skewness(d, list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#'
#' @keywords internal
S7::method(skewness, Gaussian1Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Excess Kurtosis of the Gaussian Distribution
#' @name kurtosis.Gaussian1Distrib
#'
#' @description
#' Exactly zero at every parameter value. The raw fourth standardized moment of
#' a Gaussian is 3, and [kurtosis()] reports the excess over that number, so
#' this family sits at the origin of the scale by construction. Every other
#' family's excess kurtosis is read as a comparison with it.
#'
#' @param x A `Gaussian1Distrib`, from [gaussian1_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [kurtosis()] for the convention; [kurtosis.LaplaceDistrib()], which
#'   is 3, and [kurtosis.LogisticDistrib()], which is 6/5;
#'   [gaussian1_distrib()].
#'
#' @examples
#' d <- gaussian1_distrib()
#'
#' # Zero for every Gaussian, which is what the excess convention means.
#' kurtosis(d, list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#'
#' @keywords internal
S7::method(kurtosis, Gaussian1Distrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

# --- cauchy ----------------------------------------------------------------

#' @title The Cauchy Distribution Has No Mean
#' @name mean.CauchyDistrib
#'
#' @description
#' Returns `NaN`. The Cauchy density decays like \eqn{y^{-2}}, so
#' \eqn{\int |y| f(y)\,\mathrm{d}y} diverges and the first moment does not
#' exist. `NaN` is returned directly, so that no quadrature is attempted: a
#' numerical integration over a divergent integral returns whatever its
#' truncation gives, a number that moves with the panel layout and reads like
#' an estimate.
#'
#' @details
#' What the family does have is a **median**, equal to \eqn{\mu}, and a
#' half-interquartile range, equal to \eqn{\sigma}. Both are exact at every
#' parameter value and both come from [distrib_quantile()]. A Cauchy sample
#' mean is itself Cauchy with the same scale, so it does not settle as the
#' sample grows and cannot stand in for the missing moment either.
#'
#' @section Notation:
#' \eqn{\mu} is the location and \eqn{\sigma > 0} the scale, in the
#' parametrization [cauchy_distrib()] uses. Neither is a moment.
#'
#' @param x A `CauchyDistrib`, from [cauchy_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths, and the list is still aligned and validated, so a missing or
#'   out-of-bounds component throws as it would for a family whose moments
#'   exist.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `NaN`, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [skewness.CauchyDistrib()] for the full argument;
#'   [variance.CauchyDistrib()] and [kurtosis.CauchyDistrib()], `NaN` for the
#'   same reason; [distrib_quantile()] for the median and the quartiles.
#'
#' @examples
#' d <- cauchy_distrib()
#'
#' # No moment exists, so this is NaN at every parameter value.
#' mean(d, list(mu = c(0, 1, 2), sigma = 1))
#'
#' # The median is mu and the half-interquartile range is sigma.
#' distrib_quantile(d, c(0.25, 0.5, 0.75), list(mu = 3, sigma = 2))
#'
#' @keywords internal
S7::method(mean, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

#' @title The Cauchy Distribution Has No Variance
#' @name variance.CauchyDistrib
#'
#' @description
#' Returns `NaN`. The variance is the second central moment, and it needs a
#' first moment to centre on; the Cauchy has neither. Its density decays like
#' \eqn{y^{-2}}, so \eqn{\int |y|^p f(y)\,\mathrm{d}y} diverges for every
#' \eqn{p \ge 1}, and the failure at \eqn{p = 2} is the worse of the two.
#' `NaN` is returned directly, so that no quadrature is attempted.
#'
#' @details
#' The spread a Cauchy does have is its half-interquartile range, exactly
#' \eqn{\sigma}, available from [distrib_quantile()]. Note that \eqn{\sigma} is
#' a legitimate scale parameter: it sets the width of the density in the
#' ordinary way, and only its identification with a standard deviation fails.
#'
#' @section Notation:
#' \eqn{\mu} is the location and \eqn{\sigma > 0} the scale, in the
#' parametrization [cauchy_distrib()] uses.
#'
#' @param x A `CauchyDistrib`, from [cauchy_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths, and the list is still aligned and validated.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `NaN`, of length
#'   `max(length(theta$mu), length(theta$sigma))`. [std_dev()] gives `NaN` too,
#'   being the square root of this.
#'
#' @seealso [skewness.CauchyDistrib()] for the full argument;
#'   [mean.CauchyDistrib()] and [kurtosis.CauchyDistrib()], `NaN` for the same
#'   reason; [distrib_quantile()] for the quartiles;
#'   [variance.StudentT1Distrib()], which is this family at
#'   \eqn{\nu = 1} and reports `NaN` there too.
#'
#' @examples
#' d <- cauchy_distrib()
#'
#' # NaN at every parameter value, and so is the standard deviation.
#' c(variance(d, list(mu = 0, sigma = 2)), std_dev(d, list(mu = 0, sigma = 2)))
#'
#' # The half-interquartile range is sigma, and it is exact.
#' q <- distrib_quantile(d, c(0.25, 0.75), list(mu = 3, sigma = 2))
#' diff(q) / 2
#'
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
#'
#' @description
#' Returns `NaN`. The excess kurtosis is the fourth standardized central
#' moment less three, so it needs the first four moments of \eqn{Y}, and the
#' Cauchy has none of them: its density decays like \eqn{y^{-2}}, and
#' \eqn{\int |y|^p f(y)\,\mathrm{d}y} diverges already at \eqn{p = 1}. `NaN` is
#' returned directly, so that no quadrature is attempted.
#'
#' @details
#' The reading a large positive excess kurtosis usually carries, that the tails
#' are heavy, is right here and cannot be quantified on this scale: the family
#' is heavy-tailed enough that the measure of tail weight is itself undefined.
#' A Student t with \eqn{\nu > 4} is the nearest family that reports a number,
#' and \eqn{6/(\nu-4)} grows without bound as \eqn{\nu} falls towards 4.
#'
#' @section Notation:
#' \eqn{\mu} is the location and \eqn{\sigma > 0} the scale, in the
#' parametrization [cauchy_distrib()] uses.
#'
#' @param x A `CauchyDistrib`, from [cauchy_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths, and the list is still aligned and validated.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `NaN`, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [skewness.CauchyDistrib()] for the full argument;
#'   [mean.CauchyDistrib()] and [variance.CauchyDistrib()], `NaN` for the same
#'   reason; [kurtosis.StudentT1Distrib()], which reports a number above four
#'   degrees of freedom.
#'
#' @examples
#' d <- cauchy_distrib()
#'
#' # NaN at every parameter value.
#' kurtosis(d, list(mu = 0, sigma = c(1, 2)))
#'
#' # The Student t is the family that does report a number, above nu = 4.
#' kurtosis(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(1, 5, 10)))
#'
#' @keywords internal
S7::method(kurtosis, CauchyDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, NaN)
}

# --- logistic --------------------------------------------------------------

#' @title Mean of the Logistic Distribution
#' @name mean.LogisticDistrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. The density is symmetric about \eqn{\mu}, so
#' the location is the mean, the median and the mode at once.
#'
#' @param x A `LogisticDistrib`, from [logistic_distrib()].
#' @param theta A named list with components `mu` (the location) and `sigma`
#'   (the scale, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [variance.LogisticDistrib()], where the scale does enter;
#'   [kurtosis.LogisticDistrib()]; [logistic_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#'
#' # The location is the mean, and the scale does not move it.
#' mean(d, list(mu = c(-1, 0, 4), sigma = 2))
#'
#' @keywords internal
S7::method(mean, LogisticDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Logistic Distribution
#' @name variance.LogisticDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \pi^2\sigma^2/3}. The scale is not
#' the standard deviation: it is smaller by \eqn{\pi/\sqrt3 \approx 1.8138}. A
#' logistic scale and a Gaussian one are therefore not comparable as they
#' stand, which matters whenever a logistic is used as a heavier-tailed
#' substitute for a Gaussian.
#'
#' @param x A `LogisticDistrib`, from [logistic_distrib()].
#' @param theta A named list with components `mu` (the location) and `sigma`
#'   (the scale, positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [mean.LogisticDistrib()], [kurtosis.LogisticDistrib()],
#'   [logistic_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#'
#' # pi^2 sigma^2 / 3.
#' all.equal(variance(d, list(mu = 0, sigma = 2)), pi^2 * 4 / 3)
#'
#' # The standard deviation is pi / sqrt(3) times the scale.
#' std_dev(d, list(mu = 0, sigma = 1))
#'
#' @keywords internal
S7::method(variance, LogisticDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  pi^2 * theta[[2]]^2 / 3 + moment_const(theta, 2L, 0)
}

#' @title Skewness of the Logistic Distribution
#' @name skewness.LogisticDistrib
#'
#' @description
#' Exactly zero at every parameter value. The density is symmetric about
#' \eqn{\mu}, so every odd central moment vanishes. The constant is recycled to
#' the length the parameters imply.
#'
#' @param x A `LogisticDistrib`, from [logistic_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [kurtosis.LogisticDistrib()], which is 6/5;
#'   [mean.LogisticDistrib()]; [logistic_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#'
#' # Zero at every location and scale, with one value per setting.
#' skewness(d, list(mu = c(0, 5), sigma = c(1, 2)))
#'
#' @keywords internal
S7::method(skewness, LogisticDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 2L, 0)
}

#' @title Excess Kurtosis of the Logistic Distribution
#' @name kurtosis.LogisticDistrib
#'
#' @description
#' Constant: \eqn{\gamma_2 = 6/5 = 1.2}, the excess over the Gaussian. The
#' family is location-scale with no shape parameter, so a standardized moment
#' is a number: every logistic has the same tail weight, a little heavier than
#' a Gaussian's and well short of a Laplace's 3.
#'
#' @param x A `LogisticDistrib`, from [logistic_distrib()].
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`. The values are not read, only their
#'   lengths.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of 1.2s, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [skewness.LogisticDistrib()], the other constant;
#'   [kurtosis.Gaussian1Distrib()] and [kurtosis.LaplaceDistrib()], the two
#'   families it sits between; [logistic_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#'
#' # 6 / 5 for every logistic.
#' kurtosis(d, list(mu = c(0, 5), sigma = c(1, 3)))
#'
#' # Between the Gaussian's 0 and the Laplace's 3.
#' c(gaussian = kurtosis(gaussian1_distrib(), list(mu = 0, sigma = 1)),
#'   logistic = kurtosis(d, list(mu = 0, sigma = 1)),
#'   laplace  = kurtosis(laplace_distrib(), list(mu = 0, sigma = 1)))
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu} for \eqn{\nu > 1}, and `NaN` at or below one,
#' where the first moment does not exist. At \eqn{\nu = 1} the family is the
#' Cauchy, whose divergence is the reason for the threshold; above it the
#' density is symmetric about \eqn{\mu} and the location is the mean.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale and \eqn{\nu > 0} the
#' degrees of freedom.
#'
#' @param x A `StudentT1Distrib`, from [student_t1_distrib()].
#' @param theta A named list with components `mu` (the location), `sigma` (the
#'   scale, positive) and `nu` (the degrees of freedom, positive), each a
#'   numeric vector of length 1 or `n`. Settings with \eqn{\nu \le 1} give
#'   `NaN` in their own positions and do not affect the others.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length equal to the longest of the
#'   three components, `NaN` wherever \eqn{\nu \le 1}.
#'
#' @seealso [variance.StudentT1Distrib()], whose threshold is 2;
#'   [mean.CauchyDistrib()], the \eqn{\nu = 1} case;
#'   [student_t1_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#'
#' # The location is the mean above one degree of freedom, and NaN at or below.
#' mean(d, list(mu = 3, sigma = 1, nu = c(0.5, 1, 2)))
#'
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
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \sigma^2\nu/(\nu-2)} for
#' \eqn{\nu > 2}, `Inf` for \eqn{1 < \nu \le 2}, and `NaN` at or below one. The
#' three answers are three different statements: a finite value, a second
#' moment that diverges to infinity while the first exists, and a second moment
#' that is undefined because the first is.
#'
#' @details
#' The factor \eqn{\nu/(\nu-2)} exceeds one at every finite \eqn{\nu} and falls
#' onto it as the degrees of freedom grow, so the scale is a lower bound on the
#' standard deviation and is attained only in the Gaussian limit. It is 3 at
#' \eqn{\nu = 3} and 1.25 at \eqn{\nu = 10}, so the inflation is large only
#' close to the threshold.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale and \eqn{\nu > 0} the degrees of freedom.
#'
#' @param x A `StudentT1Distrib`, from [student_t1_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length equal to the longest of the
#'   three components, `Inf` where \eqn{1 < \nu \le 2} and `NaN` where
#'   \eqn{\nu \le 1}.
#'
#' @seealso [kurtosis.StudentT1Distrib()], whose threshold is 4;
#'   [variance.CauchyDistrib()], the \eqn{\nu = 1} case;
#'   [student_t1_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#'
#' # NaN, infinite and finite, in that order, as nu crosses its two thresholds.
#' variance(d, list(mu = 0, sigma = 1, nu = c(1, 1.5, 2, 3)))
#'
#' # The scale is a lower bound, approached as the degrees of freedom grow.
#' variance(d, list(mu = 0, sigma = 2, nu = c(3, 10, 1e4)))
#'
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
#'
#' @description
#' Zero for \eqn{\nu > 3}, and `NaN` at or below three, where the third moment
#' does not exist. The density is symmetric about \eqn{\mu}, so the value is
#' zero wherever it is defined. The threshold records where the third moment
#' exists.
#'
#' @section Notation:
#' \eqn{\mu} is the location and \eqn{\nu > 0} the degrees of freedom.
#'
#' @param x A `StudentT1Distrib`, from [student_t1_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`. Only `nu` enters the
#'   value; settings with \eqn{\nu \le 3} give `NaN`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length equal to the longest of the three
#'   components, zero where \eqn{\nu > 3} and `NaN` elsewhere.
#'
#' @seealso [kurtosis.StudentT1Distrib()], whose threshold is 4;
#'   [skewness.SkewTDistrib()], the asymmetric extension;
#'   [student_t1_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#'
#' # Undefined at or below three degrees of freedom, then zero by symmetry.
#' skewness(d, list(mu = 0, sigma = 1, nu = c(2, 3, 4)))
#'
#' @keywords internal
S7::method(skewness, StudentT1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(3)]))
  nu <- rep(theta[[3]], length.out = n)
  ifelse(nu > 3, 0, NaN)
}

#' @title Excess Kurtosis of the Student t Distribution
#' @name kurtosis.StudentT1Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 6/(\nu-4)} for \eqn{\nu > 4}, `Inf` for
#' \eqn{2 < \nu \le 4}, and `NaN` at or below two. The value is the toolkit's
#' standard example of tail weight running with a shape parameter: it is 6 at
#' \eqn{\nu = 5}, 1 at \eqn{\nu = 10} and tends to the Gaussian's zero as the
#' degrees of freedom grow.
#'
#' @section Notation:
#' \eqn{\nu > 0} is the degrees of freedom. Neither the location nor the scale
#' enters a standardized moment.
#'
#' @param x A `StudentT1Distrib`, from [student_t1_distrib()].
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`. Only `nu` enters the
#'   value.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length equal to the longest
#'   of the three components, `Inf` where \eqn{2 < \nu \le 4} and `NaN` where
#'   \eqn{\nu \le 2}.
#'
#' @seealso [skewness.StudentT1Distrib()], whose threshold is 3;
#'   [kurtosis.CauchyDistrib()], the \eqn{\nu = 1} case;
#'   [kurtosis.SkewTDistrib()]; [student_t1_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#'
#' # NaN, infinite and finite, in that order, as nu crosses its two thresholds.
#' kurtosis(d, list(mu = 0, sigma = 1, nu = c(2, 3, 4, 5)))
#'
#' # 6 / (nu - 4), tending to the Gaussian's zero.
#' kurtosis(d, list(mu = 0, sigma = 1, nu = c(10, 100, 1e4)))
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. This parametrization carries the mean and the
#' variance as its two parameters, so both of the first two moments are reads,
#' and the shape and rate are recovered from them where they are needed.
#'
#' @param x A `Gamma2Distrib`, from [gamma2_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `sigma2` (the variance, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma2))`.
#'
#' @seealso [variance.Gamma2Distrib()], the other parameter;
#'   [skewness.Gamma2Distrib()]; [mean.Gamma1Distrib()] for the mean-dispersion
#'   parametrization; [gamma2_distrib()].
#'
#' @examples
#' d <- gamma2_distrib()
#'
#' # The first parameter is the mean.
#' mean(d, list(mu = c(1, 2, 3), sigma2 = 1))
#'
#' @keywords internal
S7::method(mean, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Gamma Distribution
#' @name variance.Gamma2Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \sigma^2}. The variance is the
#' second parameter of this parametrization, so the method reads it off. The
#' shape it implies is \eqn{a = \mu^2/\sigma^2}, and the higher moments are
#' functions of it.
#'
#' @param x A `Gamma2Distrib`, from [gamma2_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `sigma2` (the variance, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$sigma2))`.
#'
#' @seealso [mean.Gamma2Distrib()], the other parameter;
#'   [skewness.Gamma2Distrib()] and [kurtosis.Gamma2Distrib()], which are
#'   functions of the implied shape; [gamma2_distrib()].
#'
#' @examples
#' d <- gamma2_distrib()
#'
#' # The second parameter is the variance.
#' variance(d, list(mu = 2, sigma2 = c(0.25, 1, 4)))
#'
#' # A chi-squared with mu degrees of freedom is this family at sigma2 = 2 mu.
#' all.equal(variance(chisq_distrib(), list(mu = 5)),
#'           variance(d, list(mu = 5, sigma2 = 10)))
#'
#' @keywords internal
S7::method(variance, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[2]]
}

#' @title Skewness of the Gamma Distribution
#' @name skewness.Gamma2Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = 2/\sqrt a = 2\sigma/\mu} with shape
#' \eqn{a = \mu^2/\sigma^2}. It is positive at every parameter value, the
#' support being the positive half-line, and it is exactly twice the
#' coefficient of variation, so a gamma's asymmetry and its relative spread
#' carry the same information.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{\sigma^2 > 0} the variance and
#' \eqn{a = \mu^2/\sigma^2} the shape.
#'
#' @param x A `Gamma2Distrib`, from [gamma2_distrib()].
#' @param theta A named list with components `mu` (positive) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$sigma2))`, positive throughout.
#'
#' @seealso [kurtosis.Gamma2Distrib()], which is \eqn{3/2} times the square of
#'   this; [variance.Gamma2Distrib()]; [skewness.ExponentialDistrib()], the
#'   shape-1 case; [gamma2_distrib()].
#'
#' @examples
#' d <- gamma2_distrib()
#'
#' # Twice the coefficient of variation.
#' all.equal(skewness(d, list(mu = 2, sigma2 = 1)), 2 * sqrt(1) / 2)
#'
#' # It flattens as the implied shape grows.
#' skewness(d, list(mu = 2, sigma2 = c(4, 1, 0.25)))
#'
#' @keywords internal
S7::method(skewness, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  2 * sqrt(theta[[2]]) / theta[[1]]
}

#' @title Excess Kurtosis of the Gamma Distribution
#' @name kurtosis.Gamma2Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 6/a = 6\sigma^2/\mu^2} with shape
#' \eqn{a = \mu^2/\sigma^2}. It is positive at every parameter value and is
#' exactly \eqn{3\gamma_1^2/2}, so the family occupies one curve of the
#' skewness-kurtosis plane and cannot be tuned along the two axes separately.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{\sigma^2 > 0} the variance and
#' \eqn{a = \mu^2/\sigma^2} the shape.
#'
#' @param x A `Gamma2Distrib`, from [gamma2_distrib()].
#' @param theta A named list with components `mu` (positive) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$sigma2))`, positive throughout.
#'
#' @seealso [skewness.Gamma2Distrib()], to which this is tied;
#'   [kurtosis.ExponentialDistrib()], the shape-1 case; [gamma2_distrib()].
#'
#' @examples
#' d <- gamma2_distrib()
#'
#' # The family lies on the curve kurtosis = 1.5 skewness^2.
#' th <- list(mu = 2, sigma2 = 1)
#' all.equal(kurtosis(d, th), 1.5 * skewness(d, th)^2)
#'
#' # Six over the implied shape.
#' kurtosis(d, list(mu = 2, sigma2 = c(4, 1, 0.25)))
#'
#' @keywords internal
S7::method(kurtosis, Gamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  6 * theta[[2]] / theta[[1]]^2
}

# --- exponential -----------------------------------------------------------

#' @title Mean of the Exponential Distribution
#' @name mean.ExponentialDistrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. The family carries one parameter and it is
#' the mean, so this method reads it off. The rate is its reciprocal.
#'
#' @param x An `ExponentialDistrib`, from [exponential_distrib()].
#' @param theta A named list with one component, `mu` (the mean, positive), a
#'   numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, the length of `theta$mu`.
#'
#' @seealso [variance.ExponentialDistrib()], which is the square of this;
#'   [mean.Weibull1Distrib()] and [mean.Gamma2Distrib()], the two families that
#'   contain this one; [exponential_distrib()].
#'
#' @examples
#' d <- exponential_distrib()
#'
#' # The one parameter is the mean.
#' mean(d, list(mu = c(1, 2, 3)))
#'
#' @keywords internal
S7::method(mean, ExponentialDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Exponential Distribution
#' @name variance.ExponentialDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu^2}. The standard deviation
#' equals the mean, so the coefficient of variation is exactly one at every
#' parameter value. That single number is what identifies an exponential among
#' the gammas: a sample whose relative spread is far from one is not one.
#'
#' @param x An `ExponentialDistrib`, from [exponential_distrib()].
#' @param theta A named list with one component, `mu` (the mean, positive), a
#'   numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, the length of `theta$mu`.
#'
#' @seealso [mean.ExponentialDistrib()]; [variance.Gamma2Distrib()], where the
#'   two moments are free of each other; [exponential_distrib()].
#'
#' @examples
#' d <- exponential_distrib()
#'
#' # The square of the mean, so the coefficient of variation is one.
#' variance(d, list(mu = 3))
#' std_dev(d, list(mu = 3)) / mean(d, list(mu = 3))
#'
#' # The same law as a Weibull of shape 1.
#' all.equal(variance(d, list(mu = 3)),
#'           variance(weibull1_distrib(), list(mu = 3, sigma = 1)))
#'
#' @keywords internal
S7::method(variance, ExponentialDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]^2
}

#' @title Skewness of the Exponential Distribution
#' @name skewness.ExponentialDistrib
#'
#' @description
#' Constant: \eqn{\gamma_1 = 2}. The family is a scale family with no shape
#' parameter, so a standardized moment is a number: every exponential has the
#' same asymmetry, whatever its mean.
#'
#' @param x An `ExponentialDistrib`, from [exponential_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`. The value is not read, only its length.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of 2s, the length of `theta$mu`.
#'
#' @seealso [kurtosis.ExponentialDistrib()], the other constant;
#'   [skewness.Gamma2Distrib()], which is this at shape 1;
#'   [exponential_distrib()].
#'
#' @examples
#' d <- exponential_distrib()
#'
#' # Two, whatever the mean.
#' skewness(d, list(mu = c(0.1, 3, 100)))
#'
#' @keywords internal
S7::method(skewness, ExponentialDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 1L, 2)
}

#' @title Excess Kurtosis of the Exponential Distribution
#' @name kurtosis.ExponentialDistrib
#'
#' @description
#' Constant: \eqn{\gamma_2 = 6}, the excess over the Gaussian. Like the
#' skewness it is fixed for the whole family, there being no shape parameter
#' for a standardized moment to depend on.
#'
#' @param x An `ExponentialDistrib`, from [exponential_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`. The value is not read, only its length.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of 6s, the length of `theta$mu`.
#'
#' @seealso [skewness.ExponentialDistrib()], the other constant;
#'   [kurtosis.Gamma2Distrib()], which is this at shape 1;
#'   [kurtosis.Weibull1Distrib()], which is this at shape 1 too;
#'   [exponential_distrib()].
#'
#' @examples
#' d <- exponential_distrib()
#'
#' # Six, whatever the mean, and the two containing families agree.
#' c(exponential = kurtosis(d, list(mu = 3)),
#'   gamma       = kurtosis(gamma2_distrib(), list(mu = 3, sigma2 = 9)),
#'   weibull     = kurtosis(weibull1_distrib(), list(mu = 3, sigma = 1)))
#'
#' @keywords internal
S7::method(kurtosis, ExponentialDistrib) <- function(x, theta, ...) {
  moment_const(align_theta(x, theta), 1L, 6)
}

# --- chisq -----------------------------------------------------------------
#
# Parametrized by its mean, which for a chi-squared is its degrees of freedom.

#' @title Mean of the Chi-Squared Distribution
#' @name mean.ChisqDistrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. The family carries one parameter, which is
#' both the mean and the degrees of freedom, the two coinciding for this law.
#' It is not restricted to whole numbers here, so the family is the gamma with
#' variance tied to twice the mean.
#'
#' @param x A `ChisqDistrib`, from [chisq_distrib()].
#' @param theta A named list with one component, `mu` (the mean and the degrees
#'   of freedom, positive), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, the length of `theta$mu`.
#'
#' @seealso [variance.ChisqDistrib()], which is twice this;
#'   [mean.Gamma2Distrib()] for the containing family; [chisq_distrib()].
#'
#' @examples
#' d <- chisq_distrib()
#'
#' # The one parameter is the mean and the degrees of freedom at once.
#' mean(d, list(mu = c(1, 5, 20)))
#'
#' @keywords internal
S7::method(mean, ChisqDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Chi-Squared Distribution
#' @name variance.ChisqDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = 2\mu}. The variance is tied to the
#' mean, so the family has one parameter where a gamma has two: it is exactly
#' the gamma at \eqn{\sigma^2 = 2\mu}, and a sample whose variance is far from
#' twice its mean is not chi-squared at any degrees of freedom.
#'
#' @param x A `ChisqDistrib`, from [chisq_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, the length of `theta$mu`.
#'
#' @seealso [mean.ChisqDistrib()]; [variance.Gamma2Distrib()], where the two
#'   moments are free of each other; [chisq_distrib()].
#'
#' @examples
#' d <- chisq_distrib()
#'
#' # Twice the mean.
#' all.equal(variance(d, list(mu = 5)), 10)
#'
#' # The gamma at sigma2 = 2 mu is the same law.
#' all.equal(variance(d, list(mu = 5)),
#'           variance(gamma2_distrib(), list(mu = 5, sigma2 = 10)))
#'
#' @keywords internal
S7::method(variance, ChisqDistrib) <- function(x, theta, ...) {
  2 * align_theta(x, theta)[[1]]
}

#' @title Skewness of the Chi-Squared Distribution
#' @name skewness.ChisqDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = \sqrt{8/\mu}}. It is positive at every degrees
#' of freedom, the support being the positive half-line, and it decays like
#' \eqn{\mu^{-1/2}}, which is the rate at which the family approaches a
#' Gaussian as its degrees of freedom grow.
#'
#' @param x A `ChisqDistrib`, from [chisq_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, the length of `theta$mu`, positive throughout.
#'
#' @seealso [kurtosis.ChisqDistrib()], which decays twice as fast;
#'   [skewness.Gamma2Distrib()] for the containing family; [chisq_distrib()].
#'
#' @examples
#' d <- chisq_distrib()
#'
#' # Square root of eight over the degrees of freedom.
#' all.equal(skewness(d, list(mu = 5)), sqrt(8 / 5))
#'
#' # It vanishes as the degrees of freedom grow.
#' skewness(d, list(mu = c(1, 10, 1000)))
#'
#' @keywords internal
S7::method(skewness, ChisqDistrib) <- function(x, theta, ...) {
  sqrt(8 / align_theta(x, theta)[[1]])
}

#' @title Excess Kurtosis of the Chi-Squared Distribution
#' @name kurtosis.ChisqDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 12/\mu}, the excess over the Gaussian. It is
#' positive at every degrees of freedom and decays like \eqn{\mu^{-1}}, twice
#' as fast as the skewness, so a chi-squared looks Gaussian in its tails before
#' it looks symmetric.
#'
#' @param x A `ChisqDistrib`, from [chisq_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, the length of `theta$mu`,
#'   positive throughout.
#'
#' @seealso [skewness.ChisqDistrib()]; [kurtosis.Gamma2Distrib()] for the
#'   containing family; [chisq_distrib()].
#'
#' @examples
#' d <- chisq_distrib()
#'
#' # Twelve over the degrees of freedom.
#' all.equal(kurtosis(d, list(mu = 5)), 12 / 5)
#'
#' # It falls twice as fast as the skewness.
#' rbind(skewness = skewness(d, list(mu = c(1, 10, 1000))),
#'       kurtosis = kurtosis(d, list(mu = c(1, 10, 1000))))
#'
#' @keywords internal
S7::method(kurtosis, ChisqDistrib) <- function(x, theta, ...) {
  12 / align_theta(x, theta)[[1]]
}

# --- lognormal1 ------------------------------------------------------------
#
# The parameters live on the log scale, so none of them is a moment of Y.

#' @title Mean of the Lognormal Distribution
#' @name mean.Lognormal1Distrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \exp(\mu + \sigma^2/2)}. The parameters describe
#' \eqn{\log Y} and not \eqn{Y}, so \eqn{\mu} is the mean of the logarithm and
#' \eqn{\exp(\mu)} is the median of \eqn{Y}; the mean is larger than the
#' median by the factor \eqn{\exp(\sigma^2/2)}, which grows with the spread.
#'
#' @section Notation:
#' \eqn{\mu} is the mean of \eqn{\log Y} and \eqn{\sigma^2 > 0} its variance.
#'
#' @param x A `Lognormal1Distrib`, from [lognormal1_distrib()].
#' @param theta A named list with components `mu` (the mean of the logarithm,
#'   any real value) and `sigma2` (its variance, positive), each a numeric
#'   vector of length 1 or `n`. The mean overflows to `Inf` for large `sigma2`,
#'   the exponential of half of it being taken.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma2))`.
#'
#' @seealso [variance.Lognormal1Distrib()], [skewness.Lognormal1Distrib()];
#'   [distrib_quantile()] for the median, which is \eqn{\exp(\mu)};
#'   [lognormal1_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#'
#' # exp(mu + sigma2 / 2), which exceeds the median exp(mu).
#' all.equal(mean(d, list(mu = 0, sigma2 = 1)), exp(0.5))
#' distrib_quantile(d, 0.5, list(mu = 0, sigma2 = 1))
#'
#' @keywords internal
S7::method(mean, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  exp(theta[[1]] + theta[[2]] / 2)
}

#' @title Variance of the Lognormal Distribution
#' @name variance.Lognormal1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\operatorname{Var}(Y) = (e^{\sigma^2} - 1)\, e^{2\mu + \sigma^2}.}
#' The two factors separate the two effects: \eqn{e^{2\mu+\sigma^2}} is the
#' square of the mean and sets the scale, and \eqn{e^{\sigma^2} - 1} is the
#' squared coefficient of variation and depends on the spread of the logarithm
#' alone.
#'
#' @section Notation:
#' \eqn{\mu} is the mean of \eqn{\log Y} and \eqn{\sigma^2 > 0} its variance.
#'
#' @param x A `Lognormal1Distrib`, from [lognormal1_distrib()].
#' @param theta A named list with components `mu` (any real value) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`. The variance grows
#'   like \eqn{e^{2\sigma^2}} and overflows to `Inf` well before `sigma2`
#'   reaches 400.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$sigma2))`.
#'
#' @seealso [mean.Lognormal1Distrib()], [skewness.Lognormal1Distrib()],
#'   [lognormal1_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#'
#' # (exp(sigma2) - 1) exp(2 mu + sigma2).
#' all.equal(variance(d, list(mu = 0, sigma2 = 1)), (exp(1) - 1) * exp(1))
#'
#' # The coefficient of variation depends on sigma2 alone.
#' th <- list(mu = c(-2, 0, 5), sigma2 = 0.25)
#' std_dev(d, th) / mean(d, th)
#'
#' @keywords internal
S7::method(variance, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  expm1(theta[[2]]) * exp(2 * theta[[1]] + theta[[2]])
}

#' @title Skewness of the Lognormal Distribution
#' @name skewness.Lognormal1Distrib
#'
#' @description
#' Closed form:
#' \eqn{\gamma_1 = (e^{\sigma^2} + 2)\sqrt{e^{\sigma^2} - 1}}, free of
#' \eqn{\mu}. It is positive at every parameter value and grows without bound
#' with the spread of the logarithm: at \eqn{\sigma^2 = 1} it is already 6.18,
#' where a gamma of the same mean would need a shape near 0.1 to match it.
#'
#' @section Notation:
#' \eqn{\sigma^2 > 0} is the variance of \eqn{\log Y}. The location of the
#' logarithm does not enter a standardized moment.
#'
#' @param x A `Lognormal1Distrib`, from [lognormal1_distrib()].
#' @param theta A named list with components `mu` (any real value) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`. Only `sigma2` enters
#'   the value, and the exponentials overflow to `Inf` for large `sigma2`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$sigma2))`, positive throughout.
#'
#' @seealso [kurtosis.Lognormal1Distrib()], which grows faster still;
#'   [variance.Lognormal1Distrib()]; [lognormal1_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#'
#' # The published form, written out.
#' all.equal(skewness(d, list(mu = 0, sigma2 = 1)),
#'           (exp(1) + 2) * sqrt(exp(1) - 1))
#'
#' # The location of the logarithm does not enter it.
#' skewness(d, list(mu = c(-2, 0, 5), sigma2 = 0.25))
#'
#' @keywords internal
S7::method(skewness, Lognormal1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  s2 <- theta[[2]]
  (exp(s2) + 2) * sqrt(expm1(s2)) + moment_const(theta, 2L, 0)
}

#' @title Excess Kurtosis of the Lognormal Distribution
#' @name kurtosis.Lognormal1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\gamma_2 = e^{4\sigma^2} + 2e^{3\sigma^2} + 3e^{2\sigma^2} - 6,}
#' free of \eqn{\mu}. The leading term grows like \eqn{e^{4\sigma^2}}, so the
#' tail weight explodes with the spread of the logarithm: at
#' \eqn{\sigma^2 = 1} it is 110.9, and the numerical route would need a
#' quadrature far into the tail to see it at all.
#'
#' @section Notation:
#' \eqn{\sigma^2 > 0} is the variance of \eqn{\log Y}.
#'
#' @param x A `Lognormal1Distrib`, from [lognormal1_distrib()].
#' @param theta A named list with components `mu` (any real value) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`. Only `sigma2` enters
#'   the value; the leading exponential overflows to `Inf` past
#'   \eqn{\sigma^2 \approx 177}.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$sigma2))`, positive throughout.
#'
#' @seealso [skewness.Lognormal1Distrib()], [variance.Lognormal1Distrib()],
#'   [lognormal1_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#'
#' # The published form, written out.
#' all.equal(kurtosis(d, list(mu = 0, sigma2 = 1)),
#'           exp(4) + 2 * exp(3) + 3 * exp(2) - 6)
#'
#' # It climbs steeply with the spread of the logarithm.
#' round(kurtosis(d, list(mu = 0, sigma2 = c(0.1, 0.5, 1))), 3)
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. This parametrization carries the mean and a
#' dispersion, so the mean is a read and the higher moments are functions of
#' the product \eqn{\phi\mu}.
#'
#' @param x An `InvGauss1Distrib`, from [invgauss1_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `phi` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$phi))`.
#'
#' @seealso [variance.InvGauss1Distrib()], which is \eqn{\phi\mu^3};
#'   [skewness.InvGauss1Distrib()]; [invgauss1_distrib()].
#'
#' @examples
#' d <- invgauss1_distrib()
#'
#' # The first parameter is the mean.
#' mean(d, list(mu = c(1, 2, 3), phi = 0.5))
#'
#' @keywords internal
S7::method(mean, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Inverse Gaussian Distribution
#' @name variance.InvGauss1Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \phi\mu^3}. The cube is what
#' distinguishes this family from the gamma, whose variance function is
#' \eqn{\mu^2}: the inverse Gaussian's spread grows faster with the mean, so it
#' is the heavier-tailed of the two standard positive-response families.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean and \eqn{\phi > 0} the dispersion.
#'
#' @param x An `InvGauss1Distrib`, from [invgauss1_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `phi` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$phi))`.
#'
#' @seealso [mean.InvGauss1Distrib()], [skewness.InvGauss1Distrib()];
#'   [variance.Gamma2Distrib()], whose variance function is the square;
#'   [variance.InvGauss2Distrib()] for the other parametrization;
#'   [invgauss1_distrib()].
#'
#' @examples
#' d <- invgauss1_distrib()
#'
#' # The dispersion times the cube of the mean.
#' all.equal(variance(d, list(mu = 2, phi = 0.5)), 0.5 * 2^3)
#'
#' # The variance function is cubic, so it climbs faster than a gamma's.
#' variance(d, list(mu = c(1, 2, 4), phi = 0.5))
#'
#' @keywords internal
S7::method(variance, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]] * theta[[1]]^3
}

#' @title Skewness of the Inverse Gaussian Distribution
#' @name skewness.InvGauss1Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = 3\sqrt{\phi\mu}}. It is positive at every
#' parameter value and depends on the two parameters only through their
#' product, so the whole family is one curve indexed by \eqn{\phi\mu}.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean and \eqn{\phi > 0} the dispersion.
#'
#' @param x An `InvGauss1Distrib`, from [invgauss1_distrib()].
#' @param theta A named list with components `mu` (positive) and `phi`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$phi))`, positive throughout.
#'
#' @seealso [kurtosis.InvGauss1Distrib()], which is \eqn{5/3} times the square
#'   of this; [variance.InvGauss1Distrib()]; [invgauss1_distrib()].
#'
#' @examples
#' d <- invgauss1_distrib()
#'
#' # Three times the square root of the product.
#' all.equal(skewness(d, list(mu = 2, phi = 0.5)), 3 * sqrt(1))
#'
#' # Only the product matters.
#' skewness(d, list(mu = c(1, 2, 4), phi = c(4, 2, 1)))
#'
#' @keywords internal
S7::method(skewness, InvGauss1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  3 * sqrt(theta[[2]] * theta[[1]])
}

#' @title Excess Kurtosis of the Inverse Gaussian Distribution
#' @name kurtosis.InvGauss1Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 15\phi\mu}, the excess over the Gaussian. Like
#' the skewness it depends on the two parameters through their product alone,
#' and it is exactly \eqn{5\gamma_1^2/3}, so the family occupies one curve of
#' the skewness-kurtosis plane, as the gamma does with its own constant.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean and \eqn{\phi > 0} the dispersion.
#'
#' @param x An `InvGauss1Distrib`, from [invgauss1_distrib()].
#' @param theta A named list with components `mu` (positive) and `phi`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$phi))`, positive throughout.
#'
#' @seealso [skewness.InvGauss1Distrib()], to which this is tied;
#'   [kurtosis.Gamma2Distrib()], whose constant is \eqn{3/2};
#'   [invgauss1_distrib()].
#'
#' @examples
#' d <- invgauss1_distrib()
#'
#' # Fifteen times the product.
#' all.equal(kurtosis(d, list(mu = 2, phi = 0.5)), 15)
#'
#' # The family lies on the curve kurtosis = (5/3) skewness^2.
#' th <- list(mu = 2, phi = 0.5)
#' all.equal(kurtosis(d, th), (5 / 3) * skewness(d, th)^2)
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. This parametrization carries the mean and a
#' precision, so the mean is a read; the two shapes of the standard
#' parametrization are recovered as \eqn{a = \mu\phi} and
#' \eqn{b = (1-\mu)\phi}, and they sum to \eqn{\phi}.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the mean, \eqn{\phi > 0} the precision, and
#' \eqn{a = \mu\phi}, \eqn{b = (1-\mu)\phi} the two shapes.
#'
#' @param x A `Beta1Distrib`, from [beta1_distrib()].
#' @param theta A named list with components `mu` (the mean, strictly between
#'   0 and 1) and `phi` (the precision, positive), each a numeric vector of
#'   length 1 or `n`. Aligned and validated by name, so a mean at or outside
#'   the unit interval throws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$phi))`, strictly inside \eqn{(0,1)}.
#'
#' @seealso [variance.Beta1Distrib()], where the precision does enter;
#'   [skewness.Beta1Distrib()]; [mean.Beta2Distrib()] for the two-shape
#'   parametrization; [beta1_distrib()].
#'
#' @examples
#' d <- beta1_distrib()
#'
#' # The first parameter is the mean, and the precision does not move it.
#' mean(d, list(mu = c(0.2, 0.5, 0.8), phi = 5))
#'
#' @keywords internal
S7::method(mean, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Beta Distribution
#' @name variance.Beta1Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu(1-\mu)/(\phi+1)}. The
#' numerator is the Bernoulli variance at the same mean, the largest a
#' distribution on \eqn{(0,1)} can have, and \eqn{\phi+1} divides it down.
#' Larger \eqn{\phi} therefore means a tighter distribution, which is the
#' reason the parameter is called a precision.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the mean and \eqn{\phi > 0} the precision.
#'
#' @param x A `Beta1Distrib`, from [beta1_distrib()].
#' @param theta A named list with components `mu` (strictly between 0 and 1)
#'   and `phi` (positive), each a numeric vector of length 1 or `n`. The
#'   variance tends to \eqn{\mu(1-\mu)} as the precision goes to zero, where
#'   the mass piles up at the two endpoints.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$phi))`.
#'
#' @seealso [mean.Beta1Distrib()], [skewness.Beta1Distrib()],
#'   [kurtosis.Beta1Distrib()], [beta1_distrib()].
#'
#' @examples
#' d <- beta1_distrib()
#'
#' # mu (1 - mu) / (phi + 1).
#' all.equal(variance(d, list(mu = 0.3, phi = 5)), 0.3 * 0.7 / 6)
#'
#' # At mu = 1/2 and phi = 2 the density is uniform, of variance 1/12.
#' all.equal(variance(d, list(mu = 0.5, phi = 2)), 1 / 12)
#'
#' @keywords internal
S7::method(variance, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] * (1 - theta[[1]]) / (theta[[2]] + 1)
}

#' @title Skewness of the Beta Distribution
#' @name skewness.Beta1Distrib
#'
#' @description
#' Closed form in the two shapes \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}:
#' \deqn{\gamma_1 = \frac{2(b - a)\sqrt{a + b + 1}}
#'                       {(a + b + 2)\sqrt{ab}}.}
#' It takes the sign of \eqn{b - a}, so a beta is right-skewed below a mean of
#' one half, left-skewed above it, and exactly symmetric at \eqn{\mu = 1/2}
#' whatever the precision.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the mean, \eqn{\phi > 0} the precision, and
#' \eqn{a = \mu\phi}, \eqn{b = (1-\mu)\phi} the two shapes.
#'
#' @param x A `Beta1Distrib`, from [beta1_distrib()].
#' @param theta A named list with components `mu` (strictly between 0 and 1)
#'   and `phi` (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$phi))`.
#'
#' @seealso [kurtosis.Beta1Distrib()], written in the same two shapes;
#'   [variance.Beta1Distrib()]; [beta1_distrib()].
#'
#' @examples
#' d <- beta1_distrib()
#'
#' # The published form in the two shapes, written out.
#' a <- 0.3 * 5; b <- 0.7 * 5
#' all.equal(skewness(d, list(mu = 0.3, phi = 5)),
#'           2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b)))
#'
#' # Exactly zero at a mean of one half, at any precision.
#' skewness(d, list(mu = 0.5, phi = c(0.5, 2, 50)))
#'
#' @keywords internal
S7::method(skewness, Beta1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b))
}

#' @title Excess Kurtosis of the Beta Distribution
#' @name kurtosis.Beta1Distrib
#'
#' @description
#' Closed form in the two shapes \eqn{a = \mu\phi} and \eqn{b = (1-\mu)\phi}:
#' \deqn{\gamma_2 = \frac{6\{(a-b)^2(a+b+1) - ab(a+b+2)\}}
#'                       {ab(a+b+2)(a+b+3)}.}
#' It is the one family in the toolkit whose excess kurtosis is routinely
#' negative: a beta on a bounded support has no tails to be heavy, and at
#' \eqn{\mu = 1/2}, \eqn{\phi = 2} the density is uniform, whose excess
#' kurtosis is \eqn{-6/5}, the smallest value any distribution attains.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the mean, \eqn{\phi > 0} the precision, and
#' \eqn{a = \mu\phi}, \eqn{b = (1-\mu)\phi} the two shapes.
#'
#' @param x A `Beta1Distrib`, from [beta1_distrib()].
#' @param theta A named list with components `mu` (strictly between 0 and 1)
#'   and `phi` (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$phi))`.
#'
#' @seealso [skewness.Beta1Distrib()], written in the same two shapes;
#'   [variance.Beta1Distrib()]; [beta1_distrib()].
#'
#' @examples
#' d <- beta1_distrib()
#'
#' # The uniform case sits at the lower bound of -6/5.
#' all.equal(kurtosis(d, list(mu = 0.5, phi = 2)), -1.2)
#'
#' # Negative over the symmetric middle, positive as the mass piles at an end.
#' round(kurtosis(d, list(mu = c(0.5, 0.1, 0.02), phi = 5)), 4)
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \sigma/(1-\xi)} for \eqn{\xi < 1}, and `Inf` at or
#' above one, where the first moment diverges. The shape controls which moments
#' exist at all: the moment of order \eqn{k} is finite exactly for
#' \eqn{\xi < 1/k}, so a fitted object can report a mean and no variance.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale and \eqn{\xi} the shape, of either sign. At
#' \eqn{\xi = 0} the family is the exponential of mean \eqn{\sigma}, and at
#' \eqn{\xi < 0} the support is bounded above by \eqn{\sigma/|\xi|}.
#'
#' @param x A `GPDDistrib`, from [gpd_distrib()].
#' @param theta A named list with components `sigma` (the scale, positive) and
#'   `xi` (the shape, any real value), each a numeric vector of length 1 or
#'   `n`. Settings with \eqn{\xi \ge 1} give `Inf` in their own positions.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$sigma), length(theta$xi))`, `Inf` wherever
#'   \eqn{\xi \ge 1}.
#'
#' @seealso [variance.GPDDistrib()], whose threshold is \eqn{1/2};
#'   [mean.ExponentialDistrib()], the \eqn{\xi = 0} case; [gpd_distrib()].
#'
#' @examples
#' d <- gpd_distrib()
#'
#' # Finite below shape one and infinite at or above it.
#' mean(d, list(sigma = 1, xi = c(0, 0.5, 1, 1.5)))
#'
#' # At shape zero the family is exponential of mean sigma.
#' all.equal(mean(d, list(sigma = 2, xi = 0)),
#'           mean(exponential_distrib(), list(mu = 2)))
#'
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
#'
#' @description
#' Closed form:
#' \eqn{\operatorname{Var}(Y) = \sigma^2/\{(1-\xi)^2(1-2\xi)\}} for
#' \eqn{\xi < 1/2}, and `Inf` at or above one half. The threshold is half the
#' mean's, the second moment being the first to fail as the tail grows heavy,
#' and a fitted object at \eqn{\xi = 0.6} reports a finite mean beside an
#' infinite variance.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale and \eqn{\xi} the shape, of either sign.
#'
#' @param x A `GPDDistrib`, from [gpd_distrib()].
#' @param theta A named list with components `sigma` (positive) and `xi` (any
#'   real value), each a numeric vector of length 1 or `n`. Settings with
#'   \eqn{\xi \ge 1/2} give `Inf`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$sigma), length(theta$xi))`, `Inf` wherever
#'   \eqn{\xi \ge 1/2}.
#'
#' @seealso [mean.GPDDistrib()], whose threshold is 1;
#'   [skewness.GPDDistrib()], whose threshold is \eqn{1/3};
#'   [gpd_distrib()].
#'
#' @examples
#' d <- gpd_distrib()
#'
#' # Finite below shape one half and infinite at or above it.
#' variance(d, list(sigma = 1, xi = c(0, 0.4, 0.5, 0.6)))
#'
#' # A mean without a variance, which is what the two thresholds allow.
#' c(mean = mean(d, list(sigma = 1, xi = 0.6)),
#'   variance = variance(d, list(sigma = 1, xi = 0.6)))
#'
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
#'
#' @description
#' Closed form: \eqn{\gamma_1 = 2(1+\xi)\sqrt{1-2\xi}/(1-3\xi)} for
#' \eqn{\xi < 1/3}, and `Inf` at or above one third. The scale cancels, so the
#' value depends on the shape alone; it is 2 at \eqn{\xi = 0}, where the family
#' is exponential, and it climbs without bound as the shape approaches its
#' threshold.
#'
#' @section Notation:
#' \eqn{\xi} is the shape, of either sign. The scale does not enter a
#' standardized moment.
#'
#' @param x A `GPDDistrib`, from [gpd_distrib()].
#' @param theta A named list with components `sigma` (positive) and `xi` (any
#'   real value), each a numeric vector of length 1 or `n`. Only `xi` enters
#'   the value; settings with \eqn{\xi \ge 1/3} give `Inf`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$sigma), length(theta$xi))`, `Inf` wherever
#'   \eqn{\xi \ge 1/3}.
#'
#' @seealso [kurtosis.GPDDistrib()], whose threshold is \eqn{1/4};
#'   [skewness.ExponentialDistrib()], the \eqn{\xi = 0} case; [gpd_distrib()].
#'
#' @examples
#' d <- gpd_distrib()
#'
#' # Two at shape zero, climbing steeply, infinite at or above one third.
#' skewness(d, list(sigma = 1, xi = c(0, 0.3, 1 / 3, 0.4)))
#'
#' @keywords internal
S7::method(skewness, GPDDistrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  n <- max(lengths(theta[seq_len(2)]))
  xi <- rep(theta[[2]], length.out = n)
  ifelse(xi < 1 / 3, 2 * (1 + xi) * sqrt(1 - 2 * xi) / (1 - 3 * xi), Inf)
}

#' @title Excess Kurtosis of the Generalized Pareto Distribution
#' @name kurtosis.GPDDistrib
#'
#' @description
#' Closed form:
#' \eqn{\gamma_2 = 3(1-2\xi)(2\xi^2+\xi+3)/\{(1-3\xi)(1-4\xi)\} - 3} for
#' \eqn{\xi < 1/4}, and `Inf` at or above one quarter. The threshold is the
#' tightest of the four moments, so this is the first quantity to become
#' unreportable as a fitted shape rises.
#'
#' @section Notation:
#' \eqn{\xi} is the shape, of either sign. The scale does not enter a
#' standardized moment.
#'
#' @param x A `GPDDistrib`, from [gpd_distrib()].
#' @param theta A named list with components `sigma` (positive) and `xi` (any
#'   real value), each a numeric vector of length 1 or `n`. Only `xi` enters
#'   the value; settings with \eqn{\xi \ge 1/4} give `Inf`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$sigma), length(theta$xi))`, `Inf` wherever
#'   \eqn{\xi \ge 1/4}.
#'
#' @seealso [skewness.GPDDistrib()], whose threshold is \eqn{1/3};
#'   [kurtosis.ExponentialDistrib()], the \eqn{\xi = 0} case; [gpd_distrib()].
#'
#' @examples
#' d <- gpd_distrib()
#'
#' # Six at shape zero, and infinite at or above one quarter.
#' kurtosis(d, list(sigma = 1, xi = c(0, 0.2, 0.25, 0.3)))
#'
#' # The four thresholds, from the mean's 1 down to this one's 1/4.
#' th <- list(sigma = 1, xi = 0.3)
#' c(mean = mean(d, th), variance = variance(d, th),
#'   skewness = skewness(d, th), kurtosis = kurtosis(d, th))
#'
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

#' Raw Moments of a Generalized Gamma
#'
#' @description
#' Returns the first four raw moments
#' \eqn{E[Y^k] = a^k\,\Gamma\{(d+k)/p\}/\Gamma(d/p)} for \eqn{k = 1, \ldots, 4}.
#' The four moment methods of the family share this helper and assemble
#' different combinations of the four values.
#'
#' @details
#' The gamma ratio is formed on the log scale, as
#' `exp(k * log(a) + lgamma((d + k) / p) - lgamma(d / p))`, so it stays finite
#' at shapes where either gamma function on its own would overflow.
#'
#' @param a The scale parameter, a positive numeric vector.
#' @param d The first shape parameter, a positive numeric vector.
#' @param p The second shape parameter, a positive numeric vector. Small `p`
#'   pushes \eqn{(d+k)/p} to large arguments, where the log scale earns its
#'   keep.
#'
#' @return A list of four numeric vectors, the raw moments of order 1 to 4,
#'   each recycled to the longest of `a`, `d` and `p`.
#'
#' @seealso [mean.GenGamma1Distrib()], [variance.GenGamma1Distrib()],
#'   [skewness.GenGamma1Distrib()] and [kurtosis.GenGamma1Distrib()] for the
#'   four consumers.
#'
#' @examples
#' # At d = p the family is Weibull, and the first raw moment is its mean.
#' distributions7:::gengamma_raw_moments(2, 3, 3)[[1]]
#' mean(weibull1_distrib(), list(mu = 2, sigma = 3))
#'
#' @keywords internal
gengamma_raw_moments <- function(a, d, p) {
  k0 <- d / p
  lapply(1:4, function(k) exp(k * log(a) + lgamma(k0 + k / p) - lgamma(k0)))
}

#' @title Mean of the Generalized Gamma Distribution
#' @name mean.GenGamma1Distrib
#'
#' @description
#' Closed form: \eqn{E[Y] = a\,\Gamma\{(d+1)/p\}/\Gamma(d/p)}, the first raw
#' moment. None of the three parameters is the mean: \eqn{a} is a scale and
#' \eqn{d} and \eqn{p} are shapes, and the gamma ratio is what carries the
#' shapes into the mean.
#'
#' @details
#' The family nests four the toolkit ships separately, and each is a check on
#' this formula: the gamma at \eqn{p = 1}, the Weibull at \eqn{d = p}, the
#' exponential at \eqn{d = p = 1} and the half-normal at
#' \eqn{a = \sqrt2, d = 1, p = 2}.
#'
#' @section Notation:
#' \eqn{a > 0} is the scale, \eqn{d > 0} and \eqn{p > 0} the two shapes.
#'
#' @param x A `GenGamma1Distrib`, from [gengamma1_distrib()].
#' @param theta A named list with components `a`, `d` and `p`, all positive,
#'   each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length equal to the longest of the
#'   three components.
#'
#' @seealso [variance.GenGamma1Distrib()], [gengamma_raw_moments()] for the
#'   shared quantities, [gengamma1_distrib()].
#'
#' @examples
#' d <- gengamma1_distrib()
#'
#' # At p = 1 the family is a gamma of shape d, whose mean is a d.
#' all.equal(mean(d, list(a = 2, d = 3, p = 1)), 6)
#'
#' # At d = p it is a Weibull of scale a and shape p.
#' all.equal(mean(d, list(a = 2, d = 3, p = 3)),
#'           mean(weibull1_distrib(), list(mu = 2, sigma = 3)))
#'
#' # At a = sqrt(2), d = 1, p = 2 it is the half-normal, of mean sqrt(2/pi).
#' all.equal(mean(d, list(a = sqrt(2), d = 1, p = 2)), sqrt(2 / pi))
#'
#' @keywords internal
S7::method(mean, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])[[1]]
}

#' @title Variance of the Generalized Gamma Distribution
#' @name variance.GenGamma1Distrib
#'
#' @description
#' Closed form from the first two raw moments,
#' \eqn{\operatorname{Var}(Y) = m_2 - m_1^2} with
#' \eqn{m_k = a^k\,\Gamma\{(d+k)/p\}/\Gamma(d/p)}. The scale enters as a
#' square and the two shapes through the gamma ratios.
#'
#' @section Notation:
#' \eqn{a > 0} is the scale, \eqn{d > 0} and \eqn{p > 0} the two shapes, and
#' \eqn{m_k} the \eqn{k}-th raw moment.
#'
#' @param x A `GenGamma1Distrib`, from [gengamma1_distrib()].
#' @param theta A named list with components `a`, `d` and `p`, all positive,
#'   each a numeric vector of length 1 or `n`. Cancellation between the two
#'   raw moments costs digits where the coefficient of variation is small.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length equal to the longest of the
#'   three components.
#'
#' @seealso [mean.GenGamma1Distrib()], [skewness.GenGamma1Distrib()],
#'   [gengamma_raw_moments()], [gengamma1_distrib()].
#'
#' @examples
#' d <- gengamma1_distrib()
#'
#' # At p = 1 the family is a gamma of shape d, whose variance is a^2 d.
#' all.equal(variance(d, list(a = 2, d = 3, p = 1)), 12)
#'
#' # At d = p it agrees with the Weibull.
#' all.equal(variance(d, list(a = 2, d = 3, p = 3)),
#'           variance(weibull1_distrib(), list(mu = 2, sigma = 3)))
#'
#' @keywords internal
S7::method(variance, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])
  m[[2]] - m[[1]]^2
}

#' @title Skewness of the Generalized Gamma Distribution
#' @name skewness.GenGamma1Distrib
#'
#' @description
#' Assembled from the first three raw moments,
#' \eqn{\gamma_1 = (m_3 - 3m_1m_2 + 2m_1^3)/(m_2 - m_1^2)^{3/2}} with
#' \eqn{m_k = a^k\,\Gamma\{(d+k)/p\}/\Gamma(d/p)}. The scale cancels, so the
#' value depends on the two shapes alone; unlike a gamma's it can be negative,
#' which is part of what the second shape parameter buys.
#'
#' @section Notation:
#' \eqn{d > 0} and \eqn{p > 0} are the two shapes and \eqn{m_k} the \eqn{k}-th
#' raw moment. The scale does not enter a standardized moment.
#'
#' @param x A `GenGamma1Distrib`, from [gengamma1_distrib()].
#' @param theta A named list with components `a`, `d` and `p`, all positive,
#'   each a numeric vector of length 1 or `n`. Only `d` and `p` enter the
#'   value.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length equal to the longest of the three
#'   components.
#'
#' @seealso [kurtosis.GenGamma1Distrib()], from the same raw moments;
#'   [gengamma_raw_moments()], [gengamma1_distrib()].
#'
#' @examples
#' d <- gengamma1_distrib()
#'
#' # At d = p it agrees with the Weibull, which changes sign with its shape.
#' all.equal(skewness(d, list(a = 2, d = 3, p = 3)),
#'           skewness(weibull1_distrib(), list(mu = 2, sigma = 3)))
#'
#' # The scale does not enter a standardized moment.
#' skewness(d, list(a = c(0.1, 1, 100), d = 3, p = 1.5))
#'
#' @keywords internal
S7::method(skewness, GenGamma1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- gengamma_raw_moments(theta[[1]], theta[[2]], theta[[3]])
  v <- m[[2]] - m[[1]]^2
  (m[[3]] - 3 * m[[1]] * m[[2]] + 2 * m[[1]]^3) / v^1.5
}

#' @title Excess Kurtosis of the Generalized Gamma Distribution
#' @name kurtosis.GenGamma1Distrib
#'
#' @description
#' Assembled from the first four raw moments,
#' \eqn{\gamma_2 = (m_4 - 4m_1m_3 + 6m_1^2m_2 - 3m_1^4)/(m_2 - m_1^2)^2 - 3}
#' with \eqn{m_k = a^k\,\Gamma\{(d+k)/p\}/\Gamma(d/p)}. The scale cancels, and
#' the two shapes move the skewness and the kurtosis with some freedom, where a
#' gamma ties them by \eqn{\gamma_2 = 3\gamma_1^2/2}.
#'
#' @section Notation:
#' \eqn{d > 0} and \eqn{p > 0} are the two shapes and \eqn{m_k} the \eqn{k}-th
#' raw moment.
#'
#' @param x A `GenGamma1Distrib`, from [gengamma1_distrib()].
#' @param theta A named list with components `a`, `d` and `p`, all positive,
#'   each a numeric vector of length 1 or `n`. Only `d` and `p` enter the
#'   value; the fourth-order combination cancels heavily at small dispersions.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length equal to the longest
#'   of the three components.
#'
#' @seealso [skewness.GenGamma1Distrib()], from the same raw moments;
#'   [gengamma_raw_moments()], [gengamma1_distrib()].
#'
#' @examples
#' d <- gengamma1_distrib()
#'
#' # At d = p it agrees with the Weibull.
#' all.equal(kurtosis(d, list(a = 2, d = 3, p = 3)),
#'           kurtosis(weibull1_distrib(), list(mu = 2, sigma = 3)))
#'
#' # At p = 1 it agrees with the gamma, which ties it to the skewness.
#' th <- list(a = 2, d = 3, p = 1)
#' all.equal(kurtosis(d, th), 1.5 * skewness(d, th)^2)
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. The family carries one parameter and it is
#' the mean, so this method reads it off. It is also the variance, which is the
#' equidispersion the family is defined by.
#'
#' @param x A `PoissonDistrib`, from [poisson_distrib()].
#' @param theta A named list with one component, `mu` (the mean, positive), a
#'   numeric vector of length 1 or `n`. Aligned and validated by name, so a
#'   non-positive value throws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, the length of `theta$mu`.
#'
#' @seealso [variance.PoissonDistrib()], equal to this;
#'   [mean.NegBin2Distrib()] and [mean.NegBin1Distrib()], the overdispersed
#'   count families; [poisson_distrib()].
#'
#' @examples
#' d <- poisson_distrib()
#'
#' # The one parameter is the mean, and it is also the variance.
#' c(mean(d, list(mu = 3)), variance(d, list(mu = 3)))
#'
#' @keywords internal
S7::method(mean, PoissonDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Poisson Distribution
#' @name variance.PoissonDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu}, equal to the mean. This is
#' equidispersion, and it is the assumption every overdispersed count family in
#' the toolkit relaxes: a sample whose variance exceeds its mean is evidence
#' against a Poisson and for a negative binomial or a Poisson-inverse Gaussian.
#'
#' @param x A `PoissonDistrib`, from [poisson_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, the length of `theta$mu`.
#'
#' @seealso [mean.PoissonDistrib()], equal to this;
#'   [variance.NegBin2Distrib()], which adds \eqn{\mu^2/\theta};
#'   [poisson_distrib()].
#'
#' @examples
#' d <- poisson_distrib()
#'
#' # Equal to the mean at every setting.
#' variance(d, list(mu = c(0.5, 3, 100)))
#'
#' # The negative binomial exceeds it and falls onto it as theta grows.
#' variance(negbin2_distrib(), list(mu = 3, theta = c(1, 100, 1e8)))
#'
#' @keywords internal
S7::method(variance, PoissonDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Skewness of the Poisson Distribution
#' @name skewness.PoissonDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = 1/\sqrt\mu}. It is positive at every mean, a
#' count having a floor at zero and no ceiling, and it decays like
#' \eqn{\mu^{-1/2}}, which is the rate at which the family approaches a
#' Gaussian as the counts grow.
#'
#' @param x A `PoissonDistrib`, from [poisson_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, the length of `theta$mu`, positive throughout.
#'
#' @seealso [kurtosis.PoissonDistrib()], which decays twice as fast;
#'   [skewness.NegBin2Distrib()], which tends to this;
#'   [poisson_distrib()].
#'
#' @examples
#' d <- poisson_distrib()
#'
#' # One over the square root of the mean.
#' all.equal(skewness(d, list(mu = 4)), 0.5)
#'
#' # It vanishes as the counts grow.
#' skewness(d, list(mu = c(1, 25, 10000)))
#'
#' @keywords internal
S7::method(skewness, PoissonDistrib) <- function(x, theta, ...) {
  1 / sqrt(align_theta(x, theta)[[1]])
}

#' @title Excess Kurtosis of the Poisson Distribution
#' @name kurtosis.PoissonDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 1/\mu}, the excess over the Gaussian. It is
#' positive at every mean and decays like \eqn{\mu^{-1}}, twice as fast as the
#' skewness, so a Poisson of moderate mean looks Gaussian in its tails before
#' it looks symmetric.
#'
#' @param x A `PoissonDistrib`, from [poisson_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, the length of `theta$mu`,
#'   positive throughout.
#'
#' @seealso [skewness.PoissonDistrib()]; [kurtosis.NegBin2Distrib()], which
#'   tends to this; [poisson_distrib()].
#'
#' @examples
#' d <- poisson_distrib()
#'
#' # One over the mean, falling twice as fast as the skewness.
#' rbind(skewness = skewness(d, list(mu = c(1, 25, 10000))),
#'       kurtosis = kurtosis(d, list(mu = c(1, 25, 10000))))
#'
#' @keywords internal
S7::method(kurtosis, PoissonDistrib) <- function(x, theta, ...) {
  1 / align_theta(x, theta)[[1]]
}

# --- bernoulli -------------------------------------------------------------

#' @title Mean of the Bernoulli Distribution
#' @name mean.BernoulliDistrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}, the success probability. The response takes
#' the values 0 and 1, so the mean is the probability of the second and every
#' other moment is a function of it: a one-parameter family has nothing else
#' for them to depend on.
#'
#' @param x A `BernoulliDistrib`, from [bernoulli_distrib()].
#' @param theta A named list with one component, `mu` (the success
#'   probability, strictly between 0 and 1), a numeric vector of length 1 or
#'   `n`. Aligned and validated by name, so a value at or outside the unit
#'   interval throws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, the length of `theta$mu`, strictly inside
#'   \eqn{(0,1)}.
#'
#' @seealso [variance.BernoulliDistrib()]; [mean.BinomialDistrib()], the sum of
#'   `size` of these; [bernoulli_distrib()].
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # The one parameter is the mean.
#' mean(d, list(mu = c(0.1, 0.5, 0.9)))
#'
#' @keywords internal
S7::method(mean, BernoulliDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Bernoulli Distribution
#' @name variance.BernoulliDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu(1-\mu)}. It is largest at
#' \eqn{\mu = 1/2}, where it is \eqn{1/4}, and goes to zero at either end, the
#' response becoming deterministic there. That quarter is the largest variance
#' any distribution on \eqn{\{0,1\}} can have.
#'
#' @param x A `BernoulliDistrib`, from [bernoulli_distrib()].
#' @param theta A named list with one component, `mu` (strictly between 0 and
#'   1), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, the length of `theta$mu`, in
#'   \eqn{(0, 1/4]}.
#'
#' @seealso [mean.BernoulliDistrib()]; [variance.BinomialDistrib()], which is
#'   `size` times this; [bernoulli_distrib()].
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # Largest at one half, vanishing at either end.
#' variance(d, list(mu = c(0.01, 0.3, 0.5, 0.7, 0.99)))
#'
#' @keywords internal
S7::method(variance, BernoulliDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  p * (1 - p)
}

#' @title Skewness of the Bernoulli Distribution
#' @name skewness.BernoulliDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = (1-2\mu)/\sqrt{\mu(1-\mu)}}. It is zero at
#' \eqn{\mu = 1/2}, positive below and negative above, and it diverges at
#' either end of the unit interval, where almost all the mass sits on one of
#' the two values.
#'
#' @param x A `BernoulliDistrib`, from [bernoulli_distrib()].
#' @param theta A named list with one component, `mu` (strictly between 0 and
#'   1), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, the length of `theta$mu`.
#'
#' @seealso [kurtosis.BernoulliDistrib()]; [skewness.BinomialDistrib()], which
#'   is this divided by \eqn{\sqrt{n}}; [bernoulli_distrib()].
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # The published form, written out.
#' all.equal(skewness(d, list(mu = 0.3)), (1 - 0.6) / sqrt(0.3 * 0.7))
#'
#' # Zero at one half, and it diverges towards either end.
#' skewness(d, list(mu = c(0.01, 0.3, 0.5, 0.7, 0.99)))
#'
#' @keywords internal
S7::method(skewness, BernoulliDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  (1 - 2 * p) / sqrt(p * (1 - p))
}

#' @title Excess Kurtosis of the Bernoulli Distribution
#' @name kurtosis.BernoulliDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = \{1 - 6\mu(1-\mu)\}/\{\mu(1-\mu)\}}, the excess
#' over the Gaussian. It reaches its minimum \eqn{-2} at \eqn{\mu = 1/2}, which
#' is the smallest excess kurtosis any distribution attains, and diverges at
#' either end of the unit interval.
#'
#' @param x A `BernoulliDistrib`, from [bernoulli_distrib()].
#' @param theta A named list with one component, `mu` (strictly between 0 and
#'   1), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, the length of `theta$mu`, at or
#'   above \eqn{-2}.
#'
#' @seealso [skewness.BernoulliDistrib()]; [kurtosis.BinomialDistrib()], which
#'   is this divided by \eqn{n}; [bernoulli_distrib()].
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # The published form, written out.
#' all.equal(kurtosis(d, list(mu = 0.3)), (1 - 6 * 0.21) / 0.21)
#'
#' # Minus two at one half, the lower bound for any distribution.
#' kurtosis(d, list(mu = c(0.1, 0.5, 0.9)))
#'
#' @keywords internal
S7::method(kurtosis, BernoulliDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  v <- p * (1 - p)
  (1 - 6 * v) / v
}

# --- binomial --------------------------------------------------------------

#' @title Mean of the Binomial Distribution
#' @name mean.BinomialDistrib
#'
#' @description
#' Closed form: \eqn{E[Y] = n\mu}, with \eqn{n} the number of trials and
#' \eqn{\mu} the success probability. The trial count is a property of the
#' object and not a parameter, so it is read from `x@size` and does not appear
#' in `theta`.
#'
#' @section Notation:
#' \eqn{n} is the number of trials, held on the object, and \eqn{\mu \in (0,1)}
#' the success probability.
#'
#' @param x A `BinomialDistrib`, from [binomial_distrib()], carrying the trial
#'   count in its `size` property.
#' @param theta A named list with one component, `mu` (the success
#'   probability, strictly between 0 and 1), a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, the length of `theta$mu`, in
#'   \eqn{(0, n)}.
#'
#' @seealso [variance.BinomialDistrib()]; [mean.BernoulliDistrib()], the
#'   one-trial case; [mean.BetaBinom1Distrib()], the overdispersed extension;
#'   [binomial_distrib()].
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#'
#' # Ten trials at probability 0.3.
#' mean(d, list(mu = 0.3))
#'
#' @keywords internal
S7::method(mean, BinomialDistrib) <- function(x, theta, ...) {
  x@size * align_theta(x, theta)[[1]]
}

#' @title Variance of the Binomial Distribution
#' @name variance.BinomialDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = n\mu(1-\mu)}. The trials are
#' independent, so the variance is the trial count times the Bernoulli
#' variance; a sample of counts out of \eqn{n} whose variance exceeds this is
#' overdispersed, and a beta-binomial is the family that carries the excess.
#'
#' @param x A `BinomialDistrib`, from [binomial_distrib()], carrying the trial
#'   count in its `size` property.
#' @param theta A named list with one component, `mu` (strictly between 0 and
#'   1), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, the length of `theta$mu`.
#'
#' @seealso [mean.BinomialDistrib()]; [variance.BernoulliDistrib()], the
#'   one-trial case; [variance.BetaBinom1Distrib()], which inflates this;
#'   [binomial_distrib()].
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#'
#' # Ten times the Bernoulli variance.
#' all.equal(variance(d, list(mu = 0.3)),
#'           10 * variance(bernoulli_distrib(), list(mu = 0.3)))
#'
#' @keywords internal
S7::method(variance, BinomialDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  x@size * p * (1 - p)
}

#' @title Skewness of the Binomial Distribution
#' @name skewness.BinomialDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = (1-2\mu)/\sqrt{n\mu(1-\mu)}}. It is the
#' Bernoulli's divided by \eqn{\sqrt n}, so it vanishes as the trials
#' accumulate, at the rate the central limit theorem gives. It is zero at
#' \eqn{\mu = 1/2} at every trial count.
#'
#' @param x A `BinomialDistrib`, from [binomial_distrib()], carrying the trial
#'   count in its `size` property.
#' @param theta A named list with one component, `mu` (strictly between 0 and
#'   1), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, the length of `theta$mu`.
#'
#' @seealso [kurtosis.BinomialDistrib()], which vanishes as \eqn{1/n};
#'   [skewness.BernoulliDistrib()], the one-trial case;
#'   [binomial_distrib()].
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#'
#' # The Bernoulli's divided by the square root of the trial count.
#' all.equal(skewness(d, list(mu = 0.3)),
#'           skewness(bernoulli_distrib(), list(mu = 0.3)) / sqrt(10))
#'
#' # Zero at one half, at any trial count.
#' skewness(d, list(mu = 0.5))
#'
#' @keywords internal
S7::method(skewness, BinomialDistrib) <- function(x, theta, ...) {
  p <- align_theta(x, theta)[[1]]
  (1 - 2 * p) / sqrt(x@size * p * (1 - p))
}

#' @title Excess Kurtosis of the Binomial Distribution
#' @name kurtosis.BinomialDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = \{1 - 6\mu(1-\mu)\}/\{n\mu(1-\mu)\}}, the
#' excess over the Gaussian. It is the Bernoulli's divided by \eqn{n}, so it
#' vanishes twice as fast as the skewness and can take either sign: it is
#' negative over the middle of the unit interval, where the support is bounded
#' at both ends and the tails are lighter than a Gaussian's.
#'
#' @param x A `BinomialDistrib`, from [binomial_distrib()], carrying the trial
#'   count in its `size` property.
#' @param theta A named list with one component, `mu` (strictly between 0 and
#'   1), a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, the length of `theta$mu`.
#'
#' @seealso [skewness.BinomialDistrib()]; [kurtosis.BernoulliDistrib()], the
#'   one-trial case; [binomial_distrib()].
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#'
#' # The Bernoulli's divided by the trial count.
#' all.equal(kurtosis(d, list(mu = 0.3)),
#'           kurtosis(bernoulli_distrib(), list(mu = 0.3)) / 10)
#'
#' # Negative over the middle, positive towards either end.
#' round(kurtosis(d, list(mu = c(0.05, 0.3, 0.5, 0.95))), 4)
#'
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
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. The family is parametrized by its mean, the
#' expected number of failures before the first success, and the success
#' probability it implies is \eqn{1/(1+\mu)}.
#'
#' @param x A `GeometricDistrib`, from [geometric_distrib()].
#' @param theta A named list with one component, `mu` (the mean, positive), a
#'   numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, the length of `theta$mu`.
#'
#' @seealso [variance.GeometricDistrib()], which exceeds this;
#'   [mean.NegBin2Distrib()], of which this is the \eqn{\theta = 1} case;
#'   [geometric_distrib()].
#'
#' @examples
#' d <- geometric_distrib()
#'
#' # The one parameter is the mean.
#' mean(d, list(mu = c(0.5, 3, 20)))
#'
#' @keywords internal
S7::method(mean, GeometricDistrib) <- function(x, theta, ...) {
  align_theta(x, theta)[[1]]
}

#' @title Variance of the Geometric Distribution
#' @name variance.GeometricDistrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu(1+\mu)}. It exceeds the mean
#' at every parameter value, so the family is overdispersed relative to a
#' Poisson, and it is exactly the negative binomial's variance at
#' \eqn{\theta = 1}.
#'
#' @param x A `GeometricDistrib`, from [geometric_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, the length of `theta$mu`.
#'
#' @seealso [mean.GeometricDistrib()]; [variance.NegBin2Distrib()], of which
#'   this is the \eqn{\theta = 1} case; [variance.PoissonDistrib()], the
#'   equidispersed comparison; [geometric_distrib()].
#'
#' @examples
#' d <- geometric_distrib()
#'
#' # Mean times one plus the mean.
#' all.equal(variance(d, list(mu = 3)), 12)
#'
#' # The negative binomial at theta = 1 is the same law.
#' all.equal(variance(d, list(mu = 3)),
#'           variance(negbin2_distrib(), list(mu = 3, theta = 1)))
#'
#' @keywords internal
S7::method(variance, GeometricDistrib) <- function(x, theta, ...) {
  mu <- align_theta(x, theta)[[1]]
  mu * (1 + mu)
}

#' @title Skewness of the Geometric Distribution
#' @name skewness.GeometricDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = (1+2\mu)/\sqrt{\mu(1+\mu)}}. It is positive at
#' every mean and, unlike a Poisson's, it does not vanish as the counts grow:
#' it tends to 2 from above, the exponential's value, which is the continuous
#' limit of the family.
#'
#' @param x A `GeometricDistrib`, from [geometric_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, the length of `theta$mu`, above 2 throughout.
#'
#' @seealso [kurtosis.GeometricDistrib()], which tends to 6;
#'   [skewness.ExponentialDistrib()], the continuous limit;
#'   [geometric_distrib()].
#'
#' @examples
#' d <- geometric_distrib()
#'
#' # The published form, written out.
#' all.equal(skewness(d, list(mu = 3)), 7 / sqrt(12))
#'
#' # It tends to the exponential's 2 as the mean grows.
#' skewness(d, list(mu = c(1, 10, 1000)))
#'
#' @keywords internal
S7::method(skewness, GeometricDistrib) <- function(x, theta, ...) {
  mu <- align_theta(x, theta)[[1]]
  (1 + 2 * mu) / sqrt(mu * (1 + mu))
}

#' @title Excess Kurtosis of the Geometric Distribution
#' @name kurtosis.GeometricDistrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 6 + 1/\{\mu(1+\mu)\}}, the excess over the
#' Gaussian. It stays above 6 at every mean and tends to 6 as the counts grow,
#' the exponential's value; the family is heavy-tailed at every parameter
#' setting and does not become Gaussian as a Poisson does.
#'
#' @param x A `GeometricDistrib`, from [geometric_distrib()].
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, the length of `theta$mu`, above
#'   6 throughout.
#'
#' @seealso [skewness.GeometricDistrib()], which tends to 2;
#'   [kurtosis.ExponentialDistrib()], the continuous limit;
#'   [kurtosis.PoissonDistrib()], which does vanish; [geometric_distrib()].
#'
#' @examples
#' d <- geometric_distrib()
#'
#' # The published form, written out.
#' all.equal(kurtosis(d, list(mu = 3)), 6 + 1 / 12)
#'
#' # It stays above six where a Poisson's excess kurtosis vanishes.
#' rbind(geometric = kurtosis(d, list(mu = c(1, 10, 1000))),
#'       poisson   = kurtosis(poisson_distrib(), list(mu = c(1, 10, 1000))))
#'
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

#' @title Mean of the Negative Binomial Distribution in the NB1 Parametrization
#' @name mean.NegBin1Distrib
#'
#' @description
#' Closed form: \eqn{E[Y] = \mu}. This parametrization carries the mean as its
#' first parameter, so the method reads it off. What distinguishes NB1 from NB2
#' is the variance function and not the mean.
#'
#' @param x A `NegBin1Distrib`, from [negbin1_distrib()].
#' @param theta A named list with components `mu` (the mean, positive) and
#'   `theta` (the dispersion, positive), each a numeric vector of length 1 or
#'   `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$theta))`.
#'
#' @seealso [variance.NegBin1Distrib()], which is linear in the mean;
#'   [variance.NegBin2Distrib()], which is quadratic; [negbin1_distrib()].
#'
#' @examples
#' d <- negbin1_distrib()
#'
#' # The first parameter is the mean, and the dispersion does not move it.
#' mean(d, list(mu = c(1, 2, 3), theta = 2))
#'
#' @keywords internal
S7::method(mean, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 2L, 0) + theta[[1]]
}

#' @title Variance of the Negative Binomial Distribution in the NB1 Parametrization
#' @name variance.NegBin1Distrib
#'
#' @description
#' Closed form: \eqn{\operatorname{Var}(Y) = \mu(1+\theta)}, linear in the
#' mean. This is what NB1 means: the variance is a fixed multiple of the mean,
#' where NB2 makes it \eqn{\mu + \mu^2/\theta} and so quadratic. The two are
#' different models of overdispersion and are not reparametrizations of each
#' other.
#'
#' @details
#' The distinction follows Cameron and Trivedi's numbering, and it matters for
#' regression: under NB1 the dispersion relative to a Poisson is the same at
#' every fitted mean, and under NB2 it grows with the mean.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean and \eqn{\theta > 0} the dispersion.
#'
#' @param x A `NegBin1Distrib`, from [negbin1_distrib()].
#' @param theta A named list with components `mu` (positive) and `theta`
#'   (positive), each a numeric vector of length 1 or `n`. The Poisson limit is
#'   \eqn{\theta \to 0} here, the opposite end from NB2's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$theta))`.
#'
#' @references
#' Cameron, A. C. and Trivedi, P. K. (2013). *Regression Analysis of Count
#' Data*, 2nd edition. Cambridge University Press.
#'
#' @seealso [mean.NegBin1Distrib()]; [variance.NegBin2Distrib()], the quadratic
#'   alternative; [variance.PoissonDistrib()], the \eqn{\theta \to 0} limit;
#'   [negbin1_distrib()].
#'
#' @examples
#' d <- negbin1_distrib()
#'
#' # Linear in the mean: the ratio to the Poisson variance is constant.
#' variance(d, list(mu = c(1, 10, 100), theta = 2)) / c(1, 10, 100)
#'
#' # NB2's ratio grows with the mean instead.
#' variance(negbin2_distrib(), list(mu = c(1, 10, 100), theta = 2)) / c(1, 10, 100)
#'
#' @keywords internal
S7::method(variance, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] * (1 + theta[[2]])
}

#' @title Skewness of the Negative Binomial Distribution in the NB1 Parametrization
#' @name skewness.NegBin1Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_1 = (1+2\theta)/\sqrt{\mu(1+\theta)}}. It is
#' positive at every parameter value and tends to the Poisson's
#' \eqn{\mu^{-1/2}} as the dispersion goes to zero, which is where NB1
#' approaches a Poisson.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean and \eqn{\theta > 0} the dispersion.
#'
#' @param x A `NegBin1Distrib`, from [negbin1_distrib()].
#' @param theta A named list with components `mu` (positive) and `theta`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$theta))`, positive throughout.
#'
#' @seealso [kurtosis.NegBin1Distrib()]; [skewness.NegBin2Distrib()], the other
#'   parametrization; [skewness.PoissonDistrib()], the limit;
#'   [negbin1_distrib()].
#'
#' @examples
#' d <- negbin1_distrib()
#'
#' # The published form, written out.
#' all.equal(skewness(d, list(mu = 4, theta = 2)), 5 / sqrt(4 * 3))
#'
#' # It falls onto the Poisson's mu^(-1/2) = 0.5 as the dispersion vanishes.
#' skewness(d, list(mu = 4, theta = c(1, 0.01, 1e-8)))
#'
#' @keywords internal
S7::method(skewness, NegBin1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  mu <- theta[[1]]
  th <- theta[[2]]
  (1 + 2 * th) / sqrt(mu * (1 + th))
}

#' @title Excess Kurtosis of the Negative Binomial Distribution in the NB1 Parametrization
#' @name kurtosis.NegBin1Distrib
#'
#' @description
#' Closed form: \eqn{\gamma_2 = 6\theta/\mu + 1/\{\mu(1+\theta)\}}, the excess
#' over the Gaussian. Both terms are positive, so the family is always
#' leptokurtic, and the expression tends to the Poisson's \eqn{1/\mu} as the
#' dispersion goes to zero.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean and \eqn{\theta > 0} the dispersion.
#'
#' @param x A `NegBin1Distrib`, from [negbin1_distrib()].
#' @param theta A named list with components `mu` (positive) and `theta`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$theta))`, positive throughout.
#'
#' @seealso [skewness.NegBin1Distrib()]; [kurtosis.NegBin2Distrib()], the other
#'   parametrization; [kurtosis.PoissonDistrib()], the limit;
#'   [negbin1_distrib()].
#'
#' @examples
#' d <- negbin1_distrib()
#'
#' # The published form, written out.
#' all.equal(kurtosis(d, list(mu = 4, theta = 2)), 12 / 4 + 1 / 12)
#'
#' # It falls onto the Poisson's 1 / mu = 0.25 as the dispersion vanishes.
#' kurtosis(d, list(mu = 4, theta = c(1, 0.01, 1e-8)))
#'
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

#' Falling Factorial Moments of a Beta-Binomial
#'
#' @description
#' Returns the first four falling factorial moments
#' \deqn{E[Y^{(k)}] = n^{(k)} \prod_{j=0}^{k-1} \frac{a+j}{a+b+j},
#'       \qquad k = 1, \ldots, 4,}
#' with \eqn{x^{(k)} = x(x-1)\cdots(x-k+1)}. These are the quantities a
#' beta-binomial has in closed form; the raw and central moments follow from
#' them through [central_from_factorial()].
#'
#' @section Notation:
#' \eqn{a > 0} and \eqn{b > 0} are the two shapes of the mixing beta,
#' \eqn{n} the number of trials, and \eqn{x^{(k)}} the falling factorial.
#'
#' @param a The first shape of the mixing beta, a single positive number.
#' @param b The second shape, a single positive number.
#' @param n The number of trials, a single non-negative whole number. Factorial
#'   moments of order above `n` are zero, the falling factorial \eqn{n^{(k)}}
#'   vanishing there.
#'
#' @return An **unnamed** list of four numbers, the falling factorial moments of
#'   order 1 to 4 in that order, reached by position.
#'
#' @seealso [central_from_factorial()], which converts these;
#'   [betabinom_central()], the wrapper the moment methods call.
#'
#' @examples
#' # At a = b the mixing beta is symmetric, so the mean is half the trials.
#' distributions7:::betabinom_factorial_moments(2, 2, 10)[[1]]
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
#' Converts the first four falling factorial moments into the mean and the
#' second, third and fourth central moments. The raw moments come first,
#' through the Stirling numbers of the second kind written out:
#' \deqn{m_1 = f_1,\quad m_2 = f_2 + f_1,\quad m_3 = f_3 + 3f_2 + f_1,\quad
#'       m_4 = f_4 + 6f_3 + 7f_2 + f_1,}
#' and the central moments follow by the usual binomial expansion about
#' \eqn{m_1}.
#'
#' @param f An unnamed list of four, the falling factorial moments of order 1
#'   to 4 in that order, as [betabinom_factorial_moments()] returns. They are
#'   read by position.
#'
#' @return A named list with `mean` and the central moments `c2`, `c3` and
#'   `c4`, each a numeric vector the length of the inputs.
#'
#' @seealso [betabinom_factorial_moments()] for the input;
#'   [betabinom_central()], which chains the two.
#'
#' @examples
#' f <- distributions7:::betabinom_factorial_moments(2, 2, 10)
#' distributions7:::central_from_factorial(f)$mean
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

#' Mean and Central Moments of a Beta-Binomial
#'
#' @description
#' The mean and the second, third and fourth central moments of a
#' beta-binomial, one value per parameter setting. It maps the mean-dispersion
#' parametrization onto the two shapes, \eqn{a = \mu/\sigma} and
#' \eqn{b = (1-\mu)/\sigma}, then chains [betabinom_factorial_moments()] into
#' [central_from_factorial()] once per setting and stacks the results. The four
#' moment methods of the family share it.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the success probability, \eqn{\sigma > 0} the
#' dispersion and \eqn{n} the number of trials.
#'
#' @param mu The success probability, a numeric vector strictly inside
#'   \eqn{(0,1)}.
#' @param sigma The dispersion, a positive numeric vector. As it goes to zero
#'   the mixing beta concentrates and the family tends to a binomial.
#' @param n The number of trials, a non-negative whole number.
#'
#' @return A named list with `mean` and the central moments `c2`, `c3` and
#'   `c4`, each a numeric vector as long as the longer of `mu` and `sigma`.
#'
#' @seealso [variance.BetaBinom1Distrib()], [skewness.BetaBinom1Distrib()] and
#'   [kurtosis.BetaBinom1Distrib()] for the consumers.
#'
#' @examples
#' # The mean is n mu, whatever the dispersion.
#' distributions7:::betabinom_central(0.3, 0.5, 10)$mean
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
#'
#' @description
#' Closed form: \eqn{E[Y] = n\mu}, the same as a binomial's. Mixing the success
#' probability over a beta leaves the mean where it was, the beta's own mean
#' being \eqn{\mu}; the dispersion shows in the variance and above, not here.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the success probability, \eqn{\sigma > 0} the
#' dispersion and \eqn{n} the number of trials, held on the object.
#'
#' @param x A `BetaBinom1Distrib`, from [betabinom1_distrib()], carrying the
#'   trial count in its `size` property.
#' @param theta A named list with components `mu` (the success probability,
#'   strictly between 0 and 1) and `sigma` (the dispersion, positive), each a
#'   numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [variance.BetaBinom1Distrib()], where the dispersion does enter;
#'   [mean.BinomialDistrib()], which this equals; [betabinom1_distrib()].
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#'
#' # n mu, and the dispersion does not move it.
#' mean(d, list(mu = 0.3, sigma = c(0.01, 0.5, 5)))
#'
#' @keywords internal
S7::method(mean, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  betabinom_central(theta[[1]], theta[[2]], x@size)$mean
}

#' @title Variance of the Beta-Binomial Distribution
#' @name variance.BetaBinom1Distrib
#'
#' @description
#' Closed form:
#' \deqn{\operatorname{Var}(Y) = n\mu(1-\mu)\,\frac{1 + n\sigma}{1 + \sigma}.}
#' The first factor is the binomial variance and the second is the
#' overdispersion, at least 1 for every \eqn{n \ge 1} and growing towards
#' \eqn{n} as the dispersion rises. At \eqn{\sigma \to 0} the family is
#' binomial; at large \eqn{\sigma} the variance approaches
#' \eqn{n^2\mu(1-\mu)}, which is a Bernoulli scaled by the trial count and is
#' the largest a distribution on \eqn{\{0, \ldots, n\}} with that mean can
#' have.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the success probability, \eqn{\sigma > 0} the
#' dispersion and \eqn{n} the number of trials.
#'
#' @param x A `BetaBinom1Distrib`, from [betabinom1_distrib()], carrying the
#'   trial count in its `size` property.
#' @param theta A named list with components `mu` (strictly between 0 and 1)
#'   and `sigma` (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [mean.BetaBinom1Distrib()]; [variance.BinomialDistrib()], the
#'   \eqn{\sigma \to 0} limit; [betabinom_central()];
#'   [betabinom1_distrib()].
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#'
#' # The published form, written out.
#' all.equal(variance(d, list(mu = 0.3, sigma = 0.5)),
#'           10 * 0.3 * 0.7 * (1 + 10 * 0.5) / (1 + 0.5))
#'
#' # It falls onto the binomial variance as the dispersion vanishes.
#' c(betabinomial = variance(d, list(mu = 0.3, sigma = 1e-8)),
#'   binomial     = variance(binomial_distrib(size = 10), list(mu = 0.3)))
#'
#' @keywords internal
S7::method(variance, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  betabinom_central(theta[[1]], theta[[2]], x@size)$c2
}

#' @title Skewness of the Beta-Binomial Distribution
#' @name skewness.BetaBinom1Distrib
#'
#' @description
#' Assembled from the falling factorial moments, \eqn{\gamma_1 = c_3/c_2^{3/2}}
#' with \eqn{c_k} the central moments [betabinom_central()] returns. It is zero
#' at \eqn{\mu = 1/2} whatever the dispersion, the mixing beta being symmetric
#' there, and it takes the sign of \eqn{1 - 2\mu} elsewhere.
#'
#' @details
#' The falling factorial moments are the quantities a beta-binomial has in
#' closed form: each is a product of \eqn{k} ratios and needs no sum over the
#' support. The raw and then the central moments follow from them by two
#' written-out linear maps.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the success probability, \eqn{\sigma > 0} the
#' dispersion, \eqn{n} the number of trials and \eqn{c_k} the \eqn{k}-th
#' central moment.
#'
#' @param x A `BetaBinom1Distrib`, from [betabinom1_distrib()], carrying the
#'   trial count in its `size` property.
#' @param theta A named list with components `mu` (strictly between 0 and 1)
#'   and `sigma` (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [kurtosis.BetaBinom1Distrib()], from the same moments;
#'   [skewness.BinomialDistrib()], the \eqn{\sigma \to 0} limit;
#'   [betabinom_central()]; [betabinom1_distrib()].
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#'
#' # Zero at a success probability of one half, at any dispersion.
#' round(skewness(d, list(mu = 0.5, sigma = c(0.1, 1, 10))), 12)
#'
#' # It falls onto the binomial's as the dispersion vanishes.
#' c(betabinomial = skewness(d, list(mu = 0.3, sigma = 1e-8)),
#'   binomial     = skewness(binomial_distrib(size = 10), list(mu = 0.3)))
#'
#' @keywords internal
S7::method(skewness, BetaBinom1Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  m <- betabinom_central(theta[[1]], theta[[2]], x@size)
  m$c3 / m$c2^1.5
}

#' @title Excess Kurtosis of the Beta-Binomial Distribution
#' @name kurtosis.BetaBinom1Distrib
#'
#' @description
#' Assembled from the falling factorial moments,
#' \eqn{\gamma_2 = c_4/c_2^2 - 3} with \eqn{c_k} the central moments
#' [betabinom_central()] returns. The success probability is what sets its
#' sign. Near \eqn{\mu = 1/2} it is negative and falls
#' towards \eqn{-2} as the dispersion grows, the mass splitting between the two
#' endpoints and the law approaching a scaled Bernoulli; near either end it is
#' positive and large, the long tail towards the middle dominating.
#'
#' @section Notation:
#' \eqn{\mu \in (0,1)} is the success probability, \eqn{\sigma > 0} the
#' dispersion, \eqn{n} the number of trials and \eqn{c_k} the \eqn{k}-th
#' central moment.
#'
#' @param x A `BetaBinom1Distrib`, from [betabinom1_distrib()], carrying the
#'   trial count in its `size` property.
#' @param theta A named list with components `mu` (strictly between 0 and 1)
#'   and `sigma` (positive), each a numeric vector of length 1 or `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of length
#'   `max(length(theta$mu), length(theta$sigma))`.
#'
#' @seealso [skewness.BetaBinom1Distrib()], from the same moments;
#'   [kurtosis.BinomialDistrib()], the \eqn{\sigma \to 0} limit;
#'   [betabinom_central()]; [betabinom1_distrib()].
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#'
#' # At a central success probability it falls towards -2 with the dispersion.
#' round(kurtosis(d, list(mu = 0.5, sigma = c(0.1, 1, 100))), 4)
#'
#' # At an extreme one it is positive and large.
#' round(kurtosis(d, list(mu = 0.05, sigma = c(0.01, 0.1, 1))), 4)
#'
#' # It falls onto the binomial's as the dispersion vanishes.
#' c(betabinomial = kurtosis(d, list(mu = 0.3, sigma = 1e-8)),
#'   binomial     = kurtosis(binomial_distrib(size = 10), list(mu = 0.3)))
#'
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
