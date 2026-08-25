#' @include distrib.R generics.R skewnormal1_distrib.R reparametrize.R moments.R
NULL

# The skew normal in Azzalini's CENTERED parametrization: the mean, the standard
# deviation and the skewness itself, rather than the location, the scale and
# the shape.
#
# This is a family of its own and not a reparametrize() of skewnormal1, for two
# reasons that both come from the map passing through a cube root. It carries a
# sign, which a jet cannot take; and its derivative is unbounded as the
# skewness goes to zero, so what makes the parametrization worth having is a
# CANCELLATION between terms that individually diverge. Measured at gamma1 =
# 1e-4, d alpha / d gamma1 is 258 while the variance of the score in gamma1 is
# 0.158, the same value it has at gamma1 = 0.05: the divergence cancels.

#' @title The Constant Behind the Centered Parametrization
#'
#' @description
#' Returns \eqn{b = \sqrt{2/\pi} \approx 0.7978846}, which is
#' \eqn{E[|Z|]} for a standard Gaussian \eqn{Z}. It is the one constant the
#' skew normal's first moment introduces, and it appears in every quantity of
#' the centered parametrization: the mean is \eqn{\xi + \omega b\delta}, the
#' variance \eqn{\omega^2(1-b^2\delta^2)}, and the largest reachable skewness
#' is built from \eqn{b/\sqrt{1-b^2}}.
#'
#' @return A single number.
#'
#' @seealso [sn_max_skew()], which is written in terms of it, and
#'   [skewnormal2_distrib()] for the parametrization.
#'
#' @examples
#' distributions7:::sn_b()
#'
#' # It is the mean of the absolute value of a standard Gaussian.
#' set.seed(1)
#' c(constant = distributions7:::sn_b(), sample = mean(abs(rnorm(1e6))))
#'
#' @keywords internal
sn_b <- function() sqrt(2 / pi)

#' @title The Largest Skewness a Skew Normal Can Reach
#'
#' @description
#' Returns \eqn{\sup|\gamma_1|} over the skew normal family,
#' \deqn{\dfrac{4-\pi}{2}\left(\dfrac{b}{\sqrt{1-b^2}}\right)^3
#'       \approx 0.9952717,}
#' with \eqn{b = \sqrt{2/\pi}}. The bound is attained only in the limit
#' \eqn{\alpha \to \pm\infty}, where the family degenerates to a half-normal.
#'
#' @details
#' A skewness beyond it belongs to no skew normal, so
#' [skewnormal2_distrib()] bounds `gamma1` there and gives it a
#' [linkfunctions7::bounded_link()] over the open interval. Without the bound
#' the centered-to-direct map would return a `NaN` several frames down, where
#' the reason for it is no longer visible.
#'
#' The ceiling is approached slowly: measured, the skewness is 0.9556 at
#' \eqn{\alpha = 10} and 0.99527 at \eqn{\alpha = 10^4}. Data skewer than this
#' needs [skewt_distrib()].
#'
#' @return A single number.
#'
#' @seealso [sn_b()] for the constant it is built from,
#'   [skewnormal2_distrib()] for the parameter it bounds, and
#'   [skewt_distrib()] for the family that goes further.
#'
#' @examples
#' distributions7:::sn_max_skew()
#'
#' # The direct parametrization approaches it and does not pass it.
#' d1 <- skewnormal1_distrib()
#' vapply(c(10, 1e4, 1e8),
#'        function(a) skewness(d1, list(mu = 0, sigma = 1, alpha = a)), 0)
#'
#' # It is where the centered parametrization's link is bounded.
#' skewnormal2_distrib()@params_bounds$gamma1
#'
#' @keywords internal
sn_max_skew <- function() {
  b <- sn_b()
  (4 - pi) / 2 * (b / sqrt(1 - b^2))^3
}

#' @title From the Centered Parameters to the Direct Ones
#'
#' @description
#' Maps \eqn{(\mu, \sigma, \gamma_1)}, the mean, the standard deviation and the
#' skewness, to \eqn{(\xi, \omega, \alpha)}, the location, the scale and the
#' shape that [skewnormal1_distrib()] takes. Every probability function of the
#' centered family calls it and then delegates to the direct one.
#'
#' @details
#' With \eqn{b = \sqrt{2/\pi}},
#' \deqn{c = \mathrm{sign}(\gamma_1)
#'           \left(\dfrac{2|\gamma_1|}{4-\pi}\right)^{1/3}, \qquad
#'       \mu_z = \dfrac{c}{\sqrt{1+c^2}}, \qquad
#'       \delta = \dfrac{\mu_z}{b}, \qquad
#'       \alpha = \dfrac{\delta}{\sqrt{1-\delta^2}},}
#' and then \eqn{\omega = \sigma/\sqrt{1-\mu_z^2}} and
#' \eqn{\xi = \mu - \omega\mu_z}.
#'
#' The caller supplies the sign, so the body reads
#' \eqn{s\,(2s\gamma_1/(4-\pi))^{1/3}} with no `abs()` in it. Away from zero
#' the sign is locally constant, so the expression is exact and differentiable
#' as written, and the derivative tables of [md_skewnormal2()] differentiate it
#' directly.
#'
#' @param mu,sigma,gamma1 The centered parameters: the mean, the standard
#'   deviation and the skewness, each a numeric vector. `sigma` must be
#'   positive and `gamma1` must lie strictly inside
#'   \eqn{(-0.9952717, 0.9952717)}; nothing is validated here.
#' @param s The sign of `gamma1`, \eqn{\pm 1}, taken by the caller from its
#'   plain value.
#'
#' @return A named list with `mu`, `sigma` and `alpha`, the direct parameters,
#'   each of the length of the recycled inputs. The names are the parent's, so
#'   the result can be passed to [skewnormal1_distrib()]'s methods as they
#'   stand; `mu` there is the location \eqn{\xi} and `sigma` the scale
#'   \eqn{\omega}.
#'
#' @seealso [sn2_theta()], which supplies the sign and calls this;
#'   [md_skewnormal2()] for the derivative tables of the same map; and
#'   [skewnormal2_distrib()] for the family.
#'
#' @examples
#' # The mean and the standard deviation are not the location and the scale.
#' distributions7:::sn_cp_to_dp(0, 1, 0.5, 1)
#'
#' # At zero skewness the map is the identity on the first two, and the shape
#' # is zero: the Gaussian sits at the same point in both parametrizations.
#' distributions7:::sn_cp_to_dp(3, 2, 0, 1)
#'
#' # Round trip: the direct parameters reproduce the centered moments.
#' dp <- distributions7:::sn_cp_to_dp(0, 1, 0.5, 1)
#' d1 <- skewnormal1_distrib()
#' c(mean = mean(d1, dp), sd = sqrt(variance(d1, dp)), skew = skewness(d1, dp))
#'
#' @keywords internal
sn_cp_to_dp <- function(mu, sigma, gamma1, s) {
  b <- sn_b()
  # |gamma1| as s * gamma1: away from zero the sign is locally constant, so
  # this is exact and carries the right derivatives when the argument is a jet.
  cc <- s * (2 * (s * gamma1) / (4 - pi))^(1 / 3)
  muz <- cc / sqrt(1 + cc^2)
  del <- muz / b
  om <- sigma / sqrt(1 - muz^2)
  list(mu = mu - om * muz,
       sigma = om,
       alpha = del / sqrt(1 - del^2))
}

#' @title Skew Normal Distribution Class, Centered Parametrization
#' @name SkewNormal2Distrib
#'
#' @description
#' The S7 class of the skew normal written in its first three moments: the mean
#' \eqn{\mu}, the standard deviation \eqn{\sigma} and the skewness
#' \eqn{\gamma_1}. It is the same law as [SkewNormal1Distrib], reached through
#' the map of [sn_cp_to_dp()], and the only thing that differs is which three
#' numbers name a member of it.
#'
#' The parametrization is worth having for one property: its expected
#' information is non-singular at \eqn{\gamma_1 = 0}, where the direct
#' parametrization's loses a rank. Measured at \eqn{\mu = 0}, \eqn{\sigma = 1},
#' the information's eigenvalues tend to 2, 1 and \eqn{1/6} as the skewness
#' goes to zero.
#'
#' Build one with [skewnormal2_distrib()], which supplies the three link
#' functions and bounds the skewness at \eqn{\pm 0.9952717}. This page
#' documents the raw S7 constructor, which validates none of the relationships
#' between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `SkewNormal2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [skewnormal2_distrib()] the properties hold `"skew normal2"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma", "gamma1")`, the
#'   interpretations `c(mu = "mean", sigma = "standard deviation", gamma1 =
#'   "skewness")`, `3`, and the domains \eqn{(-\infty,\infty)},
#'   \eqn{(0,\infty)} and \eqn{(-0.9952717, 0.9952717)}.
#'
#' @section Methods:
#' The probability functions delegate to [skewnormal1_distrib()] at the implied
#' direct parameters:
#'   [`distrib_pdf()`][distrib_pdf.SkewNormal2Distrib],
#'   [`distrib_cdf()`][distrib_cdf.SkewNormal2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.SkewNormal2Distrib],
#'   [`distrib_rng()`][distrib_rng.SkewNormal2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.SkewNormal2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewNormal2Distrib].
#'
#' The parameter derivatives carry the parent's through the map by the
#' partition sum of [chain_derivatives()]:
#'   [`distrib_gradient()`][distrib_gradient.SkewNormal2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.SkewNormal2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.SkewNormal2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.SkewNormal2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.SkewNormal2Distrib].
#'
#' Three of the four moments are a parameter read back:
#'   [`mean()`][mean.SkewNormal2Distrib],
#'   [`variance()`][variance.SkewNormal2Distrib],
#'   [`skewness()`][skewness.SkewNormal2Distrib]. The fourth,
#'   [`kurtosis()`][kurtosis.SkewNormal2Distrib], follows from the other three
#'   and is the parent's at the implied direct parameters.
#'
#' @section The point at zero skewness:
#' The parameter derivatives are **rejected** at \eqn{\gamma_1 = 0} exactly.
#' The map runs through a cube root, so \eqn{\partial\alpha/\partial\gamma_1}
#' is unbounded there; the first derivatives of the log-density have a finite
#' limit but the second ones grow like \eqn{\gamma_1^{-2/3}}, so the point is
#' excluded with a message rather than approximated. The **density** and the
#' distribution function are fine there and equal the Gaussian's.
#'
#' @seealso [skewnormal2_distrib()] to build one;
#'   [skewnormal1_distrib()] for the direct parametrization;
#'   [sn_cp_to_dp()] for the map between them.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # The skewness is bounded, which the direct parametrization's shape is not.
#' d@params_bounds$gamma1
#'
#' # All three parameters are moments, which is what "centered" names.
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' c(mean = mean(d, th), sd = sqrt(variance(d, th)), skewness = skewness(d, th))
SkewNormal2Distrib <- S7::new_class("SkewNormal2Distrib",
                                    parent = continuous_distrib)

#' @title The Direct Parameters a Centered Triple Implies
#'
#' @description
#' Takes the sign of the skewness off its plain value and runs
#' [sn_cp_to_dp()], returning the location, scale and shape that
#' [skewnormal1_distrib()] takes. Every probability function of the centered
#' family calls this and then delegates.
#'
#' @param theta A list with `mu`, `sigma` and `gamma1`, in that order, each a
#'   numeric vector. It is read positionally, so it must already be aligned;
#'   the methods that call it have been through [align_theta()].
#'
#' @return A named list with `mu`, `sigma` and `alpha`: the location
#'   \eqn{\xi}, the scale \eqn{\omega} and the shape \eqn{\alpha}. The names
#'   are the parent's, so the result passes straight into
#'   [skewnormal1_distrib()]'s methods.
#'
#' @seealso [sn_cp_to_dp()] for the map itself and [sn2_chain()] for the
#'   derivative route that uses the same map.
#'
#' @examples
#' distributions7:::sn2_theta(list(mu = 0, sigma = 1, gamma1 = 0.5))
#'
#' # A negative skewness reflects the shape and moves the location the other way.
#' distributions7:::sn2_theta(list(mu = 0, sigma = 1, gamma1 = -0.5))
#'
#' @keywords internal
sn2_theta <- function(theta) {
  g <- theta[[3]]
  s <- ifelse(g >= 0, 1, -1)
  sn_cp_to_dp(theta[[1]], theta[[2]], g, s)
}

#' @title Derivatives of the Skew Normal in Its Centered Parametrization
#'
#' @description
#' Carries [skewnormal1_distrib()]'s derivatives into the centered coordinates
#' through the partition sum of [chain_derivatives()], at any order from one to
#' four, observed or expected. The parent supplies the derivatives in
#' \eqn{(\xi, \omega, \alpha)} and [md_skewnormal2()] supplies the map's
#' partial derivatives; the partition sum assembles them.
#'
#' @details
#' # The point that is excluded
#'
#' The map to the direct parametrization runs through
#' \eqn{c = \sqrt[3]{2\gamma_1/(4-\pi)}}, whose derivative grows like
#' \eqn{\gamma_1^{-2/3}}. At zero skewness the map is not differentiable and
#' the chain rule is asked for a quantity that does not exist.
#'
#' The first derivatives of the log-density survive the limit: the map's
#' divergent factor cancels and they approach a finite value from both sides.
#' The second ones do not. Measured at \eqn{y = 0.5}, \eqn{\mu = 0},
#' \eqn{\sigma = 1}, the score in \eqn{\gamma_1} runs
#' \eqn{-0.2152, -0.2257, -0.2284, -0.2290} at
#' \eqn{\gamma_1 = 10^{-2}, 10^{-4}, 10^{-6}, 10^{-8}}, while
#' \eqn{\partial^2\ell/\partial\gamma_1^2} runs
#' \eqn{0.29, 11.5, 253, 5451} over the same values, a factor of 4.642 per
#' decade against \eqn{10^{2/3} = 4.6416}.
#'
#' Zero skewness is therefore rejected here, where the map is used and the
#' reason can be named, with a message that points at
#' [skewnormal1_distrib()], whose derivatives at \eqn{\alpha = 0} are ordinary
#' numbers.
#'
#' # Where the cancellation runs out of digits
#'
#' The **expected** information stays finite as \eqn{\gamma_1 \to 0} and tends
#' to \eqn{1/6} in its own component, but it is computed as a difference of
#' terms of size \eqn{\gamma_1^{-2/3}}. Measured, it holds to seven figures
#' down to \eqn{\gamma_1 = 10^{-8}} (0.16666782 against \eqn{1/6}), loses
#' three by \eqn{10^{-10}} and is **negative at** \eqn{10^{-12}}, which no
#' information can be. That is a limit of double-precision arithmetic on a
#' parameter value no fit visits, and it is why the near-symmetric case is
#' better handled by the direct parametrization.
#'
#' @param distrib A [SkewNormal2Distrib] object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu`, `sigma` and `gamma1`. It is aligned here, so
#'   it may be given in any order and by name.
#' @param order A single integer, 1, 2, 3 or 4: the derivative order.
#' @param expected Logical of length 1. When `TRUE` the parent's **expected**
#'   derivatives are carried instead of the observed ones. Defaults to `FALSE`.
#'
#' @return A named list of numeric vectors, one per distinct component of the
#'   requested order, named in the centered parameters. At order 2 the caller
#'   subsets it by [hess_names()] to fix the ordering.
#'
#' @section Errors:
#' Signals an error when any element of `gamma1` is exactly zero, naming the
#' cube root as the cause and [skewnormal1_distrib()] as the alternative.
#'
#' @seealso [chain_derivatives()] for the partition sum,
#'   [md_skewnormal2()] for the map's derivatives, and
#'   [distrib_gradient.SkewNormal2Distrib()] for the method that calls this.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # Order one, against the method that wraps it.
#' all.equal(distributions7:::sn2_chain(d, y, th, 1L),
#'           distrib_gradient(d, y, th))
#'
#' # The score in the skewness stays of order one as the map's Jacobian
#' # diverges; the curvature does not.
#' t(vapply(10^-c(2, 4, 6, 8), function(g) {
#'   p <- list(mu = 0, sigma = 1, gamma1 = g)
#'   c(gamma1 = g,
#'     score = distrib_gradient(d, 0.5, p)$gamma1,
#'     curvature = distrib_hessian(d, 0.5, p)$gamma1_gamma1)
#' }, numeric(3)))
#'
#' # Zero skewness is rejected rather than approximated.
#' tryCatch(distrib_gradient(d, 0, list(mu = 0, sigma = 1, gamma1 = 0)),
#'          error = function(e) "rejected, as documented")
#'
#' @keywords internal
sn2_chain <- function(distrib, y, theta, order, expected = FALSE) {
  theta <- align_theta(distrib, theta)
  if (any(theta[[3L]] == 0)) {
    stop(paste0(
      "The centered parametrization has no derivatives at zero skewness:\n",
      "  the map to the direct parameters runs through the cube root of\n",
      "  gamma1, whose derivative is unbounded there. The first derivatives\n",
      "  of the log-density have a finite limit and the second ones grow\n",
      "  like gamma1^(-2/3), so the point is excluded rather than\n",
      "  approximated. skewnormal1_distrib() carries the same family in the\n",
      "  direct parametrization, whose derivatives at alpha = 0 are ordinary\n",
      "  numbers."), call. = FALSE)
  }
  chain_derivatives(
    parent = skewnormal1_distrib(),
    y = y,
    th_par = sn2_theta(theta),
    maps = md_skewnormal2(theta[1:3]),
    new_params = distrib@params,
    order = order,
    expected = expected
  )
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Skew Normal Density in the Centered Parametrization
#' @name distrib_pdf.SkewNormal2Distrib
#'
#' @description
#' Computes the skew normal density at the direct parameters that the centered
#' triple implies. With
#' \eqn{(\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)} from
#' [sn_cp_to_dp()] and \eqn{z = (y-\xi)/\omega},
#' \deqn{f(y; \mu, \sigma, \gamma_1)
#'       = \dfrac{2}{\omega}\,\phi(z)\,\Phi(\alpha z).}
#' The density is the same function of \eqn{y} as
#' [distrib_pdf.SkewNormal1Distrib()]'s; only the three numbers naming it
#' differ.
#'
#' Unlike the derivatives, the density is defined at \eqn{\gamma_1 = 0} and
#' equals the Gaussian's there: the map itself is continuous through zero, and
#' only its derivative is not.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations, anywhere on the real line.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `y`. `sigma` must be
#'   strictly positive and `gamma1` must lie in
#'   \eqn{(-0.9952717, 0.9952717)}; a skewness outside that range belongs to no
#'   skew normal and the map returns `NaN`.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma), length(gamma1))`.
#'
#' @section Notation:
#' \eqn{\mu}, \eqn{\sigma} and \eqn{\gamma_1} are the mean, the standard
#' deviation and the skewness; \eqn{\xi}, \eqn{\omega} and \eqn{\alpha} the
#' location, scale and shape they imply; \eqn{\phi} and \eqn{\Phi} the standard
#' Gaussian density and distribution function.
#'
#' @seealso [distrib_pdf.SkewNormal1Distrib()] for the same density in the
#'   direct parametrization, [sn_cp_to_dp()] for the map, and [distrib_pdf()]
#'   for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' d1 <- skewnormal1_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # The same law, reached through the map.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(d1, y, distributions7:::sn2_theta(th)))
#'
#' # It integrates to one, and its first moment is the parameter mu.
#' c(mass = integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value,
#'   mean = integrate(function(v) v * distrib_pdf(d, v, th), -Inf, Inf)$value)
#'
#' # At zero skewness the density is the Gaussian's, where the derivatives
#' # are not defined.
#' all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, gamma1 = 0)), dnorm(y))
S7::method(distrib_pdf, SkewNormal2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(skewnormal1_distrib(), y, sn2_theta(theta), log = log)
}

#' @title Skew Normal Distribution Function in the Centered Parametrization
#' @name distrib_cdf.SkewNormal2Distrib
#'
#' @description
#' Computes the skew normal distribution function at the direct parameters the
#' centered triple implies, through Azzalini's Owen's T identity
#' \eqn{F(q) = \Phi(z) - 2T(z, \alpha)} with \eqn{z = (q-\xi)/\omega}. The
#' arithmetic is [distrib_cdf.SkewNormal1Distrib()]'s; this method supplies the
#' translated parameters.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `q`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Passed to [distrib_cdf.SkewNormal1Distrib()].
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms with `log.p = TRUE`.
#'
#' @section Notation:
#' \eqn{\gamma_1} is the skewness, \eqn{(\xi, \omega, \alpha)} the implied
#' location, scale and shape, and \eqn{T} Owen's T function.
#'
#' @seealso [distrib_cdf.SkewNormal1Distrib()] for the identity,
#'   [distrib_quantile.SkewNormal2Distrib()] for its inverse, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' q <- c(-2, -0.5, 0.5, 2)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # Against a direct quadrature of the density.
#' rbind(owen = distrib_cdf(d, q, th),
#'       quadrature = vapply(q, function(u)
#'         integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#'
#' # At zero skewness it is the Gaussian's.
#' all.equal(distrib_cdf(d, q, list(mu = 0, sigma = 1, gamma1 = 0)), pnorm(q))
#'
#' # A positive skewness puts more than half the mass below the mean.
#' distrib_cdf(d, 0, th)
S7::method(distrib_cdf, SkewNormal2Distrib) <- function(distrib, q, theta,
                                                         lower.tail = TRUE,
                                                         log.p = FALSE, ...) {
  distrib_cdf(skewnormal1_distrib(), q, sn2_theta(theta),
              lower.tail = lower.tail, log.p = log.p, ...)
}

#' @title Skew Normal Quantile Function in the Centered Parametrization
#' @name distrib_quantile.SkewNormal2Distrib
#'
#' @description
#' Computes the quantiles of the skew normal at the direct parameters the
#' centered triple implies. The skew normal has no closed-form quantile
#' function, so the value comes from [continuous_distrib()]'s root finding on
#' the distribution function, reached through
#' [distrib_quantile.continuous_distrib()].
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms when `log.p = TRUE`.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `p`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Passed to [distrib_quantile.continuous_distrib()], including the
#'   root finder's tolerance.
#'
#' @return A numeric vector of quantiles, of the length of the recycled inputs.
#'
#' @seealso [distrib_cdf.SkewNormal2Distrib()], which it inverts,
#'   [distrib_rng.SkewNormal2Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # The round trip through the distribution function.
#' p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' q <- distrib_quantile(d, p, th)
#' rbind(quantile = q, back = distrib_cdf(d, q, th))
#'
#' # A right-skewed density has its median below its mean.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, SkewNormal2Distrib) <- function(distrib, p, theta,
                                                              lower.tail = TRUE,
                                                              log.p = FALSE, ...) {
  distrib_quantile(skewnormal1_distrib(), p, sn2_theta(theta),
                   lower.tail = lower.tail, log.p = log.p, ...)
}

#' @title Skew Normal Random Generation in the Centered Parametrization
#' @name distrib_rng.SkewNormal2Distrib
#'
#' @description
#' Draws from the skew normal at the direct parameters the centered triple
#' implies, through [distrib_rng.SkewNormal1Distrib()]'s stochastic
#' representation. The draws are exact and cost two `rnorm` calls, whatever the
#' skewness is.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of length `n`.
#'
#' @return A numeric vector of `n` draws, whose first three sample moments
#'   estimate `mu`, `sigma^2` and `gamma1`.
#'
#' @seealso [distrib_rng.SkewNormal1Distrib()] for the representation,
#'   [mean.SkewNormal2Distrib()] for the moments the draws reproduce, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 3, sigma = 2, gamma1 = 0.6)
#'
#' set.seed(2)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # The three parameters are the three sample moments, which is the whole
#' # point of this parametrization.
#' rbind(sample = c(mean(x), sd(x), mean((x - mean(x))^3) / sd(x)^3),
#'       parameter = c(th$mu, th$sigma, th$gamma1))
S7::method(distrib_rng, SkewNormal2Distrib) <- function(distrib, n, theta) {
  distrib_rng(skewnormal1_distrib(), n, sn2_theta(theta))
}

#' @title Skew Normal Score in the Centered Parametrization
#' @name distrib_gradient.SkewNormal2Distrib
#'
#' @description
#' Computes the three first derivatives of the log-density in the centered
#' parameters, by carrying [distrib_gradient.SkewNormal1Distrib()]'s through
#' the Jacobian of [sn_cp_to_dp()]:
#' \deqn{\dfrac{\partial \ell}{\partial \psi_j}
#'       = \sum_{k} \dfrac{\partial \ell}{\partial \theta_k}
#'                  \dfrac{\partial \theta_k}{\partial \psi_j},
#'       \qquad \psi = (\mu, \sigma, \gamma_1),\;
#'              \theta = (\xi, \omega, \alpha).}
#'
#' The component in \eqn{\gamma_1} stays of order one however small
#' \eqn{\gamma_1} is, although the Jacobian itself grows without bound.
#' Measured at \eqn{y = 0.5}, \eqn{\mu = 0}, \eqn{\sigma = 1}, the score reads
#' \eqn{-0.2152, -0.2257, -0.2284, -0.2290} at
#' \eqn{\gamma_1 = 10^{-2}, 10^{-4}, 10^{-6}, 10^{-8}} while
#' \eqn{\partial\alpha/\partial\gamma_1} reads \eqn{5.1, 258} at the first two.
#' The divergent parts cancel, and that cancellation is the reason the
#' parametrization exists.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`. The
#'   skewness must not be exactly zero; see the error below.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body, so this method always
#'   returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu`, `sigma` and `gamma1`,
#'   each of the length of the recycled inputs.
#'
#' @section Errors:
#' Signals an error when any element of `gamma1` is exactly zero: the map runs
#' through a cube root and is not differentiable there. The density is defined
#' at that point, and [skewnormal1_distrib()] carries the same family with
#' ordinary derivatives at symmetry.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\psi} the centered
#' parameters and \eqn{\theta} the direct ones.
#'
#' @seealso [sn2_chain()] for the partition sum this calls,
#'   [distrib_hessian.SkewNormal2Distrib()] for the next order, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # Against numerical differentiation of the log-density itself.
#' f <- function(p) sum(distrib_pdf(d, y, as.list(setNames(p, names(th))),
#'                                  log = TRUE))
#' rbind(analytic = vapply(g, sum, 0),
#'       numeric = numDeriv::grad(f, unlist(th)))
#'
#' # The score in the skewness stays bounded as the map's Jacobian diverges.
#' vapply(10^-c(2, 4, 6, 8),
#'        function(v) distrib_gradient(d, 0.5,
#'                      list(mu = 0, sigma = 1, gamma1 = v))$gamma1, 0)
S7::method(distrib_gradient, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"), ...) {
  sn2_chain(distrib, y, theta, 1L)
}

#' @title Skew Normal Observed Hessian in the Centered Parametrization
#' @name distrib_hessian.SkewNormal2Distrib
#'
#' @description
#' Computes the six second derivatives of the log-density in the centered
#' parameters, by the second-order chain rule through [sn_cp_to_dp()]:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \psi_i \partial \psi_j}
#'       = \sum_{k,l} \dfrac{\partial^2 \ell}{\partial\theta_k\partial\theta_l}
#'         \dfrac{\partial\theta_k}{\partial\psi_i}
#'         \dfrac{\partial\theta_l}{\partial\psi_j}
#'       + \sum_{k} \dfrac{\partial \ell}{\partial\theta_k}
#'         \dfrac{\partial^2\theta_k}{\partial\psi_i\partial\psi_j}.}
#' Both terms come from [chain_derivatives()], with the map's partial
#' derivatives supplied by [md_skewnormal2()] as a written-out table.
#'
#' The **observed** curvature in \eqn{\gamma_1} diverges as the skewness goes
#' to zero, at the rate \eqn{\gamma_1^{-2/3}} the cube root sets. Measured at
#' \eqn{y = 0.5}, \eqn{\mu = 0}, \eqn{\sigma = 1}, it is 0.29, 11.5, 253 and
#' 5451 at \eqn{\gamma_1 = 10^{-2}, 10^{-4}, 10^{-6}, 10^{-8}}, a factor of
#' 4.642 per decade against \eqn{10^{2/3} = 4.6416}. The **expected**
#' curvature does not: see [distrib_expected_hessian.SkewNormal2Distrib()].
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`. The
#'   skewness must not be exactly zero.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors, in [hess_names()]'s order:
#'   `mu_mu`, `sigma_sigma`, `gamma1_gamma1`, `mu_sigma`, `mu_gamma1`,
#'   `sigma_gamma1`.
#'
#' @section Errors:
#' Signals an error when any element of `gamma1` is exactly zero.
#'
#' @seealso [distrib_gradient.SkewNormal2Distrib()] for the order below,
#'   [distrib_expected_hessian.SkewNormal2Distrib()] for the expectation,
#'   [sn2_chain()] for the partition sum, and [distrib_hessian()] for the
#'   generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # Against numerical differentiation of the log-density.
#' f <- function(p) sum(distrib_pdf(d, y, as.list(setNames(p, names(th))),
#'                                  log = TRUE))
#' H <- numDeriv::hessian(f, unlist(th))
#' rbind(analytic = c(sum(h$mu_gamma1), sum(h$gamma1_gamma1)),
#'       numeric = c(H[1, 3], H[3, 3]))
#'
#' # The observed curvature in the skewness diverges as gamma1^(-2/3). The
#' # ratio over two decades converges to 10^(2/3) = 4.6416.
#' cv <- vapply(10^-c(2, 4, 6, 8),
#'              function(v) distrib_hessian(d, 0.5,
#'                            list(mu = 0, sigma = 1, gamma1 = v))$gamma1_gamma1, 0)
#' rbind(curvature = cv, per_decade = c(NA, (cv[-1] / cv[-4])^(1 / 2)))
S7::method(distrib_hessian, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  sn2_chain(distrib, y, theta, 2L)[hess_names(distrib@params)]
}

#' @title Skew Normal Expected Information in the Centered Parametrization
#' @name distrib_expected_hessian.SkewNormal2Distrib
#'
#' @description
#' Computes the expected second derivatives by carrying the parent's expected
#' information through the same congruence the observed Hessian uses,
#' \eqn{J^\top E[\ell''] J} with \eqn{J} the Jacobian of [sn_cp_to_dp()]. The
#' first-order term of the chain rule drops out under expectation, the score
#' having mean zero.
#'
#' The matrix is **non-singular at zero skewness**, which the direct
#' parametrization's is not: there the score for \eqn{\alpha} is exactly
#' proportional to the score for the location and the information loses a rank.
#' Measured at \eqn{\mu = 0}, \eqn{\sigma = 1}, the eigenvalues here tend to
#' 2, 1 and \eqn{1/6} as \eqn{\gamma_1 \to 0}. Removing that singularity is
#' what the centered parametrization is for.
#'
#' @details
#' # Cost, and where the digits run out
#'
#' The parent's own expected information is the base class's quadrature, so
#' this method is a chain on top of a numerical quantity: measured at 100
#' observations it costs about 5.2 seconds against the parent's 2.2, where a
#' family that writes its information out answers in a median of 0.18
#' milliseconds. [expected_hessian_exact()] therefore returns `FALSE` here, and
#' `approx` is read.
#'
#' The congruence is a difference of terms of size \eqn{\gamma_1^{-2/3}}, so
#' the limit is approached and then lost. Measured, the \eqn{\gamma_1}
#' component is 0.16666782 against \eqn{1/6 = 0.16666667} at
#' \eqn{\gamma_1 = 10^{-8}}, 0.1655 at \eqn{10^{-10}}, and **negative** at
#' \eqn{10^{-12}}. A fit does not visit those values, and a genuinely
#' symmetric problem is better posed in [skewnormal1_distrib()].
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector. Its values do not enter the result, which is an
#'   expectation; only its length does, through recycling.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`. The
#'   skewness must not be exactly zero.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"bartlett"`, `"integrate"`, `"mc"` or `"opg"`, the
#'   strategy the parent uses for its own expectation. Defaults to
#'   `"bartlett"`, the variance of the score.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors, in [hess_names()]'s order.
#'   Every entry is an expectation, so it does not depend on `y`.
#'
#' @section Errors:
#' Signals an error when any element of `gamma1` is exactly zero.
#'
#' @seealso [distrib_hessian.SkewNormal2Distrib()] for the observed curvature,
#'   [expected_hessian_exact.SkewNormal2Distrib()] for why this counts as
#'   approximated, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' e <- distrib_expected_hessian(d, 0, th)
#' names(e)
#'
#' # The information is positive definite, and stays so into symmetry, where
#' # the direct parametrization loses a rank.
#' info <- function(g) {
#'   e <- distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, gamma1 = g))
#'   M <- matrix(c(e$mu_mu, e$mu_sigma, e$mu_gamma1,
#'                 e$mu_sigma, e$sigma_sigma, e$sigma_gamma1,
#'                 e$mu_gamma1, e$sigma_gamma1, e$gamma1_gamma1), 3, 3)
#'   eigen(-M, only.values = TRUE)$values
#' }
#' rbind(gamma1_0.5 = info(0.5), gamma1_1e_6 = info(1e-6))
#'
#' # Its own component tends to 1/6.
#' c(limit = 1 / 6,
#'   at_1e_6 = -distrib_expected_hessian(d, 0,
#'               list(mu = 0, sigma = 1, gamma1 = 1e-6))$gamma1_gamma1)
S7::method(distrib_expected_hessian, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                                      scale = c("parameter", "link"),
                                                                      approx = c("bartlett", "integrate", "mc", "opg"),
                                                                      nsim = 10000, ...) {
  sn2_chain(distrib, y, theta, 2L, expected = TRUE)[hess_names(distrib@params)]
}

#' @title The Centered Skew Normal Does Not Write Its Expected Information Out
#' @name expected_hessian_exact.SkewNormal2Distrib
#'
#' @description
#' Returns `FALSE`, by asking [skewnormal1_distrib()] the same question. The
#' registration of [distrib_expected_hessian.SkewNormal2Distrib()] says where
#' the arithmetic is assembled, not that it is closed form: it is a chain onto
#' the parent, whose own expected information is the base class's quadrature.
#'
#' @details
#' Reading the owning class would answer `TRUE` here and be wrong, which is
#' this predicate's whole reason for existing as a generic. The cost separates
#' the two cases by four orders of magnitude: measured at 100 observations,
#' this family costs 5220 milliseconds and the parent it chains onto 2230,
#' where the families that do write their information out answer in a median of
#' 0.183 milliseconds.
#'
#' Reported as exact, it made [fit_distrib()] reject a legitimate
#' `fisher_scoring(approx = )` here with a message saying the family computes
#' its expected information in closed form, which is untrue.
#'
#' @param x A `SkewNormal2Distrib` object.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `FALSE`, a logical of length 1.
#'
#' @seealso [expected_hessian_exact()] for the generic and the rule it
#'   encodes, and [distrib_expected_hessian.SkewNormal2Distrib()] for the
#'   method it describes.
#'
#' @examples
#' # Neither parametrization of the skew normal writes its information out.
#' c(centered = expected_hessian_exact(skewnormal2_distrib()),
#'   direct = expected_hessian_exact(skewnormal1_distrib()),
#'   gaussian = expected_hessian_exact(gaussian1_distrib()))
#'
#' @keywords internal
S7::method(expected_hessian_exact, SkewNormal2Distrib) <- function(x, ...) {
  expected_hessian_exact(skewnormal1_distrib())
}

#' @title Skew Normal Third Derivatives in the Centered Parametrization
#' @name distrib_deriv3.SkewNormal2Distrib
#'
#' @description
#' Computes the ten third derivatives of the log-density in the centered
#' parameters, by the third-order partition sum of [chain_derivatives()] over
#' the map of [sn_cp_to_dp()]. The parent's third derivatives are closed form
#' in a compiled kernel and the map's partial derivatives are a written-out
#' table in [md_skewnormal2()], so nothing here is a finite difference.
#'
#' With `expected = TRUE` the parent's expected derivatives are carried
#' instead. Those are numerical, so `approx` and `nsim` are read.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`. The
#'   skewness must not be exactly zero.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, the
#'   strategy the parent uses when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of ten numeric vectors, one per distinct third-order
#'   component in the centered parameters, from `mu_mu_mu` to
#'   `gamma1_gamma1_gamma1` as [deriv_names()] names them.
#'
#' @section Errors:
#' Signals an error when any element of `gamma1` is exactly zero.
#'
#' @seealso [distrib_hessian.SkewNormal2Distrib()] for the order below,
#'   [distrib_deriv4.SkewNormal2Distrib()] for the order above,
#'   [chain_derivatives()] for the partition sum, and [distrib_deriv3()] for the
#'   generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Against a central difference of the analytic Hessian.
#' eps <- 1e-5
#' rbind(analytic = d3$mu_mu_gamma1,
#'       numeric = (distrib_hessian(d, y, list(mu = 0, sigma = 1,
#'                                             gamma1 = 0.5 + eps))$mu_mu -
#'                  distrib_hessian(d, y, list(mu = 0, sigma = 1,
#'                                             gamma1 = 0.5 - eps))$mu_mu) /
#'                 (2 * eps))
S7::method(distrib_deriv3, SkewNormal2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  sn2_chain(distrib, y, theta, 3L, expected = expected)
}

#' @title Skew Normal Fourth Derivatives in the Centered Parametrization
#' @name distrib_deriv4.SkewNormal2Distrib
#'
#' @description
#' Computes the fifteen fourth derivatives of the log-density in the centered
#' parameters, by the fourth-order partition sum of [chain_derivatives()] over
#' the map of [sn_cp_to_dp()]. The map's fourth partial derivatives are the
#' last entries of [md_skewnormal2()]'s table, so the order costs one more term
#' of the same enumeration.
#'
#' With `expected = TRUE` the parent's expected derivatives are carried
#' instead, and those are numerical.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`. The
#'   skewness must not be exactly zero.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of fifteen numeric vectors, one per distinct
#'   fourth-order component in the centered parameters.
#'
#' @section Errors:
#' Signals an error when any element of `gamma1` is exactly zero.
#'
#' @seealso [distrib_deriv3.SkewNormal2Distrib()] for the order below,
#'   [chain_derivatives()] for the partition sum, and [distrib_deriv4()] for
#'   the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' length(distrib_deriv4(d, y, th))
#'
#' # Against a central difference of the third order.
#' eps <- 1e-5
#' rbind(analytic = distrib_deriv4(d, y, th)$mu_mu_gamma1_gamma1,
#'       numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1,
#'                                            gamma1 = 0.5 + eps))$mu_mu_gamma1 -
#'                  distrib_deriv3(d, y, list(mu = 0, sigma = 1,
#'                                            gamma1 = 0.5 - eps))$mu_mu_gamma1) /
#'                 (2 * eps))
S7::method(distrib_deriv4, SkewNormal2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  sn2_chain(distrib, y, theta, 4L, expected = expected)
}

#' @title Skew Normal Response Derivative in the Centered Parametrization
#' @name distrib_grad_y.SkewNormal2Distrib
#'
#' @description
#' Computes \eqn{\partial\ell/\partial y} by delegating to
#' [distrib_grad_y.SkewNormal1Distrib()] at the implied direct parameters. The
#' value is the parent's unchanged: a change of parameters does not touch a
#' derivative in the response, the two variables being separate arguments of
#' the same log-density.
#'
#' It is therefore defined at \eqn{\gamma_1 = 0}, where the parameter
#' derivatives are not: nothing here differentiates the map.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param ... Passed to [distrib_grad_y.SkewNormal1Distrib()].
#'
#' @return A numeric vector of the length of the recycled inputs.
#'
#' @seealso [distrib_hess_y.SkewNormal2Distrib()] for the second derivative,
#'   [distrib_grad_y.SkewNormal1Distrib()] for the closed form it delegates to,
#'   and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # Against a central difference of the log-density in the response.
#' eps <- 1e-6
#' rbind(analytic = distrib_grad_y(d, y, th),
#'       numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
#'                  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#'
#' # Defined at zero skewness, where the parameter derivatives are not.
#' all.equal(distrib_grad_y(d, y, list(mu = 0, sigma = 1, gamma1 = 0)), -y)
S7::method(distrib_grad_y, SkewNormal2Distrib) <- function(distrib, y, theta, ...) {
  distrib_grad_y(skewnormal1_distrib(), y, sn2_theta(theta), ...)
}

#' @title Skew Normal Second Response Derivative in the Centered Parametrization
#' @name distrib_hess_y.SkewNormal2Distrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y^2} by delegating to
#' [distrib_hess_y.SkewNormal1Distrib()] at the implied direct parameters. Like
#' the first response derivative it is the parent's unchanged, and is defined
#' at \eqn{\gamma_1 = 0}.
#'
#' The value is strictly negative at every observation and every skewness: the
#' skew normal log-density is concave in the response.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param ... Passed to [distrib_hess_y.SkewNormal1Distrib()].
#'
#' @return A numeric vector of the length of the recycled inputs, negative
#'   throughout.
#'
#' @seealso [distrib_grad_y.SkewNormal2Distrib()] for the first derivative,
#'   [distrib_hess_y.SkewNormal1Distrib()] for the closed form, and
#'   [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # Against a central difference of the response derivative.
#' eps <- 1e-5
#' rbind(analytic = distrib_hess_y(d, y, th),
#'       numeric = (distrib_grad_y(d, y + eps, th) -
#'                  distrib_grad_y(d, y - eps, th)) / (2 * eps))
#'
#' # Concave in the response at every skewness the family reaches.
#' vapply(c(-0.9, -0.3, 0.3, 0.9), function(g)
#'   max(distrib_hess_y(d, seq(-6, 6, by = 0.5),
#'                      list(mu = 0, sigma = 1, gamma1 = g))), 0)
S7::method(distrib_hess_y, SkewNormal2Distrib) <- function(distrib, y, theta, ...) {
  distrib_hess_y(skewnormal1_distrib(), y, sn2_theta(theta), ...)
}

#' @title Mean of the Skew Normal in the Centered Parametrization
#' @name mean.SkewNormal2Distrib
#'
#' @description
#' Returns \eqn{\mu}, the first parameter, which in this parametrization is the
#' mean by construction. The addition of [moment_const()] recycles the value to
#' the length the three parameters imply, so a `theta` whose components vary by
#' observation gets one mean per observation.
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of the length the recycled parameters
#'   imply.
#'
#' @seealso [variance.SkewNormal2Distrib()] and
#'   [skewness.SkewNormal2Distrib()], the other two parameters read back;
#'   [kurtosis.SkewNormal2Distrib()], which is not a parameter; and
#'   [mean.SkewNormal1Distrib()] for the same quantity in the direct
#'   parametrization, where it is not.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' mean(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#'
#' # One value per observation when a parameter varies.
#' mean(d, list(mu = c(0, 3, 7), sigma = 2, gamma1 = 0.6))
#'
#' @keywords internal
S7::method(mean, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 3L, 0) + theta[[1]]
}

#' @title Variance of the Skew Normal in the Centered Parametrization
#' @name variance.SkewNormal2Distrib
#'
#' @description
#' Returns \eqn{\sigma^2}, the square of the second parameter, which in this
#' parametrization is the variance by construction. [moment_const()] recycles
#' the result to the length the three parameters imply.
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of the length the recycled parameters
#'   imply.
#'
#' @seealso [mean.SkewNormal2Distrib()] and [skewness.SkewNormal2Distrib()],
#'   and [variance.SkewNormal1Distrib()], where the scale is not the standard
#'   deviation.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' variance(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#'
#' # It agrees with a quadrature of the second central moment.
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' c(parameter = variance(d, th),
#'   quadrature = integrate(function(v) v^2 * distrib_pdf(d, v, th),
#'                          -Inf, Inf)$value)
#'
#' @keywords internal
S7::method(variance, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 3L, 0)
}

#' @title Skewness of the Skew Normal in the Centered Parametrization
#' @name skewness.SkewNormal2Distrib
#'
#' @description
#' Returns \eqn{\gamma_1}, the third parameter, which in this parametrization
#' is the standardized third central moment by construction. [moment_const()]
#' recycles the result to the length the three parameters imply.
#'
#' The value cannot leave \eqn{(-0.9952717, 0.9952717)}, the constructor having
#' bounded the parameter at the supremum the family reaches; see
#' [sn_max_skew()].
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of skewnesses, of the length the recycled
#'   parameters imply.
#'
#' @seealso [sn_max_skew()] for the bound,
#'   [skewness.SkewNormal1Distrib()] for the same quantity computed from a
#'   shape, and [kurtosis.SkewNormal2Distrib()], which is not free.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' skewness(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#'
#' # It agrees with a quadrature of the standardized third central moment.
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' c(parameter = skewness(d, th),
#'   quadrature = integrate(function(v) v^3 * distrib_pdf(d, v, th),
#'                          -Inf, Inf)$value)
#'
#' @keywords internal
S7::method(skewness, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 3L, 0) + theta[[3]]
}

#' @title Kurtosis of the Skew Normal in the Centered Parametrization
#' @name kurtosis.SkewNormal2Distrib
#'
#' @description
#' Returns the excess kurtosis, computed from
#' [kurtosis.SkewNormal1Distrib()] at the implied direct parameters. It is the
#' one moment this parametrization does not name: fixing the first three uses
#' up all three parameters, so the fourth follows from them.
#'
#' The consequence for use is that a skew normal cannot match an arbitrary
#' first four moments. At \eqn{\gamma_1 = 0.5} the excess kurtosis is 0.347,
#' and it is determined by the skewness alone.
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of the length the recycled
#'   parameters imply.
#'
#' @seealso [skewness.SkewNormal2Distrib()], which does name a parameter;
#'   [kurtosis.SkewNormal1Distrib()] for the closed form; and
#'   [skewt_distrib()], whose degrees of freedom give the fourth moment a
#'   parameter of its own.
#'
#' @examples
#' d <- skewnormal2_distrib()
#'
#' # The kurtosis follows from the skewness and from nothing else: changing
#' # the mean and the standard deviation leaves it alone.
#' c(a = kurtosis(d, list(mu = 0, sigma = 1, gamma1 = 0.5)),
#'   b = kurtosis(d, list(mu = 9, sigma = 4, gamma1 = 0.5)))
#'
#' # It rises with the skewness, and is zero at symmetry.
#' vapply(c(0.001, 0.3, 0.6, 0.9),
#'        function(g) kurtosis(d, list(mu = 0, sigma = 1, gamma1 = g)), 0)
#'
#' @keywords internal
S7::method(kurtosis, SkewNormal2Distrib) <- function(x, theta, ...) {
  kurtosis(skewnormal1_distrib(), sn2_theta(align_theta(x, theta)))
}


#' @title Skew Normal Distribution Object, Centered Parametrization
#'
#' @description
#' Builds a skew normal distribution object parametrized by its mean
#' \eqn{\mu}, its standard deviation \eqn{\sigma} and its skewness
#' \eqn{\gamma_1}. It is the same family as [skewnormal1_distrib()], named by
#' three moments instead of by a location, a scale and a shape.
#'
#' The parametrization is Azzalini's, and the property it exists for is that
#' its expected information stays non-singular at \eqn{\gamma_1 = 0}, where the
#' direct parametrization's loses a rank. The price is that the map between the
#' two runs through a cube root, so no parameter derivative exists at exactly
#' zero skewness.
#'
#' @param link_mu A `linkfunctions7` link object for the mean, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_sigma A link object for the standard deviation, which must be
#'   strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_gamma1 A link object for the skewness, which must lie strictly
#'   inside \eqn{(-0.9952717, 0.9952717)}. Defaults to
#'   [linkfunctions7::bounded_link()] over exactly that interval, so any real
#'   linear predictor maps to a skewness the family can reach.
#'
#' @details
#' # This is a family, not a reparametrize()
#'
#' The map passes through
#' \eqn{c = \mathrm{sign}(\gamma_1)(2|\gamma_1|/(4-\pi))^{1/3}}, and two things
#' follow from it. It carries a sign, so [sn2_theta()] reads that off the plain
#' value and hands it to [sn_cp_to_dp()] as an argument, leaving a body with no
#' `abs()` in it to differentiate. And its derivatives are written out by hand
#' in [md_skewnormal2()], as a keyed table of the map's partials, which
#' [chain_derivatives()] consumes. The toolkit assembles higher-order
#' derivatives from written-out tables throughout.
#'
#' # What the map costs, and what it buys
#'
#' \eqn{\partial\alpha/\partial\gamma_1} grows without bound as
#' \eqn{\gamma_1 \to 0}: measured, 3.9 at \eqn{\gamma_1 = 0.5}, 12.8 at 0.01
#' and 258 at \eqn{10^{-4}}. The score does not follow it: the divergent
#' contributions cancel, so the information in \eqn{\gamma_1} tends to
#' \eqn{1/6} and the whole matrix stays positive definite at symmetry.
#'
#' The **observed** curvature does diverge, at the rate the cube root sets:
#' \eqn{\gamma_1^{-2/3}}, measured at 4.642 per decade against
#' \eqn{10^{2/3} = 4.6416}. The parameter derivatives are therefore rejected at
#' \eqn{\gamma_1 = 0} exactly, with a message naming the cause. The density,
#' the distribution function, the quantile function, the generator and both
#' response derivatives are fine there and equal the Gaussian's.
#'
#' The cancellation is between terms of size \eqn{\gamma_1^{-2/3}}, so it
#' eventually runs out of digits. Measured, the expected information's
#' \eqn{\gamma_1} component holds seven figures at \eqn{\gamma_1 = 10^{-8}},
#' loses three by \eqn{10^{-10}} and is negative at \eqn{10^{-12}}. Those are
#' values no fit visits, and a genuinely symmetric problem is better posed in
#' [skewnormal1_distrib()].
#'
#' # The bound on the skewness
#'
#' A skew normal cannot reach \eqn{|\gamma_1| > 0.9952717} whatever its shape,
#' so `gamma1` is bounded there and carries a bounded link by default. That
#' ceiling is why [skewt_distrib()] exists.
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (-\infty, \infty)}
#' - \eqn{\sigma \in (0, \infty)}
#' - \eqn{\gamma_1 \in (-0.9952717, 0.9952717)}
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{2}{\omega}\,\phi\!\left(\frac{y-\xi}{\omega}\right)\Phi\!\left(\alpha\,\frac{y-\xi}{\omega}\right), \qquad (\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)}
#' on \eqn{y \in \mathbb{R}}, with \eqn{\mathrm{DP}} the map of
#' [sn_cp_to_dp()].
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2},
#'       \qquad \gamma_1(Y) = \gamma_1.}
#'
#' @return An S7 object of class [SkewNormal2Distrib], inheriting from
#'   `continuous_distrib`. Its `params` are `c("mu", "sigma", "gamma1")`, its
#'   `bounds` `c(-Inf, Inf)`, and its `link_params` the three links given here.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
#' Families*. Cambridge University Press. The centered parametrization is
#' section 3.1.4.
#'
#' @seealso [skewnormal1_distrib()] for the direct parametrization,
#'   [sn_cp_to_dp()] for the map between them,
#'   [skewt_distrib()] for a family that reaches a larger skewness, and
#'   [SkewNormal2Distrib] for the class and its method list.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # All three parameters are moments, which is what "centered" names.
#' c(mean = mean(d, th), sd = sqrt(variance(d, th)),
#'   skewness = skewness(d, th))
#'
#' # The fourth moment is not free: it follows from the skewness.
#' vapply(c(0.001, 0.5, 0.9),
#'        function(g) kurtosis(d, list(mu = 0, sigma = 1, gamma1 = g)), 0)
#'
#' # The information stays invertible at symmetry, where the direct
#' # parametrization's does not.
#' smallest_eigenvalue <- function(dd, p) {
#'   e <- distrib_expected_hessian(dd, 0, p)
#'   M <- matrix(c(e[[1]], e[[4]], e[[5]],
#'                 e[[4]], e[[2]], e[[6]],
#'                 e[[5]], e[[6]], e[[3]]), 3, 3)   # hess_names() order
#'   min(eigen(-M, only.values = TRUE)$values)
#' }
#' c(centered = smallest_eigenvalue(d, list(mu = 0, sigma = 1, gamma1 = 1e-6)),
#'   direct = smallest_eigenvalue(skewnormal1_distrib(),
#'                                list(mu = 0, sigma = 1, alpha = 0)))
#'
#' # A fit recovers all three moments.
#' set.seed(11)
#' x <- distrib_rng(d, 4000, list(mu = 3, sigma = 2, gamma1 = 0.6))
#' coef(fit_distrib(d, x))
#'
#' @export
skewnormal2_distrib <- function(link_mu = identity_link(),
                                link_sigma = log_link(),
                                link_gamma1 = bounded_link(
                                  lwr = -sn_max_skew(), upr = sn_max_skew()
                                )) {
  g <- sn_max_skew()
  SkewNormal2Distrib(
    distrib_name = "skew normal2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma", "gamma1"),
    params_interpretation = c(mu = "mean", sigma = "standard deviation",
                              gamma1 = "skewness"),
    n_params = 3,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf),
                         gamma1 = c(-g, g)),
    link_params = list(mu = link_mu, sigma = link_sigma,
                       gamma1 = link_gamma1)
  )
}
