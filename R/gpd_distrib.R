#' @include distrib.R generics.R
NULL

#' @title Generalized Pareto Distribution Class
#' @name GPDDistrib
#'
#' @description
#' The S7 class of the generalized Pareto family, the law of exceedances over a
#' high threshold, parametrized by a scale \eqn{\sigma > 0} and a shape
#' \eqn{\xi}:
#' \eqn{f(y) = \sigma^{-1}(1 + \xi y/\sigma)^{-1/\xi-1}} on \eqn{y \ge 0}. At
#' \eqn{\xi = 0} it is the exponential, reached by a series, so a fit may move
#' the shape through zero.
#'
#' It is the first family in the package whose **support moves with its
#' parameters**: for \eqn{\xi < 0} it ends at \eqn{-\sigma/\xi}. That is what
#' makes it non-regular, and it is why its expected information exists only
#' above \eqn{\xi = -1/2}.
#'
#' Build one with [gpd_distrib()], which supplies the two link functions. This
#' page documents the raw S7 constructor, which validates none of the
#' relationships between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `GPDDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [gpd_distrib()] the properties hold `"generalized pareto"`,
#'   `"univariate"`, `c(0, Inf)`, `c("sigma", "xi")`, the interpretations
#'   `c(sigma = "scale", xi = "shape")`, `2`, and the domains
#'   \eqn{(0,\infty)} and \eqn{(-\infty,\infty)}.
#'
#' @section Methods:
#' Registered in this file, all compiled apart from the first four:
#'   [`distrib_pdf()`][distrib_pdf.GPDDistrib],
#'   [`distrib_cdf()`][distrib_cdf.GPDDistrib],
#'   [`distrib_quantile()`][distrib_quantile.GPDDistrib],
#'   [`distrib_rng()`][distrib_rng.GPDDistrib],
#'   [`distrib_gradient()`][distrib_gradient.GPDDistrib],
#'   [`distrib_hessian()`][distrib_hessian.GPDDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GPDDistrib].
#'
#' Registered elsewhere: the third and fourth orders in `gpd_higher.R`
#' ([`distrib_deriv3()`][distrib_deriv3.GPDDistrib],
#' [`distrib_deriv4()`][distrib_deriv4.GPDDistrib]); the response derivatives
#' and the mixed one in `cross_derivatives_families.R`
#' ([`distrib_grad_y()`][distrib_grad_y], [`distrib_hess_y()`][distrib_hess_y],
#' [`distrib_cross_y()`][distrib_cross_y]); the four moments in `moments.R`;
#' and the second-order response derivatives in `theta2_more.R`.
#'
#' @section What the moving endpoint costs:
#' The derivatives are correct as derivatives of the log-density at every
#' admissible point, whatever the sign of \eqn{\xi}. What is not automatic is
#' the license to differentiate under the integral sign, which the Bartlett
#' identities rest on. The expected information exists for \eqn{\xi > -1/2}
#' and is Smith's closed form; at or below that it does not exist and
#' [distrib_expected_hessian()] returns `NA`, along with the classical
#' asymptotics of the maximum likelihood estimator.
#'
#' @seealso [gpd_distrib()] to build one;
#'   [exponential_distrib()] for the \eqn{\xi = 0} case;
#'   [gumbel_distrib()] for the companion family in an analysis of extremes;
#'   [gpd_endpoint()] for the upper limit of the support.
#'
#' @examples
#' d <- gpd_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # The declared bounds are fixed at construction; the true endpoint moves
#' # with the parameters and is finite whenever the shape is negative.
#' c(declared = d@bounds,
#'   actual = c(0, distributions7:::gpd_endpoint(2, -0.4)))
#'
#' # The scale rides a log link and the shape may take either sign.
#' vapply(d@link_params, function(l) l@link_name, character(1))
GPDDistrib <- S7::new_class("GPDDistrib", parent = continuous_distrib)

#' @title The Upper Endpoint of a Generalized Pareto
#'
#' @description
#' Returns \eqn{-\sigma/\xi} where \eqn{\xi < 0} and `Inf` elsewhere: the point
#' beyond which a generalized Pareto puts no mass.
#'
#' @details
#' The endpoint depends on the parameters, which is the whole reason this
#' family needs care. For \eqn{\xi < 0} the support is \eqn{[0, -\sigma/\xi]}
#' and both ends of the interval move as the parameters do, so
#' differentiating an expectation under the integral sign is not automatically
#' licensed. The consequence a user meets is that the expected information
#' exists only for \eqn{\xi > -1/2}; see [gpd_distrib()].
#'
#' `ifelse()` is correct here and is deliberate: both arguments are of the
#' length of the test after recycling, so a vector of shapes gives a vector of
#' endpoints.
#'
#' @param sigma A numeric vector of scales, strictly positive. Nothing is
#'   validated.
#' @param xi A numeric vector of shapes, of any sign.
#'
#' @return A numeric vector of endpoints, of the length of the recycled
#'   inputs. Entries are `Inf` wherever `xi >= 0`.
#'
#' @seealso [gpd_distrib()] for the family and what the moving endpoint costs,
#'   and [distrib_expected_hessian.GPDDistrib()] for where it stops existing.
#'
#' @examples
#' # Finite only for a negative shape.
#' distributions7:::gpd_endpoint(2, c(-0.4, -0.1, 0, 0.3))
#'
#' # The density really is zero at and beyond it.
#' d <- gpd_distrib()
#' distrib_pdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
#'
#' @keywords internal
gpd_endpoint <- function(sigma, xi) ifelse(xi < 0, -sigma / xi, Inf)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Generalized Pareto Density
#' @name distrib_pdf.GPDDistrib
#'
#' @description
#' Computes the generalized Pareto density
#' \deqn{f(y; \sigma, \xi) = \dfrac{1}{\sigma}
#'       \left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi - 1},
#'       \qquad y \ge 0,\; 1 + \xi y/\sigma > 0,}
#' and 0 outside that region. At \eqn{\xi = 0} the limit is the exponential
#' density \eqn{e^{-y/\sigma}/\sigma}.
#'
#' The compiled kernel does not branch on \eqn{\xi = 0}. It writes the
#' log-survival as \eqn{-(y/\sigma)\Lambda(u)} with \eqn{u = \xi y/\sigma} and
#' \eqn{\Lambda(u) = \log(1+u)/u}, which is analytic with \eqn{\Lambda(0) = 1},
#' so every division by the shape disappears and the exponential limit is an
#' ordinary point of the formula. Measured, the density agrees with `dexp` to
#' \eqn{1.4\times10^{-17}} at \eqn{\xi = 0} exactly.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector of observations. A negative value, or one beyond
#'   [gpd_endpoint()], gives a density of 0.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`. `sigma` must be strictly
#'   positive; `xi` may take either sign.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A numeric vector of densities, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale and \eqn{\xi} the shape; \eqn{\sigma} is not
#' the mean, which is \eqn{\sigma/(1-\xi)} and exists only for \eqn{\xi < 1}.
#'
#' @seealso [distrib_cdf.GPDDistrib()] for the distribution function,
#'   [distrib_gradient.GPDDistrib()] for the score,
#'   [gpd_endpoint()] for where the support ends, and [distrib_pdf()] for the
#'   generic.
#'
#' @examples
#' d <- gpd_distrib()
#' y <- c(0.2, 1, 4)
#' th <- list(sigma = 1.5, xi = 0.3)
#'
#' # The formula written out.
#' all.equal(distrib_pdf(d, y, th), (1 + 0.3 * y / 1.5)^(-1 / 0.3 - 1) / 1.5)
#'
#' # It integrates to one.
#' integrate(function(v) distrib_pdf(d, v, th), 0, Inf)$value
#'
#' # Shape zero is the exponential, and the limit is reached by a series
#' # rather than by a branch: the error is linear in xi.
#' vapply(c(0, 1e-10, 1e-6, 1e-3), function(x)
#'   max(abs(distrib_pdf(d, y, list(sigma = 1.5, xi = x)) -
#'           dexp(y, rate = 1 / 1.5))), 0)
#'
#' # A negative shape bounds the support at -sigma/xi.
#' distrib_pdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
S7::method(distrib_pdf, GPDDistrib) <- function(distrib, y, theta, log = FALSE, ..., threads = 1L) {
  out <- gpd_logpdf_cpp(y, theta[[1]], theta[[2]], threads)
  if (log) out else exp(out)
}

#' @title Generalized Pareto Distribution Function
#' @name distrib_cdf.GPDDistrib
#'
#' @description
#' Computes the generalized Pareto distribution function
#' \deqn{F(q; \sigma, \xi) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi},}
#' with \eqn{1 - e^{-q/\sigma}} at \eqn{\xi = 0}. The survival function is
#' formed first and the tail and the logarithm applied to it, so the far tail
#' keeps its digits.
#'
#' The value is 0 below \eqn{q = 0} and 1 at and beyond [gpd_endpoint()], and
#' is clamped to \eqn{[0, 1]} so that rounding cannot return a probability
#' outside the unit interval.
#'
#' @details
#' The branch on \eqn{|\xi| < 10^{-8}} is taken by indexing, and that is not a
#' stylistic choice. `ifelse()` returns a result the length of its **test**, so
#' with it a scalar shape beside a vector of quantiles would collapse the
#' answer to one number. The test is recycled to the answer's length first.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `q`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value
#'   is \eqn{P(Y \le q)}; when `FALSE` it is the survival function, which is
#'   the quantity actually computed.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms with `log.p = TRUE`, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale and \eqn{\xi} the shape.
#'
#' @seealso [distrib_quantile.GPDDistrib()], which inverts this in closed
#'   form, [distrib_pdf.GPDDistrib()] for the density, and [distrib_cdf()] for
#'   the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' q <- c(0.2, 1, 4)
#' th <- list(sigma = 1.5, xi = 0.3)
#'
#' # The formula written out.
#' all.equal(distrib_cdf(d, q, th), 1 - (1 + 0.3 * q / 1.5)^(-1 / 0.3))
#'
#' # A scalar shape with a vector of quantiles returns one value per
#' # quantile, which an ifelse() branch would not.
#' length(distrib_cdf(d, seq(0, 5, by = 1), th))
#'
#' # Shape zero is the exponential.
#' all.equal(distrib_cdf(d, q, list(sigma = 1.5, xi = 0)),
#'           pexp(q, rate = 1 / 1.5))
#'
#' # A negative shape reaches one at the endpoint and stays there.
#' distrib_cdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
S7::method(distrib_cdf, GPDDistrib) <- function(distrib, q, theta,
                                                 lower.tail = TRUE,
                                                 log.p = FALSE, ...) {
  # ifelse() returns a result the length of its TEST, so a scalar shape would
  # collapse a vector of quantiles to one number. The branch is taken with a
  # test recycled to the answer's length instead.
  s <- theta[[1]]
  x <- theta[[2]]
  z <- q / s
  n <- length(z)
  small <- rep_len(abs(x) < 1e-8, n)
  t <- 1 + rep_len(x, n) * z
  surv <- numeric(n)
  surv[small] <- exp(-z[small])
  surv[!small] <- t[!small]^(-1 / rep_len(x, n)[!small])
  surv[q < 0] <- 1
  surv[t <= 0] <- 0
  surv <- pmin(pmax(surv, 0), 1)
  p <- if (lower.tail) 1 - surv else surv
  if (log.p) log(p) else p
}

#' @title Generalized Pareto Quantile Function
#' @name distrib_quantile.GPDDistrib
#'
#' @description
#' Computes the quantiles in closed form,
#' \deqn{Q(p; \sigma, \xi) = \dfrac{\sigma}{\xi}\left\{(1-p)^{-\xi} - 1\right\},}
#' with \eqn{-\sigma\log(1-p)} at \eqn{\xi = 0}. The generalized Pareto is one
#' of the few families here whose quantile function is elementary, so nothing
#' is inverted numerically and [distrib_rng.GPDDistrib()] can use the inverse
#' transform.
#'
#' Near zero the shape branch is taken on \eqn{|\xi| < 10^{-8}} and the value
#' comes from `log1p(-p)`, which keeps its digits for a probability close to
#' zero where `log(1 - p)` would not.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms when `log.p = TRUE`. At `p = 1` the value is the endpoint,
#'   which is `Inf` for a non-negative shape.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `p`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is the survival probability.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm
#'   and is exponentiated first. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of the length of the recycled
#'   inputs.
#'
#' @seealso [distrib_cdf.GPDDistrib()], which it inverts,
#'   [distrib_rng.GPDDistrib()], which draws from it, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' th <- list(sigma = 1.5, xi = 0.3)
#'
#' # The round trip is exact, both functions being closed form.
#' y <- c(0.2, 1, 4)
#' all.equal(distrib_quantile(d, distrib_cdf(d, y, th), th), y)
#'
#' # Shape zero is the exponential quantile.
#' p <- c(0.1, 0.5, 0.9)
#' all.equal(distrib_quantile(d, p, list(sigma = 1.5, xi = 0)),
#'           qexp(p, rate = 1 / 1.5))
#'
#' # A negative shape has a finite upper quantile: the endpoint itself.
#' c(q999 = distrib_quantile(d, 0.999, list(sigma = 2, xi = -0.4)),
#'   endpoint = distributions7:::gpd_endpoint(2, -0.4))
S7::method(distrib_quantile, GPDDistrib) <- function(distrib, p, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  s <- theta[[1]]
  x <- theta[[2]]
  n <- length(p)
  sv <- rep_len(s, n)
  xv <- rep_len(x, n)
  small <- abs(xv) < 1e-8
  out <- numeric(n)
  out[small] <- -sv[small] * log1p(-p[small])
  out[!small] <- sv[!small] * ((1 - p[!small])^(-xv[!small]) - 1) / xv[!small]
  out
}

#' @title Generalized Pareto Random Generation
#' @name distrib_rng.GPDDistrib
#'
#' @description
#' Draws by the inverse transform, `distrib_quantile(distrib, runif(n),
#' theta)`. The quantile function of this family is elementary, so the
#' transform is exact and costs one uniform per draw; the base class's
#' ratio-of-uniforms fallback is not needed.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of length `n`; a component of length 1 is recycled,
#'   so a parameter varying by observation draws one value per observation from
#'   its own member of the family.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws, non-negative and bounded above by
#'   [gpd_endpoint()] when the shape is negative.
#'
#' @seealso [distrib_quantile.GPDDistrib()], which it inverts through,
#'   [distrib_pdf.GPDDistrib()] for the density the draws follow, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' th <- list(sigma = 1.5, xi = 0.3)
#'
#' set.seed(41)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # The mean exists here (xi < 1) and matches sigma/(1 - xi).
#' c(sample = mean(x), theory = mean(d, th))
#'
#' # A negative shape draws inside a bounded support.
#' set.seed(42)
#' xn <- distrib_rng(d, 1e4, list(sigma = 2, xi = -0.4))
#' c(max_draw = max(xn), endpoint = distributions7:::gpd_endpoint(2, -0.4))
#'
#' # Shape zero draws an exponential sample.
#' set.seed(43)
#' x0 <- distrib_rng(d, 1e5, list(sigma = 1.5, xi = 0))
#' c(sample_mean = mean(x0), theory = 1.5)
S7::method(distrib_rng, GPDDistrib) <- function(distrib, n, theta, ...) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Generalized Pareto Score
#' @name distrib_gradient.GPDDistrib
#'
#' @description
#' Computes the two first derivatives of the log-density in closed form, in a
#' compiled kernel. With \eqn{z = y/\sigma}, \eqn{t = 1 + \xi z} and
#' \eqn{u = z/t},
#' \deqn{\dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{(\xi+1)u - 1}{\sigma}, \qquad
#'       \dfrac{\partial \ell}{\partial \xi}
#'         = \dfrac{\log t}{\xi^2} - \left(1 + \dfrac{1}{\xi}\right)u.}
#'
#' The shape component is written this way only away from zero. Both of its
#' terms blow up as \eqn{\xi \to 0} and their difference has the finite limit
#' \eqn{z^2/2 - z}, so the kernel evaluates it through the analytic function
#' \eqn{\Lambda(u) = \log(1+u)/u} instead, whose derivatives come from a
#' recursion above \eqn{|u| = 1/2} and from a Taylor series below it. Measured
#' at \eqn{z = 1}, the component reads \eqn{-0.4967}, \eqn{-0.49997},
#' \eqn{-0.5000000} at \eqn{\xi = 10^{-2}, 10^{-4}, 10^{-8}} against the limit
#' \eqn{-1/2}.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of two numeric vectors, `sigma` and `xi`, each of the
#'   length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\sigma > 0} the
#' scale, \eqn{\xi} the shape, \eqn{z = y/\sigma} and \eqn{t = 1 + \xi z}.
#'
#' @seealso [distrib_hessian.GPDDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.GPDDistrib()] for Smith's closed form, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' y <- c(0.2, 1, 4)
#' th <- list(sigma = 1.5, xi = 0.3)
#' g <- distrib_gradient(d, y, th)
#'
#' # Against a central difference of the log-density.
#' eps <- 1e-6
#' rbind(analytic = g$xi,
#'       numeric = (distrib_pdf(d, y, list(sigma = 1.5, xi = 0.3 + eps),
#'                              log = TRUE) -
#'                  distrib_pdf(d, y, list(sigma = 1.5, xi = 0.3 - eps),
#'                              log = TRUE)) / (2 * eps))
#'
#' # The shape component has a finite limit at zero that its written form
#' # does not: both of its terms diverge and the difference does not.
#' vapply(c(1e-2, 1e-4, 1e-8, 1e-14),
#'        function(x) distrib_gradient(d, 1, list(sigma = 1, xi = x))$xi, 0)
#'
#' # ...and the limit is z^2/2 - z, which at z = 1 is -0.5.
#' distrib_gradient(d, 1, list(sigma = 1, xi = 0))$xi
S7::method(distrib_gradient, GPDDistrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ..., threads = 1L) {
  gpd_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Generalized Pareto Observed Hessian
#' @name distrib_hessian.GPDDistrib
#'
#' @description
#' Computes the three second derivatives of the log-density in closed form, in
#' a compiled kernel. In the notation of [distrib_gradient.GPDDistrib()] the
#' expressions stay short because \eqn{t - \xi z = 1}, which makes
#' \eqn{\partial u/\partial\sigma} equal to \eqn{-z/(\sigma t^2)}.
#'
#' The pure-\eqn{\xi} component carries the same removable singularity the
#' score does, one order worse, and goes through the same analytic function
#' \eqn{\Lambda(u) = \log(1+u)/u}. The two terms that individually diverge as
#' \eqn{\xi \to 0} cancel, and the limit is an ordinary point of the formula.
#'
#' This is the **observed** curvature at the data. The expected one is closed
#' form as well, but only above \eqn{\xi = -1/2}; see
#' [distrib_expected_hessian.GPDDistrib()].
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of three numeric vectors in [hess_names()]'s order:
#'   `sigma_sigma`, `xi_xi`, `sigma_xi`.
#'
#' @section Notation:
#' \eqn{z = y/\sigma}, \eqn{t = 1 + \xi z}, \eqn{u = z/t}, and \eqn{\ell} is
#' the log-density of one observation.
#'
#' @seealso [distrib_gradient.GPDDistrib()] for the order below,
#'   [distrib_deriv3.GPDDistrib()] for the order above,
#'   [distrib_expected_hessian.GPDDistrib()] for the expectation, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' y <- c(0.2, 1, 4)
#' th <- list(sigma = 1.5, xi = 0.3)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # Against a central difference of the score.
#' eps <- 1e-5
#' rbind(analytic = h$sigma_xi,
#'       numeric = (distrib_gradient(d, y, list(sigma = 1.5,
#'                                              xi = 0.3 + eps))$sigma -
#'                  distrib_gradient(d, y, list(sigma = 1.5,
#'                                              xi = 0.3 - eps))$sigma) /
#'                 (2 * eps))
#'
#' # The pure-shape component passes through zero without a branch.
#' vapply(c(1e-2, 1e-6, 0, -1e-6),
#'        function(x) distrib_hessian(d, 1, list(sigma = 1, xi = x))$xi_xi, 0)
S7::method(distrib_hessian, GPDDistrib) <- function(distrib, y, theta,
                                                     scale = c("parameter", "link"), ..., threads = 1L) {
  gpd_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Generalized Pareto Expected Information
#' @name distrib_expected_hessian.GPDDistrib
#'
#' @description
#' Returns Smith's (1985) closed form, valid for \eqn{\xi > -1/2}:
#' \deqn{E\!\left[\dfrac{\partial^2\ell}{\partial\sigma^2}\right]
#'         = \dfrac{-1}{(1+2\xi)\sigma^2}, \qquad
#'       E\!\left[\dfrac{\partial^2\ell}{\partial\sigma\partial\xi}\right]
#'         = \dfrac{-1}{(1+2\xi)\sigma(1+\xi)}, \qquad
#'       E\!\left[\dfrac{\partial^2\ell}{\partial\xi^2}\right]
#'         = \dfrac{-2}{(1+2\xi)(1+\xi)}.}
#' Every entry is free of the data, so `approx` and `nsim` are ignored and
#' [expected_hessian_exact()] answers `TRUE`.
#'
#' @details
#' # Where it stops existing
#'
#' At \eqn{\xi \le -1/2} the information does **not exist** and every component
#' is `NA`. The condition is exactly that the integrand be integrable: written
#' on the probability scale, the second derivative grows like
#' \eqn{(1-u)^{-2|\xi|}} near the upper endpoint, which is integrable if and
#' only if \eqn{|\xi| < 1/2}. Below that point the classical asymptotics of the
#' maximum likelihood estimator do not hold either, so the `NA` reports a
#' quantity that does not exist.
#'
#' Approaching the boundary the information diverges: measured at
#' \eqn{\sigma = 1.5}, the \eqn{\xi} component is \eqn{-2.00} at
#' \eqn{\xi = 0}, \eqn{-7.14} at \eqn{-0.3} and \eqn{-196.1} at \eqn{-0.49}.
#'
#' # Checking it
#'
#' A Monte Carlo average of the observed Hessian is the weaker reference here,
#' and disagreed with this formula by 9 per cent at \eqn{\xi = -0.3} while the
#' formula was right: the second derivative blows up at the upper endpoint, so
#' the sample mean converges slowly. Integrating on the **probability** scale
#' instead, \eqn{E[h] = \int_0^1 h(Q(u))\,du}, turns the endpoint into an
#' ordinary point and agrees to \eqn{10^{-11}}.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector. Its values do not enter the result, which is an
#'   expectation; only its length does, through recycling.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx Ignored: the expectation is closed form. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors in [hess_names()]'s order:
#'   `sigma_sigma`, `xi_xi`, `sigma_xi`. Every entry is `NA` where
#'   `xi <= -0.5`.
#'
#' @references
#' Smith, R. L. (1985). Maximum likelihood estimation in a class of nonregular
#' cases. *Biometrika* 72, 67-90.
#'
#' @seealso [distrib_hessian.GPDDistrib()] for the observed curvature,
#'   [gpd_endpoint()] for the endpoint that causes the condition, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' th <- list(sigma = 1.5, xi = 0.3)
#' e <- distrib_expected_hessian(d, 0, th)
#'
#' # Smith's formula written out.
#' all.equal(unlist(e[c("sigma_sigma", "sigma_xi", "xi_xi")], use.names = FALSE),
#'           c(-1 / (1.6 * 1.5^2), -1 / (1.6 * 1.5 * 1.3), -2 / (1.6 * 1.3)))
#'
#' # It diverges as the shape approaches -1/2, and does not exist at or
#' # below it.
#' t(vapply(c(-0.7, -0.5, -0.49, -0.3, 0),
#'          function(x) {
#'            u <- distrib_expected_hessian(d, 0, list(sigma = 1.5, xi = x))
#'            c(xi = x, sigma_sigma = u$sigma_sigma, xi_xi = u$xi_xi)
#'          }, numeric(3)))
#'
#' # The strategy argument is ignored, the expectation being closed form.
#' identical(distrib_expected_hessian(d, 0, th),
#'           distrib_expected_hessian(d, 0, th, approx = "mc", nsim = 10))
S7::method(distrib_expected_hessian, GPDDistrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"),
                                                              approx = c("opg", "bartlett", "integrate", "mc"),
                                                              nsim = 10000, ..., threads = 1L) {
  gpd_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Generalized Pareto Distribution Object
#'
#' @description
#' Builds a generalized Pareto distribution object with scale \eqn{\sigma > 0}
#' and shape \eqn{\xi}. It is the family of exceedances over a high threshold,
#' and the natural companion of [gumbel_distrib()] in an analysis of extremes.
#'
#' At \eqn{\xi = 0} it is the exponential, at \eqn{\xi > 0} it has a
#' polynomial tail, and at \eqn{\xi < 0} it has a **finite upper endpoint** at
#' \eqn{-\sigma/\xi}. That last case makes the family non-regular, and it is
#' the one thing to know before using it.
#'
#' @param link_sigma A `linkfunctions7` link object for the scale
#'   \eqn{\sigma}, which must be strictly positive. Defaults to
#'   [linkfunctions7::log_link()].
#' @param link_xi A link object for the shape \eqn{\xi}, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()], the shape
#'   being free to take either sign.
#'
#' @details
#' # Density and distribution function
#'
#' \deqn{f(y) = \dfrac{1}{\sigma}\left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi-1},
#'       \qquad
#'       F(q) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi}.}
#' At \eqn{\xi = 0} both reduce to the exponential. The implementation reaches
#' that limit through a series, so the parameter may pass through zero during
#' a fit. The quantile function is elementary,
#' \eqn{Q(p) = \sigma\{(1-p)^{-\xi}-1\}/\xi}, so the generator is an exact
#' inverse transform.
#'
#' # It is not parametrized by its mean
#'
#' The mean is \eqn{\sigma/(1-\xi)} and exists only for \eqn{\xi < 1}, so a
#' mean parametrization would leave the family undescribable exactly where it
#' is most used, in the heavy-tailed regime. The same argument keeps
#' [mv_sigma()] and [variance()] apart for the multivariate \eqn{t}: a
#' parametrization must not depend on a moment that need not exist.
#'
#' # The support depends on the parameters
#'
#' For \eqn{\xi < 0} the support is \eqn{[0, -\sigma/\xi]}, and this is the
#' first family here of which that is true. What it costs is the automatic
#' license to differentiate under the integral sign, on which the Bartlett
#' identities rest. Two things survive and one does not:
#'
#' - the derivatives returned are correct as derivatives of the log-density at
#'   every admissible point, whatever the sign of \eqn{\xi};
#' - the expected information exists and is Smith's closed form for
#'   \eqn{\xi > -1/2}. The condition is exactly that the integrand be
#'   integrable: near the upper endpoint the second derivative grows like
#'   \eqn{(1-u)^{-2|\xi|}} on the probability scale, which is integrable if and
#'   only if \eqn{|\xi| < 1/2};
#' - below \eqn{\xi = -1/2} the information does not exist,
#'   [distrib_expected_hessian()] returns `NA`, and the classical asymptotics
#'   of the maximum likelihood estimator do not hold (Smith, 1985).
#'
#' The object's `bounds` are `c(0, Inf)` because they are fixed at
#' construction while the true endpoint moves with the parameters. The density
#' is zero beyond it, so nothing computes a wrong number, and a caller reading
#' `bounds` should ask [gpd_endpoint()] instead.
#'
#' # Moments
#'
#' The mean is \eqn{\sigma/(1-\xi)} for \eqn{\xi < 1}, the variance
#' \eqn{\sigma^2/\{(1-\xi)^2(1-2\xi)\}} for \eqn{\xi < 1/2}, the skewness needs
#' \eqn{\xi < 1/3} and the kurtosis \eqn{\xi < 1/4. } Each returns `Inf` above
#' its threshold, so a fitted object can report a mean where it has no
#' variance.
#'
#' # Parameter domains
#'
#' - \eqn{\sigma \in (0, \infty)}
#' - \eqn{\xi \in (-\infty, \infty)}
#'
#' @return An S7 object of class [GPDDistrib], inheriting from
#'   `continuous_distrib`. Its `params` are `c("sigma", "xi")`, its `bounds`
#'   `c(0, Inf)`, and its `link_params` the two links given here.
#'
#' @references
#' Smith, R. L. (1985). Maximum likelihood estimation in a class of nonregular
#' cases. *Biometrika* 72, 67-90.
#'
#' Davison, A. C. and Smith, R. L. (1990). Models for exceedances over high
#' thresholds. *Journal of the Royal Statistical Society B* 52, 393-442.
#'
#' @seealso [exponential_distrib()] for the \eqn{\xi = 0} case,
#'   [gumbel_distrib()] and [weibull1_distrib()] for the other extreme-value
#'   families, [gpd_endpoint()] for the moving support, and [GPDDistrib] for
#'   the class and its method list.
#'
#' @importFrom linkfunctions7 log_link identity_link
#' @importFrom stats runif
#'
#' @examples
#' d <- gpd_distrib()
#' d@params
#' th <- list(sigma = 1.5, xi = 0.3)
#'
#' distrib_pdf(d, c(0.2, 1, 4), th)
#'
#' # Shape zero is the exponential, reached by a series.
#' all.equal(distrib_pdf(d, c(0.2, 1, 4), list(sigma = 1.5, xi = 0)),
#'           dexp(c(0.2, 1, 4), rate = 1 / 1.5))
#'
#' # Which moments exist depends on the shape, and the mean can be finite
#' # where the variance is not.
#' t(vapply(c(-0.3, 0.2, 0.4, 0.6), function(x) {
#'   p <- list(sigma = 1.5, xi = x)
#'   c(xi = x, mean = mean(d, p), variance = variance(d, p))
#' }, numeric(3)))
#'
#' # The information exists only above -1/2.
#' distrib_expected_hessian(d, 0, list(sigma = 1.5, xi = -0.7))
#'
#' # A fit recovers both parameters.
#' set.seed(44)
#' x <- distrib_rng(d, 4000, th)
#' coef(fit_distrib(d, x))
#'
#' @export
gpd_distrib <- function(link_sigma = log_link(), link_xi = identity_link()) {
  GPDDistrib(
    distrib_name = "generalized pareto", dimension = "univariate",
    bounds = c(0, Inf),
    params = c("sigma", "xi"),
    params_interpretation = c(sigma = "scale", xi = "shape"),
    n_params = 2,
    params_bounds = list(sigma = c(0, Inf), xi = c(-Inf, Inf)),
    link_params = list(sigma = link_sigma, xi = link_xi)
  )
}
