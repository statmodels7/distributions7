#' @include distrib.R generics.R
NULL

#' @title Generalized Gamma Distribution Class
#' @name GenGamma1Distrib
#'
#' @description
#' The S7 class of the generalized gamma in Stacy's three-parameter form: a
#' scale \eqn{a > 0} and two shapes \eqn{d > 0} and \eqn{p > 0}, with density
#' \eqn{f(y) = p\,y^{d-1}e^{-(y/a)^p}/\{a^d\Gamma(d/p)\}} on \eqn{y > 0}.
#'
#' It is the flexible family for a positive response, and it nests four others
#' exactly: the gamma at \eqn{p = 1}, the Weibull at \eqn{d = p}, the
#' exponential at \eqn{d = p = 1} and the half-normal at
#' \eqn{a = \sqrt2, d = 1, p = 2}. Choosing between them becomes something to
#' estimate.
#'
#' Build one with [gengamma1_distrib()], which supplies the three link
#' functions. This page documents the raw S7 constructor, which validates none
#' of the relationships between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `GenGamma1Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [gengamma1_distrib()] the properties hold `"gengamma1"`, `"univariate"`,
#'   `c(0, Inf)`, `c("a", "d", "p")`, the interpretations
#'   `c(a = "scale", d = "shape", p = "power")`, `3`, and three domains
#'   \eqn{(0, \infty)}.
#'
#' @section Methods:
#' Registered in this file, the last three compiled:
#'   [`distrib_pdf()`][distrib_pdf.GenGamma1Distrib],
#'   [`distrib_cdf()`][distrib_cdf.GenGamma1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.GenGamma1Distrib],
#'   [`distrib_rng()`][distrib_rng.GenGamma1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.GenGamma1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.GenGamma1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GenGamma1Distrib].
#'
#' Registered elsewhere: the third and fourth orders in `gengamma1_higher.R`;
#' the response derivatives in `cross_derivatives_families.R`; the mixed one in
#' `cross_derivatives_simple.R`; and the four moments in `moments.R`.
#'
#' @section The one representation everything rests on:
#' \eqn{u = (Y/a)^p} is Gamma with shape \eqn{k = d/p} and unit rate, exactly.
#' The distribution function is that Gamma's, the quantile function inverts it,
#' the generator raises a Gamma draw to the power \eqn{1/p}, and every
#' expectation the information needs is a moment of \eqn{u}. It is the same
#' device the Weibull and the Gumbel use, where the corresponding variable is
#' standard exponential.
#'
#' @seealso [gengamma1_distrib()] to build one;
#'   [gamma2_distrib()], [weibull1_distrib()] and [exponential_distrib()] for
#'   the families it nests; [lognormal1_distrib()] for the limit it approaches.
#'
#' @examples
#' d <- gengamma1_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # All three parameters are positive, so all three ride a log link.
#' vapply(d@link_params, function(l) l@link_name, character(1))
#'
#' # Three of the four exact special cases, at one point each.
#' y <- c(0.5, 1.5, 4)
#' rbind(gamma = distrib_pdf(d, y, list(a = 2, d = 3, p = 1)) -
#'               dgamma(y, shape = 3, scale = 2),
#'       weibull = distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)) -
#'                 dweibull(y, shape = 1.5, scale = 2),
#'       half_normal = distrib_pdf(d, y, list(a = sqrt(2), d = 1, p = 2)) -
#'                     2 * dnorm(y))
GenGamma1Distrib <- S7::new_class("GenGamma1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Generalized Gamma Density
#' @name distrib_pdf.GenGamma1Distrib
#'
#' @description
#' Computes the generalized gamma density in Stacy's form,
#' \deqn{f(y; a, d, p) = \dfrac{p}{a^{d}\,\Gamma(d/p)}\, y^{d-1}
#'       e^{-(y/a)^{p}}, \qquad y > 0,}
#' and 0 at or below zero. The compiled kernel works on the log scale
#' throughout, so a large \eqn{d} does not form \eqn{a^d} or \eqn{y^{d-1}}
#' separately.
#'
#' The three parameters do three separate things: \eqn{a} sets the scale,
#' \eqn{d} the behavior near the origin (the density vanishes there for
#' \eqn{d > 1} and diverges for \eqn{d < 1}), and \eqn{p} the weight of the
#' upper tail.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of observations. A non-positive value gives a
#'   density of 0.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`. All three must be strictly
#'   positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A numeric vector of densities, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{a > 0} is the scale, \eqn{d > 0} and \eqn{p > 0} the two shapes, and
#' \eqn{\Gamma} the gamma function. \eqn{a} is not the mean, which is
#' \eqn{a\,\Gamma\{(d+1)/p\}/\Gamma(d/p)}.
#'
#' @seealso [distrib_cdf.GenGamma1Distrib()] for the distribution function,
#'   [distrib_gradient.GenGamma1Distrib()] for the score, and [distrib_pdf()]
#'   for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' y <- c(0.5, 1.5, 4)
#' th <- list(a = 2, d = 3, p = 1.5)
#'
#' # The formula written out.
#' all.equal(distrib_pdf(d, y, th),
#'           1.5 / (2^3 * gamma(3 / 1.5)) * y^(3 - 1) * exp(-(y / 2)^1.5))
#'
#' # It integrates to one.
#' integrate(function(v) distrib_pdf(d, v, th), 0, Inf)$value
#'
#' # The four exact special cases.
#' c(gamma = max(abs(distrib_pdf(d, y, list(a = 2, d = 3, p = 1)) -
#'                   dgamma(y, shape = 3, scale = 2))),
#'   weibull = max(abs(distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)) -
#'                     dweibull(y, shape = 1.5, scale = 2))),
#'   exponential = max(abs(distrib_pdf(d, y, list(a = 2, d = 1, p = 1)) -
#'                         dexp(y, rate = 1 / 2))),
#'   half_normal = max(abs(distrib_pdf(d, y, list(a = sqrt(2), d = 1, p = 2)) -
#'                         2 * dnorm(y))))
#'
#' # d decides what happens at the origin.
#' distrib_pdf(d, 1e-8, list(a = 2, d = c(0.5, 1, 2), p = 1.5))
S7::method(distrib_pdf, GenGamma1Distrib) <- function(distrib, y, theta, log = FALSE, ..., threads = 1L) {
  out <- gengamma_logpdf_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  if (log) out else exp(out)
}

#' @title Generalized Gamma Distribution Function
#' @name distrib_cdf.GenGamma1Distrib
#'
#' @description
#' Computes \eqn{F(q) = P(d/p,\, (q/a)^{p})}, the regularized lower incomplete
#' gamma function. The identity is exact:
#' \eqn{u = (Y/a)^p} is Gamma with shape \eqn{d/p} and unit rate, so the
#' distribution function is that Gamma's evaluated at \eqn{(q/a)^p} and the
#' whole computation is one `pgamma` call.
#'
#' `lower.tail` and `log.p` are passed straight to `pgamma`, so the survival
#' function and the log-probability keep the accuracy R's own routine gives
#' them in the far tail.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param q A numeric vector of quantiles. A negative value is clamped to zero
#'   before the transformation, so the value there is 0.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `q`. All three must be strictly
#'   positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value
#'   is \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}, computed by
#'   `pgamma` rather than as one minus the other.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms with `log.p = TRUE`, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{a} is the scale, \eqn{d} and \eqn{p} the shapes, and
#' \eqn{P(s, x) = \gamma(s,x)/\Gamma(s)} the regularized lower incomplete gamma
#' function, which is `pgamma(x, shape = s)` in R.
#'
#' @seealso [distrib_quantile.GenGamma1Distrib()], which inverts this,
#'   [distrib_pdf.GenGamma1Distrib()] for the density, and [distrib_cdf()] for
#'   the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' q <- c(0.5, 1.5, 4)
#' th <- list(a = 2, d = 3, p = 1.5)
#'
#' # The identity written out.
#' all.equal(distrib_cdf(d, q, th), pgamma((q / 2)^1.5, shape = 3 / 1.5))
#'
#' # Against a direct quadrature of the density.
#' rbind(gamma_identity = distrib_cdf(d, q, th),
#'       quadrature = vapply(q, function(u)
#'         integrate(function(v) distrib_pdf(d, v, th), 0, u)$value, 0))
#'
#' # The upper tail keeps its digits where one minus the lower one would not.
#' c(upper = distrib_cdf(d, 20, th, lower.tail = FALSE),
#'   one_minus_lower = 1 - distrib_cdf(d, 20, th))
S7::method(distrib_cdf, GenGamma1Distrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  a <- theta[[1]]; d <- theta[[2]]; p <- theta[[3]]
  w <- pmax(q, 0)^p / a^p
  stats::pgamma(w, shape = d / p, rate = 1,
                lower.tail = lower.tail, log.p = log.p)
}

#' @title Generalized Gamma Quantile Function
#' @name distrib_quantile.GenGamma1Distrib
#'
#' @description
#' Computes \eqn{Q(u) = a\,\{Q_\Gamma(u;\, d/p)\}^{1/p}}, inverting the same
#' representation the distribution function uses: \eqn{(Y/a)^p} is Gamma with
#' shape \eqn{d/p}, so its quantile raised to the power \eqn{1/p} and scaled by
#' \eqn{a} is the generalized gamma's. The whole computation is one `qgamma`
#' call, and nothing is inverted by root finding.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms when `log.p = TRUE`. Note that the family's third parameter is
#'   also called `p`; here the argument is the probability, and the parameter
#'   is read from `theta`.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of the probabilities. All three must
#'   be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is the survival probability, passed
#'   through to `qgamma`.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of the length of the recycled
#'   inputs, strictly positive for a probability in \eqn{(0, 1)}.
#'
#' @seealso [distrib_cdf.GenGamma1Distrib()], which it inverts,
#'   [distrib_rng.GenGamma1Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' th <- list(a = 2, d = 3, p = 1.5)
#'
#' # The round trip through the distribution function.
#' y <- c(0.5, 1.5, 4)
#' all.equal(distrib_quantile(d, distrib_cdf(d, y, th), th), y)
#'
#' # The identity written out.
#' u <- c(0.1, 0.5, 0.9)
#' all.equal(distrib_quantile(d, u, th),
#'           2 * qgamma(u, shape = 3 / 1.5)^(1 / 1.5))
#'
#' # The median sits below the mean, the family being right-skewed here.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, GenGamma1Distrib) <- function(distrib, p, theta,
                                                           lower.tail = TRUE,
                                                           log.p = FALSE, ...) {
  a <- theta[[1]]; d <- theta[[2]]; pw <- theta[[3]]
  g <- stats::qgamma(p, shape = d / pw, rate = 1,
                     lower.tail = lower.tail, log.p = log.p)
  a * g^(1 / pw)
}

#' @title Generalized Gamma Random Generation
#' @name distrib_rng.GenGamma1Distrib
#'
#' @description
#' Draws by the family's own defining representation: a Gamma variate with
#' shape \eqn{d/p} and unit rate, raised to the power \eqn{1/p} and multiplied
#' by \eqn{a}. It is exact and costs one `rgamma` call, so neither inversion
#' nor the base class's ratio-of-uniforms fallback is involved.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of length `n`; a component of length 1 is recycled,
#'   so a parameter varying by observation draws one value per observation from
#'   its own member of the family.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_pdf.GenGamma1Distrib()] for the density the draws follow,
#'   [mean.GenGamma1Distrib()] for the moments they reproduce, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' th <- list(a = 2, d = 3, p = 1.5)
#'
#' set.seed(51)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # Two moments against their closed forms.
#' rbind(sample = c(mean(x), var(x)),
#'       theory = c(mean(d, th), variance(d, th)))
#'
#' # The representation, written out at the same seed.
#' set.seed(51)
#' all.equal(x, 2 * rgamma(2e5, shape = 3 / 1.5)^(1 / 1.5))
S7::method(distrib_rng, GenGamma1Distrib) <- function(distrib, n, theta, ...) {
  a <- theta[[1]]; d <- theta[[2]]; p <- theta[[3]]
  a * stats::rgamma(n, shape = d / p, rate = 1)^(1 / p)
}

#' @title Generalized Gamma Score
#' @name distrib_gradient.GenGamma1Distrib
#'
#' @description
#' Computes the three first derivatives of the log-density in closed form, in a
#' compiled kernel. With \eqn{w = (y/a)^{p}}, \eqn{L = \log(y/a)} and
#' \eqn{k = d/p},
#' \deqn{\dfrac{\partial\ell}{\partial a} = \dfrac{pw - d}{a}, \qquad
#'       \dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
#'       \dfrac{\partial\ell}{\partial p} = \dfrac{1}{p}
#'         + \dfrac{d\,\psi(k)}{p^{2}} - wL.}
#' Only the digamma function appears beyond elementary operations, and it
#' appears at one argument, \eqn{k = d/p}.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of observations, strictly positive.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`. All three must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of three numeric vectors, `a`, `d` and `p`, each of the
#'   length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\psi} the digamma
#' function `digamma()`, and \eqn{a}, \eqn{d}, \eqn{p} the scale and the two
#' shapes.
#'
#' @seealso [distrib_hessian.GenGamma1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.GenGamma1Distrib()] for the closed-form
#'   information, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' y <- c(0.5, 1.5, 4)
#' th <- list(a = 2, d = 3, p = 1.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The scale component written out.
#' all.equal(g$a, (1.5 * (y / 2)^1.5 - 3) / 2)
#'
#' # All three against numerical differentiation of the log-likelihood.
#' f <- function(v) sum(distrib_pdf(d, y, list(a = v[1], d = v[2], p = v[3]),
#'                                  log = TRUE))
#' rbind(analytic = vapply(g, sum, 0),
#'       numeric = numDeriv::grad(f, c(2, 3, 1.5)))
#'
#' # The score sums to nearly zero at the maximum likelihood estimate.
#' set.seed(52)
#' x <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, x)
#' vapply(distrib_gradient(d, x, as.list(coef(fit))), sum, 0) / 2000
S7::method(distrib_gradient, GenGamma1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  gengamma_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Generalized Gamma Observed Hessian
#' @name distrib_hessian.GenGamma1Distrib
#'
#' @description
#' Computes the six second derivatives of the log-density in closed form, in a
#' compiled kernel, by differentiating the expressions of
#' [distrib_gradient.GenGamma1Distrib()] again.
#'
#' One component is worth naming: the mixed \eqn{a}-\eqn{d} entry is
#' \eqn{-1/a}, free of the data. The scale and the first shape meet in the
#' log-density only through the term \eqn{-d\log a}, which is bilinear, so
#' its second derivative carries no observation.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector of observations, strictly positive.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`. All three must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of six numeric vectors in [hess_names()]'s order:
#'   `a_a`, `d_d`, `p_p`, `a_d`, `a_p`, `d_p`.
#'
#' @section Notation:
#' \eqn{w = (y/a)^p}, \eqn{L = \log(y/a)}, \eqn{k = d/p}, and \eqn{\psi},
#' \eqn{\psi'} are the digamma and trigamma functions.
#'
#' @seealso [distrib_gradient.GenGamma1Distrib()] for the order below,
#'   [distrib_deriv3.GenGamma1Distrib()] for the order above,
#'   [distrib_expected_hessian.GenGamma1Distrib()] for its expectation, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' y <- c(0.5, 1.5, 4)
#' th <- list(a = 2, d = 3, p = 1.5)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # The scale-shape entry is -1/a at every observation.
#' h$a_d
#'
#' # Against a central difference of the score.
#' eps <- 1e-5
#' rbind(analytic = h$d_p,
#'       numeric = (distrib_gradient(d, y, list(a = 2, d = 3,
#'                                              p = 1.5 + eps))$d -
#'                  distrib_gradient(d, y, list(a = 2, d = 3,
#'                                              p = 1.5 - eps))$d) / (2 * eps))
S7::method(distrib_hessian, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ..., threads = 1L) {
  gengamma_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Generalized Gamma Expected Information
#' @name distrib_expected_hessian.GenGamma1Distrib
#'
#' @description
#' Returns the expected second derivatives in closed form. Every expectation
#' the observed Hessian needs is a moment of \eqn{u = (Y/a)^{p}}, which is
#' Gamma with shape \eqn{k = d/p} and unit rate:
#' \deqn{E[u] = k, \qquad E[u\log u] = k\,\psi(k+1), \qquad
#'       E[u(\log u)^{2}] = k\{\psi(k+1)^{2} + \psi'(k+1)\}.}
#' So `approx` and `nsim` are ignored and [expected_hessian_exact()] answers
#' `TRUE`. It is the same device the Weibull and the Gumbel use, where the
#' corresponding variable is standard exponential and the expectations are
#' derivatives of \eqn{\Gamma} at 2.
#'
#' @param distrib A `GenGamma1Distrib` object, from [gengamma1_distrib()].
#' @param y A numeric vector. Its values do not enter the result, which is an
#'   expectation; only its length does, through recycling.
#' @param theta A named list with components `a`, `d` and `p`, each a numeric
#'   vector of length 1 or of the length of `y`. All three must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx Ignored: the expectation is closed form. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of six numeric vectors in [hess_names()]'s order:
#'   `a_a`, `d_d`, `p_p`, `a_d`, `a_p`, `d_p`. Every entry is negative
#'   definite as a matrix and free of the data.
#'
#' @section Identification:
#' The three parameters are weakly identified together. \eqn{d} and \eqn{p}
#' enter the density largely through their ratio \eqn{k = d/p}, and the
#' information reflects that: measured at \eqn{a = 2, d = 3, p = 1.5} its
#' eigenvalues are 4.72, 0.138 and 0.00321, a condition number of 1470, and the
#' flat direction is \eqn{(0.79, -0.52, 0.33)} in \eqn{(a, d, p)}. A fit of all
#' three wants several hundred observations; holding one with [fixed()] is
#' often the better model.
#'
#' @section Notation:
#' \eqn{k = d/p}, \eqn{\psi} and \eqn{\psi'} the digamma and trigamma
#' functions, and \eqn{u = (Y/a)^p} the Gamma-distributed transformation the
#' family is built on.
#'
#' @seealso [distrib_hessian.GenGamma1Distrib()] for the observed curvature,
#'   [expected_hessian_exact()] for the predicate, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- gengamma1_distrib()
#' th <- list(a = 2, d = 3, p = 1.5)
#' e <- distrib_expected_hessian(d, 0, th)
#' names(e)
#'
#' # The scale-shape entry is the observed one, -1/a, the term being bilinear.
#' c(expected = e$a_d, observed = -1 / 2)
#'
#' # The strategy argument is ignored, the expectation being closed form.
#' identical(e, distrib_expected_hessian(d, 0, th, approx = "mc", nsim = 5))
#' distributions7:::expected_hessian_exact(d)
#'
#' # The information is ill conditioned: d and p act largely through d/p.
#' nm <- c("a", "d", "p")
#' M <- matrix(c(e$a_a, e$a_d, e$a_p, e$a_d, e$d_d, e$d_p,
#'               e$a_p, e$d_p, e$p_p), 3, 3, dimnames = list(nm, nm))
#' ev <- eigen(-M)
#' c(eigenvalues = ev$values, condition = max(ev$values) / min(ev$values))
#' round(ev$vectors[, 3], 3)   # the flat direction
S7::method(distrib_expected_hessian, GenGamma1Distrib) <- function(distrib, y, theta,
                                                                   scale = c("parameter", "link"),
                                                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                                                   nsim = 10000, ..., threads = 1L) {
  gengamma_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Generalized Gamma Distribution Object
#'
#' @description
#' Builds a generalized gamma distribution object in Stacy's form, with a scale
#' \eqn{a > 0} and two shapes \eqn{d > 0} and \eqn{p > 0}. It is the flexible
#' family for a positive response, and it makes the choice between the gamma,
#' the Weibull and the exponential something to estimate.
#'
#' @param link_a A `linkfunctions7` link object for the scale \eqn{a}, which
#'   must be strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_d A link object for the first shape \eqn{d}, also strictly
#'   positive. Defaults to [linkfunctions7::log_link()].
#' @param link_p A link object for the second shape \eqn{p}, also strictly
#'   positive. Defaults to [linkfunctions7::log_link()].
#'
#' @details
#' # Density
#'
#' \deqn{f(y) = \dfrac{p}{a^{d}\,\Gamma(d/p)}\,y^{d-1}e^{-(y/a)^{p}},
#'       \qquad y > 0.}
#' The three parameters do three separate things: \eqn{a} sets the scale,
#' \eqn{d} the behavior at the origin, and \eqn{p} the weight of the upper
#' tail.
#'
#' # What it nests
#'
#' Stacy's parametrization is chosen to make these visible:
#'
#' - \eqn{p = 1} is the [gamma][gamma2_distrib] with shape \eqn{d} and scale
#'   \eqn{a};
#' - \eqn{d = p} is the [Weibull][weibull1_distrib] with shape \eqn{p} and
#'   scale \eqn{a};
#' - \eqn{d = p = 1} is the [exponential][exponential_distrib];
#' - \eqn{a = \sqrt2}, \eqn{d = 1}, \eqn{p = 2} is the half-normal;
#' - \eqn{p \to 0} with \eqn{d/p} held large approaches the
#'   [lognormal][lognormal1_distrib].
#'
#' The first four are exact and are checked in the examples. The fifth is a
#' limit, reached at no admissible value.
#'
#' # Score and information
#'
#' With \eqn{w = (y/a)^{p}}, \eqn{L = \log(y/a)} and \eqn{k = d/p},
#' \deqn{\dfrac{\partial\ell}{\partial a} = \dfrac{pw-d}{a}, \qquad
#'       \dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
#'       \dfrac{\partial\ell}{\partial p} = \dfrac{1}{p}
#'         + \dfrac{d\psi(k)}{p^{2}} - wL,}
#' and the expected information is **closed form**. The reason is the one
#' representation the whole family rests on: \eqn{u = (Y/a)^p} is Gamma with
#' shape \eqn{k} and unit rate, so every expectation the Hessian needs is one
#' of \eqn{E[u] = k}, \eqn{E[u\log u] = k\psi(k+1)} and
#' \eqn{E[u(\log u)^2] = k\{\psi(k+1)^2 + \psi'(k+1)\}}. The same
#' representation gives the distribution function, the quantile function and
#' the generator.
#'
#' # Moments
#'
#' \eqn{E[Y^{r}] = a^{r}\Gamma\{(d+r)/p\}/\Gamma(d/p)}, finite for every
#' \eqn{r > -d}. Unlike the generalized Pareto's, they exist at every parameter
#' value.
#'
#' # Identification
#'
#' The three parameters are **weakly identified together**: \eqn{d} and \eqn{p}
#' enter the density largely through their ratio, and the profile likelihood in
#' that direction is flat. Measured at \eqn{a = 2, d = 3, p = 1.5}, the
#' information's eigenvalues are 4.72, 0.138 and 0.00321, a condition number of
#' 1470, with the flat direction \eqn{(0.79, -0.52, 0.33)}. A fit of all three
#' wants several hundred observations: on 200 the standard errors are 0.88,
#' 0.59 and 0.44 against estimates of 2.20, 2.89 and 1.68, and on 2000 they are
#' 0.31, 0.21 and 0.14. Holding one with [fixed()] is often the better model.
#'
#' # Parameter domains
#'
#' - \eqn{a \in (0, \infty)}
#' - \eqn{d \in (0, \infty)}
#' - \eqn{p \in (0, \infty)}
#'
#' @return An S7 object of class [GenGamma1Distrib], inheriting from
#'   `continuous_distrib`. Its `params` are `c("a", "d", "p")`, its `bounds`
#'   `c(0, Inf)`, and its `link_params` the three links given here.
#'
#' @references
#' Stacy, E. W. (1962). A generalization of the gamma distribution.
#' *Annals of Mathematical Statistics* 33, 1187-1192.
#'
#' @seealso [gamma2_distrib()], [weibull1_distrib()],
#'   [exponential_distrib()] and [lognormal1_distrib()] for the families it
#'   nests or approaches, [fixed()] for holding a shape, and
#'   [GenGamma1Distrib] for the class and its method list.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats pgamma qgamma rgamma
#'
#' @examples
#' d <- gengamma1_distrib()
#' d@params
#' th <- list(a = 2, d = 3, p = 1.5)
#'
#' distrib_pdf(d, c(0.5, 2, 5), th)
#'
#' # The four exact special cases, at three points each.
#' y <- c(0.5, 2, 5)
#' c(gamma = max(abs(distrib_pdf(d, y, list(a = 2, d = 3, p = 1)) -
#'                   dgamma(y, shape = 3, scale = 2))),
#'   weibull = max(abs(distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)) -
#'                     dweibull(y, shape = 1.5, scale = 2))),
#'   exponential = max(abs(distrib_pdf(d, y, list(a = 2, d = 1, p = 1)) -
#'                         dexp(y, rate = 1 / 2))),
#'   half_normal = max(abs(distrib_pdf(d, y, list(a = sqrt(2), d = 1, p = 2)) -
#'                         2 * dnorm(y))))
#'
#' # The moment formula.
#' c(closed = 2 * gamma((3 + 1) / 1.5) / gamma(3 / 1.5), ours = mean(d, th))
#'
#' # Two fits at two sizes: the standard errors say how much data the three
#' # parameters want.
#' for (n in c(200, 2000)) {
#'   set.seed(52)
#'   f <- fit_distrib(d, distrib_rng(d, n, th))
#'   print(rbind(estimate = coef(f), se = sqrt(diag(vcov(f)))))
#' }
#'
#' @export
gengamma1_distrib <- function(link_a = log_link(), link_d = log_link(),
                             link_p = log_link()) {
  GenGamma1Distrib(
    distrib_name = "gengamma1", dimension = "univariate",
    bounds = c(0, Inf),
    params = c("a", "d", "p"),
    params_interpretation = c(a = "scale", d = "shape", p = "power"),
    n_params = 3,
    params_bounds = list(a = c(0, Inf), d = c(0, Inf), p = c(0, Inf)),
    link_params = list(a = link_a, d = link_d, p = link_p)
  )
}
