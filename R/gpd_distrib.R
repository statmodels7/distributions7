#' @include distrib.R generics.R
NULL

#' @title Generalized Pareto Distribution Class
#' @name GPDDistrib
#'
#' @description
#' This is the S7 class of the generalized Pareto family, the law of exceedances
#' over a high threshold, parametrized by a scale \eqn{\sigma > 0} and a shape
#' \eqn{\xi}: \eqn{f(y) = \sigma^{-1}(1 + \xi y/\sigma)^{-1/\xi-1}} on \eqn{y
#' \ge 0}. At \eqn{\xi = 0} it is the exponential, reached by a series, so a fit
#' may move the shape through zero.
#'
#' For \eqn{\xi < 0} the support ends at \eqn{-\sigma/\xi}, so it moves with
#' the parameters. Because of this the family is non-regular, and for the
#' same reason its expected information exists only above \eqn{\xi = -1/2}.
#'
#' Build one with [gpd_distrib()], which supplies the two link functions. This
#' page documents the raw S7 constructor, which takes the parent's
#' properties; the validator of `distrib` checks them, and the constructor
#' adds no check of its own.
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
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.GPDDistrib],
#'   [`distrib_cdf()`][distrib_cdf.GPDDistrib],
#'   [`distrib_quantile()`][distrib_quantile.GPDDistrib] and
#'   [`distrib_rng()`][distrib_rng.GPDDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.GPDDistrib],
#'   [`distrib_hessian()`][distrib_hessian.GPDDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.GPDDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.GPDDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.GPDDistrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.GPDDistrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.GPDDistrib] and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.GPDDistrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.GPDDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.GPDDistrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.GPDDistrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.GPDDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.GPDDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.GPDDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.GPDDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.GPDDistrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.GPDDistrib],
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.GPDDistrib],
#'   [`distrib_deriv3_cdf()`][distrib_deriv3_cdf.GPDDistrib] and
#'   [`distrib_deriv4_cdf()`][distrib_deriv4_cdf.GPDDistrib];
#' - the moments: [`mean()`][mean.GPDDistrib],
#'   [`variance()`][variance.GPDDistrib], [`skewness()`][skewness.GPDDistrib]
#'   and [`kurtosis()`][kurtosis.GPDDistrib].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @section Support depending on the parameters:
#' The derivatives are correct as derivatives of the log-density at every
#' admissible point, whatever the sign of \eqn{\xi}. The license to
#' differentiate under the integral sign, on which the Bartlett identities
#' rest, is not automatic. The expected information exists for \eqn{\xi > -1/2}
#' and is Smith's closed form. At or below that value it does not exist, and
#' [distrib_expected_hessian()] returns `NA`; the classical asymptotics of the
#' maximum likelihood estimator do not hold either.
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
#' licensed. As a consequence, the expected information exists only for
#' \eqn{\xi > -1/2}; see [gpd_distrib()].
#'
#' The two arguments are recycled to a common length before the branch on
#' the sign of the shape, so a vector of scales beside a scalar shape gives a
#' vector of endpoints.
#'
#' @param sigma A numeric vector of scales, strictly positive. The function
#'   does not validate its arguments.
#' @param xi A numeric vector of shapes, of any sign.
#'
#' @return A numeric vector of endpoints, of the length of the recycled
#'   inputs. Entries are `Inf` wherever `xi >= 0`.
#'
#' @seealso [gpd_distrib()] for the family and what the moving endpoint costs,
#'   and [distrib_expected_hessian.GPDDistrib()] for where it stops existing.
#'
#' @examples
#' # The endpoint is finite only for a negative shape.
#' distributions7:::gpd_endpoint(2, c(-0.4, -0.1, 0, 0.3))
#'
#' # The density really is zero at and beyond it.
#' d <- gpd_distrib()
#' distrib_pdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
#'
#' @keywords internal
gpd_endpoint <- function(sigma, xi) {
  n <- max(length(sigma), length(xi))
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)
  ifelse(xi < 0, -sigma / xi, Inf)
}

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
#' ordinary point of the formula. At \eqn{\xi = 0} the density agrees with
#' `dexp` to rounding error.
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
#' # The density is compared with the formula written out.
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
#' with \eqn{1 - e^{-q/\sigma}} at \eqn{\xi = 0}. A compiled kernel computes
#' the logarithm of the survival function,
#' \eqn{\log S = -(q/\sigma)\Lambda(u)} with \eqn{u = \xi q/\sigma} and
#' \eqn{\Lambda(u) = \log(1+u)/u}, as the density does, so there is no branch
#' at \eqn{\xi = 0}. The upper tail is \eqn{e^{\log S}}, the lower tail
#' \eqn{-\mathrm{expm1}(\log S)}, and their logarithms are \eqn{\log S} and
#' \eqn{\log(1 - e^{\log S})} taken with `expm1()` or `log1p()`, so both
#' tails and their logarithms keep their digits far from the median.
#'
#' The value is 0 at and below \eqn{q = 0} and 1 at and beyond
#' [gpd_endpoint()].
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `q`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value
#'   is \eqn{P(Y \le q)}; when `FALSE` it is the survival function.
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
#' # The distribution function is compared with the formula written out.
#' all.equal(distrib_cdf(d, q, th), 1 - (1 + 0.3 * q / 1.5)^(-1 / 0.3))
#'
#' # A scalar shape with a vector of quantiles returns one value per
#' # quantile.
#' length(distrib_cdf(d, seq(0, 5, by = 1), th))
#'
#' # Both tails keep their digits: a tiny lower tail and the logarithm of a
#' # far upper tail.
#' distrib_cdf(d, 1e-12, list(sigma = 1, xi = 0))
#' distrib_cdf(d, 1000, list(sigma = 1, xi = 0), lower.tail = FALSE,
#'             log.p = TRUE)
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
  # log S = -log(t)/xi from the compiled kernel, which forms t from the
  # exact product; each tail and its logarithm are taken from log S.
  n <- max(length(q), length(theta[[1]]), length(theta[[2]]))
  ls <- gpd_logsurv_cpp(rep_len(as.numeric(q), n), as.numeric(theta[[1]]),
                        as.numeric(theta[[2]]))
  if (!lower.tail) return(if (log.p) ls else exp(ls))
  if (log.p) log1mexp_neg(ls) else -expm1(ls)
}

#' @title Generalized Pareto Quantile Function
#' @name distrib_quantile.GPDDistrib
#'
#' @description
#' Computes the quantiles in closed form,
#' \deqn{Q(p; \sigma, \xi) = \dfrac{\sigma}{\xi}\left\{(1-p)^{-\xi} - 1\right\},}
#' with \eqn{-\sigma\log(1-p)} at \eqn{\xi = 0}. The quantile function is
#' elementary, so the method does not invert numerically and
#' [distrib_rng.GPDDistrib()] uses the inverse transform.
#'
#' The method first forms the logarithm of the survival probability,
#' \eqn{\log S}, from `p` without computing \eqn{1 - p} or exponentiating:
#' `log1p(-p)` for a lower-tail probability, \eqn{\log p} for an upper-tail
#' one, and the corresponding forms when `log.p = TRUE`. The quantile is then
#' \eqn{\sigma\,\mathrm{expm1}(-\xi\log S)/\xi}, and \eqn{-\sigma\log S} at
#' \eqn{\xi = 0}, so a probability close to zero or close to one keeps its
#' digits in whichever tail it is given.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms when `log.p = TRUE`. At `p = 1` the value is the endpoint,
#'   which is `Inf` for a non-negative shape; a value outside the range gives
#'   `NaN`.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `p`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is the survival probability.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
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
#' # The round trip is exact because both functions are closed form.
#' y <- c(0.2, 1, 4)
#' all.equal(distrib_quantile(d, distrib_cdf(d, y, th), th), y)
#'
#' # Shape zero is the exponential quantile.
#' p <- c(0.1, 0.5, 0.9)
#' all.equal(distrib_quantile(d, p, list(sigma = 1.5, xi = 0)),
#'           qexp(p, rate = 1 / 1.5))
#'
#' # A negative shape bounds the quantiles by the endpoint, which p = 1
#' # reaches.
#' c(q999 = distrib_quantile(d, 0.999, list(sigma = 2, xi = -0.4)),
#'   q1 = distrib_quantile(d, 1, list(sigma = 2, xi = -0.4)),
#'   endpoint = distributions7:::gpd_endpoint(2, -0.4))
S7::method(distrib_quantile, GPDDistrib) <- function(distrib, p, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  # ls is the logarithm of the survival probability, formed from the input
  # without 1 - p or exp(); Q = sigma expm1(-xi ls)/xi, and -sigma ls at 0.
  n <- max(length(p), length(theta[[1]]), length(theta[[2]]))
  p <- rep_len(as.numeric(p), n)
  ls <- if (lower.tail) {
    if (log.p) log1mexp_neg(p) else log1p(-p)
  } else {
    if (log.p) p else log(p)
  }
  # a log-survival above zero comes from a probability outside [0, 1]
  ls[!is.na(ls) & ls > 0] <- NaN
  sv <- rep_len(theta[[1]], n)
  xv <- rep_len(theta[[2]], n)
  out <- -sv * ls
  nz <- !is.na(xv) & xv != 0
  out[nz] <- sv[nz] * expm1(-xv[nz] * ls[nz]) / xv[nz]
  out
}

#' @title Generalized Pareto Random Generation
#' @name distrib_rng.GPDDistrib
#'
#' @description
#' Draws by the inverse transform, `distrib_quantile(distrib, runif(n),
#' theta)`. The quantile function of this family is elementary, so the
#' transform is exact and costs one uniform per draw. The base class draws
#' the same way for a family with its own quantile method.
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
#' @seealso [distrib_quantile.GPDDistrib()], through which it inverts the
#'   distribution function, [distrib_pdf.GPDDistrib()] for the density that the
#'   draws follow, and [distrib_rng()] for the generic.
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
#' \eqn{z^2/2 - z}, so the kernel evaluates it as
#' \eqn{-z/t - z^2\phi'(u)}, \eqn{u = \xi z}, with
#' \eqn{\phi(u) = \log(1+u)/u} computed as [distrib_deriv3.GPDDistrib()]
#' describes; \eqn{\xi = 0} is an ordinary point of that form. \eqn{t} is
#' formed from the exact product \eqn{\xi y}.
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
#' # The analytic score is compared with a central difference of the
#' # log-density.
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
#' The pure-\eqn{\xi} component has the same removable singularity as the
#' score, to one order higher, and is computed from the second derivative of
#' the analytic function \eqn{\Lambda(u) = \log(1+u)/u} (written
#' \eqn{\phi} on the pages of the score and of the higher orders). The two
#' terms that individually diverge as \eqn{\xi \to 0} cancel, and the limit
#' is an ordinary point of the formula.
#'
#' This is the observed curvature at the data. The expected one is closed
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
#' # The analytic Hessian is compared with a central difference of the score.
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
#' Every entry is free of the data, so `approx` is ignored and
#' [expected_hessian_exact()] returns `TRUE`.
#'
#' @details
#' # Where it stops existing
#'
#' At \eqn{\xi \le -1/2} the information does not exist and every component
#' is `NA`. The condition is exactly that the integrand be integrable: written
#' on the probability scale, the second derivative grows like
#' \eqn{(1-u)^{-2|\xi|}} near the upper endpoint, which is integrable if and
#' only if \eqn{|\xi| < 1/2}. Below that point the classical asymptotics of the
#' maximum likelihood estimator do not hold either, so the `NA` means that the
#' quantity does not exist.
#'
#' Approaching the boundary the information diverges: the \eqn{\xi}
#' component \eqn{-2/\{(1+2\xi)(1+\xi)\}}, which does not depend on
#' \eqn{\sigma}, equals \eqn{-2} at \eqn{\xi = 0}, about \eqn{-7.14} at
#' \eqn{-0.3} and about \eqn{-196.1} at \eqn{-0.49}.
#'
#' # Checking it
#'
#' A Monte Carlo average of the observed Hessian converges slowly here,
#' because the second derivative is unbounded at the upper endpoint when
#' \eqn{\xi < 0}. If the integral is taken on the probability scale instead,
#' \eqn{E[h] = \int_0^1 h(Q(u))\,du}, the endpoint becomes an ordinary point
#' and the result agrees to \eqn{10^{-11}}.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector. Its values do not enter the result, which is an
#'   expectation; only its length does, through recycling.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx An object inheriting from [ExpectApprox]. Ignored: the
#'   expectation is closed form. It is accepted so that the signature matches
#'   the generic's. See [expected_derivative_methods].
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
#' # The result is compared with Smith's formula written out.
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
#' # The approx argument is ignored because the expectation is closed form.
#' identical(distrib_expected_hessian(d, 0, th),
#'           distrib_expected_hessian(d, 0, th, approx = expect_mc(nsim = 10)))
S7::method(distrib_expected_hessian, GPDDistrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"),
                                                              approx = expect_opg(), ..., threads = 1L) {
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
#' polynomial tail, and at \eqn{\xi < 0} it has a finite upper endpoint at
#' \eqn{-\sigma/\xi}. The last case makes the family non-regular, so the
#' expected information exists only for \eqn{\xi > -1/2}.
#'
#' @param link_sigma A `linkfunctions7` link object for the scale
#'   \eqn{\sigma}, which must be strictly positive. Defaults to
#'   [linkfunctions7::log_link()].
#' @param link_xi A link object for the shape \eqn{\xi}, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()], because the
#'   shape is free to take either sign.
#'
#' @details
#' # Density and distribution function
#'
#' \deqn{f(y) = \dfrac{1}{\sigma}\left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi-1},
#'       \qquad
#'       F(q) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi}.}
#' At \eqn{\xi = 0} both reduce to the exponential. The density, the
#' distribution function and the derivatives are evaluated through the
#' analytic function \eqn{\log(1+u)/u} of \eqn{u = \xi y/\sigma}, so the shape
#' may pass through zero during a fit. The quantile function is elementary,
#' \eqn{Q(p) = \sigma\{(1-p)^{-\xi}-1\}/\xi}, so the generator is an exact
#' inverse transform.
#'
#' # Parametrization
#'
#' The family is parametrized by the scale \eqn{\sigma} and the shape
#' \eqn{\xi}. The mean is \eqn{\sigma/(1-\xi)} and exists only for
#' \eqn{\xi < 1}, so a mean parametrization would not describe the
#' heavy-tailed members of the family. For the same reason [mv_sigma()] and
#' [variance()] are distinct for the multivariate \eqn{t}.
#'
#' # The support depends on the parameters
#'
#' For \eqn{\xi < 0} the support is \eqn{[0, -\sigma/\xi]}, so it depends on
#' the parameters. This removes the automatic
#' license to differentiate under the integral sign, on which the Bartlett
#' identities rest. The first two results below remain valid and the third
#' fails:
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
#' is zero beyond the endpoint and the distribution function is one, and a
#' caller that needs the upper limit of the support should use
#' [gpd_endpoint()] in place of `bounds`.
#'
#' # Moments
#'
#' The mean is \eqn{\sigma/(1-\xi)} for \eqn{\xi < 1} and the variance
#' \eqn{\sigma^2/\{(1-\xi)^2(1-2\xi)\}} for \eqn{\xi < 1/2}. The skewness
#' requires \eqn{\xi < 1/3} and the kurtosis \eqn{\xi < 1/4}. Each of the four
#' moments is returned as `Inf` at and above its threshold, so a fitted
#' object can report a mean where the variance is infinite.
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
