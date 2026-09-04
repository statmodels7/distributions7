#' @include distrib.R generics.R pig1_distrib.R
NULL

#' @title Poisson-Inverse Gaussian Distribution Class, Orthogonal Parametrization
#' @name Pig2Distrib
#'
#' @description
#' The S7 class of the Poisson-inverse Gaussian in the parametrization whose
#' two parameters are orthogonal, gamlss's `PIG2`. It is the same law as
#' [Pig1Distrib]; the mean stays \eqn{\mu} and the dispersion is replaced by
#' \eqn{\alpha}, the argument the mass function's Bessel function is evaluated
#' at.
#'
#' Orthogonal means the expected information is diagonal. Measured at
#' \eqn{\mu = 3}, its mixed entry summed over the support is
#' \eqn{-8.8\times10^{-15}} here against 7.39 in the mean-dispersion
#' parametrization.
#'
#' Build one with [pig2_distrib()], which supplies the two link functions. This
#' page documents the raw S7 constructor, which validates none of the
#' relationships between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Pig2Distrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. For an object built by
#'   [pig2_distrib()] the properties hold
#'   `"poisson-inverse gaussian (orthogonal)"`, `"univariate"`, `c(0, Inf)`,
#'   `c("mu", "alpha")`, the interpretations
#'   `c(mu = "mean", alpha = "bessel argument")`, `2`, and two domains
#'   \eqn{(0, \infty)}.
#'
#' @section Methods:
#' Registered in this file, the middle five reading one compiled kernel:
#'   [`distrib_pdf()`][distrib_pdf.Pig2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Pig2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Pig2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Pig2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Pig2Distrib],
#'   [`distrib_rng()`][distrib_rng.Pig2Distrib].
#'
#' Registered elsewhere: the four moments in `moments.R` and the data-based
#' starting value in `starting_values.R`. The distribution function and the
#' quantile come from [discrete_distrib()], both exact sums over the support.
#'
#' @seealso [pig2_distrib()] to build one;
#'   [pig1_distrib()] for the mean-dispersion parametrization;
#'   [pig2_sigma()] for the map between them;
#'   the compiled kernels for the kernel all five derivative methods read.
#'
#' @examples
#' d <- pig2_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # The same law as pig1, at the dispersion alpha implies.
#' al <- 3.010399
#' d1 <- pig1_distrib()
#' rbind(pig2 = distrib_pdf(d, 0:5, list(mu = 3, alpha = al)),
#'       pig1 = distrib_pdf(d1, 0:5,
#'                          list(mu = 3,
#'                               sigma = distributions7:::pig2_sigma(3, al))))
Pig2Distrib <- S7::new_class("Pig2Distrib", parent = discrete_distrib)

#' @title The Dispersion a Poisson-Inverse Gaussian Alpha Implies
#'
#' @description
#' Converts the orthogonal parametrization's \eqn{\alpha} into the dispersion
#' \eqn{\sigma} of [pig1_distrib()]:
#' \deqn{\sigma = \dfrac{\mu + \sqrt{\mu^2 + \alpha^2}}{\alpha^2},}
#' the positive root of \eqn{\alpha^2\sigma^2 - 2\mu\sigma - 1 = 0}.
#'
#' @details
#' The relation is the inverse of \eqn{\alpha = \sqrt{1 + 2\sigma\mu}/\sigma}.
#' Written this way there is no cancellation anywhere in the domain: both
#' \eqn{\mu} and \eqn{\sqrt{\mu^2 + \alpha^2}} are positive, so the numerator
#' is a sum of positive terms. The other root of the quadratic is negative and
#' is discarded.
#'
#' Only [distrib_rng.Pig2Distrib()] calls this. The derivative methods do not:
#' the compiled kernel takes \eqn{\alpha} as a variable of its own, so the map
#' is never differentiated.
#'
#' @param mu A numeric vector of means, strictly positive. Nothing is
#'   validated.
#' @param alpha A numeric vector of Bessel arguments, strictly positive.
#'
#' @return A numeric vector of dispersions, of the length of the recycled
#'   inputs, strictly positive.
#'
#' @seealso [pig2_distrib()] for the parametrization, [pig1_distrib()] for the
#'   one it maps onto, and [distrib_rng.Pig2Distrib()] for its only caller.
#'
#' @examples
#' # The map, and the round trip back through alpha = sqrt(1 + 2 sigma mu)/sigma.
#' sg <- distributions7:::pig2_sigma(3, 3.010399)
#' c(sigma = sg, alpha_back = sqrt(1 + 2 * sg * 3) / sg)
#'
#' # A large alpha is a small dispersion: the family tends to the Poisson.
#' distributions7:::pig2_sigma(3, c(1, 3, 30, 300))
#'
#' @keywords internal
pig2_sigma <- function(mu, alpha) (mu + sqrt(mu^2 + alpha^2)) / alpha^2

# --- S7 METHODS IMPLEMENTATION ---

#' @title Orthogonal Poisson-Inverse Gaussian Probability Mass Function
#' @name distrib_pdf.Pig2Distrib
#'
#' @description
#' Computes the same mass as [distrib_pdf.Pig1Distrib()], at the dispersion
#' \eqn{\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2} that \eqn{\alpha}
#' implies. The parameter \eqn{\alpha} of this parametrization **is** the
#' argument \eqn{\sqrt{1 + 2\sigma\mu}/\sigma} at which the Bessel function is
#' evaluated, so in this coordinate the mass function's own argument is a
#' parameter.
#'
#' The compiled kernel takes \eqn{\alpha} directly, so [pig2_sigma()] is not
#' called here and no map is composed.
#'
#' @param distrib A `Pig2Distrib` object, from [pig2_distrib()].
#' @param y A numeric vector of counts. A negative, non-integer or non-finite
#'   value is off the support and gives a probability of 0, or `-Inf` with
#'   `log = TRUE`.
#' @param theta A named list with components `mu` and `alpha`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param log Logical of length 1. When `TRUE` the log-probability is
#'   returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A numeric vector of probabilities, of the length of the recycled
#'   inputs.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\alpha} the Bessel argument, and \eqn{\sigma}
#' the dispersion of [pig1_distrib()], which \eqn{\alpha} determines.
#'
#' @seealso [distrib_pdf.Pig1Distrib()] for the same mass in mean and
#'   dispersion, [pig2_sigma()] for the map, and [distrib_pdf()] for the
#'   generic.
#'
#' @examples
#' d <- pig2_distrib()
#' d1 <- pig1_distrib()
#' y <- 0:6
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#'
#' # The same law, reached through the map.
#' all.equal(distrib_pdf(d, y, list(mu = 3, alpha = al)),
#'           distrib_pdf(d1, y, list(mu = 3, sigma = 0.8)))
#'
#' # The mass sums to one.
#' sum(distrib_pdf(d, 0:300, list(mu = 3, alpha = al)))
#'
#' # A large alpha is a small dispersion: the mass tends to the Poisson's.
#' rbind(pig2 = distrib_pdf(d, y, list(mu = 3, alpha = 1e4)),
#'       poisson = dpois(y, 3))
S7::method(distrib_pdf, Pig2Distrib) <- function(distrib, y, theta, log = FALSE, ..., threads = 1L) {
  out <- pig2_pdf_cpp(y, theta[[1]], theta[[2]], threads)
  if (log) out else exp(out)
}

#' @title Orthogonal Poisson-Inverse Gaussian Score
#' @name distrib_gradient.Pig2Distrib
#'
#' @description
#' Returns the exact first derivatives of the log-mass in \eqn{(\mu, \alpha)},
#' read off columns `d10` and `d01` of the compiled kernel of
#' `pig2_gradient_cpp`. The kernel takes \eqn{\alpha} as a variable of its own, so
#' these are derivatives in the orthogonal coordinates directly and no chain
#' rule through [pig2_sigma()] is composed.
#'
#' The property that gives this parametrization its name lives one order up:
#' the **expected** value of \eqn{\partial^2\ell/\partial\mu\partial\alpha} is
#' zero, so the two scores are uncorrelated and their maximum likelihood
#' estimates asymptotically independent.
#'
#' @param distrib A `Pig2Distrib` object, from [pig2_distrib()].
#' @param y A numeric vector of counts. A value off the support gives `NaN`.
#' @param theta A named list with components `mu` and `alpha`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `alpha`, each of the
#'   length of the recycled inputs.
#'
#' @seealso [distrib_hessian.Pig2Distrib()] for the second derivatives,
#'   [distrib_gradient.Pig1Distrib()] for the same quantity in mean and
#'   dispersion, `pig2_gradient_cpp` for the kernel, and [distrib_gradient()] for
#'   the generic.
#'
#' @examples
#' d <- pig2_distrib()
#' y <- 0:6
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' th <- list(mu = 3, alpha = al)
#' g <- distrib_gradient(d, y, th)
#'
#' # Against numerical differentiation of the log-likelihood.
#' f <- function(p) sum(distrib_pdf(d, y, list(mu = p[1], alpha = p[2]),
#'                                  log = TRUE))
#' rbind(analytic = vapply(g, sum, 0),
#'       numeric = numDeriv::grad(f, c(3, al)))
#'
#' # The two scores are uncorrelated under the model, which is what
#' # orthogonality means and what pig1 does not have.
#' c(pig2 = sum(distrib_expected_hessian(d, 0:200, th,
#'                                       approx = "bartlett")$mu_alpha),
#'   pig1 = sum(distrib_expected_hessian(pig1_distrib(), 0:200,
#'                                       list(mu = 3, sigma = 0.8),
#'                                       approx = "bartlett")$mu_sigma))
S7::method(distrib_gradient, Pig2Distrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ..., threads = 1L) {
  pig2_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Orthogonal Poisson-Inverse Gaussian Observed Hessian
#' @name distrib_hessian.Pig2Distrib
#'
#' @description
#' Returns the exact second derivatives of the log-mass in
#' \eqn{(\mu, \alpha)}, read off columns `d20`, `d02` and `d11` of the compiled
#' kernel `pig2_hessian_cpp`.
#'
#' The **observed** mixed entry is not zero at any single observation; what
#' vanishes is its expectation. Measured at \eqn{\mu = 3},
#' \eqn{\alpha = 3.0104}, the expectation summed over the support is
#' \eqn{-8.8\times10^{-15}} while the individual entries are of order one.
#'
#' @param distrib A `Pig2Distrib` object, from [pig2_distrib()].
#' @param y A numeric vector of counts. A value off the support gives `NaN`.
#' @param theta A named list with components `mu` and `alpha`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors in [hess_names()]'s order:
#'   `mu_mu`, `alpha_alpha`, `mu_alpha`.
#'
#' @seealso [distrib_gradient.Pig2Distrib()] for the order below,
#'   [distrib_deriv3.Pig2Distrib()] for the order above, `pig2_hessian_cpp` for
#'   the kernel, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- pig2_distrib()
#' y <- 0:6
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' th <- list(mu = 3, alpha = al)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # Against a central difference of the score.
#' eps <- 1e-5
#' rbind(analytic = h$mu_alpha,
#'       numeric = (distrib_gradient(d, y, list(mu = 3, alpha = al + eps))$mu -
#'                  distrib_gradient(d, y, list(mu = 3, alpha = al - eps))$mu) /
#'                 (2 * eps))
#'
#' # The mixed entry is not zero observation by observation; its expectation
#' # is.
#' c(observed = h$mu_alpha[1],
#'   expected = sum(distrib_expected_hessian(d, 0:200, th,
#'                                           approx = "bartlett")$mu_alpha))
S7::method(distrib_hessian, Pig2Distrib) <- function(distrib, y, theta,
                                                     scale = c("parameter", "link"), ..., threads = 1L) {
  pig2_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Orthogonal Poisson-Inverse Gaussian Third Derivatives
#' @name distrib_deriv3.Pig2Distrib
#'
#' @description
#' Returns the exact third derivatives of the log-mass in
#' \eqn{(\mu, \alpha)}, read off columns `d30`, `d21`, `d12` and `d03` of the
#' compiled kernel `pig2_deriv3_cpp`. All four orders come out
#' of one pass of the kernel, so this order costs what the score does.
#'
#' With `expected = TRUE` the value is an expectation instead, and there it is
#' **not** closed form: the call routes to `expected_derivative()`, so `approx`
#' and `nsim` are read.
#'
#' @param distrib A `Pig2Distrib` object, from [pig2_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` its values are
#'   the support points the expectation is summed over.
#' @param theta A named list with components `mu` and `alpha`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned from the compiled kernel.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`; for a discrete family `"integrate"` is an
#'   exact sum over the support.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors: `mu_mu_mu`, `mu_mu_alpha`,
#'   `mu_alpha_alpha` and `alpha_alpha_alpha`.
#'
#' @seealso [distrib_hessian.Pig2Distrib()] for the order below,
#'   [distrib_deriv4.Pig2Distrib()] for the order above, `pig2_deriv3_cpp` for
#'   the kernel, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- pig2_distrib()
#' y <- 0:6
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' th <- list(mu = 3, alpha = al)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Against a central difference of the analytic Hessian.
#' eps <- 1e-5
#' rbind(analytic = d3$mu_mu_alpha,
#'       numeric = (distrib_hessian(d, y, list(mu = 3, alpha = al + eps))$mu_mu -
#'                  distrib_hessian(d, y, list(mu = 3, alpha = al - eps))$mu_mu) /
#'                 (2 * eps))
S7::method(distrib_deriv3, Pig2Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ..., threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig2_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Orthogonal Poisson-Inverse Gaussian Fourth Derivatives
#' @name distrib_deriv4.Pig2Distrib
#'
#' @description
#' Returns the exact fourth derivatives of the log-mass in
#' \eqn{(\mu, \alpha)}, read off the last five columns of the compiled kernel
#' of `pig2_deriv4_cpp`. The kernel is asked for this order alone, so this is
#' the order it was written for.
#'
#' With `expected = TRUE` the value is an expectation and is not closed form,
#' as at third order.
#'
#' @param distrib A `Pig2Distrib` object, from [pig2_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` its values are
#'   the support points the expectation is summed over.
#' @param theta A named list with components `mu` and `alpha`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned from the compiled kernel.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors: `mu_mu_mu_mu`,
#'   `mu_mu_mu_alpha`, `mu_mu_alpha_alpha`, `mu_alpha_alpha_alpha` and
#'   `alpha_alpha_alpha_alpha`.
#'
#' @seealso [distrib_deriv3.Pig2Distrib()] for the order below,
#'   `pig2_deriv4_cpp` for the kernel, and [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- pig2_distrib()
#' y <- 0:6
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' th <- list(mu = 3, alpha = al)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Against a central difference of the third order.
#' eps <- 1e-5
#' rbind(analytic = d4$mu_mu_mu_alpha,
#'       numeric = (distrib_deriv3(d, y,
#'                    list(mu = 3, alpha = al + eps))$mu_mu_mu -
#'                  distrib_deriv3(d, y,
#'                    list(mu = 3, alpha = al - eps))$mu_mu_mu) / (2 * eps))
S7::method(distrib_deriv4, Pig2Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ..., threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig2_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Orthogonal Poisson-Inverse Gaussian Random Generation
#' @name distrib_rng.Pig2Distrib
#'
#' @description
#' Draws from the mixture representation of [distrib_rng.Pig1Distrib()], at the
#' dispersion [pig2_sigma()] implies: \eqn{\lambda} from the inverse Gaussian
#' with mean \eqn{\mu} and shape \eqn{\mu/\sigma}, then \eqn{Y \mid \lambda}
#' from the Poisson. This is the one method of the family that composes the
#' map; the derivative methods take \eqn{\alpha} directly.
#'
#' @param distrib A `Pig2Distrib` object, from [pig2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `alpha`, each a numeric
#'   vector of length 1 or of length `n`; a component of length 1 is recycled,
#'   so a parameter varying by observation draws one value per observation from
#'   its own member of the family.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An integer-valued numeric vector of `n` draws.
#'
#' @seealso [distrib_rng.Pig1Distrib()] for the same sampler in mean and
#'   dispersion, [pig2_sigma()] for the map it composes, and [distrib_rng()]
#'   for the generic.
#'
#' @examples
#' d <- pig2_distrib()
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' th <- list(mu = 3, alpha = al)
#'
#' set.seed(64)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # Two moments against their closed forms.
#' rbind(sample = c(mean(x), var(x)),
#'       theory = c(mean(d, th), variance(d, th)))
#'
#' # The empirical mass against the exact one, over the head of the support.
#' rbind(sample = as.numeric(table(factor(x, levels = 0:6))) / 2e5,
#'       exact = distrib_pdf(d, 0:6, th))
S7::method(distrib_rng, Pig2Distrib) <- function(distrib, n, theta, ...) {
  mu <- rep_len(theta[[1]], n)
  sg <- pig2_sigma(mu, rep_len(theta[[2]], n))
  lam <- statmod::rinvgauss(n, mean = mu, shape = mu / sg)
  stats::rpois(n, lam)
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Poisson-Inverse Gaussian Distribution Object, Orthogonal Parametrization
#'
#' @description
#' Builds a Poisson-inverse Gaussian distribution in the parametrization whose
#' two parameters are orthogonal, which is gamlss's `PIG2`. The mean stays
#' \eqn{\mu}; the second parameter \eqn{\alpha} is exactly the argument at
#' which the mass function's Bessel function is evaluated, related to the
#' dispersion of [pig1_distrib()] by
#' \eqn{\alpha = \sqrt{1 + 2\sigma\mu}/\sigma}.
#'
#' Orthogonal means the expected information is diagonal. Measured at
#' \eqn{\mu = 3}, its mixed entry summed over the support is
#' \eqn{-8.8\times10^{-15}} here against 7.39 in the mean-dispersion
#' parametrization.
#'
#' @details
#' # What orthogonality buys, and what it costs
#'
#' The maximum likelihood estimates of \eqn{\mu} and \eqn{\alpha} are
#' asymptotically independent, so a Fisher scoring step in one parameter does
#' not disturb the other and the two can be modeled with separate linear
#' predictors without their estimates fighting. The cost is that \eqn{\alpha}
#' has no moment reading of its own: a reader who wants to interpret the
#' overdispersion directly wants [pig1_distrib()]'s \eqn{\sigma}, which enters
#' the variance as \eqn{\mu + \sigma\mu^2}.
#'
#' The two are one law and the map between them is closed both ways:
#' \eqn{\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2} through
#' [pig2_sigma()], and \eqn{\alpha = \sqrt{1 + 2\sigma\mu}/\sigma} back. gamlss
#' states the first as \eqn{\alpha = 1/(\sqrt{\mu^2+\sigma_2^2} - \mu)} in its
#' own \eqn{\sigma_2}, which coincides with this \eqn{\alpha}.
#'
#' # How the derivatives are computed
#'
#' By the same compiled kernel as [pig1_distrib()], with \eqn{\alpha} as a
#' variable of its own, so the Bessel argument needs no chain rule at all.
#' Every partial to fourth order is exact and comes out of one pass; see
#' the compiled kernels. The expected information has no closed form and goes
#' through [expected_derivative_methods()].
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (0, \infty)}
#' - \eqn{\alpha \in (0, \infty)}
#'
#' A large \eqn{\alpha} is a small dispersion, and the family tends to the
#' Poisson: measured at \eqn{\alpha = 10^4} the mass agrees with `dpois` over
#' the head of the support.
#'
#' @section The distribution:
#' \deqn{P(Y=y) = \sqrt{\frac{2\alpha}{\pi}}\,\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\,y!}\,K_{y-1/2}(\alpha), \qquad \sigma = \frac{1}{\sqrt{\mu^{2}+\alpha^{2}} - \mu}}
#' on \eqn{y \in \{0, 1, \dots\}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \mu + \sigma\mu^{2}}
#'
#' @param link_mu A `linkfunctions7` link object for the mean \eqn{\mu}, which
#'   must be strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_alpha A link object for \eqn{\alpha}, also strictly positive.
#'   Defaults to [linkfunctions7::log_link()].
#'
#' @return An S7 object of class [Pig2Distrib], inheriting from
#'   `discrete_distrib`. Its `params` are `c("mu", "alpha")`, its `bounds`
#'   `c(0, Inf)`, and its `link_params` the two links given here.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. *Applied Statistics* 54(3), 507--554.
#'
#' Heller, G. Z., Couturier, D.-L., and Heritier, S. R. (2019). Beyond mean
#' modelling: bias due to misspecification of dispersion in Poisson-inverse
#' Gaussian regression. *Biometrical Journal* 61(2), 333--342.
#'
#' @seealso [pig1_distrib()] for the mean-dispersion parametrization,
#'   [pig2_sigma()] for the map, [negbin2_distrib()] for the other
#'   overdispersed count family, and [Pig2Distrib] for the class.
#'
#' @examples
#' d <- pig2_distrib()
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' th <- list(mu = 3, alpha = al)
#'
#' distrib_pdf(d, 0:5, th)
#' c(mean = mean(d, th), variance = variance(d, th))
#'
#' # The same law as pig1 at the dispersion alpha implies.
#' all.equal(distrib_pdf(d, 0:5, th),
#'           distrib_pdf(pig1_distrib(), 0:5, list(mu = 3, sigma = 0.8)))
#'
#' # The property the parametrization exists for.
#' c(pig2 = sum(distrib_expected_hessian(d, 0:200, th,
#'                                       approx = "bartlett")$mu_alpha),
#'   pig1 = sum(distrib_expected_hessian(pig1_distrib(), 0:200,
#'                                       list(mu = 3, sigma = 0.8),
#'                                       approx = "bartlett")$mu_sigma))
#'
#' # A large alpha is a small dispersion, and the family tends to the Poisson.
#' rbind(pig2 = distrib_pdf(d, 0:5, list(mu = 3, alpha = 1e4)),
#'       poisson = dpois(0:5, 3))
#'
#' # A fit recovers both parameters.
#' set.seed(65)
#' x <- distrib_rng(d, 4000, th)
#' coef(fit_distrib(d, x))
#'
#' @export
pig2_distrib <- function(link_mu = log_link(), link_alpha = log_link()) {
  Pig2Distrib(
    distrib_name = "poisson-inverse gaussian (orthogonal)",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "alpha"),
    params_interpretation = c(mu = "mean", alpha = "bessel argument"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), alpha = c(0, Inf)),
    link_params = list(mu = link_mu, alpha = link_alpha)
  )
}
