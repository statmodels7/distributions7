#' @include distrib.R generics.R
NULL

#' @title Poisson-Inverse Gaussian Distribution Class
#' @name Pig1Distrib
#'
#' @description
#' The S7 class of the Poisson-inverse Gaussian on \eqn{\{0, 1, 2, \dots\}} in
#' its mean-dispersion parametrization, gamlss's `PIG`: the mean is \eqn{\mu}
#' and the variance \eqn{\mu + \sigma\mu^2}. The family is a Poisson mixed
#' over an inverse Gaussian rate, an overdispersed count model with a heavier
#' tail than the negative binomial at the same variance.
#'
#' Build one with [pig1_distrib()], which supplies the two link functions. This
#' page documents the raw S7 constructor, which validates none of the
#' relationships between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Pig1Distrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. For an object built by
#'   [pig1_distrib()] the properties hold `"poisson-inverse gaussian"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "sigma")`, the interpretations
#'   `c(mu = "mean", sigma = "dispersion")`, `2`, and two domains
#'   \eqn{(0, \infty)}.
#'
#' @section Methods:
#' Registered in this file, the middle five reading one compiled kernel:
#'   [`distrib_pdf()`][distrib_pdf.Pig1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Pig1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Pig1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Pig1Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Pig1Distrib],
#'   [`distrib_rng()`][distrib_rng.Pig1Distrib].
#'
#' Registered elsewhere: the four moments in `moments.R` and the data-based
#' starting value in `starting_values.R`.
#'
#' The **distribution function** and the **quantile** come from
#' [discrete_distrib()], where both are exact sums over the support. The
#' **expected information** has no closed form and goes
#' through [expected_derivative_methods()].
#'
#' @seealso [pig1_distrib()] to build one;
#'   [pig2_distrib()] for the parametrization whose two parameters are
#'   orthogonal; [negbin2_distrib()] for the other overdispersed count family;
#'   [pig_hd_block()] for the kernel all five derivative methods read.
#'
#' @examples
#' d <- pig1_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # The mean is a parameter and the variance follows from both.
#' th <- list(mu = 3, sigma = 0.8)
#' c(mean = mean(d, th), variance = variance(d, th),
#'   formula = 3 + 0.8 * 3^2)
#'
#' # The mass sums to one over the support.
#' sum(distrib_pdf(d, 0:300, th))
Pig1Distrib <- S7::new_class("Pig1Distrib", parent = discrete_distrib)

#' @title Log-Likelihood Derivatives of the Poisson-Inverse Gaussian
#'
#' @description
#' Evaluates the log-likelihood and its fourteen partial derivatives to fourth
#' order at once, through one compiled kernel, and returns the block of
#' components an order asks for. Both parametrizations use it, with a different
#' `kernel` argument; every derivative method of both families is one call of
#' this with a different `cols`.
#'
#' @details
#' # The closed form the kernel evaluates
#'
#' With \eqn{c = 1 + 2\sigma\mu} and \eqn{\alpha = \sqrt{c}/\sigma}, the
#' half-integer order collapses the Bessel function to a finite sum and the
#' log-likelihood to
#' \deqn{\ell(y) = y\log\mu - \tfrac{y}{2}\log c + \tfrac{1}{\sigma}
#'   + \psi(\alpha) - \log y!,\qquad
#'   \psi(\alpha) = -\alpha + \log S_y(\alpha),}
#' \deqn{S_y(\alpha) = \sum_{k=0}^{y-1}
#'   \dfrac{\Gamma(y+k)}{\Gamma(k+1)\Gamma(y-k)}\,(2\alpha)^{-k},
#'   \qquad S_0 = 1.}
#' \eqn{S_y} sums \eqn{y} positive terms on the log scale, so nothing cancels.
#' The derivatives of \eqn{\psi} in \eqn{\alpha} are the weighted
#' rising-factorial moments of \eqn{k} under those terms; everything else is
#' elementary.
#'
#' # Which kernel runs
#'
#' `pig1_hd_cpp` and `pig2_hd_cpp` are **explicit closed-form kernels**, every
#' partial written out by hand, and they are what the package methods run. A
#' second pair, `pig1_hd_jet_cpp` and `pig2_hd_jet_cpp`, carries a bivariate
#' jet truncated at total order four through the same expression; it shares no
#' algebra with the explicit route, so the tests compare the two with no
#' tolerance to hide behind, and it is not used in production. Measured over
#' \eqn{2\times10^5} observations, the explicit kernel takes 0.24 seconds
#' against the jet's 1.30, and the two agree to \eqn{5\times10^{-14}}.
#'
#' # Rows that are not admissible
#'
#' A `y` that is negative, non-integer or non-finite never reaches the kernel
#' and gets a row of `NaN`. [distrib_pdf.Pig1Distrib()] turns that into a
#' log-probability of `-Inf`, the right answer for a point off the support.
#'
#' @param y A numeric vector of observations. Recycled with `theta` to the
#'   common length.
#' @param theta The parameter list, already aligned and read positionally: the
#'   mean first and the dispersion or \eqn{\alpha} second.
#' @param cols A named character vector selecting the kernel columns wanted.
#'   The names become the names of the result and the values are among
#'   `"l"`, `"d10"`, `"d01"`, `"d20"`, `"d11"`, `"d02"`, `"d30"`, `"d21"`,
#'   `"d12"`, `"d03"`, `"d40"`, `"d31"`, `"d22"`, `"d13"`, `"d04"`, where
#'   `dij` is \eqn{\partial^{i+j}\ell/\partial\theta_1^i\partial\theta_2^j}.
#' @param kernel The compiled kernel: `pig1_hd_cpp` for the mean-dispersion
#'   parametrization or `pig2_hd_cpp` for the orthogonal one.
#' @param threads A single positive integer, how many threads that kernel may
#'   use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of numeric vectors, one per entry of `cols`, each of
#'   the recycled length, named by `names(cols)`.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\sigma} the dispersion, \eqn{\alpha} the Bessel
#' argument, \eqn{K_\nu} the modified Bessel function of the second kind, and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [pig1_distrib()] and [pig2_distrib()] for the two families, and
#'   [distrib_gradient.Pig1Distrib()] for a method that calls this.
#'
#' @examples
#' # The score, taken directly.
#' distributions7:::pig_hd_block(0:4, list(mu = 3, sigma = 0.8),
#'                               c(mu = "d10", sigma = "d01"),
#'                               distributions7:::pig1_hd_cpp)
#'
#' # ...which is what the method returns.
#' all.equal(distributions7:::pig_hd_block(0:4, list(mu = 3, sigma = 0.8),
#'                                         c(mu = "d10", sigma = "d01"),
#'                                         distributions7:::pig1_hd_cpp),
#'           distrib_gradient(pig1_distrib(), 0:4,
#'                            list(mu = 3, sigma = 0.8)))
#'
#' # The explicit kernel and the jet twin agree, sharing no algebra.
#' y <- as.numeric(0:6)
#' max(abs(distributions7:::pig1_hd_cpp(y, rep(3, 7), rep(0.8, 7), 1L) -
#'         distributions7:::pig1_hd_jet_cpp(y, rep(3, 7), rep(0.8, 7))))
#'
#' # A point off the support gives NaN here and -Inf in the density.
#' distributions7:::pig_hd_block(c(-1, 1.5, 2), list(mu = 3, sigma = 0.8),
#'                               c(l = "l"), distributions7:::pig1_hd_cpp)
#'
#' @keywords internal
pig_hd_block <- function(y, theta, cols, kernel, threads = 1L) {
  n <- max(length(y), length(theta[[1]]), length(theta[[2]]))
  y <- rep_len(y, n)
  p1 <- rep_len(theta[[1]], n)
  p2 <- rep_len(theta[[2]], n)
  ok <- is.finite(y) & y >= 0 & y == floor(y)
  m <- matrix(NaN, n, 15)
  if (any(ok)) m[ok, ] <- kernel(y[ok], p1[ok], p2[ok], threads)
  colnames(m) <- c("l", "d10", "d01", "d20", "d11", "d02",
                   "d30", "d21", "d12", "d03",
                   "d40", "d31", "d22", "d13", "d04")
  stats::setNames(lapply(cols, function(cc) m[, cc]), names(cols))
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Poisson-Inverse Gaussian Probability Mass Function
#' @name distrib_pdf.Pig1Distrib
#'
#' @description
#' Computes the Poisson-inverse Gaussian mass. With \eqn{c = 1 + 2\sigma\mu}
#' and \eqn{\alpha = \sqrt{c}/\sigma},
#' \deqn{P(Y = y) = \sqrt{\dfrac{2\alpha}{\pi}}\,
#'   \dfrac{\mu^y e^{1/\sigma}}{(\alpha\sigma)^y\, y!}\, K_{y-1/2}(\alpha).}
#'
#' The value is not computed that way. At half-integer order the Bessel
#' function is a finite sum, the prefactors cancel, and what the compiled
#' kernel evaluates is
#' \eqn{\ell(y) = y\log\mu - (y/2)\log c + (1-\sqrt c)/\sigma
#'   + \log S_y(\alpha) - \log y!} with \eqn{S_y} a sum of \eqn{y} positive
#' terms taken on the log scale. Nothing cancels and no Bessel routine is
#' called.
#'
#' @param distrib A `Pig1Distrib` object, from [pig1_distrib()].
#' @param y A numeric vector of counts. A negative, non-integer or non-finite
#'   value is off the support and gives a probability of 0, or `-Inf` with
#'   `log = TRUE`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
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
#' \eqn{\mu} is the mean, \eqn{\sigma} the dispersion, \eqn{K_\nu} the modified
#' Bessel function of the second kind, and \eqn{S_y} the finite sum defined in
#' [pig_hd_block()].
#'
#' @seealso [distrib_gradient.Pig1Distrib()] for the score,
#'   [pig_hd_block()] for the kernel, [pig2_distrib()] for the same law in
#'   orthogonal coordinates, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- pig1_distrib()
#' y <- 0:6
#' th <- list(mu = 3, sigma = 0.8)
#'
#' # Against the Bessel formula written out, which the kernel does not use.
#' al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
#' all.equal(distrib_pdf(d, y, th),
#'           sqrt(2 * al / pi) * 3^y * exp(1 / 0.8) /
#'             ((al * 0.8)^y * factorial(y)) * besselK(al, y - 0.5))
#'
#' # The mass sums to one.
#' sum(distrib_pdf(d, 0:300, th))
#'
#' # The upper tail is heavier than a negative binomial's at the same
#' # variance, which is the reason to prefer this family.
#' nb <- negbin2_distrib()
#' rbind(pig = distrib_pdf(d, c(20, 40, 60), th),
#'       negbin = distrib_pdf(nb, c(20, 40, 60), list(mu = 3, theta = 1 / 0.8)))
#'
#' # Off the support.
#' distrib_pdf(d, c(-1, 1.5), th)
S7::method(distrib_pdf, Pig1Distrib) <- function(distrib, y, theta, log = FALSE, ..., threads = 1L) {
  out <- pig_hd_block(y, theta, c(l = "l"), pig1_hd_cpp, threads)$l
  out[is.nan(out)] <- -Inf
  if (log) out else exp(out)
}

#' @title Poisson-Inverse Gaussian Score
#' @name distrib_gradient.Pig1Distrib
#'
#' @description
#' Returns the exact first derivatives of the log-mass in \eqn{(\mu, \sigma)},
#' read off columns `d10` and `d01` of the compiled fourth-order kernel of
#' [pig_hd_block()]. Nothing is differenced: the kernel writes every partial
#' out in closed form from the finite Bessel sum.
#'
#' @param distrib A `Pig1Distrib` object, from [pig1_distrib()].
#' @param y A numeric vector of counts. A value off the support gives `NaN`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`. The result does not depend on the count.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each of the
#'   length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu} the mean and
#' \eqn{\sigma} the dispersion.
#'
#' @seealso [pig_hd_block()] for the kernel and the closed form it evaluates,
#'   [distrib_hessian.Pig1Distrib()] for the second derivatives, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- pig1_distrib()
#' y <- 0:6
#' th <- list(mu = 3, sigma = 0.8)
#' g <- distrib_gradient(d, y, th)
#'
#' # Against numerical differentiation of the log-likelihood.
#' f <- function(p) sum(distrib_pdf(d, y, list(mu = p[1], sigma = p[2]),
#'                                  log = TRUE))
#' rbind(analytic = vapply(g, sum, 0),
#'       numeric = numDeriv::grad(f, c(3, 0.8)))
#'
#' # The mean and the dispersion are not orthogonal here: the mixed entry of
#' # the expected information is far from zero. pig2_distrib() removes that.
#' sum(distrib_expected_hessian(d, 0:200, th, approx = "bartlett")$mu_sigma)
S7::method(distrib_gradient, Pig1Distrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ..., threads = 1L) {
  pig_hd_block(y, theta, c(mu = "d10", sigma = "d01"), pig1_hd_cpp, threads)
}

#' @title Poisson-Inverse Gaussian Observed Hessian
#' @name distrib_hessian.Pig1Distrib
#'
#' @description
#' Returns the exact second derivatives of the log-mass in
#' \eqn{(\mu, \sigma)}, read off columns `d20`, `d02` and `d11` of the
#' compiled fourth-order kernel of [pig_hd_block()].
#'
#' This is the **observed** curvature at the data. The expected information has
#' no closed form for this family and comes from
#' [expected_derivative_methods()], whose default here is the exact sum over
#' the support; see [distrib_expected_hessian()].
#'
#' @param distrib A `Pig1Distrib` object, from [pig1_distrib()].
#' @param y A numeric vector of counts. A value off the support gives `NaN`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors in [hess_names()]'s order:
#'   `mu_mu`, `sigma_sigma`, `mu_sigma`.
#'
#' @seealso [distrib_gradient.Pig1Distrib()] for the order below,
#'   [distrib_deriv3.Pig1Distrib()] for the order above, [pig_hd_block()] for
#'   the kernel, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- pig1_distrib()
#' y <- 0:6
#' th <- list(mu = 3, sigma = 0.8)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # Against a central difference of the score.
#' eps <- 1e-5
#' rbind(analytic = h$mu_sigma,
#'       numeric = (distrib_gradient(d, y, list(mu = 3, sigma = 0.8 + eps))$mu -
#'                  distrib_gradient(d, y, list(mu = 3, sigma = 0.8 - eps))$mu) /
#'                 (2 * eps))
#'
#' # The curvature in the mean is not of one sign: a count of zero and a
#' # count in the tail pull it in opposite directions.
#' distrib_hessian(d, c(0, 3, 30), th)$mu_mu
S7::method(distrib_hessian, Pig1Distrib) <- function(distrib, y, theta,
                                                     scale = c("parameter", "link"), ..., threads = 1L) {
  pig_hd_block(y, theta,
               c(mu_mu = "d20", sigma_sigma = "d02", mu_sigma = "d11"),
               pig1_hd_cpp, threads)
}

#' @title Poisson-Inverse Gaussian Third Derivatives
#' @name distrib_deriv3.Pig1Distrib
#'
#' @description
#' Returns the exact third derivatives of the log-mass in \eqn{(\mu, \sigma)},
#' read off columns `d30`, `d21`, `d12` and `d03` of the compiled fourth-order
#' kernel of [pig_hd_block()]. The kernel computes all four orders in one pass,
#' so this order costs no more than the score does.
#'
#' With `expected = TRUE` the value is an expectation instead, and there it is
#' **not** closed form: the call routes to `expected_derivative()`, so `approx`
#' and `nsim` are read.
#'
#' @param distrib A `Pig1Distrib` object, from [pig1_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` its values are
#'   the support points the expectation is summed over.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned from the compiled kernel.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, the
#'   strategy for the expectation. Read only when `expected = TRUE`; for a
#'   discrete family `"integrate"` is an exact sum over the support, not a
#'   quadrature.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the compiled
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors: `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_sigma_sigma` and `sigma_sigma_sigma`.
#'
#' @seealso [distrib_hessian.Pig1Distrib()] for the order below,
#'   [distrib_deriv4.Pig1Distrib()] for the order above, [pig_hd_block()] for
#'   the kernel, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- pig1_distrib()
#' y <- 0:6
#' th <- list(mu = 3, sigma = 0.8)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Against a central difference of the analytic Hessian.
#' eps <- 1e-5
#' rbind(analytic = d3$mu_mu_sigma,
#'       numeric = (distrib_hessian(d, y, list(mu = 3, sigma = 0.8 + eps))$mu_mu -
#'                  distrib_hessian(d, y, list(mu = 3, sigma = 0.8 - eps))$mu_mu) /
#'                 (2 * eps))
#'
#' # The expected version is not closed form and reads the strategy.
#' distrib_deriv3(d, 0:200, th, expected = TRUE)$mu_mu_mu[1]
S7::method(distrib_deriv3, Pig1Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ..., threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig_hd_block(y, theta,
               c(mu_mu_mu = "d30", mu_mu_sigma = "d21",
                 mu_sigma_sigma = "d12", sigma_sigma_sigma = "d03"),
               pig1_hd_cpp, threads)
}

#' @title Poisson-Inverse Gaussian Fourth Derivatives
#' @name distrib_deriv4.Pig1Distrib
#'
#' @description
#' Returns the exact fourth derivatives of the log-mass in
#' \eqn{(\mu, \sigma)}, read off the last five columns of the compiled kernel
#' of [pig_hd_block()]. The kernel is a fourth-order one throughout, so this is
#' the order it was written for and it costs the same as the score.
#'
#' With `expected = TRUE` the value is an expectation and is not closed form,
#' as at third order.
#'
#' @param distrib A `Pig1Distrib` object, from [pig1_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` its values are
#'   the support points the expectation is summed over.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
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
#'   `mu_mu_mu_sigma`, `mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
#'   `sigma_sigma_sigma_sigma`.
#'
#' @seealso [distrib_deriv3.Pig1Distrib()] for the order below,
#'   [pig_hd_block()] for the kernel, and [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- pig1_distrib()
#' y <- 0:6
#' th <- list(mu = 3, sigma = 0.8)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Against a central difference of the third order.
#' eps <- 1e-5
#' rbind(analytic = d4$mu_mu_mu_sigma,
#'       numeric = (distrib_deriv3(d, y, list(mu = 3, sigma = 0.8 + eps))$mu_mu_mu -
#'                  distrib_deriv3(d, y, list(mu = 3, sigma = 0.8 - eps))$mu_mu_mu) /
#'                 (2 * eps))
#'
#' # All four orders come from one pass of the kernel, so the fourth costs
#' # what the first does.
#' n <- 2e4
#' set.seed(62)
#' x <- distrib_rng(d, n, th)
#' rbind(score = system.time(distrib_gradient(d, x, th))[["elapsed"]],
#'       fourth = system.time(distrib_deriv4(d, x, th))[["elapsed"]])
S7::method(distrib_deriv4, Pig1Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ..., threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig_hd_block(y, theta,
               c(mu_mu_mu_mu = "d40", mu_mu_mu_sigma = "d31",
                 mu_mu_sigma_sigma = "d22", mu_sigma_sigma_sigma = "d13",
                 sigma_sigma_sigma_sigma = "d04"),
               pig1_hd_cpp, threads)
}

#' @title Poisson-Inverse Gaussian Random Generation
#' @name distrib_rng.Pig1Distrib
#'
#' @description
#' Draws from the family's own mixture representation, exactly: \eqn{\lambda}
#' from the inverse Gaussian with mean \eqn{\mu} and shape \eqn{\mu/\sigma},
#' then \eqn{Y \mid \lambda} from the Poisson with that rate. The inverse
#' Gaussian's variance is then \eqn{\sigma\mu^2}, the value the mixing
#' construction requires, and the marginal variance comes out as
#' \eqn{\mu + \sigma\mu^2}.
#'
#' The inverse Gaussian draw comes from `statmod::rinvgauss`, so the whole
#' generator is two vectorized calls and no rejection step.
#'
#' @param distrib A `Pig1Distrib` object, from [pig1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`; a component of length 1 is recycled,
#'   so a parameter varying by observation draws one value per observation from
#'   its own member of the family.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An integer-valued numeric vector of `n` draws.
#'
#' @seealso [distrib_pdf.Pig1Distrib()] for the mass the draws follow,
#'   [distrib_rng.Pig2Distrib()] for the same sampler in orthogonal
#'   coordinates, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- pig1_distrib()
#' th <- list(mu = 3, sigma = 0.8)
#'
#' set.seed(61)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # Two moments against their closed forms.
#' rbind(sample = c(mean(x), var(x)),
#'       theory = c(mean(d, th), variance(d, th)))
#'
#' # The empirical mass against the exact one, over the head of the support.
#' rbind(sample = as.numeric(table(factor(x, levels = 0:6))) / 2e5,
#'       exact = distrib_pdf(d, 0:6, th))
S7::method(distrib_rng, Pig1Distrib) <- function(distrib, n, theta, ...) {
  lam <- statmod::rinvgauss(n, mean = rep_len(theta[[1]], n),
                            shape = rep_len(theta[[1]], n) / rep_len(theta[[2]], n))
  stats::rpois(n, lam)
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Poisson-Inverse Gaussian Distribution Object
#'
#' @description
#' Builds a Poisson-inverse Gaussian distribution in its mean-dispersion
#' parametrization: mean \eqn{\mu} and variance \eqn{\mu + \sigma\mu^2}, which
#' is gamlss's `PIG`. The family is the Poisson mixed over an inverse Gaussian
#' rate, an overdispersed count model with a **heavier upper tail than the
#' negative binomial at the same variance**.
#'
#' @details
#' # The mass function, and how it is computed
#'
#' The mass carries the modified Bessel function \eqn{K_{y-1/2}}, which at
#' half-integer order is a finite sum. What the kernel evaluates is that sum,
#' on the log scale, after the prefactors have canceled; no Bessel routine is
#' called. The log-likelihood and its fourteen partial derivatives to fourth
#' order come out of one pass, all exact, so the fourth order costs what the
#' score does. See [pig_hd_block()] for the expression and for the jet twin the
#' tests hold it to.
#'
#' The **expected information** has no closed form and goes through the
#' summation strategies of [expected_derivative_methods()].
#'
#' # The tail
#'
#' Measured at \eqn{\mu = 3}, \eqn{\sigma = 0.8}, against a negative binomial
#' matched on the variance (\eqn{\theta = 1/\sigma}, both at 10.2): the mass at
#' \eqn{y = 20} is \eqn{6.8\times10^{-4}} against \eqn{4.8\times10^{-4}}, at 40
#' it is \eqn{5.8\times10^{-6}} against \eqn{5.4\times10^{-7}}, and at 60
#' \eqn{7.2\times10^{-8}} against \eqn{5.6\times10^{-10}}. The gap widens with
#' the count, and that widening is what a heavier tail means here.
#'
#' # Which parametrization to use
#'
#' In this one \eqn{\mu} and \eqn{\sigma} are **not** orthogonal: measured at
#' the same setting, the mixed entry of the expected information summed over
#' the support is 7.39. [pig2_distrib()] replaces \eqn{\sigma} by the Bessel
#' argument \eqn{\alpha} and makes that entry zero, at the cost of a second
#' parameter with no moment reading. Use this one to model the dispersion
#' directly and that one to fit the two parameters independently.
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (0, \infty)}
#' - \eqn{\sigma \in (0, \infty)}
#'
#' @section The distribution:
#' \deqn{P(Y=y) = \sqrt{\frac{2\alpha}{\pi}}\,\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\,y!}\,K_{y-1/2}(\alpha), \qquad \alpha = \sqrt{\frac{1}{\sigma^{2}} + \frac{2\mu}{\sigma}}}
#' on \eqn{y \in \{0, 1, \dots\}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \mu + \sigma\mu^{2}}
#'
#' @param link_mu A `linkfunctions7` link object for the mean \eqn{\mu}, which
#'   must be strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_sigma A link object for the dispersion \eqn{\sigma}, also
#'   strictly positive. Defaults to [linkfunctions7::log_link()].
#'
#' @return An S7 object of class [Pig1Distrib], inheriting from
#'   `discrete_distrib`. Its `params` are `c("mu", "sigma")`, its `bounds`
#'   `c(0, Inf)`, and its `link_params` the two links given here.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. *Applied Statistics* 54(3), 507--554.
#'
#' Dean, C., Lawless, J. F., and Willmot, G. E. (1989). A mixed
#' Poisson-inverse-Gaussian regression model. *Canadian Journal of
#' Statistics* 17(2), 171--181.
#'
#' @seealso [pig2_distrib()] for the orthogonal parametrization,
#'   [negbin2_distrib()] for the other overdispersed count family,
#'   [poisson_distrib()] for the limit at \eqn{\sigma \to 0},
#'   [pig_hd_block()] for the kernel, and [Pig1Distrib] for the class.
#'
#' @examples
#' d <- pig1_distrib()
#' th <- list(mu = 3, sigma = 0.8)
#'
#' distrib_pdf(d, 0:5, th)
#' c(mean = mean(d, th), variance = variance(d, th),
#'   skewness = skewness(d, th))
#'
#' # Heavier in the tail than a negative binomial of the same variance.
#' nb <- negbin2_distrib()
#' nbth <- list(mu = 3, theta = 1 / 0.8)
#' c(variance_pig = variance(d, th), variance_negbin = variance(nb, nbth))
#' rbind(pig = distrib_pdf(d, c(20, 40, 60), th),
#'       negbin = distrib_pdf(nb, c(20, 40, 60), nbth))
#'
#' # The two parameters are not orthogonal, which pig2_distrib() repairs.
#' c(pig1 = sum(distrib_expected_hessian(d, 0:200, th,
#'                                       approx = "bartlett")$mu_sigma),
#'   pig2 = sum(distrib_expected_hessian(pig2_distrib(), 0:200,
#'                                       list(mu = 3, alpha = 3.010399),
#'                                       approx = "bartlett")$mu_alpha))
#'
#' # A fit recovers both parameters.
#' set.seed(63)
#' x <- distrib_rng(d, 4000, th)
#' coef(fit_distrib(d, x))
#'
#' @export
pig1_distrib <- function(link_mu = log_link(), link_sigma = log_link()) {
  Pig1Distrib(
    distrib_name = "poisson-inverse gaussian", dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean", sigma = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), sigma = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma)
  )
}
