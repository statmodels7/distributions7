#' @include distrib.R generics.R
NULL

#' @title Inverse Gaussian Distribution Class, Mean and Shape
#' @name InvGauss2Distrib
#'
#' @description
#' This is the S7 class of the inverse Gaussian family on \eqn{(0, \infty)} in
#' its classical parametrization, the mean \eqn{\mu > 0} and the shape
#' \eqn{\lambda > 0}, so that \eqn{\operatorname{Var}(Y) = \mu^3/\lambda}. It
#' inherits from `continuous_distrib`, so every generic of the `distrib` class
#' applies to it; the methods listed below are registered on it in this file and
#' everything else comes from the parent.
#'
#' Build one with [invgauss2_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties; the validator of `distrib` checks that
#' they are consistent (the parameter count, a bound and a link for every
#' parameter, ordered support bounds).
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `InvGauss2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [invgauss2_distrib()] they hold `"invgauss2"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "lambda")`, the interpretations
#'   `c(mu = "mean", lambda = "shape")`, `2`, the domain \eqn{(0, \infty)} for
#'   both parameters, and the two links.
#'
#' @seealso [invgauss2_distrib()] to build one;
#'   [invgauss1_distrib()] for the same law in mean and dispersion
#'   \eqn{\phi = 1/\lambda};
#'   [distrib_pdf.InvGauss2Distrib()] and
#'   [distrib_gradient.InvGauss2Distrib()] for the closed forms of this class.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.InvGauss2Distrib],
#'   [`distrib_cdf()`][distrib_cdf.InvGauss2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.InvGauss2Distrib] and
#'   [`distrib_rng()`][distrib_rng.InvGauss2Distrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.InvGauss2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.InvGauss2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.InvGauss2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.InvGauss2Distrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.InvGauss2Distrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.InvGauss2Distrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.InvGauss2Distrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.InvGauss2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.InvGauss2Distrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.InvGauss2Distrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.InvGauss2Distrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.InvGauss2Distrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.InvGauss2Distrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.InvGauss2Distrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.InvGauss2Distrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.InvGauss2Distrib],
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.InvGauss2Distrib],
#'   [`distrib_deriv3_cdf()`][distrib_deriv3_cdf.InvGauss2Distrib] and
#'   [`distrib_deriv4_cdf()`][distrib_deriv4_cdf.InvGauss2Distrib];
#' - the moments: [`mean()`][mean.InvGauss2Distrib],
#'   [`variance()`][variance.InvGauss2Distrib],
#'   [`skewness()`][skewness.InvGauss2Distrib] and
#'   [`kurtosis()`][kurtosis.InvGauss2Distrib].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- invgauss2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # These properties let a consumer drive the family with no code specific to
#' # it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # lambda is a shape, so a larger value is a tighter law at a fixed mean.
#' vapply(c(1, 3, 30), function(l) variance(d, list(mu = 2, lambda = l)),
#'        numeric(1))
InvGauss2Distrib <- S7::new_class("InvGauss2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Inverse Gaussian Probability Density Function in Mean and Shape
#' @name distrib_pdf.InvGauss2Distrib
#' @description
#' Computes the inverse Gaussian density
#' \deqn{f(y; \mu, \lambda) = \sqrt{\dfrac{\lambda}{2\pi y^3}}
#'       \exp\left\{-\dfrac{\lambda(y-\mu)^2}{2\mu^2 y}\right\},
#'       \qquad y > 0,}
#' by calling [statmod::dinvgauss()] at `mean = mu` and
#' `dispersion = 1/lambda`, because that package writes the family in the
#' dispersion.
#' With `log = TRUE` the logarithm is formed inside `dinvgauss()` and stays
#' finite where the density itself underflows.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param y A numeric vector of observations. The support is
#'   \eqn{(0, \infty)}; a value at or below zero gives 0.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(lambda))`, one value per observation.
#'
#' @seealso [distrib_cdf.InvGauss2Distrib()] for the distribution function,
#'   [distrib_gradient.InvGauss2Distrib()] for the derivatives of the
#'   log-density, [distrib_pdf.InvGauss1Distrib()] for the same density in the
#'   dispersion, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' y <- c(1, 2, 3)
#' th <- list(mu = 2, lambda = 3)
#'
#' # The method is statmod::dinvgauss at dispersion 1/lambda.
#' all.equal(distrib_pdf(d, y, th),
#'           statmod::dinvgauss(y, mean = 2, dispersion = 1 / 3))
#'
#' # The density is the same as that of invgauss1 at phi = 1/lambda.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(invgauss1_distrib(), y, list(mu = 2, phi = 1 / 3)))
#'
#' # A parameter may vary by observation, one value each.
#' distrib_pdf(d, y, list(mu = c(1, 2, 4), lambda = 3))
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e4, th)
#' distrib_pdf(d, 1e4, th, log = TRUE)
S7::method(distrib_pdf, InvGauss2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  statmod::dinvgauss(y, mean = theta[[1]], dispersion = 1 / theta[[2]], log = log)
}

#' @title Inverse Gaussian Cumulative Distribution Function in Mean and Shape
#' @name distrib_cdf.InvGauss2Distrib
#' @description
#' Computes the inverse Gaussian distribution function, which is elementary in
#' the standard normal distribution function \eqn{\Phi}:
#' \deqn{F(q; \mu, \lambda) = \Phi\left\{\sqrt{\dfrac{\lambda}{q}}
#'         \left(\dfrac{q}{\mu} - 1\right)\right\}
#'       + e^{2\lambda/\mu}\,
#'       \Phi\left\{-\sqrt{\dfrac{\lambda}{q}}
#'         \left(\dfrac{q}{\mu} + 1\right)\right\},}
#' by calling [statmod::pinvgauss()] at `mean = mu` and
#' `dispersion = 1/lambda`. Both tails are available exactly, and
#' `log.p = TRUE` returns a logarithm that stays finite where the probability
#' itself underflows. The factor \eqn{e^{2\lambda/\mu}} overflows in double
#' precision when \eqn{2\lambda/\mu} exceeds about 709, so that function
#' evaluates the expression on the log scale.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param q A numeric vector of quantiles. A value at or below zero gives a
#'   lower-tail probability of 0.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(lambda))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.InvGauss2Distrib()] for the inverse,
#'   [distrib_pdf.InvGauss2Distrib()] for the density,
#'   [distrib_grad_cdf.InvGauss2Distrib()] for the derivatives of this function
#'   in the parameters, which are closed form here, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' th <- list(mu = 2, lambda = 3)
#'
#' # The method is statmod::pinvgauss at dispersion 1/lambda.
#' all.equal(distrib_cdf(d, c(1, 2, 3), th),
#'           statmod::pinvgauss(c(1, 2, 3), mean = 2, dispersion = 1 / 3))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 2, th) + distrib_cdf(d, 2, th, lower.tail = FALSE)
#'
#' # The law is right skewed, so most of the mass sits below the mean.
#' distrib_cdf(d, 2, th)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 1e4, th, lower.tail = FALSE)
#' distrib_cdf(d, 1e4, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, InvGauss2Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  statmod::pinvgauss(q, mean = theta[[1]], dispersion = 1 / theta[[2]],
                     lower.tail = lower.tail, log.p = log.p)
}

#' @title Inverse Gaussian Quantile Function in Mean and Shape
#' @name distrib_quantile.InvGauss2Distrib
#' @description
#' Computes the inverse Gaussian quantile function by calling
#' [statmod::qinvgauss()] at `mean = mu` and `dispersion = 1/lambda`. The
#' quantile function has no closed form: the distribution function is
#' elementary but cannot be inverted in elementary terms, so `qinvgauss()`
#' inverts it numerically, with the
#' far-tail correction of [invgauss_quantile()]. The
#' distribution function is strictly increasing on \eqn{(0, \infty)}, so the
#' round trip through [distrib_cdf.InvGauss2Distrib()] returns `p`.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A probability outside \eqn{[0, 1]} gives
#'   `NA` without a warning; `p = 0` gives 0 and `p = 1` gives `Inf`. With
#'   `log.p = TRUE` a positive value is passed to [statmod::qinvgauss()],
#'   which signals an error when `p` has more than one element.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles in \eqn{[0, \infty]}, of length
#'   `max(length(p), length(mu), length(lambda))`.
#'
#' @seealso [distrib_cdf.InvGauss2Distrib()], which this inverts;
#'   [distrib_rng.InvGauss2Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' th <- list(mu = 2, lambda = 3)
#'
#' # These quantiles bound a central 95 percent interval, which is strongly
#' # asymmetric about the mean of 2.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities that it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # The median falls below the mean.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, InvGauss2Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  invgauss_quantile(distrib, p, theta, theta[[1]], 1 / theta[[2]], lower.tail,
                    log.p)
}

#' @title Inverse Gaussian Random Number Generator in Mean and Shape
#' @name distrib_rng.InvGauss2Distrib
#' @description
#' Draws `n` independent inverse Gaussian variates by calling
#' [statmod::rinvgauss()] at `mean = mu` and `dispersion = 1/lambda`, so the
#' draws come from that package's generator and depend on `.Random.seed` in the
#' usual way. The inverse-transform fallback of the base class is not used.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. Both
#'   must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` strictly positive draws.
#'
#' @seealso [distrib_quantile.InvGauss2Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#'
#' # The method uses the same generator as statmod::rinvgauss at dispersion
#' # 1/lambda.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(mu = 2, lambda = 3))
#' set.seed(2)
#' identical(a, statmod::rinvgauss(3, mean = 2, dispersion = 1 / 3))
#'
#' # Both maximum likelihood estimates are closed form: the mean is the sample
#' # mean and the shape the reciprocal of mean(1/y) - 1/ybar.
#' set.seed(3)
#' z <- distrib_rng(d, 2e4, list(mu = 2, lambda = 3))
#' c(mu = mean(z), lambda = 1 / (mean(1 / z) - 1 / mean(z)))
S7::method(distrib_rng, InvGauss2Distrib) <- function(distrib, n, theta, ...) {
  statmod::rinvgauss(n, mean = theta[[1]], dispersion = 1 / theta[[2]])
}

#' @title Inverse Gaussian Score in Mean and Shape
#' @name distrib_gradient.InvGauss2Distrib
#' @description
#' Computes the first derivatives of the inverse Gaussian log-density with
#' respect to \eqn{\mu} and \eqn{\lambda}, one value per observation, in closed
#' form:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\lambda(y-\mu)}{\mu^3},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{2\lambda}
#'         - \dfrac{(y-\mu)^2}{2\mu^2 y}.}
#' The log-density is linear in the shape apart from
#' \eqn{\tfrac12\log\lambda}. Every derivative of this family is elementary for
#' that reason, and every expectation reduces to
#' \eqn{\mathbb{E}[Y] = \mu} and \eqn{\mathbb{E}[(Y-\mu)^2/Y] = \mu^2/\lambda}.
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' With `scale = "link"` the generic applies the chain rule for the links of the
#' family before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below a fixed number of observations, set per kernel, the kernel
#'   runs sequentially. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `lambda`, each of
#'   length `max(length(y), length(mu), length(lambda))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\lambda > 0} the shape, with
#' \eqn{\operatorname{Var}(Y) = \mu^3/\lambda}. Here \eqn{\lambda} names this
#' family's shape parameter throughout.
#'
#' @seealso [distrib_hessian.InvGauss2Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.InvGauss2Distrib()] for their expectation,
#'   [distrib_gradient.InvGauss1Distrib()] for the same score in the
#'   dispersion, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' y <- c(1, 2, 3)
#' th <- list(mu = 2, lambda = 3)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms are written out.
#' all.equal(g$mu, 3 * (y - 2) / 2^3)
#' all.equal(g$lambda, 1 / (2 * 3) - (y - 2)^2 / (2 * 2^2 * y))
#'
#' # The mean component vanishes at y = mu, whatever the shape.
#' distrib_gradient(d, 2, list(mu = 2, lambda = c(0.5, 3, 30)))$mu
#'
#' # Both estimating equations solve in closed form, so the summed score
#' # vanishes at the sample mean and 1/(mean(1/y) - 1/mean(y)).
#' set.seed(3)
#' z <- distrib_rng(d, 2000, th)
#' mle <- list(mu = mean(z), lambda = 1 / (mean(1 / z) - 1 / mean(z)))
#' vapply(distrib_gradient(d, z, mle), sum, numeric(1))
S7::method(distrib_gradient, InvGauss2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  invgauss2_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Observed Hessian in Mean and Shape
#' @name distrib_hessian.InvGauss2Distrib
#' @description
#' Computes the three distinct second derivatives of the inverse Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\lambda}, one value per
#' observation, in closed form:
#' \deqn{\ell^{(\mu\mu)} = \dfrac{\lambda(2\mu - 3y)}{\mu^4}, \qquad
#'       \ell^{(\mu\lambda)} = \dfrac{y - \mu}{\mu^3}, \qquad
#'       \ell^{(\lambda\lambda)} = -\dfrac{1}{2\lambda^2}.}
#'
#' The shape parametrization is better behaved than the dispersion one at
#' second order. The pure shape entry is a negative constant, free of the data,
#' where the corresponding entry of [distrib_hessian.InvGauss1Distrib()] can be
#' positive; only the curvature in \eqn{\mu} can turn positive here, and it
#' does so wherever \eqn{y < 2\mu/3}.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below a fixed number of observations, set per kernel, the kernel
#'   runs sequentially. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_lambda` and
#'   `lambda_lambda`, in that order, each of length
#'   `max(length(y), length(mu), length(lambda))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#' \eqn{\lambda} names this family's shape parameter.
#'
#' @seealso [distrib_gradient.InvGauss2Distrib()] for the score,
#'   [distrib_expected_hessian.InvGauss2Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.InvGauss2Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' y <- c(1, 2, 3)
#' th <- list(mu = 2, lambda = 3)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # The three closed forms are written out.
#' all.equal(h$mu_mu, 3 * (2 * 2 - 3 * y) / 2^4)
#' all.equal(h$mu_lambda, (y - 2) / 2^3)
#' all.equal(h$lambda_lambda, rep(-1 / (2 * 3^2), 3))
#'
#' # The curvature in mu is positive below 2 mu/3 and negative above.
#' c(at_1 = distrib_hessian(d, 1, th)$mu_mu,
#'   at_3 = distrib_hessian(d, 3, th)$mu_mu)
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 2 + eps, lambda = 3))$mu
#' dn <- distrib_gradient(d, y, list(mu = 2 - eps, lambda = 3))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
S7::method(distrib_hessian, InvGauss2Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  invgauss2_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Expected Hessian in Mean and Shape
#' @name distrib_expected_hessian.InvGauss2Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. Every second derivative is at
#' most linear in the response, so the expectations need only
#' \eqn{\mathbb{E}[Y] = \mu}:
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{\lambda}{\mu^3},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\lambda)}\right] = 0, \qquad
#'       \mathbb{E}\left[\ell^{(\lambda\lambda)}\right] =
#'         -\dfrac{1}{2\lambda^2}.}
#' The pure shape entry equals the observed value itself, because it is free of
#' the data.
#'
#' Both diagonal entries are negative at every parameter setting, so the
#' information is positive definite everywhere. The zero off-diagonal entry
#' shows that the mean and the shape are orthogonal, so the mean equation can
#' be fitted without bias when the shape is held at any value.
#'
#' Because the values do not depend on the data, `y` is read only for its
#' length.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL` or an object inheriting from [ExpectMethod]. The generic
#'   applies an object before dispatching on the family, so this method
#'   receives `NULL`, and it returns the closed form of the family. See
#'   [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below a fixed number of observations, set per kernel, the kernel
#'   runs sequentially. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_lambda` and
#'   `lambda_lambda`, in that order, each of length
#'   `max(length(y), length(mu), length(lambda))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The expected information is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the observed information under the model. The inverse
#' Gaussian is a regular family, so the second Bartlett identity holds and this
#' equals the variance of the score. \eqn{\lambda} names this family's shape
#' parameter.
#'
#' @seealso [distrib_hessian.InvGauss2Distrib()] for the observed Hessian, of
#'   which this is the expectation, [fisher_scoring()], which inverts it at
#'   each step,
#'   and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' th <- list(mu = 2, lambda = 3)
#'
#' # Each of the three constants has one value per observation.
#' lapply(distrib_expected_hessian(d, c(1, 2, 3), th), unique)
#' c(-3 / 2^3, 0, -1 / (2 * 3^2))
#'
#' # The pure shape entry equals the observed one at every observation.
#' identical(distrib_expected_hessian(d, c(1, 2, 3), th)$lambda_lambda,
#'           distrib_hessian(d, c(1, 2, 3), th)$lambda_lambda)
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_expected_hessian, InvGauss2Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    method = NULL, ..., threads = 1L) {
  invgauss2_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Inverse Gaussian Third-Order Derivatives in Mean and Shape
#' @name distrib_deriv3.InvGauss2Distrib
#' @description
#' Computes the four distinct third derivatives of the inverse Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\lambda}, in closed form:
#' \deqn{\ell^{(\mu\mu\mu)} = \dfrac{\lambda(12y - 6\mu)}{\mu^5}, \qquad
#' \ell^{(\mu\mu\lambda)} = \dfrac{2\mu - 3y}{\mu^4}, \qquad
#' \ell^{(\mu\lambda\lambda)} = 0, \qquad \ell^{(\lambda\lambda\lambda)} =
#' \dfrac{1}{\lambda^3}.} The log-density is linear in \eqn{\lambda} apart from
#' \eqn{\tfrac12\log\lambda}, so any component naming the shape twice and the
#' mean once is exactly zero and the pure shape component carries only the
#' logarithm. `expected` is passed straight to the kernel, which substitutes
#' \eqn{\mathbb{E}[Y] = \mu}; both routes are closed form.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL` or an object inheriting from [ExpectMethod]. The generic
#'   applies an object before dispatching on the family, so this method
#'   receives `NULL`, and it returns the closed form of the family. See
#'   [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below a fixed number of observations, set per kernel, the kernel
#'   runs sequentially. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_lambda`,
#'   `mu_lambda_lambda` and `lambda_lambda_lambda`, each of length
#'   `max(length(y), length(mu), length(lambda))`. The names enumerate the
#'   distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\lambda} names this family's shape
#' parameter.
#'
#' @seealso [distrib_hessian.InvGauss2Distrib()] for the order below and
#'   [distrib_deriv4.InvGauss2Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- invgauss2_distrib()
#' y <- c(1, 2, 3)
#' th <- list(mu = 2, lambda = 3)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # The mixed mu-lambda-lambda component is exactly zero.
#' d3$mu_lambda_lambda
#'
#' # The pure shape component is 1/lambda^3, free of the data.
#' unique(d3$lambda_lambda_lambda)
#' 1 / 3^3
#'
#' # Expected values: y is replaced by mu.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 2 + eps, lambda = 3))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 2 - eps, lambda = 3))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
S7::method(distrib_deriv3, InvGauss2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          method = NULL, ..., threads = 1L) {
  invgauss2_deriv3_cpp(y, theta[[1]], theta[[2]], expected, threads)
}

#' @title Inverse Gaussian Fourth-Order Derivatives in Mean and Shape
#' @name distrib_deriv4.InvGauss2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the inverse Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\lambda}, in closed form. As
#' at third order the log-density is linear in \eqn{\lambda} apart from
#' \eqn{\tfrac12\log\lambda}, so every component naming the shape twice or more
#' alongside the mean is exactly zero, and the pure shape component is
#' \eqn{-3/\lambda^4}.
#'
#' `expected` is passed straight to the kernel, which substitutes
#' \eqn{\mathbb{E}[Y] = \mu}; both routes are closed form.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param y A numeric vector of strictly positive observations. With
#'   `expected = TRUE` only its length is used.
#' @param theta A named list with components `mu` and `lambda`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL` or an object inheriting from [ExpectMethod]. The generic
#'   applies an object before dispatching on the family, so this method
#'   receives `NULL`, and it returns the closed form of the family. See
#'   [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below a fixed number of observations, set per kernel, the kernel
#'   runs sequentially. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_lambda`, `mu_mu_lambda_lambda`, `mu_lambda_lambda_lambda` and
#'   `lambda_lambda_lambda_lambda`, each of length
#'   `max(length(y), length(mu), length(lambda))`.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\lambda} names this
#' family's shape parameter.
#'
#' @seealso [distrib_deriv3.InvGauss2Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.InvGauss2Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- invgauss2_distrib()
#' y <- c(1, 2, 3)
#' th <- list(mu = 2, lambda = 3)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # Two components are exactly zero and the pure shape one is -3/lambda^4.
#' c(unique(d4$mu_mu_lambda_lambda),
#'   unique(d4$mu_lambda_lambda_lambda),
#'   unique(d4$lambda_lambda_lambda_lambda))
#' -3 / 3^4
#'
#' # Expected values: y is replaced by mu.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 2 + eps, lambda = 3))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 2 - eps, lambda = 3))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
S7::method(distrib_deriv4, InvGauss2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          method = NULL, ..., threads = 1L) {
  invgauss2_deriv4_cpp(y, theta[[1]], theta[[2]], expected, threads)
}


#' Inverse Gaussian Distribution, Mean and Shape
#'
#' @description
#' Builds the distribution object for the inverse Gaussian family on
#' \eqn{(0, \infty)} in its classical parametrization, the mean \eqn{\mu > 0}
#' and the shape \eqn{\lambda > 0}, so that
#' \eqn{\operatorname{Var}(Y) = \mu^3/\lambda}. The returned object carries
#' closed-form derivatives of the log-density to fourth order, in the
#' parameters and in the response, and closed-form moments. The fifth order and
#' [distrib_cross3_y()] come from the numerical fallbacks of the base class.
#'
#' The two arguments set the links that carry each parameter to the
#' unconstrained scale on which an optimizer works. Both default to the
#' logarithm, because both parameters are positive.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted mean positive.
#' @param link_lambda A `link` object from `linkfunctions7` for the shape
#'   \eqn{\lambda}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, \infty)} is
#' \deqn{f(y; \mu, \lambda) = \sqrt{\dfrac{\lambda}{2\pi y^{3}}}
#'       \exp\left\{-\dfrac{\lambda(y-\mu)^{2}}{2\mu^{2}y}\right\},}
#' the mean is \eqn{\mu}, the variance \eqn{\mu^3/\lambda}, the skewness
#' \eqn{3\sqrt{\mu/\lambda}} and the excess kurtosis \eqn{15\mu/\lambda}. The
#' distribution function is elementary in the standard normal one and the
#' quantile function is its numerical inverse.
#'
#' This is the same law as [invgauss1_distrib()], which carries a dispersion
#' \eqn{\phi = 1/\lambda}. The map between the two moves one coordinate at a
#' time and leaves the mean untouched, so both sets of derivatives stay
#' elementary.
#' The family is usually written in the shape parametrization. In it the law is
#' the distribution of the first-passage time to a level \eqn{a > 0} of a
#' Brownian motion with drift \eqn{\nu > 0} and variance parameter
#' \eqn{s^2}, with \eqn{\mu = a/\nu} and \eqn{\lambda = a^2/s^2}.
#'
#' # Derivatives
#'
#' The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\lambda(y-\mu)}{\mu^3},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{2\lambda}
#'         - \dfrac{(y-\mu)^2}{2\mu^2 y},}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{\lambda}{\mu^3},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\lambda)}\right] = 0, \qquad
#'       \mathbb{E}\left[\ell^{(\lambda\lambda)}\right] =
#'         -\dfrac{1}{2\lambda^2}.}
#' The zero off-diagonal makes the mean and the shape orthogonal.
#'
#' The log-density is linear in \eqn{\lambda} apart from
#' \eqn{\tfrac12\log\lambda}, and that shapes every higher order. The pure
#' shape derivatives are those of the logarithm alone,
#' \eqn{-1/(2\lambda^2)}, \eqn{1/\lambda^3} and \eqn{-3/\lambda^4}; every
#' component naming the shape twice or more alongside the mean is exactly zero;
#' and each remaining component is at most linear in the response, so every
#' expectation needs only \eqn{\mathbb{E}[Y] = \mu} and all four orders are
#' closed form. It also makes the pure shape entry of the observed Hessian a
#' negative constant, where the corresponding entry in the dispersion
#' parametrization can be positive.
#'
#' Third and fourth orders are in [distrib_deriv3.InvGauss2Distrib()] and
#' [distrib_deriv4.InvGauss2Distrib()]. The derivatives in the response and the
#' mixed block are registered elsewhere in the package and are closed form as
#' well; see [distrib_grad_y.InvGauss2Distrib()].
#'
#' # Estimation
#'
#' Both maximum likelihood estimates are available in closed form:
#' \deqn{\hat\mu = \bar y, \qquad
#'       \dfrac{1}{\hat\lambda} = \dfrac{1}{n}\sum_i \dfrac{1}{y_i}
#'         - \dfrac{1}{\bar y}.}
#' [fit_distrib()] reaches them numerically on the link scale, and the example
#' below checks both against the sample.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\lambda > 0} the shape, with
#' \eqn{\operatorname{Var}(Y) = \mu^3/\lambda}. Here \eqn{\lambda} is this
#' family's shape and not a penalty parameter or an eigenvalue. \eqn{\eta} is a
#' parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `InvGauss2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"invgauss2"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "lambda")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [invgauss1_distrib()] for the same law in the mean and a dispersion, which
#' is the generalized linear model parametrization; [gamma1_distrib()] and
#' [lognormal1_distrib()] for other positive families;
#' [pig1_distrib()] for the Poisson mixed over this law;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' user-defined family against the same battery of checks that this family
#' passes;
#' [InvGauss2Distrib] for the class.
#'
#' @references
#' Chhikara, R. S. and Folks, J. L. (1989). *The Inverse Gaussian
#' Distribution: Theory, Methodology, and Applications*. Marcel Dekker,
#' New York.
#'
#' @examples
#' d <- invgauss2_distrib()
#' d
#'
#' # The density is the same as that of invgauss1 at phi = 1/lambda.
#' y <- c(1, 2, 3)
#' th <- list(mu = 2, lambda = 3)
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(invgauss1_distrib(), y, list(mu = 2, phi = 1 / 3)))
#'
#' # The variance is mu^3/lambda and the skewness is 3 sqrt(mu/lambda).
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(2^3 / 3, 3 * sqrt(2 / 3), 15 * 2 / 3)
#'
#' # Fitting recovers the closed-form maximum likelihood estimates.
#' set.seed(3)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted = coef(fit),
#'       closed = c(mu = mean(z),
#'                  lambda = 1 / (mean(1 / z) - 1 / mean(z))))
#'
#' # The mean and the shape are orthogonal: the mixed entry is 0.
#' distrib_expected_hessian(d, 1, th)$mu_lambda
#'
#' @export
invgauss2_distrib <- function(link_mu = log_link(), link_lambda = log_link()) {
  InvGauss2Distrib(
    distrib_name = "invgauss2",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "lambda"),
    params_interpretation = c(mu = "mean", lambda = "shape"),
    n_params = 2,
    params_bounds = list(mu = c(0, Inf), lambda = c(0, Inf)),
    link_params = list(mu = link_mu, lambda = link_lambda)
  )
}
