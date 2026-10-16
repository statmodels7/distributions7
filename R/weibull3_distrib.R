#' @include distrib.R generics.R numerical_functions.R cdf_derivatives.R cdf_higher.R dexpected_families.R moments.R
NULL

#' @title Weibull Distribution Class, Mean and Shape
#' @name Weibull3Distrib
#'
#' @description
#' This is the S7 class of the Weibull family parametrized by its mean \eqn{m}
#' and its shape \eqn{\sigma}. It inherits from `continuous_distrib`. Build one
#' with [weibull3_distrib()]; this page documents the raw S7 constructor; the
#' validator of `distrib` checks that the properties are consistent (the
#' parameter count, a bound and a link for every parameter, ordered support
#' bounds).
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Weibull3Distrib`, inheriting from
#'   `continuous_distrib`. For an object built by [weibull3_distrib()] the
#'   properties hold `"weibull3"`, `"univariate"`, `c(0, Inf)`,
#'   `c("mean", "sigma")`, the interpretations `c(mean = "mean",
#'   sigma = "shape")`, `2` and the domains \eqn{(0, \infty)} for both.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.Weibull3Distrib],
#'   [`distrib_cdf()`][distrib_cdf.Weibull3Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Weibull3Distrib] and
#'   [`distrib_rng()`][distrib_rng.Weibull3Distrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.Weibull3Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Weibull3Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Weibull3Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Weibull3Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Weibull3Distrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.Weibull3Distrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.Weibull3Distrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.Weibull3Distrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.Weibull3Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Weibull3Distrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.Weibull3Distrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.Weibull3Distrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.Weibull3Distrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.Weibull3Distrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.Weibull3Distrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.Weibull3Distrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.Weibull3Distrib],
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.Weibull3Distrib],
#'   [`distrib_deriv3_cdf()`][distrib_deriv3_cdf.Weibull3Distrib] and
#'   [`distrib_deriv4_cdf()`][distrib_deriv4_cdf.Weibull3Distrib];
#' - the moments: [`mean()`][mean.Weibull3Distrib],
#'   [`variance()`][variance.Weibull3Distrib],
#'   [`skewness()`][skewness.Weibull3Distrib] and
#'   [`kurtosis()`][kurtosis.Weibull3Distrib].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @seealso [weibull3_distrib()] to build one; [Weibull1Distrib] for the
#'   parametrization by the scale.
#'
#' @examples
#' d <- weibull3_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params
Weibull3Distrib <- S7::new_class("Weibull3Distrib", parent = continuous_distrib)

#' @title The Scale of a Weibull with a Given Mean
#' @name weibull3_scale
#' @description `weibull3_scale()` returns \eqn{b = m / \Gamma(1 + 1/\sigma)},
#'   formed on the log scale with `lgamma()`, and `weibull3_log_scale()`
#'   returns \eqn{\log b}. `weibull3_scale_ok()` reports where \eqn{b} itself
#'   is a finite positive number; at a small shape it underflows to 0, and
#'   the density, distribution, quantile and generator methods then work with
#'   \eqn{\log b}.
#' @param theta A list with the mean and the shape, in that order.
#' @param lb A numeric vector of log-scales, for `weibull3_scale_ok()`.
#' @return A numeric vector of scales, of log-scales, or of logicals.
#' @seealso [distrib_pdf.Weibull3Distrib()]
#' @keywords internal
weibull3_scale <- function(theta) {
  exp(weibull3_log_scale(theta))
}

#' @rdname weibull3_scale
#' @keywords internal
weibull3_log_scale <- function(theta) {
  log(theta[[1]]) - lgamma(1 + 1 / theta[[2]])
}

#' @rdname weibull3_scale
#' @keywords internal
weibull3_scale_ok <- function(lb) {
  b <- exp(lb)
  is.finite(b) & b > 0
}


#' @title Weibull Density, Distribution and Generator in the Mean and Shape
#' @name distrib_pdf.Weibull3Distrib
#' @aliases distrib_cdf.Weibull3Distrib distrib_quantile.Weibull3Distrib
#'   distrib_rng.Weibull3Distrib
#'
#' @description
#' The methods of this page describe the Weibull with shape \eqn{\sigma} and
#' scale \eqn{b = m / \Gamma(1 + 1/\sigma)}, so that \eqn{E[Y] = m}. The density,
#' the distribution function, the quantile function and the generator are
#' [stats::dweibull()], [stats::pweibull()], [stats::qweibull()] and
#' [stats::rweibull()] at that scale. At a small shape \eqn{b} underflows to
#' 0; there the four methods evaluate the same formulas through
#' \eqn{\log b} and \eqn{\log(y/b)}. On the log scale of the lower tail the
#' distribution function and the quantile go through the logarithm of the
#' cumulative hazard ([weibull_log_lower()], [log_neg_log1mexp()]), so that
#' they do not round where the hazard or \eqn{1 - p} would.
#'
#' @param distrib A `Weibull3Distrib` object, from [weibull3_distrib()].
#' @param y,q A numeric vector of observations or quantiles.
#' @param p A numeric vector of probabilities.
#' @param n The number of draws.
#' @param theta A named list with components `mean` and `sigma`.
#' @param log,log.p,lower.tail As in [stats::dweibull()] and its siblings.
#' @param ... Unused.
#'
#' @return A numeric vector.
#'
#' @seealso [weibull3_distrib()].
#'
#' @examples
#' d <- weibull3_distrib()
#' th <- list(mean = 4, sigma = 1.7)
#' integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf)$value
S7::method(distrib_pdf, Weibull3Distrib) <- function(distrib, y, theta,
                                                     log = FALSE, ...) {
  lb <- weibull3_log_scale(theta)
  if (all(weibull3_scale_ok(lb))) {
    return(stats::dweibull(y, shape = theta[[2]], scale = exp(lb), log = log))
  }
  n <- max(length(y), length(lb))
  y <- rep_len(y, n); s <- rep_len(theta[[2]], n); lb <- rep_len(lb, n)
  lz <- suppressWarnings(log(y)) - lb
  out <- log(s) - lb + (s - 1) * lz - exp(s * lz)
  out[y < 0 | y == Inf] <- -Inf
  z <- y == 0
  out[z] <- ifelse(s[z] < 1, Inf, ifelse(s[z] > 1, -Inf, -lb[z]))
  if (log) out else exp(out)
}

S7::method(distrib_cdf, Weibull3Distrib) <- function(distrib, q, theta,
                                                     lower.tail = TRUE,
                                                     log.p = FALSE, ...) {
  lb <- weibull3_log_scale(theta)
  # log F from the logarithm of the cumulative hazard, which pweibull() loses
  # once the hazard itself underflows
  if (lower.tail && log.p) {
    return(weibull_log_lower(theta[[2]] * (suppressWarnings(log(pmax(q, 0))) - lb)))
  }
  if (all(weibull3_scale_ok(lb))) {
    return(stats::pweibull(q, shape = theta[[2]], scale = exp(lb),
                           lower.tail = lower.tail, log.p = log.p))
  }
  # H = (q / b)^sigma is the cumulative hazard; F = 1 - exp(-H)
  H <- exp(theta[[2]] * (suppressWarnings(log(pmax(q, 0))) - lb))
  if (lower.tail) {
    if (log.p) log(-expm1(-H)) else -expm1(-H)
  } else {
    if (log.p) -H else exp(-H)
  }
}

S7::method(distrib_quantile, Weibull3Distrib) <- function(distrib, p, theta,
                                                          lower.tail = TRUE,
                                                          log.p = FALSE, ...) {
  lb <- weibull3_log_scale(theta)
  # the cumulative hazard from log F without forming 1 - F, which qweibull()
  # rounds to 1 once F is below eps
  if (lower.tail && log.p) return(exp(lb + log_neg_log1mexp(p) / theta[[2]]))
  if (all(weibull3_scale_ok(lb))) {
    return(stats::qweibull(p, shape = theta[[2]], scale = exp(lb),
                           lower.tail = lower.tail, log.p = log.p))
  }
  # the cumulative hazard at the quantile, then q = b H^(1/sigma)
  H <- if (lower.tail) {
    if (log.p) -log(-expm1(p)) else -log1p(-p)
  } else {
    if (log.p) -p else -log(p)
  }
  exp(lb + log(H) / theta[[2]])
}

S7::method(distrib_rng, Weibull3Distrib) <- function(distrib, n, theta, ...) {
  lb <- weibull3_log_scale(theta)
  if (all(weibull3_scale_ok(lb))) {
    return(stats::rweibull(n, shape = theta[[2]], scale = exp(lb)))
  }
  # the inversion stats::rweibull() performs, b (-log U)^(1/sigma), on the log
  # scale
  exp(lb + log(-log(stats::runif(n))) / theta[[2]])
}


#' @title Weibull Derivatives in the Mean and Shape
#' @name distrib_gradient.Weibull3Distrib
#' @aliases distrib_hessian.Weibull3Distrib distrib_deriv3.Weibull3Distrib
#'   distrib_deriv4.Weibull3Distrib distrib_deriv5.Weibull3Distrib
#'   distrib_expected_hessian.Weibull3Distrib
#'   distrib_dexpected_hessian.Weibull3Distrib
#'   distrib_d2expected_hessian.Weibull3Distrib
#'
#' @description
#' Return the derivatives of the log-density in \eqn{(m, \sigma)} of orders
#' one to five, the expected information and its expected third and fourth
#' derivatives, and the first two derivatives of the expected information,
#' and each of them is computed by its own compiled kernel.
#'
#' @details
#' With \eqn{w = \log(y/m) + \log\Gamma(1 + 1/\sigma) = \log(y/b)} the
#' log-density is
#' \deqn{\ell = \log\sigma - \log y + \sigma w - e^{\sigma w}.}
#' Every component is a closed form in \eqn{\log(y/m)},
#' \eqn{\log\Gamma(1 + 1/\sigma)} and the polygamma functions at
#' \eqn{1 + 1/\sigma}, derived offline and written out per component, one
#' kernel per order. \eqn{T = e^{\sigma w}} is a standard exponential, so an
#' expected derivative replaces each \eqn{T^j (\log T)^k} by
#' \eqn{\Gamma^{(k)}(j + 1)}.
#'
#' @param distrib A `Weibull3Distrib` object.
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mean` and `sigma`.
#' @param scale `"parameter"` or `"link"`. For orders one to four the
#'   generic applies the link scale after dispatch; the method for
#'   [distrib_deriv5()] applies it itself.
#' @param expected Logical; for orders three and four, whether the expected
#'   derivative is returned.
#' @param approx,nsim Accepted for the generic's signature; the expectations
#'   are closed.
#' @param ...,threads Unused.
#'
#' @return A named list with one numeric vector per component: 2, 3, 4, 5 and
#'   6 components at orders one to five. At order two the components are keyed
#'   `mean_mean`, `sigma_sigma` and `mean_sigma`.
#'
#' @seealso [distrib_pdf.Weibull3Distrib()].
#'
#' @examples
#' d <- weibull3_distrib()
#' th <- list(mean = 4, sigma = 1.7)
#' y <- c(0.7, 2.5, 6)
#' g <- distrib_gradient(d, y, th)
#' h <- 1e-6
#' up <- distrib_pdf(d, y, list(mean = 4, sigma = 1.7 + h), log = TRUE)
#' dn <- distrib_pdf(d, y, list(mean = 4, sigma = 1.7 - h), log = TRUE)
#' all.equal(g$sigma, (up - dn) / (2 * h), tolerance = 1e-6)
S7::method(distrib_gradient, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  weibull3_gradient_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_hessian, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  weibull3_hessian_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_deriv3, Weibull3Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) weibull3_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  else weibull3_deriv3_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_deriv4, Weibull3Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) weibull3_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  else weibull3_deriv4_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_deriv5, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, weibull3_deriv5_cpp(y, theta[[1]], theta[[2]]),
               match.arg(scale))
}

S7::method(distrib_expected_hessian, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  weibull3_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

register_dexpected(Weibull3Distrib, function(d, y, th, k, t) {
  if (k == 1L) weibull3_dexpected1_cpp(y, th[[1]], th[[2]])
  else weibull3_dexpected2_cpp(y, th[[1]], th[[2]])
})


#' Response Derivatives of the Weibull at Zero
#'
#' @description
#' Replaces the entries of `out` at `y = 0` by the derivative of order `k` of
#' \eqn{(\sigma - 1)\log y - b^{-\sigma} y^{\sigma}} there, from [dy_log()] and
#' [dy_pow()]. The compiled response kernels of [weibull3_distrib()] divide by
#' \eqn{y} and give `NaN` at zero, where the derivative exists when
#' \eqn{\sigma = 1} and is infinite otherwise.
#'
#' @param out The kernel's result, a numeric vector.
#' @param y The response, of the length of `out`.
#' @param theta A list with the mean and the shape, in that order.
#' @param k The derivative order, 1 to 4.
#'
#' @return `out` with its entries at `y = 0` replaced.
#'
#' @keywords internal
weibull3_dy_at_zero <- function(out, y, theta, k) {
  z <- which(y == 0)
  if (length(z) == 0L) return(out)
  n <- length(out)
  s <- rep_len(theta[[2]], n)[z]
  lb <- rep_len(weibull3_log_scale(theta), n)[z]
  out[z] <- dy_log(s - 1, 0, k) + dy_pow(-exp(-s * lb), s, 0, k)
  out
}

#' @title Weibull Derivatives in the Response, Mean and Shape
#' @name distrib_grad_y.Weibull3Distrib
#' @aliases distrib_hess_y.Weibull3Distrib distrib_deriv3_y.Weibull3Distrib
#'   distrib_deriv4_y.Weibull3Distrib distrib_cross_y.Weibull3Distrib
#'   distrib_cross2_y.Weibull3Distrib distrib_grad_y_hess.Weibull3Distrib
#'   distrib_hess_y_hess.Weibull3Distrib
#'
#' @description
#' Computes the derivatives of the log-density in the response to order four,
#' and the mixed derivatives of orders one and two in the response and one and
#' two in \eqn{(m, \sigma)}. Each is a closed form computed by its own
#' compiled kernel.
#'
#' @param distrib A `Weibull3Distrib` object.
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mean` and `sigma`.
#' @param scale `"parameter"` or `"link"`, for the mixed derivatives.
#' @param ... Unused.
#'
#' @return A numeric vector for the derivatives in the response; a named list,
#'   keyed by parameter or by parameter pair, for the mixed derivatives.
#'
#' @seealso [distrib_gradient.Weibull3Distrib()].
#'
#' @examples
#' d <- weibull3_distrib()
#' distrib_grad_y(d, c(0.7, 2.5), list(mean = 4, sigma = 1.7))
S7::method(distrib_grad_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy_at_zero(weibull3_dy1_cpp(y, theta[[1]], theta[[2]])$y, y, theta, 1L)
}
S7::method(distrib_hess_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy_at_zero(weibull3_dy2_cpp(y, theta[[1]], theta[[2]])$y, y, theta, 2L)
}
S7::method(distrib_deriv3_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy_at_zero(weibull3_dy3_cpp(y, theta[[1]], theta[[2]])$y, y, theta, 3L)
}
S7::method(distrib_deriv4_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy_at_zero(weibull3_dy4_cpp(y, theta[[1]], theta[[2]])$y, y, theta, 4L)
}
S7::method(distrib_cross_y, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  weibull3_cross_y_cpp(y, theta[[1]], theta[[2]])
}
S7::method(distrib_cross2_y, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  weibull3_cross2_y_cpp(y, theta[[1]], theta[[2]])
}
S7::method(distrib_grad_y_hess, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  weibull3_grad_y_hess_cpp(y, theta[[1]], theta[[2]])
}
S7::method(distrib_hess_y_hess, Weibull3Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  weibull3_hess_y_hess_cpp(y, theta[[1]], theta[[2]])
}


#' @title Weibull Distribution-Function Derivatives, Mean and Shape
#' @name distrib_grad_cdf.Weibull3Distrib
#' @aliases distrib_hess_cdf.Weibull3Distrib
#'   distrib_deriv3_cdf.Weibull3Distrib distrib_deriv4_cdf.Weibull3Distrib
#'
#' @description
#' Computes the derivatives of \eqn{F(q) = 1 - \exp\{-(q/b)^\sigma\}} in
#' \eqn{(m, \sigma)} to order four. Each is a closed form computed by its own
#' compiled kernel and carried to the upper tail and the log scale by
#' [cdf_tail_scale()] and [cdf_scale_k()]. The derivatives of the log
#' survival function (`lower.tail = FALSE, log = TRUE`) are those of
#' \eqn{-(q/b)^\sigma} itself, from [weibull3_log_surv_derivs()], so that
#' they are finite where \eqn{1 - F} underflows.
#'
#' @param distrib A `Weibull3Distrib` object.
#' @param q A numeric vector of strictly positive quantiles.
#' @param theta A named list with components `mean` and `sigma`.
#' @param lower.tail,log As in [distrib_grad_cdf()].
#' @param ... Unused.
#'
#' @return A named list, one component per parameter or per multi-index.
#'
#' @seealso [distrib_cdf.Weibull3Distrib()].
#'
#' @examples
#' d <- weibull3_distrib()
#' distrib_grad_cdf(d, c(1, 3), list(mean = 4, sigma = 1.7), log = FALSE)
S7::method(distrib_grad_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  if (!lower.tail && log) return(weibull3_log_surv_derivs(distrib, q, theta, 1L))
  d1 <- weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]])
  cdf_tail_scale(distrib, tail_cdf(distrib, q, theta, lower.tail), d1, NULL,
                 lower.tail, log)
}

S7::method(distrib_hess_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  if (!lower.tail && log) return(weibull3_log_surv_derivs(distrib, q, theta, 2L))
  d1 <- weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]])
  d2 <- weibull3_dcdf2_cpp(q, theta[[1]], theta[[2]])
  cdf_tail_scale(distrib, tail_cdf(distrib, q, theta, lower.tail), d1, d2,
                 lower.tail, log)
}

S7::method(distrib_deriv3_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  if (!lower.tail && log) return(weibull3_log_surv_derivs(distrib, q, theta, 3L))
  dF <- list(weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf2_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf3_cpp(q, theta[[1]], theta[[2]]))
  cdf_scale_k(distrib, tail_cdf(distrib, q, theta, lower.tail), dF, 3L, lower.tail, log)
}

S7::method(distrib_deriv4_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  if (!lower.tail && log) return(weibull3_log_surv_derivs(distrib, q, theta, 4L))
  dF <- list(weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf2_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf3_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf4_cpp(q, theta[[1]], theta[[2]]))
  cdf_scale_k(distrib, tail_cdf(distrib, q, theta, lower.tail), dF, 4L, lower.tail, log)
}


#' Log-Survival Derivatives of the Weibull, Mean and Shape
#'
#' @description
#' Computes the derivatives of order `o` of
#' \eqn{\log S(q) = -e^{h}}, \eqn{h = \sigma\{\log q - \log m +
#' \log\Gamma(1 + 1/\sigma)\}}, in \eqn{(m, \sigma)}, as
#' \eqn{-e^{h}} times the complete Bell polynomial in the partial derivatives
#' of \eqn{h} ([bell_f_ratio()]). The function \eqn{h} is \eqn{\sigma} times a
#' sum of a function of \eqn{m} and a function of \eqn{\sigma}, so a mixed
#' partial derivative of \eqn{h} is nonzero only when it names \eqn{\sigma}
#' at most once; the pure derivatives in \eqn{\sigma} are
#' \eqn{\sigma g_k + k g_{k-1}}, where \eqn{g_k} is the derivative of order
#' \eqn{k} of \eqn{\log\Gamma(1 + 1/\sigma)} by the chain rule over the
#' polygamma functions and \eqn{g_0} is the whole sum. At \eqn{q \le 0} every
#' derivative is zero.
#'
#' @param distrib A `Weibull3Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mean` and `sigma`.
#' @param o The derivative order, 1 to 4.
#'
#' @return A named list, one component per multi-index of order `o`, in the
#'   order of [hess_names()] at order two.
#'
#' @examples
#' d <- weibull3_distrib()
#' distributions7:::weibull3_log_surv_derivs(d, c(18, 480),
#'                                           list(mean = 60, sigma = 6), 2L)
#'
#' @keywords internal
weibull3_log_surv_derivs <- function(distrib, q, theta, o) {
  params <- distrib@params
  n <- max(length(q), length(theta[[1]]), length(theta[[2]]))
  q <- rep_len(q, n)
  m <- rep_len(theta[[1]], n)
  sigma <- rep_len(theta[[2]], n)
  inside <- q > 0
  # log q is taken at 1 outside the support, where the result is masked
  lq <- log(ifelse(inside, q, 1))
  x <- 1 / sigma
  # derivatives of x = 1/sigma and of lgamma(1 + x) in sigma
  xd <- list(-x^2, 2 * x^3, -6 * x^4, 24 * x^5)
  p <- lapply(0:3, function(j) psigamma(1 + x, j))
  g <- list(
    p[[1]] * xd[[1]],
    p[[2]] * xd[[1]]^2 + p[[1]] * xd[[2]],
    p[[3]] * xd[[1]]^3 + 3 * p[[2]] * xd[[1]] * xd[[2]] + p[[1]] * xd[[3]],
    p[[4]] * xd[[1]]^4 + 6 * p[[3]] * xd[[1]]^2 * xd[[2]] +
      p[[2]] * (3 * xd[[2]]^2 + 4 * xd[[1]] * xd[[3]]) + p[[1]] * xd[[4]]
  )
  A <- lq - log(m) + lgamma(1 + x)
  h <- sigma * A
  hderiv <- function(block) {
    ns <- sum(block == params[2L])
    j <- sum(block == params[1L])
    if (j == 0L) {
      return(sigma * g[[ns]] + ns * (if (ns == 1L) A else g[[ns - 1L]]))
    }
    if (ns >= 2L) return(0 * q)
    base <- (-1)^j * factorial(j - 1L) / m^j
    if (ns == 1L) base else sigma * base
  }
  eh <- exp(h)
  out <- lapply(deriv_indices(params, o), function(I) {
    v <- -eh * bell_f_ratio(params[I], hderiv)
    ifelse(inside, v, 0)
  })
  out <- stats::setNames(out, deriv_names(params, o))
  if (o == 2L) out[hess_names(params)] else out
}


#' @title Weibull Moments, Mean and Shape
#' @name mean.Weibull3Distrib
#' @aliases variance.Weibull3Distrib skewness.Weibull3Distrib
#'   kurtosis.Weibull3Distrib
#'
#' @description
#' With \eqn{g_k = \Gamma(1 + k/\sigma)}, the mean is \eqn{m}, the variance
#' \eqn{m^2 (g_2/g_1^2 - 1)}, the skewness
#' \eqn{(g_3 - 3 g_1 g_2 + 2 g_1^3)/(g_2 - g_1^2)^{3/2}} and the excess
#' kurtosis \eqn{(g_4 - 4 g_1 g_3 + 6 g_1^2 g_2 - 3 g_1^4)/(g_2 - g_1^2)^2 - 3}.
#' The ratios are formed on the log scale with `lgamma()`.
#'
#' @param x A `Weibull3Distrib` object.
#' @param theta A named list with components `mean` and `sigma`.
#' @param ... Unused.
#'
#' @return A numeric vector, one value per parameter row.
#'
#' @seealso [weibull3_distrib()].
#'
#' @examples
#' d <- weibull3_distrib()
#' th <- list(mean = 4, sigma = 1.7)
#' c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
S7::method(mean, Weibull3Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + moment_const(theta, 2L, 0)
}

S7::method(variance, Weibull3Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  s <- theta[[2]]
  theta[[1]]^2 * expm1(lgamma(1 + 2 / s) - 2 * lgamma(1 + 1 / s)) +
    moment_const(theta, 2L, 0)
}

#' @title The Skewness and Excess Kurtosis of a Weibull Shape
#' @name weibull3_shape_moments
#' @description Returns the skewness and the excess kurtosis of a Weibull of
#'   shape \eqn{s}, computed from the ratios
#'   \eqn{\Gamma(1 + k/s)/\Gamma(1 + 1/s)^k} formed with `lgamma()`. Both are
#'   independent of the scale.
#' @param s A numeric vector of shapes.
#' @return A list with components `skew` and `kurt`.
#' @seealso [skewness.Weibull3Distrib()]
#' @keywords internal
weibull3_shape_moments <- function(s) {
  g <- function(k) exp(lgamma(1 + k / s) - k * lgamma(1 + 1 / s))
  v <- g(2) - 1
  list(skew = (g(3) - 3 * g(2) + 2) / v^1.5,
       kurt = (g(4) - 4 * g(3) + 6 * g(2) - 3) / v^2 - 3)
}

S7::method(skewness, Weibull3Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  weibull3_shape_moments(theta[[2]])$skew + moment_const(theta, 2L, 0)
}

S7::method(kurtosis, Weibull3Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  weibull3_shape_moments(theta[[2]])$kurt + moment_const(theta, 2L, 0)
}


#' Weibull Distribution in the Mean
#'
#' @description
#' Creates a Weibull distribution object parametrized by its mean and its
#' shape.
#'
#' @details
#' The first parameter of [weibull1_distrib()] is the scale and not the mean:
#' the mean is \eqn{b\,\Gamma(1 + 1/\sigma)}. Here the scale is
#' \deqn{b = \dfrac{m}{\Gamma(1 + 1/\sigma)},}
#' and every derivative in \eqn{(m, \sigma)} to order five, the expected
#' information and its derivatives, the derivatives in the response and the
#' derivatives of the distribution function are closed forms, and each order
#' is computed in its own compiled kernel (see
#' [distrib_gradient.Weibull3Distrib()]).
#'
#' The number follows gamlss, where the Weibull in the mean is `WEI3`.
#' This package does not use the name `weibull2`, because that name denotes a
#' different parametrization in gamlss.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{\sigma}{b}\left(\frac{y}{b}\right)^{\sigma-1}
#'   e^{-(y/b)^{\sigma}}, \qquad b = \frac{m}{\Gamma(1+1/\sigma)}}
#' on \eqn{y \in (0, \infty)}, with \eqn{E[Y] = m}.
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_sigma Link function for the shape. Defaults to the log.
#'
#' @return An S7 object of class `Weibull3Distrib`, inheriting from
#'   `continuous_distrib`, with `params` `c("mean", "sigma")` and
#'   `link_params` the two links given here.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. *Journal of the Royal Statistical
#' Society, Series C* 54, 507-554.
#'
#' @seealso [weibull1_distrib()]; [Weibull3Distrib] for the class.
#'
#' @examples
#' d <- weibull3_distrib()
#' theta <- list(mean = 4, sigma = 1.7)
#' mean(d, theta)
#'
#' @importFrom linkfunctions7 log_link
#' @export
weibull3_distrib <- function(link_mean = log_link(), link_sigma = log_link()) {
  Weibull3Distrib(
    distrib_name = "weibull3",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mean", "sigma"),
    params_interpretation = c(mean = "mean", sigma = "shape"),
    n_params = 2,
    params_bounds = list(mean = c(0, Inf), sigma = c(0, Inf)),
    link_params = list(mean = link_mean, sigma = link_sigma)
  )
}
