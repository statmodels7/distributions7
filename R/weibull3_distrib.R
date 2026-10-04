#' @include distrib.R generics.R numerical_functions.R cdf_derivatives.R cdf_higher.R dexpected_families.R moments.R
NULL

#' @title Weibull Distribution Class, Mean and Shape
#' @name Weibull3Distrib
#'
#' @description
#' The S7 class of the Weibull family parametrized by its mean \eqn{m} and its
#' shape \eqn{\sigma}. It inherits from `continuous_distrib`. Build one with
#' [weibull3_distrib()]; this page documents the raw S7 constructor, which
#' validates none of the relationships between the properties.
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
#' Registered on this class in this file: [distrib_pdf()], [distrib_cdf()],
#' [distrib_quantile()], [distrib_rng()], [distrib_gradient()],
#' [distrib_hessian()], [distrib_deriv3()], [distrib_deriv4()],
#' [distrib_deriv5()], [distrib_expected_hessian()], the derivatives of the
#' expected information through [register_dexpected()], [distrib_grad_y()],
#' [distrib_hess_y()], [distrib_deriv3_y()], [distrib_deriv4_y()],
#' [distrib_cross_y()], [distrib_cross2_y()], [distrib_grad_y_hess()],
#' [distrib_hess_y_hess()], [distrib_grad_cdf()], [distrib_hess_cdf()],
#' [distrib_deriv3_cdf()], [distrib_deriv4_cdf()], and the four moments
#' [mean()], [variance()], [skewness()] and [kurtosis()].
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
#' @description Returns \eqn{b = m / \Gamma(1 + 1/\sigma)}, formed on the log
#'   scale with `lgamma()`.
#' @param theta A list with the mean and the shape, in that order.
#' @return A numeric vector of scales.
#' @seealso [distrib_pdf.Weibull3Distrib()]
#' @keywords internal
weibull3_scale <- function(theta) {
  exp(log(theta[[1]]) - lgamma(1 + 1 / theta[[2]]))
}


#' @title Weibull Density, Distribution and Generator in the Mean and Shape
#' @name distrib_pdf.Weibull3Distrib
#' @aliases distrib_cdf.Weibull3Distrib distrib_quantile.Weibull3Distrib
#'   distrib_rng.Weibull3Distrib
#'
#' @description
#' The Weibull with shape \eqn{\sigma} and scale
#' \eqn{b = m / \Gamma(1 + 1/\sigma)}, so that \eqn{E[Y] = m}. The density,
#' the distribution function, the quantile function and the generator are
#' [stats::dweibull()], [stats::pweibull()], [stats::qweibull()] and
#' [stats::rweibull()] at that scale.
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
  stats::dweibull(y, shape = theta[[2]], scale = weibull3_scale(theta), log = log)
}

S7::method(distrib_cdf, Weibull3Distrib) <- function(distrib, q, theta,
                                                     lower.tail = TRUE,
                                                     log.p = FALSE, ...) {
  stats::pweibull(q, shape = theta[[2]], scale = weibull3_scale(theta),
                  lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_quantile, Weibull3Distrib) <- function(distrib, p, theta,
                                                          lower.tail = TRUE,
                                                          log.p = FALSE, ...) {
  stats::qweibull(p, shape = theta[[2]], scale = weibull3_scale(theta),
                  lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_rng, Weibull3Distrib) <- function(distrib, n, theta, ...) {
  stats::rweibull(n, shape = theta[[2]], scale = weibull3_scale(theta))
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
#' each from its own compiled kernel.
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
#' @param scale `"parameter"` or `"link"`; the link scale is applied by the
#'   generic.
#' @param expected Logical; for orders three and four, whether the expected
#'   derivative is returned.
#' @param approx,nsim Accepted for the generic's signature; the expectations
#'   are closed.
#' @param ...,threads Unused.
#'
#' @return A named list with one numeric vector per component: 2, 3, 4, 5 and
#'   6 components at orders one to five, the second order keyed
#'   `mean_mean`, `sigma_sigma`, `mean_sigma`.
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


#' @title Weibull Derivatives in the Response, Mean and Shape
#' @name distrib_grad_y.Weibull3Distrib
#' @aliases distrib_hess_y.Weibull3Distrib distrib_deriv3_y.Weibull3Distrib
#'   distrib_deriv4_y.Weibull3Distrib distrib_cross_y.Weibull3Distrib
#'   distrib_cross2_y.Weibull3Distrib distrib_grad_y_hess.Weibull3Distrib
#'   distrib_hess_y_hess.Weibull3Distrib
#'
#' @description
#' The derivatives of the log-density in the response to order four, and the
#' mixed derivatives of orders one and two in the response and one and two in
#' \eqn{(m, \sigma)}, each a closed form from its own compiled kernel.
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
  weibull3_dy1_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_hess_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy2_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_deriv3_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy3_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_deriv4_y, Weibull3Distrib) <- function(distrib, y, theta, ...) {
  weibull3_dy4_cpp(y, theta[[1]], theta[[2]])$y
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
#' The derivatives of \eqn{F(q) = 1 - \exp\{-(q/b)^\sigma\}} in
#' \eqn{(m, \sigma)} to order four, each a closed form from its own compiled
#' kernel, carried to the upper tail and the log scale by
#' [cdf_tail_scale()] and [cdf_scale_k()].
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
  d1 <- weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]])
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL,
                 lower.tail, log)
}

S7::method(distrib_hess_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  d1 <- weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]])
  d2 <- weibull3_dcdf2_cpp(q, theta[[1]], theta[[2]])
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, d2,
                 lower.tail, log)
}

S7::method(distrib_deriv3_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  dF <- list(weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf2_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf3_cpp(q, theta[[1]], theta[[2]]))
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), dF, 3L, lower.tail, log)
}

S7::method(distrib_deriv4_cdf, Weibull3Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  dF <- list(weibull3_dcdf1_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf2_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf3_cpp(q, theta[[1]], theta[[2]]),
             weibull3_dcdf4_cpp(q, theta[[1]], theta[[2]]))
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), dF, 4L, lower.tail, log)
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
#'   shape \eqn{s}, which do not depend on the scale, from the ratios
#'   \eqn{\Gamma(1 + k/s)/\Gamma(1 + 1/s)^k} formed with `lgamma()`.
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
#' derivatives of the distribution function are closed forms, each order in
#' its own compiled kernel (see [distrib_gradient.Weibull3Distrib()]).
#'
#' The number follows gamlss, where the Weibull in the mean is `WEI3`.
#' Leaving `weibull2` unused is deliberate: it names a different
#' parametrization there.
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
