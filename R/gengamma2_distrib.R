#' @include distrib.R generics.R numerical_functions.R dexpected_families.R moments.R gengamma1_distrib.R cdf_compiled.R y_higher.R cross_derivatives.R cross2_derivatives.R cross_theta2_derivatives.R
NULL

#' @title Generalized Gamma Distribution Class, Mean
#' @name GenGamma2Distrib
#'
#' @description
#' The S7 class of the generalized gamma family parametrized by its mean
#' \eqn{m}, its shape \eqn{d} and its power \eqn{p}. It inherits from
#' `continuous_distrib`. Build one with [gengamma2_distrib()]; this page
#' documents the raw S7 constructor, which validates none of the relationships
#' between the properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `GenGamma2Distrib`, inheriting from
#'   `continuous_distrib`. For an object built by [gengamma2_distrib()] the
#'   properties hold `"gengamma2"`, `"univariate"`, `c(0, Inf)`,
#'   `c("mean", "d", "p")`, the interpretations `c(mean = "mean",
#'   d = "shape", p = "power")`, `3` and the domain \eqn{(0, \infty)} for
#'   each parameter.
#'
#' @section Methods:
#' Registered on this class in this file: [distrib_pdf()], [distrib_cdf()],
#' [distrib_quantile()], [distrib_rng()], [distrib_gradient()],
#' [distrib_hessian()], [distrib_deriv3()], [distrib_deriv4()],
#' [distrib_deriv5()], [distrib_expected_hessian()], the derivatives of the
#' expected information through [register_dexpected()], [distrib_grad_y()],
#' [distrib_hess_y()], [distrib_deriv3_y()], [distrib_deriv4_y()],
#' [distrib_cross_y()], [distrib_cross2_y()], [distrib_grad_y_hess()],
#' [distrib_hess_y_hess()], and the four moments [mean()], [variance()],
#' [skewness()] and [kurtosis()]. The derivatives of the distribution function
#' are the numerical ones of the base class, as for [gengamma1_distrib()].
#'
#' @seealso [gengamma2_distrib()] to build one; [GenGamma1Distrib] for the
#'   parametrization by the scale.
#'
#' @examples
#' d <- gengamma2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params
GenGamma2Distrib <- S7::new_class("GenGamma2Distrib", parent = continuous_distrib)

#' @title The Scale of a Generalized Gamma with a Given Mean
#' @name gengamma2_scale
#' @description Returns \eqn{a = m\,\Gamma(d/p)/\Gamma((d+1)/p)}, formed on the
#'   log scale with `lgamma()`.
#' @param theta A list with the mean, the shape and the power, in that order.
#' @return A numeric vector of scales.
#' @seealso [distrib_pdf.GenGamma2Distrib()]
#' @keywords internal
gengamma2_scale <- function(theta) {
  d <- theta[[2]]; p <- theta[[3]]
  exp(log(theta[[1]]) + lgamma(d / p) - lgamma((d + 1) / p))
}


#' @title Generalized Gamma Density, Distribution and Generator in the Mean
#' @name distrib_pdf.GenGamma2Distrib
#' @aliases distrib_cdf.GenGamma2Distrib distrib_quantile.GenGamma2Distrib
#'   distrib_rng.GenGamma2Distrib
#'
#' @description
#' The generalized gamma with shape \eqn{d}, power \eqn{p} and scale
#' \eqn{a = m\,\Gamma(d/p)/\Gamma((d+1)/p)}, so that \eqn{E[Y] = m}. The
#' density, the distribution function, the quantile function and the generator
#' are those of [gengamma1_distrib()] at that scale: \eqn{(Y/a)^p} is a gamma
#' variable with shape \eqn{d/p} and unit rate.
#'
#' @param distrib A `GenGamma2Distrib` object, from [gengamma2_distrib()].
#' @param y,q A numeric vector of observations or quantiles.
#' @param p A numeric vector of probabilities.
#' @param n The number of draws.
#' @param theta A named list with components `mean`, `d` and `p`.
#' @param log,log.p,lower.tail As in [stats::pgamma()] and its siblings.
#' @param ... Unused.
#'
#' @return A numeric vector.
#'
#' @seealso [gengamma2_distrib()].
#'
#' @examples
#' d <- gengamma2_distrib()
#' th <- list(mean = 5, d = 3, p = 1.5)
#' integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf)$value
S7::method(distrib_pdf, GenGamma2Distrib) <- function(distrib, y, theta,
                                                      log = FALSE, ...) {
  out <- gengamma_logpdf_cpp(y, gengamma2_scale(theta), theta[[2]], theta[[3]], 1L)
  if (log) out else exp(out)
}

S7::method(distrib_cdf, GenGamma2Distrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  pw <- theta[[3]]
  stats::pgamma((pmax(q, 0) / gengamma2_scale(theta))^pw,
                shape = theta[[2]] / pw, rate = 1,
                lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_quantile, GenGamma2Distrib) <- function(distrib, p, theta,
                                                           lower.tail = TRUE,
                                                           log.p = FALSE, ...) {
  pw <- theta[[3]]
  gengamma2_scale(theta) *
    stats::qgamma(p, shape = theta[[2]] / pw, rate = 1,
                  lower.tail = lower.tail, log.p = log.p)^(1 / pw)
}

S7::method(distrib_rng, GenGamma2Distrib) <- function(distrib, n, theta, ...) {
  pw <- theta[[3]]
  gengamma2_scale(theta) * stats::rgamma(n, shape = theta[[2]] / pw, rate = 1)^(1 / pw)
}


#' @title Generalized Gamma Derivatives in the Mean
#' @name distrib_gradient.GenGamma2Distrib
#' @aliases distrib_hessian.GenGamma2Distrib distrib_deriv3.GenGamma2Distrib
#'   distrib_deriv4.GenGamma2Distrib distrib_deriv5.GenGamma2Distrib
#'   distrib_expected_hessian.GenGamma2Distrib
#'   distrib_dexpected_hessian.GenGamma2Distrib
#'   distrib_d2expected_hessian.GenGamma2Distrib
#'
#' @description
#' Return the derivatives of the log-density in \eqn{(m, d, p)} of orders one
#' to five, the expected information and its expected third and fourth
#' derivatives, and the first two derivatives of the expected information,
#' each from its own compiled kernel.
#'
#' @details
#' With \eqn{k = d/p}, \eqn{k_1 = (d+1)/p} and
#' \eqn{w = \log(y/m) - \log\Gamma(k) + \log\Gamma(k_1) = \log(y/a)}, the
#' log-density is
#' \deqn{\ell = \log p - \log y + d\,w - \log\Gamma(k) - e^{p w}.}
#' Every component is a closed form derived offline and written out per
#' component, one kernel per order. With \eqn{U = p\,w}, the data enter
#' through \eqn{X = U - \psi(k)} and \eqn{Q = e^U/k - 1}, and each component
#' is a polynomial in them whose coefficients depend on \eqn{(d, p)} alone,
#' through \eqn{\psi^{(n)}(k)} and the differences
#' \eqn{\psi^{(n)}(k_1) - \psi^{(n)}(k)}. As \eqn{k} and \eqn{d} grow the
#' terms of these coefficients agree to several orders, so above
#' \eqn{k = 10} or \eqn{d = 10} the coefficients are formed in double-double
#' arithmetic (about 32 digits) and rounded once; below, in double. Where
#' \eqn{|X + \psi(k) - \log k| < 1} the polynomial is written in that variable,
#' so that the cancellation between \eqn{X} and \eqn{Q} is also carried out on
#' the coefficients. \eqn{e^U} is a gamma variable with shape \eqn{k} and unit
#' rate, so the expected derivatives are polynomials in the central moments of
#' its logarithm, whose cumulants are the polygamma functions at \eqn{k}.
#'
#' @param distrib A `GenGamma2Distrib` object.
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mean`, `d` and `p`.
#' @param scale `"parameter"` or `"link"`; the link scale is applied by the
#'   generic.
#' @param expected Logical; for orders three and four, whether the expected
#'   derivative is returned.
#' @param approx,nsim Accepted for the generic's signature; the expectations
#'   are closed.
#' @param ...,threads Unused.
#'
#' @return A named list with one numeric vector per component: 3, 6, 10, 15
#'   and 21 components at orders one to five, the second order keyed
#'   `mean_mean`, `d_d`, `p_p`, `mean_d`, `mean_p`, `d_p`.
#'
#' @seealso [distrib_pdf.GenGamma2Distrib()].
#'
#' @examples
#' d <- gengamma2_distrib()
#' th <- list(mean = 5, d = 3, p = 1.5)
#' y <- c(1.2, 4, 9)
#' g <- distrib_gradient(d, y, th)
#' h <- 1e-6
#' up <- distrib_pdf(d, y, list(mean = 5, d = 3, p = 1.5 + h), log = TRUE)
#' dn <- distrib_pdf(d, y, list(mean = 5, d = 3, p = 1.5 - h), log = TRUE)
#' all.equal(g$p, (up - dn) / (2 * h), tolerance = 1e-6)
S7::method(distrib_gradient, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma2_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_hessian, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma2_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_deriv3, GenGamma2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) gengamma2_deriv3_expected_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  else gengamma2_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_deriv4, GenGamma2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) gengamma2_deriv4_expected_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  else gengamma2_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_deriv5, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, gengamma2_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]]),
               match.arg(scale))
}

S7::method(distrib_expected_hessian, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  gengamma2_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

register_dexpected(GenGamma2Distrib, function(d, y, th, k, t) {
  if (k == 1L) gengamma2_dexpected1_cpp(y, th[[1]], th[[2]], th[[3]])
  else gengamma2_dexpected2_cpp(y, th[[1]], th[[2]], th[[3]])
})


#' @title Generalized Gamma Derivatives in the Response, Mean
#' @name distrib_grad_y.GenGamma2Distrib
#' @aliases distrib_hess_y.GenGamma2Distrib distrib_deriv3_y.GenGamma2Distrib
#'   distrib_deriv4_y.GenGamma2Distrib distrib_cross_y.GenGamma2Distrib
#'   distrib_cross2_y.GenGamma2Distrib distrib_grad_y_hess.GenGamma2Distrib
#'   distrib_hess_y_hess.GenGamma2Distrib
#'
#' @description
#' The derivatives of the log-density in the response to order four, and the
#' mixed derivatives of orders one and two in the response and one and two in
#' \eqn{(m, d, p)}, each a closed form from its own compiled kernel.
#'
#' @param distrib A `GenGamma2Distrib` object.
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mean`, `d` and `p`.
#' @param scale `"parameter"` or `"link"`, for the mixed derivatives.
#' @param ... Unused.
#'
#' @return A numeric vector for the derivatives in the response; a named list,
#'   keyed by parameter or by parameter pair, for the mixed derivatives.
#'
#' @seealso [distrib_gradient.GenGamma2Distrib()].
#'
#' @examples
#' d <- gengamma2_distrib()
#' distrib_grad_y(d, c(1.2, 4), list(mean = 5, d = 3, p = 1.5))
S7::method(distrib_grad_y, GenGamma2Distrib) <- function(distrib, y, theta, ...) {
  gengamma2_dy1_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}
S7::method(distrib_hess_y, GenGamma2Distrib) <- function(distrib, y, theta, ...) {
  gengamma2_dy2_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}
S7::method(distrib_deriv3_y, GenGamma2Distrib) <- function(distrib, y, theta, ...) {
  gengamma2_dy3_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}
S7::method(distrib_deriv4_y, GenGamma2Distrib) <- function(distrib, y, theta, ...) {
  gengamma2_dy4_cpp(y, theta[[1]], theta[[2]], theta[[3]])$y
}
S7::method(distrib_cross_y, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma2_cross_y_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}
S7::method(distrib_cross2_y, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma2_cross2_y_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}
S7::method(distrib_grad_y_hess, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma2_grad_y_hess_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}
S7::method(distrib_hess_y_hess, GenGamma2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  gengamma2_hess_y_hess_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}


#' @title Generalized Gamma Moments, Mean
#' @name mean.GenGamma2Distrib
#' @aliases variance.GenGamma2Distrib skewness.GenGamma2Distrib
#'   kurtosis.GenGamma2Distrib
#'
#' @description
#' With \eqn{g_j = \Gamma((d+j)/p)\,\Gamma(d/p)^{j-1}/\Gamma((d+1)/p)^j}, the
#' ratio \eqn{E[Y^j]/E[Y]^j}, the mean is \eqn{m}, the variance
#' \eqn{m^2 (g_2 - 1)}, the skewness
#' \eqn{(g_3 - 3 g_2 + 2)/(g_2 - 1)^{3/2}} and the excess kurtosis
#' \eqn{(g_4 - 4 g_3 + 6 g_2 - 3)/(g_2 - 1)^2 - 3}.
#'
#' @details
#' These combinations of the ratios are not formed: towards the lognormal they
#' cancel to several orders. The central moments of \eqn{Y/m} and the fourth
#' cumulant come from the series of the cumulant generating function of
#' \eqn{\log Y}, in double-double arithmetic, as
#' [variance.GenGamma1Distrib()] describes; the skewness is
#' \eqn{\mu_3/\mu_2^{3/2}} and the excess kurtosis
#' \eqn{(\mu_4 - 3\mu_2^2)/\mu_2^2}.
#'
#' @param x A `GenGamma2Distrib` object.
#' @param theta A named list with components `mean`, `d` and `p`.
#' @param ... Unused.
#'
#' @return A numeric vector, one value per parameter row.
#'
#' @seealso [gengamma2_distrib()].
#'
#' @examples
#' d <- gengamma2_distrib()
#' th <- list(mean = 5, d = 3, p = 1.5)
#' c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
S7::method(mean, GenGamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + moment_const(theta, 3L, 0)
}

S7::method(variance, GenGamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]]^2 * gengamma_mu2_cpp(theta[[2]] / theta[[3]], 1 / theta[[3]]) +
    moment_const(theta, 3L, 0)
}

S7::method(skewness, GenGamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  k <- theta[[2]] / theta[[3]]
  h <- 1 / theta[[3]]
  gengamma_mu3_cpp(k, h) / gengamma_mu2_cpp(k, h)^1.5 + moment_const(theta, 3L, 0)
}

S7::method(kurtosis, GenGamma2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  k <- theta[[2]] / theta[[3]]
  h <- 1 / theta[[3]]
  gengamma_kappa4_cpp(k, h) / gengamma_mu2_cpp(k, h)^2 + moment_const(theta, 3L, 0)
}


#' Generalized Gamma Distribution in the Mean
#'
#' @description
#' Creates a generalized gamma distribution object whose first parameter is the
#' mean.
#'
#' @details
#' The Stacy parametrization of [gengamma1_distrib()] carries a scale, a shape
#' and a power, and exposes no mean at all, which is awkward for a family a
#' regression would put a linear predictor on. Since
#' \eqn{\mathbb{E}[Y] = a\,\Gamma((d+1)/p)/\Gamma(d/p)}, the scale here is
#' \deqn{a = m\,\dfrac{\Gamma(d/p)}{\Gamma((d+1)/p)},}
#' and every derivative in \eqn{(m, d, p)} to order five, the expected
#' information and its derivatives, and the derivatives in the response are
#' closed forms, each order in its own compiled kernel (see
#' [distrib_gradient.GenGamma2Distrib()]).
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{p\,y^{d-1}}{a^{d}\,\Gamma(d/p)}\,e^{-(y/a)^{p}}, \qquad a = m\,\frac{\Gamma(d/p)}{\Gamma((d+1)/p)}}
#' on \eqn{y \in (0, \infty)}, with \eqn{E[Y] = m}.
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_d Link function for the shape. Defaults to the log.
#' @param link_p Link function for the power. Defaults to the log.
#'
#' @return An S7 object of class `GenGamma2Distrib`, inheriting from
#'   `continuous_distrib`, with `params` `c("mean", "d", "p")` and
#'   `link_params` the three links given here.
#'
#' @seealso [gengamma1_distrib()]; [GenGamma2Distrib] for the class.
#'
#' @examples
#' d <- gengamma2_distrib()
#' theta <- list(mean = 5, d = 3, p = 1.5)
#' mean(d, theta)
#'
#' @importFrom linkfunctions7 log_link
#' @export
gengamma2_distrib <- function(link_mean = log_link(), link_d = log_link(),
                              link_p = log_link()) {
  GenGamma2Distrib(
    distrib_name = "gengamma2",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mean", "d", "p"),
    params_interpretation = c(mean = "mean", d = "shape", p = "power"),
    n_params = 3,
    params_bounds = list(mean = c(0, Inf), d = c(0, Inf), p = c(0, Inf)),
    link_params = list(mean = link_mean, d = link_d, p = link_p)
  )
}

# the distribution function's derivatives from the compiled kernels
# (compiled_cdf())
S7::method(distrib_grad_cdf, GenGamma2Distrib) <- compiled_grad_cdf
S7::method(distrib_hess_cdf, GenGamma2Distrib) <- compiled_hess_cdf
