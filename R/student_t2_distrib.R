#' @include distrib.R generics.R numerical_functions.R cdf_derivatives.R cdf_derivatives_families.R cdf_mapped_higher.R dexpected_families.R moments.R y_higher.R
NULL

#' @title Student t Distribution Class, Standard Deviation
#' @name StudentT2Distrib
#'
#' @description
#' The S7 class of the Student t family parametrized by its location
#' \eqn{\mu}, its standard deviation \eqn{\sigma} and its degrees of freedom
#' \eqn{\nu > 2}. It inherits from `continuous_distrib`. Build one with
#' [student_t2_distrib()]; this page documents the raw S7 constructor, which
#' validates none of the relationships between the properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `StudentT2Distrib`, inheriting from
#'   `continuous_distrib`. For an object built by [student_t2_distrib()] the
#'   properties hold `"student t2"`, `"univariate"`, `c(-Inf, Inf)`,
#'   `c("mu", "sigma", "nu")`, the interpretations `c(mu = "location",
#'   sigma = "standard deviation", nu = "degrees of freedom")`, `3` and the
#'   domains \eqn{\mathbb{R}}, \eqn{(0, \infty)} and \eqn{(2, \infty)}.
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
#' @seealso [student_t2_distrib()] to build one; [StudentT1Distrib] for the
#'   parametrization by the scale.
#'
#' @examples
#' d <- student_t2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params
StudentT2Distrib <- S7::new_class("StudentT2Distrib", parent = continuous_distrib)

#' @title The Scale of a Student t with a Given Standard Deviation
#' @name student_t2_scale
#' @description Returns \eqn{s_0 = \sigma\sqrt{1 - 2/\nu}}, the scale of the
#'   Student t whose standard deviation is \eqn{\sigma}.
#' @param theta A list with the location, the standard deviation and the
#'   degrees of freedom, in that order.
#' @return A numeric vector of scales.
#' @seealso [distrib_pdf.StudentT2Distrib()]
#' @keywords internal
student_t2_scale <- function(theta) {
  theta[[2]] * sqrt(1 - 2 / theta[[3]])
}


#' @title Student t Density, Distribution and Generator in the Standard Deviation
#' @name distrib_pdf.StudentT2Distrib
#' @aliases distrib_cdf.StudentT2Distrib distrib_quantile.StudentT2Distrib
#'   distrib_rng.StudentT2Distrib
#'
#' @description
#' The Student t with \eqn{\nu} degrees of freedom, location \eqn{\mu} and
#' scale \eqn{s_0 = \sigma\sqrt{1 - 2/\nu}}, so that the standard deviation is
#' \eqn{\sigma}. The density, the distribution function, the quantile function
#' and the generator are [stats::dt()], [stats::pt()], [stats::qt()] and
#' [stats::rt()] at that scale.
#'
#' @param distrib A `StudentT2Distrib` object, from [student_t2_distrib()].
#' @param y,q A numeric vector of observations or quantiles.
#' @param p A numeric vector of probabilities.
#' @param n The number of draws.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param log,log.p,lower.tail As in [stats::dt()] and its siblings.
#' @param ... Unused.
#'
#' @return A numeric vector.
#'
#' @seealso [student_t2_distrib()].
#'
#' @examples
#' d <- student_t2_distrib()
#' th <- list(mu = 1, sigma = 2, nu = 6)
#' integrate(function(t) (t - 1)^2 * distrib_pdf(d, t, th), -Inf, Inf)$value
S7::method(distrib_pdf, StudentT2Distrib) <- function(distrib, y, theta,
                                                      log = FALSE, ...) {
  s0 <- student_t2_scale(theta)
  val <- stats::dt((y - theta[[1]]) / s0, df = theta[[3]], log = log)
  if (log) val - log(s0) else val / s0
}

S7::method(distrib_cdf, StudentT2Distrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  stats::pt((q - theta[[1]]) / student_t2_scale(theta), df = theta[[3]],
            lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_quantile, StudentT2Distrib) <- function(distrib, p, theta,
                                                           lower.tail = TRUE,
                                                           log.p = FALSE, ...) {
  theta[[1]] + student_t2_scale(theta) *
    stats::qt(p, df = theta[[3]], lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_rng, StudentT2Distrib) <- function(distrib, n, theta, ...) {
  theta[[1]] + student_t2_scale(theta) * stats::rt(n, df = theta[[3]])
}


#' @title Student t Derivatives in the Standard Deviation
#' @name distrib_gradient.StudentT2Distrib
#' @aliases distrib_hessian.StudentT2Distrib distrib_deriv3.StudentT2Distrib
#'   distrib_deriv4.StudentT2Distrib distrib_deriv5.StudentT2Distrib
#'   distrib_expected_hessian.StudentT2Distrib
#'   distrib_dexpected_hessian.StudentT2Distrib
#'   distrib_d2expected_hessian.StudentT2Distrib
#'
#' @description
#' Return the derivatives of the log-density in \eqn{(\mu, \sigma, \nu)} of
#' orders one to five, the expected information and its expected third and
#' fourth derivatives, and the first two derivatives of the expected
#' information, each from its own compiled kernel.
#'
#' @details
#' With \eqn{z = (y - \mu)/\sigma}, \eqn{k = \nu - 2} and \eqn{q = z^2/k}, the
#' log-density is
#' \deqn{\ell = c(\nu) - \log\sigma - \tfrac12\log\pi - \frac{k + 3}{2}\log(1 + q),
#'   \qquad c(\nu) = \log\Gamma\!\left(\frac{\nu+1}{2}\right) -
#'   \log\Gamma\!\left(\frac{\nu}{2}\right) - \tfrac12\log(\nu - 2).}
#' Every component is a closed form derived offline and written out per
#' component, one kernel per order. The part in the data is a rational
#' function of \eqn{z} and \eqn{k}, reduced symbolically and written in
#' \eqn{q}, \eqn{1/k} and \eqn{t = 1/(1 + q)}; the score in \eqn{\nu} carries
#' the logarithm through \eqn{q/(1+q) - \log(1 + q)}. The derivatives of
#' \eqn{c(\nu)}, and the expected derivatives in \eqn{\nu}, are evaluated from
#' the polygamma functions below \eqn{\nu = 20} and from their asymptotic
#' series in \eqn{1/\nu} above it, because their terms cancel to leading
#' order as \eqn{\nu} grows. The expectations are closed forms: under the
#' model \eqn{t} follows a beta distribution with parameters \eqn{\nu/2} and
#' \eqn{1/2}.
#'
#' @param distrib A `StudentT2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param scale `"parameter"` or `"link"`; the link scale is applied by the
#'   generic.
#' @param expected Logical; for orders three and four, whether the expected
#'   derivative is returned.
#' @param approx,nsim Accepted for the generic's signature; the expectations
#'   are closed.
#' @param ... Unused.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list with one numeric vector per component: 3, 6, 10, 15
#'   and 21 components at orders one to five, the second order keyed
#'   `mu_mu`, `sigma_sigma`, `nu_nu`, `mu_sigma`, `mu_nu`, `sigma_nu`.
#'
#' @seealso [distrib_pdf.StudentT2Distrib()].
#'
#' @examples
#' d <- student_t2_distrib()
#' th <- list(mu = 1, sigma = 2, nu = 6)
#' y <- c(-2, 0.5, 4)
#' g <- distrib_gradient(d, y, th)
#' h <- 1e-6
#' up <- distrib_pdf(d, y, list(mu = 1, sigma = 2, nu = 6 + h), log = TRUE)
#' dn <- distrib_pdf(d, y, list(mu = 1, sigma = 2, nu = 6 - h), log = TRUE)
#' all.equal(g$nu, (up - dn) / (2 * h), tolerance = 1e-6)
S7::method(distrib_gradient, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  student_t2_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

S7::method(distrib_hessian, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  student_t2_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

S7::method(distrib_deriv3, StudentT2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) student_t2_deriv3_expected_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  else student_t2_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

S7::method(distrib_deriv4, StudentT2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) student_t2_deriv4_expected_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  else student_t2_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

S7::method(distrib_deriv5, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, student_t2_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads),
               match.arg(scale))
}

S7::method(distrib_expected_hessian, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  student_t2_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

register_dexpected(StudentT2Distrib, function(d, y, th, k, t) {
  if (k == 1L) student_t2_dexpected1_cpp(y, th[[1]], th[[2]], th[[3]], t)
  else student_t2_dexpected2_cpp(y, th[[1]], th[[2]], th[[3]], t)
})


#' @title Student t Derivatives in the Response, Standard Deviation
#' @name distrib_grad_y.StudentT2Distrib
#' @aliases distrib_hess_y.StudentT2Distrib distrib_deriv3_y.StudentT2Distrib
#'   distrib_deriv4_y.StudentT2Distrib distrib_cross_y.StudentT2Distrib
#'   distrib_cross2_y.StudentT2Distrib distrib_grad_y_hess.StudentT2Distrib
#'   distrib_hess_y_hess.StudentT2Distrib
#'
#' @description
#' The derivatives of the log-density in the response to order four, and the
#' mixed derivatives of orders one and two in the response and one and two in
#' \eqn{(\mu, \sigma, \nu)}, each a closed form from its own compiled kernel.
#'
#' @param distrib A `StudentT2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param scale `"parameter"` or `"link"`, for the mixed derivatives.
#' @param ... Unused.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A numeric vector for the derivatives in the response; a named list,
#'   keyed by parameter or by parameter pair, for the mixed derivatives.
#'
#' @seealso [distrib_gradient.StudentT2Distrib()].
#'
#' @examples
#' d <- student_t2_distrib()
#' distrib_grad_y(d, c(-2, 0.5), list(mu = 1, sigma = 2, nu = 6))
S7::method(distrib_grad_y, StudentT2Distrib) <- function(distrib, y, theta, ..., threads = 1L) {
  student_t2_dy1_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)$y
}
S7::method(distrib_hess_y, StudentT2Distrib) <- function(distrib, y, theta, ..., threads = 1L) {
  student_t2_dy2_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)$y
}
S7::method(distrib_deriv3_y, StudentT2Distrib) <- function(distrib, y, theta, ..., threads = 1L) {
  student_t2_dy3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)$y
}
S7::method(distrib_deriv4_y, StudentT2Distrib) <- function(distrib, y, theta, ..., threads = 1L) {
  student_t2_dy4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)$y
}
S7::method(distrib_cross_y, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  student_t2_cross_y_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}
S7::method(distrib_cross2_y, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  student_t2_cross2_y_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}
S7::method(distrib_grad_y_hess, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  student_t2_grad_y_hess_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}
S7::method(distrib_hess_y_hess, StudentT2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  student_t2_hess_y_hess_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}


#' @title Student t Distribution-Function Derivatives, Standard Deviation
#' @name distrib_grad_cdf.StudentT2Distrib
#' @aliases distrib_hess_cdf.StudentT2Distrib
#'   distrib_deriv3_cdf.StudentT2Distrib distrib_deriv4_cdf.StudentT2Distrib
#'
#' @description
#' The derivatives of \eqn{F(q)} in \eqn{(\mu, \sigma, \nu)} to order four.
#' The family is location-scale in \eqn{(\mu, \sigma)} at fixed \eqn{\nu}, so
#' those components are closed forms in the density and its derivatives; the
#' components in \eqn{\nu}, derivatives of an incomplete beta function in its
#' parameter, have no elementary form and are differenced, as for
#' [student_t1_distrib()].
#'
#' @param distrib A `StudentT2Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param lower.tail,log As in [distrib_grad_cdf()].
#' @param ... Unused.
#'
#' @return A named list, one component per parameter or per multi-index.
#'
#' @seealso [partial_loc_scale_grad_cdf()], [partial_loc_scale_hess_cdf()]
#'   and [partial_loc_scale_deriv_cdf_k()], the bodies registered here.
#'
#' @examples
#' d <- student_t2_distrib()
#' distrib_grad_cdf(d, c(-2, 0.5), list(mu = 1, sigma = 2, nu = 6), log = FALSE)
S7::method(distrib_grad_cdf, StudentT2Distrib) <- partial_loc_scale_grad_cdf
S7::method(distrib_hess_cdf, StudentT2Distrib) <- partial_loc_scale_hess_cdf
S7::method(distrib_deriv3_cdf, StudentT2Distrib) <- partial_loc_scale_deriv_cdf_k(3L)
S7::method(distrib_deriv4_cdf, StudentT2Distrib) <- partial_loc_scale_deriv_cdf_k(4L)


#' @title Student t Moments, Standard Deviation
#' @name mean.StudentT2Distrib
#' @aliases variance.StudentT2Distrib skewness.StudentT2Distrib
#'   kurtosis.StudentT2Distrib
#'
#' @description
#' The mean is \eqn{\mu} and the variance \eqn{\sigma^2}, both finite on the
#' whole domain \eqn{\nu > 2}. The skewness is zero for \eqn{\nu > 3} and
#' `NaN` otherwise; the excess kurtosis is \eqn{6/(\nu - 4)} for
#' \eqn{\nu > 4} and `Inf` otherwise.
#'
#' @param x A `StudentT2Distrib` object.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param ... Unused.
#'
#' @return A numeric vector, one value per parameter row.
#'
#' @seealso [student_t2_distrib()].
#'
#' @examples
#' d <- student_t2_distrib()
#' th <- list(mu = 1, sigma = 2, nu = 6)
#' c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
S7::method(mean, StudentT2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + moment_const(theta, 3L, 0)
}

S7::method(variance, StudentT2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 3L, 0)
}

S7::method(skewness, StudentT2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  nu <- theta[[3]] + moment_const(theta, 3L, 0)
  ifelse(nu > 3, 0, NaN)
}

S7::method(kurtosis, StudentT2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  nu <- theta[[3]] + moment_const(theta, 3L, 0)
  ifelse(nu > 4, 6 / (nu - 4), Inf)
}


#' Student t Distribution in the Standard Deviation
#'
#' @description
#' Creates a Student t distribution object whose second parameter is the
#' standard deviation rather than the scale.
#'
#' @details
#' The scale of [student_t1_distrib()] is not the standard deviation: the two
#' differ by \eqn{\sqrt{\nu/(\nu-2)}}. Here the scale is
#' \deqn{s_0 = \sigma\sqrt{\dfrac{\nu-2}{\nu}},}
#' which exists only for \eqn{\nu > 2}, so the degrees of freedom are bounded
#' below at two. Every derivative in \eqn{(\mu, \sigma, \nu)} to order five,
#' the expected information and its derivatives, and the derivatives in the
#' response are closed forms, each order in its own compiled kernel (see
#' [distrib_gradient.StudentT2Distrib()]). This is `TF2` in gamlss.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{s_0}\,t_{\nu}\!\left(\frac{y-\mu}{s_0}\right), \qquad s_0 = \sigma\sqrt{\frac{\nu-2}{\nu}}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}}
#'
#' @param link_mu Link function for the location. Defaults to the identity.
#' @param link_sigma Link function for the standard deviation. Defaults to the
#'   log.
#' @param link_nu Link function for the degrees of freedom. Defaults to a link
#'   bounded below at two.
#'
#' @return An S7 object of class `StudentT2Distrib`, inheriting from
#'   `continuous_distrib`, with `params` `c("mu", "sigma", "nu")` and
#'   `link_params` the three links given here.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. *Journal of the Royal Statistical
#' Society, Series C* 54, 507-554.
#'
#' @seealso [student_t1_distrib()]; [StudentT2Distrib] for the class.
#'
#' @examples
#' d <- student_t2_distrib()
#' theta <- list(mu = 0, sigma = 2, nu = 8)
#' variance(d, theta)
#'
#' @importFrom linkfunctions7 identity_link log_link bounded_link
#' @export
student_t2_distrib <- function(link_mu = identity_link(),
                               link_sigma = log_link(),
                               link_nu = bounded_link(lwr = 2)) {
  StudentT2Distrib(
    distrib_name = "student t2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma", "nu"),
    params_interpretation = c(mu = "location", sigma = "standard deviation",
                              nu = "degrees of freedom"),
    n_params = 3,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf), nu = c(2, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma, nu = link_nu)
  )
}
