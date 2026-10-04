#' @include distrib.R generics.R numerical_functions.R cdf_derivatives.R cdf_higher.R dexpected_families.R moments.R
NULL

#' @title Lognormal Distribution Class, Mean and Variance of Y
#' @name Lognormal2Distrib
#'
#' @description
#' The S7 class of the lognormal family parametrized by the mean \eqn{m} and
#' the variance \eqn{v} of \eqn{Y}. It inherits from `continuous_distrib`.
#' Build one with [lognormal2_distrib()]; this page documents the raw S7
#' constructor, which validates none of the relationships between the
#' properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Lognormal2Distrib`, inheriting from
#'   `continuous_distrib`. For an object built by [lognormal2_distrib()] the
#'   properties hold `"lognormal2"`, `"univariate"`, `c(0, Inf)`,
#'   `c("mean", "var")`, the interpretations `c(mean = "mean",
#'   var = "variance")`, `2` and the domains \eqn{(0, \infty)} for both.
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
#' @seealso [lognormal2_distrib()] to build one; [Lognormal1Distrib] for the
#'   parametrization on the log scale.
#'
#' @examples
#' d <- lognormal2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params
Lognormal2Distrib <- S7::new_class("Lognormal2Distrib",
                                   parent = continuous_distrib)

#' @title The Log-Scale Parameters of a Mean and a Variance
#' @name lognormal2_log_params
#' @description Returns the mean and the standard deviation of \eqn{\log Y}
#'   for the mean \eqn{m} and the variance \eqn{v} of \eqn{Y}:
#'   \eqn{S = \log(1 + v/m^2)}, formed by `log1p()`, `sdlog`
#'   \eqn{\sqrt{S}} and `meanlog` \eqn{\log m - S/2}.
#' @param theta A list with the mean and the variance, in that order.
#' @return A list with components `meanlog` and `sdlog`.
#' @seealso [distrib_pdf.Lognormal2Distrib()]
#' @keywords internal
lognormal2_log_params <- function(theta) {
  m <- theta[[1]]
  s <- log1p(theta[[2]] / m^2)
  list(meanlog = log(m) - s / 2, sdlog = sqrt(s))
}


#' @title Lognormal Density, Distribution and Generator in the Mean and
#'   Variance of Y
#' @name distrib_pdf.Lognormal2Distrib
#' @aliases distrib_cdf.Lognormal2Distrib distrib_quantile.Lognormal2Distrib
#'   distrib_rng.Lognormal2Distrib
#'
#' @description
#' The lognormal with \eqn{\log Y \sim N(\mu_l, S)}, where
#' \eqn{S = \log(1 + v/m^2)} and \eqn{\mu_l = \log m - S/2}, so that
#' \eqn{E[Y] = m} and \eqn{\operatorname{Var}(Y) = v}. The density, the
#' distribution function, the quantile function and the generator are
#' [stats::dlnorm()], [stats::plnorm()], [stats::qlnorm()] and
#' [stats::rlnorm()] at `meanlog` \eqn{\mu_l} and `sdlog` \eqn{\sqrt{S}}, with
#' \eqn{S} formed by `log1p()`.
#'
#' @param distrib A `Lognormal2Distrib` object, from [lognormal2_distrib()].
#' @param y,q A numeric vector of observations or quantiles.
#' @param p A numeric vector of probabilities.
#' @param n The number of draws.
#' @param theta A named list with components `mean` and `var`.
#' @param log,log.p,lower.tail As in [stats::dlnorm()] and its siblings.
#' @param ... Unused.
#'
#' @return A numeric vector.
#'
#' @seealso [lognormal2_distrib()].
#'
#' @examples
#' d <- lognormal2_distrib()
#' th <- list(mean = 3, var = 2)
#' integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf)$value
#' distrib_quantile(d, distrib_cdf(d, c(1, 3), th), th)
S7::method(distrib_pdf, Lognormal2Distrib) <- function(distrib, y, theta,
                                                       log = FALSE, ...) {
  lp <- lognormal2_log_params(theta)
  stats::dlnorm(y, meanlog = lp$meanlog, sdlog = lp$sdlog, log = log)
}

S7::method(distrib_cdf, Lognormal2Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  lp <- lognormal2_log_params(theta)
  stats::plnorm(q, meanlog = lp$meanlog, sdlog = lp$sdlog,
                lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_quantile, Lognormal2Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  lp <- lognormal2_log_params(theta)
  stats::qlnorm(p, meanlog = lp$meanlog, sdlog = lp$sdlog,
                lower.tail = lower.tail, log.p = log.p)
}

S7::method(distrib_rng, Lognormal2Distrib) <- function(distrib, n, theta, ...) {
  lp <- lognormal2_log_params(theta)
  stats::rlnorm(n, meanlog = lp$meanlog, sdlog = lp$sdlog)
}


#' @title Lognormal Derivatives in the Mean and Variance of Y
#' @name distrib_gradient.Lognormal2Distrib
#' @aliases distrib_hessian.Lognormal2Distrib distrib_deriv3.Lognormal2Distrib
#'   distrib_deriv4.Lognormal2Distrib distrib_deriv5.Lognormal2Distrib
#'   distrib_expected_hessian.Lognormal2Distrib
#'   distrib_dexpected_hessian.Lognormal2Distrib
#'   distrib_d2expected_hessian.Lognormal2Distrib
#'
#' @description
#' Return the derivatives of the log-density in \eqn{(m, v)} of orders one to
#' five, the expected information and its expected third and fourth
#' derivatives, and the first two derivatives of the expected information,
#' each from its own compiled kernel.
#'
#' @details
#' With \eqn{L = \log y}, \eqn{S = \log(1 + v/m^2)} and
#' \eqn{\mu_l = \log m - S/2} the log-density is
#' \deqn{\ell = -L - \tfrac12 \log 2\pi - \tfrac12 \log S -
#'       \frac{(L - \mu_l)^2}{2S}.}
#' Every component is a closed form in \eqn{(L, S, m, v)}, derived offline
#' and written out per component, one kernel per order. \eqn{\ell} is a
#' polynomial of degree two in \eqn{L}, so each expected derivative replaces
#' \eqn{L} and \eqn{L^2} by \eqn{E[L] = \mu_l} and \eqn{E[L^2] = \mu_l^2 + S}.
#'
#' @param distrib A `Lognormal2Distrib` object.
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mean` and `var`.
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
#'   `mean_mean`, `var_var`, `mean_var`.
#'
#' @seealso [distrib_pdf.Lognormal2Distrib()].
#'
#' @examples
#' d <- lognormal2_distrib()
#' th <- list(mean = 3, var = 2)
#' y <- c(0.7, 2.5, 6)
#' g <- distrib_gradient(d, y, th)
#' h <- 1e-6
#' up <- distrib_pdf(d, y, list(mean = 3 + h, var = 2), log = TRUE)
#' dn <- distrib_pdf(d, y, list(mean = 3 - h, var = 2), log = TRUE)
#' all.equal(g$mean, (up - dn) / (2 * h), tolerance = 1e-6)
S7::method(distrib_gradient, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal2_gradient_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_hessian, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal2_hessian_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_deriv3, Lognormal2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) lognormal2_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  else lognormal2_deriv3_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_deriv4, Lognormal2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) lognormal2_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  else lognormal2_deriv4_cpp(y, theta[[1]], theta[[2]])
}

S7::method(distrib_deriv5, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, lognormal2_deriv5_cpp(y, theta[[1]], theta[[2]]),
               match.arg(scale))
}

S7::method(distrib_expected_hessian, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  lognormal2_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

register_dexpected(Lognormal2Distrib, function(d, y, th, k, t) {
  if (k == 1L) lognormal2_dexpected1_cpp(y, th[[1]], th[[2]])
  else lognormal2_dexpected2_cpp(y, th[[1]], th[[2]])
})


#' @title Lognormal Derivatives in the Response, Mean and Variance of Y
#' @name distrib_grad_y.Lognormal2Distrib
#' @aliases distrib_hess_y.Lognormal2Distrib distrib_deriv3_y.Lognormal2Distrib
#'   distrib_deriv4_y.Lognormal2Distrib distrib_cross_y.Lognormal2Distrib
#'   distrib_cross2_y.Lognormal2Distrib distrib_grad_y_hess.Lognormal2Distrib
#'   distrib_hess_y_hess.Lognormal2Distrib
#'
#' @description
#' The derivatives of the log-density in the response to order four, and the
#' mixed derivatives of orders one and two in the response and one and two in
#' \eqn{(m, v)}, each a closed form from its own compiled kernel.
#'
#' @param distrib A `Lognormal2Distrib` object.
#' @param y A numeric vector of strictly positive observations.
#' @param theta A named list with components `mean` and `var`.
#' @param scale `"parameter"` or `"link"`, for the mixed derivatives.
#' @param ... Unused.
#'
#' @return A numeric vector for the derivatives in the response; a named list,
#'   keyed by parameter or by parameter pair, for the mixed derivatives.
#'
#' @seealso [distrib_gradient.Lognormal2Distrib()].
#'
#' @examples
#' d <- lognormal2_distrib()
#' th <- list(mean = 3, var = 2)
#' distrib_grad_y(d, c(0.7, 2.5), th)
S7::method(distrib_grad_y, Lognormal2Distrib) <- function(distrib, y, theta, ...) {
  lognormal2_dy1_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_hess_y, Lognormal2Distrib) <- function(distrib, y, theta, ...) {
  lognormal2_dy2_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_deriv3_y, Lognormal2Distrib) <- function(distrib, y, theta, ...) {
  lognormal2_dy3_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_deriv4_y, Lognormal2Distrib) <- function(distrib, y, theta, ...) {
  lognormal2_dy4_cpp(y, theta[[1]], theta[[2]])$y
}
S7::method(distrib_cross_y, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal2_cross_y_cpp(y, theta[[1]], theta[[2]])
}
S7::method(distrib_cross2_y, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal2_cross2_y_cpp(y, theta[[1]], theta[[2]])
}
S7::method(distrib_grad_y_hess, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal2_grad_y_hess_cpp(y, theta[[1]], theta[[2]])
}
S7::method(distrib_hess_y_hess, Lognormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal2_hess_y_hess_cpp(y, theta[[1]], theta[[2]])
}


#' @title Lognormal Distribution-Function Derivatives, Mean and Variance of Y
#' @name distrib_grad_cdf.Lognormal2Distrib
#' @aliases distrib_hess_cdf.Lognormal2Distrib
#'   distrib_deriv3_cdf.Lognormal2Distrib distrib_deriv4_cdf.Lognormal2Distrib
#'
#' @description
#' The derivatives of \eqn{F(q) = \Phi\{(\log q - \mu_l)/\sqrt{S}\}} in
#' \eqn{(m, v)} to order four, each a closed form from its own compiled
#' kernel, carried to the upper tail and the log scale by
#' [cdf_tail_scale()] and [cdf_scale_k()].
#'
#' @param distrib A `Lognormal2Distrib` object.
#' @param q A numeric vector of strictly positive quantiles.
#' @param theta A named list with components `mean` and `var`.
#' @param lower.tail,log As in [distrib_grad_cdf()].
#' @param ... Unused.
#'
#' @return A named list, one component per parameter or per multi-index.
#'
#' @seealso [distrib_cdf.Lognormal2Distrib()].
#'
#' @examples
#' d <- lognormal2_distrib()
#' distrib_grad_cdf(d, c(1, 3), list(mean = 3, var = 2), log = FALSE)
S7::method(distrib_grad_cdf, Lognormal2Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  d1 <- lognormal2_dcdf1_cpp(q, theta[[1]], theta[[2]])
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL,
                 lower.tail, log)
}

S7::method(distrib_hess_cdf, Lognormal2Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  d1 <- lognormal2_dcdf1_cpp(q, theta[[1]], theta[[2]])
  d2 <- lognormal2_dcdf2_cpp(q, theta[[1]], theta[[2]])
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, d2,
                 lower.tail, log)
}

S7::method(distrib_deriv3_cdf, Lognormal2Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  dF <- list(lognormal2_dcdf1_cpp(q, theta[[1]], theta[[2]]),
             lognormal2_dcdf2_cpp(q, theta[[1]], theta[[2]]),
             lognormal2_dcdf3_cpp(q, theta[[1]], theta[[2]]))
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), dF, 3L, lower.tail, log)
}

S7::method(distrib_deriv4_cdf, Lognormal2Distrib) <- function(
    distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
  dF <- list(lognormal2_dcdf1_cpp(q, theta[[1]], theta[[2]]),
             lognormal2_dcdf2_cpp(q, theta[[1]], theta[[2]]),
             lognormal2_dcdf3_cpp(q, theta[[1]], theta[[2]]),
             lognormal2_dcdf4_cpp(q, theta[[1]], theta[[2]]))
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), dF, 4L, lower.tail, log)
}


#' @title Lognormal Moments, Mean and Variance of Y
#' @name mean.Lognormal2Distrib
#' @aliases variance.Lognormal2Distrib skewness.Lognormal2Distrib
#'   kurtosis.Lognormal2Distrib
#'
#' @description
#' The mean is \eqn{m} and the variance \eqn{v}. With \eqn{w = 1 + v/m^2} the
#' skewness is \eqn{(w + 2)\sqrt{w - 1}} and the excess kurtosis
#' \eqn{w^4 + 2w^3 + 3w^2 - 6}.
#'
#' @param x A `Lognormal2Distrib` object.
#' @param theta A named list with components `mean` and `var`.
#' @param ... Unused.
#'
#' @return A numeric vector, one value per parameter row.
#'
#' @seealso [lognormal2_distrib()].
#'
#' @examples
#' d <- lognormal2_distrib()
#' th <- list(mean = 3, var = 2)
#' c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
S7::method(mean, Lognormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[1]] + moment_const(theta, 2L, 0)
}

S7::method(variance, Lognormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]] + moment_const(theta, 2L, 0)
}

S7::method(skewness, Lognormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  w <- 1 + theta[[2]] / theta[[1]]^2
  (w + 2) * sqrt(w - 1) + moment_const(theta, 2L, 0)
}

S7::method(kurtosis, Lognormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  w <- 1 + theta[[2]] / theta[[1]]^2
  w^4 + 2 * w^3 + 3 * w^2 - 6 + moment_const(theta, 2L, 0)
}


#' Lognormal Distribution in the Mean and Variance of Y
#'
#' @description
#' Creates a lognormal distribution object parametrized by the mean and the
#' variance of \eqn{Y} itself, rather than of \eqn{\log Y}.
#'
#' @details
#' The parameters of [lognormal1_distrib()] describe \eqn{\log Y}, so neither
#' of them is a moment of \eqn{Y}. Here they are, through
#' \deqn{\mu_{\log} = \log\dfrac{m^2}{\sqrt{v + m^2}}, \qquad
#'       \sigma^2_{\log} = \log\left(1 + \dfrac{v}{m^2}\right).}
#' Every derivative of the log-density in \eqn{(m, v)} to order five, the
#' expected information and its derivatives, the derivatives in the response
#' and the derivatives of the distribution function are closed forms, each
#' order in its own compiled kernel (see
#' [distrib_gradient.Lognormal2Distrib()]).
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{y\sqrt{2\pi s^{2}}}\exp\!\left\{-\frac{(\log y -
#'   \mu_l)^{2}}{2s^{2}}\right\}, \quad s^{2} = \log\!\left(1+\frac{v}{m^{2}}
#'   \right)\!, \; \mu_l = \log m - \frac{s^{2}}{2}}
#' on \eqn{y \in (0, \infty)}, with \eqn{E[Y] = m} and
#' \eqn{\operatorname{Var}(Y) = v}.
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_var Link function for the variance. Defaults to the log.
#'
#' @return An S7 object of class `Lognormal2Distrib`, inheriting from
#'   `continuous_distrib`, with `params` `c("mean", "var")` and `link_params`
#'   the two links given here.
#'
#' @seealso [lognormal1_distrib()]; [Lognormal2Distrib] for the class.
#'
#' @examples
#' d <- lognormal2_distrib()
#' theta <- list(mean = 3, var = 2)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' @importFrom linkfunctions7 log_link
#' @export
lognormal2_distrib <- function(link_mean = log_link(), link_var = log_link()) {
  Lognormal2Distrib(
    distrib_name = "lognormal2",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mean", "var"),
    params_interpretation = c(mean = "mean", var = "variance"),
    n_params = 2,
    params_bounds = list(mean = c(0, Inf), var = c(0, Inf)),
    link_params = list(mean = link_mean, var = link_var)
  )
}
