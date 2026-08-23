#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Gaussian Distribution in Mean and Variance
#' @name Gaussian2Distrib
#'
#' @description A subclass of `continuous_distrib` for the Gaussian
#'   written in its mean and its **variance**.
#' @inheritParams distrib
#' @return An object of class `Gaussian2Distrib`.
#' @seealso [gaussian2_distrib()], [gaussian1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.Gaussian2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gaussian2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gaussian2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gaussian2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gaussian2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gaussian2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gaussian2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gaussian2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gaussian2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gaussian2Distrib],
#'   [`distrib_rng()`][distrib_rng.Gaussian2Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
Gaussian2Distrib <- S7::new_class("Gaussian2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Density in Mean and Variance
#' @name distrib_pdf.Gaussian2Distrib
#' @description
#' \deqn{f(y) = \dfrac{1}{\sqrt{2\pi\sigma^2}}
#'       \exp\left\{-\dfrac{(y-\mu)^2}{2\sigma^2}\right\}}
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_pdf, Gaussian2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(y, mean = theta[[1]], sd = sqrt(theta[[2]]), log = log)
}

#' @title Gaussian Distribution Function in Mean and Variance
#' @name distrib_cdf.Gaussian2Distrib
#' @description The normal distribution function, with the standard deviation
#'   taken as the square root of the variance parameter.
#' @param distrib A `Gaussian2Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with `mu` and `sigma2`.
#' @param lower.tail Logical; if `TRUE`, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, returns log-probabilities.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_cdf, Gaussian2Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  stats::pnorm(q, mean = theta[[1]], sd = sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Quantile Function in Mean and Variance
#' @name distrib_quantile.Gaussian2Distrib
#' @description The normal quantile function.
#' @param distrib A `Gaussian2Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with `mu` and `sigma2`.
#' @param lower.tail Logical; if `TRUE`, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, `p` is given as a log-probability.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_quantile, Gaussian2Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  stats::qnorm(p, mean = theta[[1]], sd = sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Random Generation in Mean and Variance
#' @name distrib_rng.Gaussian2Distrib
#' @description Delegates to [stats::rnorm()].
#' @param distrib A `Gaussian2Distrib` object.
#' @param n The number of draws.
#' @param theta A list with `mu` and `sigma2`.
#' @return A numeric vector.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_rng, Gaussian2Distrib) <- function(distrib, n, theta) {
  stats::rnorm(n, mean = theta[[1]], sd = sqrt(theta[[2]]))
}

#' @title Gaussian Analytical Gradient in Mean and Variance
#' @name distrib_gradient.Gaussian2Distrib
#' @description
#' With \eqn{r = y - \mu} and \eqn{v = \sigma^2},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \dfrac{r}{v}, \qquad
#'       \dfrac{\partial\ell}{\partial v} = \dfrac{r^2 - v}{2v^2}}
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_gradient, Gaussian2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian2_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Observed Hessian in Mean and Variance
#' @name distrib_hessian.Gaussian2Distrib
#' @description
#' \deqn{\ell^{(\mu\mu)} = -\dfrac{1}{v}, \qquad
#'       \ell^{(\mu v)} = -\dfrac{r}{v^2}, \qquad
#'       \ell^{(vv)} = \dfrac{1}{2v^2} - \dfrac{r^2}{v^3}}
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_hessian, Gaussian2Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian2_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Expected Hessian in Mean and Variance
#' @name distrib_expected_hessian.Gaussian2Distrib
#' @description
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\dfrac{1}{v}, \qquad
#'       \mathbb{E}[\ell^{(\mu v)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(vv)}] = -\dfrac{1}{2v^2}}
#' The two parameters are orthogonal, as they are in every parametrization of
#' this family.
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_expected_hessian, Gaussian2Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ..., threads = 1L) {
  gaussian2_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Third-Order Derivatives in Mean and Variance
#' @name distrib_deriv3.Gaussian2Distrib
#' @description Closed form, observed or expected.
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param expected Logical; if `TRUE`, returns the expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_deriv3, Gaussian2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian2_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian2_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian Fourth-Order Derivatives in Mean and Variance
#' @name distrib_deriv4.Gaussian2Distrib
#' @description Closed form, observed or expected.
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param expected Logical; if `TRUE`, returns the expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_deriv4, Gaussian2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  if (expected) {
    gaussian2_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  } else {
    gaussian2_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
  }
}

#' @title Gaussian Response Derivatives in Mean and Variance
#' @name distrib_grad_y.Gaussian2Distrib
#' @description \eqn{\partial\ell/\partial y = -(y-\mu)/v}.
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_grad_y, Gaussian2Distrib) <- function(distrib, y, theta, ...) {
  -(y - theta[[1]]) / theta[[2]]
}

#' @title Gaussian Second Response Derivative in Mean and Variance
#' @name distrib_hess_y.Gaussian2Distrib
#' @description \eqn{\partial^2\ell/\partial y^2 = -1/v}, free of the response.
#' @param distrib A `Gaussian2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `sigma2`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian2_distrib()]
S7::method(distrib_hess_y, Gaussian2Distrib) <- function(distrib, y, theta, ...) {
  rep(-1 / theta[[2]], length.out = length(y))
}


#' Gaussian Distribution in Mean and Variance
#'
#' @description
#' Creates a Gaussian distribution object parametrized by its mean and its
#' **variance**.
#'
#' @details
#' This is the same law as [gaussian1_distrib()] in different
#' coordinates: \eqn{\sigma^2} here is the square of the \eqn{\sigma} there.
#' The two are separate families rather than one family under a link, because a
#' link changes the scale a parameter is *modeled* on while leaving the
#' parameter what it was, whereas here the parameter, its interpretation, its
#' standard error and its confidence interval are all about the variance.
#'
#' The numbering follows the literature where it has one: this parametrization
#' is `NO2` in \pkg{gamlss}.
#'
#' Derivatives are closed form to fourth order, observed and expected, and the
#' two parameters are orthogonal.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{\sqrt{2\pi\sigma^{2}}}\exp\!\left\{-\frac{(y-\mu)^{2}}{2\sigma^{2}}\right\}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}}
#'
#' @param link_mu Link function for \eqn{\mu}. Defaults to the identity.
#' @param link_sigma2 Link function for \eqn{\sigma^2}. Defaults to the log.
#'
#' @return An S7 object of class [Gaussian2Distrib()].
#'
#' @seealso [gaussian1_distrib()], [gaussian3_distrib()]
#'
#' @examples
#' d <- gaussian2_distrib()
#' theta <- list(mu = 1, sigma2 = 4)
#' distrib_pdf(d, c(0, 1, 2), theta)
#' variance(d, theta)
#'
#' # the same law as gaussian1 with sigma = 2
#' distrib_pdf(gaussian1_distrib(), 0.5, list(mu = 1, sigma = 2))
#'
#' @export
gaussian2_distrib <- function(link_mu = identity_link(),
                              link_sigma2 = log_link()) {
  Gaussian2Distrib(
    distrib_name = "gaussian2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma2"),
    params_interpretation = c(mu = "mean", sigma2 = "variance"),
    n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
    link_params = list(mu = link_mu, sigma2 = link_sigma2)
  )
}
