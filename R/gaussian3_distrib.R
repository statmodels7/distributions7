#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Gaussian Distribution in Mean and Precision
#' @name Gaussian3Distrib
#'
#' @description A subclass of `continuous_distrib` for the Gaussian
#'   written in its mean and its **precision**.
#' @inheritParams distrib
#' @return An object of class `Gaussian3Distrib`.
#' @seealso [gaussian3_distrib()], [gaussian1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.Gaussian3Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Gaussian3Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Gaussian3Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Gaussian3Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Gaussian3Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Gaussian3Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Gaussian3Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Gaussian3Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Gaussian3Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Gaussian3Distrib],
#'   [`distrib_rng()`][distrib_rng.Gaussian3Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
Gaussian3Distrib <- S7::new_class("Gaussian3Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Density in Mean and Precision
#' @name distrib_pdf.Gaussian3Distrib
#' @description
#' \deqn{f(y) = \sqrt{\dfrac{\tau}{2\pi}}
#'       \exp\left\{-\dfrac{\tau(y-\mu)^2}{2}\right\}}
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_pdf, Gaussian3Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(y, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]), log = log)
}

#' @title Gaussian Distribution Function in Mean and Precision
#' @name distrib_cdf.Gaussian3Distrib
#' @description The normal distribution function at standard deviation
#'   \eqn{1/\sqrt{\tau}}.
#' @param distrib A `Gaussian3Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with `mu` and `tau`.
#' @param lower.tail Logical; if `TRUE`, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, returns log-probabilities.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_cdf, Gaussian3Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  stats::pnorm(q, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Quantile Function in Mean and Precision
#' @name distrib_quantile.Gaussian3Distrib
#' @description The normal quantile function at standard deviation
#'   \eqn{1/\sqrt{\tau}}.
#' @param distrib A `Gaussian3Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with `mu` and `tau`.
#' @param lower.tail Logical; if `TRUE`, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, `p` is given as a log-probability.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_quantile, Gaussian3Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  stats::qnorm(p, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Random Generation in Mean and Precision
#' @name distrib_rng.Gaussian3Distrib
#' @description Delegates to [stats::rnorm()].
#' @param distrib A `Gaussian3Distrib` object.
#' @param n The number of draws.
#' @param theta A list with `mu` and `tau`.
#' @return A numeric vector.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_rng, Gaussian3Distrib) <- function(distrib, n, theta) {
  stats::rnorm(n, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]))
}

#' @title Gaussian Analytical Gradient in Mean and Precision
#' @name distrib_gradient.Gaussian3Distrib
#' @description
#' With \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \tau r, \qquad
#'       \dfrac{\partial\ell}{\partial\tau} = \dfrac{1}{2\tau} - \dfrac{r^2}{2}}
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_gradient, Gaussian3Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian3_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Observed Hessian in Mean and Precision
#' @name distrib_hessian.Gaussian3Distrib
#' @description
#' \deqn{\ell^{(\mu\mu)} = -\tau, \qquad \ell^{(\mu\tau)} = r, \qquad
#'       \ell^{(\tau\tau)} = -\dfrac{1}{2\tau^2}}
#' The mean block is free of the data here, as it is in every parametrization
#' of this family, and so is the precision block.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_hessian, Gaussian3Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  gaussian3_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Analytical Expected Hessian in Mean and Precision
#' @name distrib_expected_hessian.Gaussian3Distrib
#' @description
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\tau, \qquad
#'       \mathbb{E}[\ell^{(\mu\tau)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\tau\tau)}] = -\dfrac{1}{2\tau^2}}
#' Only the mixed entry differs from the observed Hessian, which is what makes
#' Fisher scoring and Newton's method take the same step on the parameter scale
#' here.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_expected_hessian, Gaussian3Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ..., threads = 1L) {
  gaussian3_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Third-Order Derivatives in Mean and Precision
#' @name distrib_deriv3.Gaussian3Distrib
#' @description
#' Closed form. Every third derivative is free of the response, so the observed
#' and the expected ones coincide and `expected` changes nothing.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param expected Logical; makes no difference here.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_deriv3, Gaussian3Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  gaussian3_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Fourth-Order Derivatives in Mean and Precision
#' @name distrib_deriv4.Gaussian3Distrib
#' @description
#' Closed form, and free of the response: the only non-zero component is
#' \eqn{\ell^{(\tau\tau\tau\tau)} = -3/\tau^4}.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param expected Logical; makes no difference here.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_deriv4, Gaussian3Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ..., threads = 1L) {
  gaussian3_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gaussian Response Derivatives in Mean and Precision
#' @name distrib_grad_y.Gaussian3Distrib
#' @description \eqn{\partial\ell/\partial y = -\tau(y-\mu)}.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_grad_y, Gaussian3Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * (y - theta[[1]])
}

#' @title Gaussian Second Response Derivative in Mean and Precision
#' @name distrib_hess_y.Gaussian3Distrib
#' @description \eqn{\partial^2\ell/\partial y^2 = -\tau}, free of the response.
#' @param distrib A `Gaussian3Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `mu` and `tau`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [gaussian3_distrib()]
S7::method(distrib_hess_y, Gaussian3Distrib) <- function(distrib, y, theta, ...) {
  rep(-theta[[2]], length.out = length(y))
}


#' Gaussian Distribution in Mean and Precision
#'
#' @description
#' Creates a Gaussian distribution object parametrized by its mean and its
#' **precision** \eqn{\tau = 1/\sigma^2}.
#'
#' @details
#' This is the same law as [gaussian1_distrib()] in different
#' coordinates, and a separate family for the same reason
#' [gaussian2_distrib()] is: the parameter here *is* the
#' precision, and that is what the estimate, the standard error and the
#' interval describe.
#'
#' It is the flattest of the three parametrizations. Every third derivative is
#' free of the response, so the observed and the expected ones coincide, and
#' the only non-zero fourth derivative is
#' \eqn{\ell^{(\tau\tau\tau\tau)} = -3/\tau^4}.
#'
#' The precision is the parametrization a Bayesian conjugate analysis uses,
#' the gamma being conjugate for \eqn{\tau} at known \eqn{\mu}.
#'
#' @section The distribution:
#' \deqn{f(y) = \sqrt{\frac{\tau}{2\pi}}\exp\!\left\{-\frac{\tau(y-\mu)^{2}}{2}\right\}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = 1/\tau}
#'
#' @param link_mu Link function for \eqn{\mu}. Defaults to the identity.
#' @param link_tau Link function for \eqn{\tau}. Defaults to the log.
#'
#' @return An S7 object of class [Gaussian3Distrib()].
#'
#' @seealso [gaussian1_distrib()], [gaussian2_distrib()]
#'
#' @examples
#' d <- gaussian3_distrib()
#' theta <- list(mu = 1, tau = 0.25)
#' distrib_pdf(d, c(0, 1, 2), theta)
#' variance(d, theta)
#'
#' # the same law as gaussian1 with sigma = 2
#' distrib_pdf(gaussian1_distrib(), 0.5, list(mu = 1, sigma = 2))
#'
#' @export
gaussian3_distrib <- function(link_mu = identity_link(),
                              link_tau = log_link()) {
  Gaussian3Distrib(
    distrib_name = "gaussian3",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "tau"),
    params_interpretation = c(mu = "mean", tau = "precision"),
    n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), tau = c(0, Inf)),
    link_params = list(mu = link_mu, tau = link_tau)
  )
}
