#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Inverse Gaussian in Mean and Shape
#' @name InvGauss2Distrib
#'
#' @description A subclass of \code{continuous_distrib} for the inverse
#'   gaussian in its classical parametrisation, the mean and the shape
#'   \eqn{\lambda}.
#' @inheritParams distrib
#' @return An object of class \code{InvGauss2Distrib}.
#' @seealso \code{\link{invgauss2_distrib}}, \code{\link{invgauss1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.InvGauss2Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.InvGauss2Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.InvGauss2Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.InvGauss2Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.InvGauss2Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.InvGauss2Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.InvGauss2Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.InvGauss2Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.InvGauss2Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
InvGauss2Distrib <- S7::new_class("InvGauss2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Inverse Gaussian Density in Mean and Shape
#' @name distrib_pdf.InvGauss2Distrib
#' @description
#' \deqn{f(y) = \sqrt{\dfrac{\lambda}{2\pi y^3}}
#'       \exp\left\{-\dfrac{\lambda(y-\mu)^2}{2\mu^2 y}\right\}}
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_pdf, InvGauss2Distrib) <- function(distrib, y, theta, log = FALSE) {
  statmod::dinvgauss(y, mean = theta[[1]], dispersion = 1 / theta[[2]], log = log)
}

#' @title Inverse Gaussian Distribution Function in Mean and Shape
#' @name distrib_cdf.InvGauss2Distrib
#' @description The inverse gaussian distribution function at dispersion
#'   \eqn{1/\lambda}.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns log-probabilities.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_cdf, InvGauss2Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  statmod::pinvgauss(q, mean = theta[[1]], dispersion = 1 / theta[[2]],
                     lower.tail = lower.tail, log.p = log.p)
}

#' @title Inverse Gaussian Quantile Function in Mean and Shape
#' @name distrib_quantile.InvGauss2Distrib
#' @description The inverse gaussian quantile function at dispersion
#'   \eqn{1/\lambda}.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is a log-probability.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_quantile, InvGauss2Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  statmod::qinvgauss(p, mean = theta[[1]], dispersion = 1 / theta[[2]],
                     lower.tail = lower.tail, log.p = log.p)
}

#' @title Inverse Gaussian Random Generation in Mean and Shape
#' @name distrib_rng.InvGauss2Distrib
#' @description Delegates to \code{\link[statmod]{rinvgauss}}.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param n The number of draws.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @return A numeric vector.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_rng, InvGauss2Distrib) <- function(distrib, n, theta) {
  statmod::rinvgauss(n, mean = theta[[1]], dispersion = 1 / theta[[2]])
}

#' @title Inverse Gaussian Analytical Gradient in Mean and Shape
#' @name distrib_gradient.InvGauss2Distrib
#' @description
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \dfrac{\lambda(y-\mu)}{\mu^3},
#'       \qquad
#'       \dfrac{\partial\ell}{\partial\lambda}
#'         = \dfrac{1}{2\lambda} - \dfrac{(y-\mu)^2}{2\mu^2 y}}
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_gradient, InvGauss2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  invgauss2_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Inverse Gaussian Analytical Observed Hessian in Mean and Shape
#' @name distrib_hessian.InvGauss2Distrib
#' @description
#' \deqn{\ell^{(\mu\mu)} = \dfrac{\lambda(2\mu-3y)}{\mu^4}, \qquad
#'       \ell^{(\mu\lambda)} = \dfrac{y-\mu}{\mu^3}, \qquad
#'       \ell^{(\lambda\lambda)} = -\dfrac{1}{2\lambda^2}}
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_hessian, InvGauss2Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  invgauss2_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Inverse Gaussian Analytical Expected Hessian in Mean and Shape
#' @name distrib_expected_hessian.InvGauss2Distrib
#' @description
#' The response enters every derivative linearly, so the expectations need only
#' \eqn{\mathbb{E}[Y] = \mu}:
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\dfrac{\lambda}{\mu^3}, \qquad
#'       \mathbb{E}[\ell^{(\mu\lambda)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\lambda\lambda)}] = -\dfrac{1}{2\lambda^2}}
#' The mean and the shape are orthogonal.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_expected_hessian, InvGauss2Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  invgauss2_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Inverse Gaussian Third-Order Derivatives in Mean and Shape
#' @name distrib_deriv3.InvGauss2Distrib
#' @description Closed form, observed or expected.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param expected Logical; if \code{TRUE}, returns the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_deriv3, InvGauss2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ...) {
  invgauss2_deriv3_cpp(y, theta[[1]], theta[[2]], expected)
}

#' @title Inverse Gaussian Fourth-Order Derivatives in Mean and Shape
#' @name distrib_deriv4.InvGauss2Distrib
#' @description Closed form, observed or expected.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{lambda}.
#' @param expected Logical; if \code{TRUE}, returns the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{invgauss2_distrib}}
S7::method(distrib_deriv4, InvGauss2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ...) {
  invgauss2_deriv4_cpp(y, theta[[1]], theta[[2]], expected)
}


#' Inverse Gaussian Distribution in Mean and Shape
#'
#' @description
#' Creates an inverse gaussian distribution object in the classical
#' parametrisation, the mean \eqn{\mu} and the shape \eqn{\lambda}, with
#' \eqn{\operatorname{Var}(Y) = \mu^3/\lambda}.
#'
#' @details
#' The same law as \code{\link{invgauss1_distrib}}, which carries a dispersion
#' \eqn{\phi = 1/\lambda}. The map between them is one coordinate at a time,
#' the mean being untouched, which is what makes both sets of derivatives
#' elementary.
#'
#' The response enters every derivative linearly, so every expectation needs
#' only \eqn{\mathbb{E}[Y] = \mu} and all four orders are closed form. The mean
#' and the shape are orthogonal.
#'
#' @param link_mu Link function for \eqn{\mu}. Defaults to the log.
#' @param link_lambda Link function for \eqn{\lambda}. Defaults to the log.
#'
#' @return An S7 object of class \code{\link{InvGauss2Distrib}}.
#'
#' @seealso \code{\link{invgauss1_distrib}}
#'
#' @examples
#' d <- invgauss2_distrib()
#' theta <- list(mu = 2, lambda = 3)
#' distrib_pdf(d, c(1, 2, 3), theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' # the same law as invgauss1 with phi = 1 / lambda
#' distrib_pdf(invgauss1_distrib(), 2, list(mu = 2, phi = 1 / 3))
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
