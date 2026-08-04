#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Exponential Distribution
#' @name ExponentialDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the
#'   exponential distribution in its mean parametrisation.
#' @inheritParams distrib
#' @return An object of class \code{ExponentialDistrib}.
#' @seealso \code{\link{exponential_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.ExponentialDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.ExponentialDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.ExponentialDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.ExponentialDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.ExponentialDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_grad_y.ExponentialDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_hessian.ExponentialDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_hess_y.ExponentialDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_pdf.ExponentialDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.ExponentialDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.ExponentialDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
ExponentialDistrib <- S7::new_class("ExponentialDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Exponential Density
#' @name distrib_pdf.ExponentialDistrib
#' @description
#' \deqn{f(y; \mu) = \dfrac{1}{\mu} e^{-y/\mu}, \qquad y > 0}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_pdf, ExponentialDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dexp(y, rate = 1 / theta[[1]], log = log)
}

#' @title Exponential Distribution Function
#' @name distrib_cdf.ExponentialDistrib
#' @description
#' \deqn{F(q; \mu) = 1 - e^{-q/\mu}}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_cdf, ExponentialDistrib) <- function(distrib, q, theta,
                                                        lower.tail = TRUE,
                                                        log.p = FALSE) {
  stats::pexp(q, rate = 1 / theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Exponential Quantile Function
#' @name distrib_quantile.ExponentialDistrib
#' @description
#' \deqn{Q(p; \mu) = -\mu \log(1 - p)}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), \code{p} is \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_quantile, ExponentialDistrib) <- function(distrib, p, theta,
                                                             lower.tail = TRUE,
                                                             log.p = FALSE) {
  stats::qexp(p, rate = 1 / theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Exponential Random Generation
#' @name distrib_rng.ExponentialDistrib
#' @description Draws from the exponential distribution through
#'   \code{\link[stats]{rexp}}.
#' @param distrib An \code{ExponentialDistrib} object.
#' @param n The number of draws.
#' @param theta A list containing the parameter \code{mu}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_rng, ExponentialDistrib) <- function(distrib, n, theta) {
  stats::rexp(n, rate = 1 / theta[[1]])
}

#' @title Exponential Analytical Gradient
#' @name distrib_gradient.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}}
#' the score of a one-parameter family written as the deviation from the mean
#' over the variance.
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu} component.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_gradient, ExponentialDistrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  exponential_gradient_cpp(y, theta[[1]])
}

#' @title Exponential Analytical Observed Hessian
#' @name distrib_hessian.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} - \dfrac{2y}{\mu^3}}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu} component.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_hessian, ExponentialDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  exponential_hessian_cpp(y, theta[[1]])
}

#' @title Exponential Analytical Expected Hessian
#' @name distrib_expected_hessian.ExponentialDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu^2}}
#' obtained from the observed form by \eqn{\mathbb{E}[y] = \mu}.
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu} component.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_expected_hessian, ExponentialDistrib) <- function(distrib, y, theta,
                                                                     scale = c("parameter", "link"),
                                                                     approx = c("bartlett", "integrate", "mc", "opg"),
                                                                     nsim = 10000, ...) {
  exponential_expected_hessian_cpp(y, theta[[1]])
}

#' @title Exponential Analytical Third-Order Derivative
#' @name distrib_deriv3.ExponentialDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu)} = -\dfrac{2}{\mu^3} + \dfrac{6y}{\mu^4},
#'       \qquad \mathbb{E}[\ell^{(\mu\mu\mu)}] = \dfrac{4}{\mu^3}}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; if \code{TRUE}, returns the expected derivative.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu_mu} component.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_deriv3, ExponentialDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) exponential_deriv3_expected_cpp(y, theta[[1]])
  else exponential_deriv3_cpp(y, theta[[1]])
}

#' @title Exponential Analytical Fourth-Order Derivative
#' @name distrib_deriv4.ExponentialDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu\mu)} = \dfrac{6}{\mu^4} - \dfrac{24y}{\mu^5},
#'       \qquad \mathbb{E}[\ell^{(\mu\mu\mu\mu)}] = -\dfrac{18}{\mu^4}}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; if \code{TRUE}, returns the expected derivative.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu_mu_mu} component.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_deriv4, ExponentialDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) exponential_deriv4_expected_cpp(y, theta[[1]])
  else exponential_deriv4_cpp(y, theta[[1]])
}

#' @title Exponential Response Gradient
#' @name distrib_grad_y.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{1}{\mu}}
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_grad_y, ExponentialDistrib) <- function(distrib, y, theta, ...) {
  rep(-1 / theta[[1]], length.out = length(y))
}

#' @title Exponential Response Hessian
#' @name distrib_hess_y.ExponentialDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = 0}
#' the log-density being linear in the response.
#' @param distrib An \code{ExponentialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param ... Unused.
#' @return A numeric vector of zeros.
#' @seealso \code{\link{exponential_distrib}}
S7::method(distrib_hess_y, ExponentialDistrib) <- function(distrib, y, theta, ...) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Exponential Distribution Object
#'
#' @description
#' Creates a distribution object for the exponential distribution parametrised
#' by its mean \eqn{\mu}.
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The exponential distribution models waiting times on \eqn{y > 0}, with mean
#' \eqn{\mu} and variance \eqn{\mu^2}: the coefficient of variation is one, and
#' fixing it is what distinguishes the family from the Gamma.
#'
#' \strong{Density:} \deqn{f(y; \mu) = \dfrac{1}{\mu} e^{-y/\mu}}
#'
#' \strong{Distribution function:} \deqn{F(q; \mu) = 1 - e^{-q/\mu}}
#'
#' \strong{Score, observed and expected Hessian:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} - \dfrac{2y}{\mu^3},
#'       \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu^2}}
#'
#' Every order follows the same pattern, the log-density being a logarithm plus
#' a reciprocal:
#' \deqn{\ell^{(k)} = \dfrac{(-1)^k (k-1)!}{\mu^k}
#'       + \dfrac{(-1)^{k+1} k!\, y}{\mu^{k+1}}, \qquad
#'       \mathbb{E}[\ell^{(k)}] = \dfrac{(-1)^k (k-1)! (1-k)}{\mu^k}}
#' so the expected orders are closed form as well, and vanish at \eqn{k = 1}.
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\mu^2}, skewness 2, excess
#' kurtosis 6.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#' }
#'
#' The family is the Weibull with unit shape, so
#' \code{fixed(weibull_distrib(), sigma = 1)} describes the same law and is
#' used in the tests as an independent implementation. It is \strong{not} a
#' Gamma with a fixed parameter: this package writes the Gamma in
#' \eqn{(\mu, \sigma^2)}, whose shape is \eqn{\mu^2/\sigma^2}, so unit shape is
#' the relation \eqn{\sigma^2 = \mu^2} between two parameters rather than a
#' value one of them can be held at.
#'
#' @return An S7 object of class \code{ExponentialDistrib}.
#'
#' @seealso \code{\link{gamma_distrib}}, \code{\link{weibull_distrib}},
#'   \code{\link{geometric_distrib}}
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dexp pexp qexp rexp
#' @examples
#' d <- exponential_distrib()
#' d@params
#'
#' theta <- list(mu = 2)
#' distrib_pdf(d, c(0.5, 1, 3), theta)
#' distrib_gradient(d, c(0.5, 1, 3), theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' @export
exponential_distrib <- function(link_mu = log_link()) {
  ExponentialDistrib(
    distrib_name = "exponential", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = link_mu)
  )
}
