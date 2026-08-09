#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Gamma Distribution in Mean and Dispersion
#' @name Gamma1Distrib
#'
#' @description A subclass of \code{continuous_distrib} for the gamma written
#'   in its mean and a dispersion, with \eqn{\operatorname{Var}(Y) = \phi\mu^2}.
#' @inheritParams distrib
#' @return An object of class \code{Gamma1Distrib}.
#' @seealso \code{\link{gamma1_distrib}}, \code{\link{gamma2_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.Gamma1Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.Gamma1Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Gamma1Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.Gamma1Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.Gamma1Distrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.Gamma1Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.Gamma1Distrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.Gamma1Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.Gamma1Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.Gamma1Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.Gamma1Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
Gamma1Distrib <- S7::new_class("Gamma1Distrib", parent = continuous_distrib)

#' The Shape and Rate a Mean and Dispersion Imply
#'
#' @description
#' \eqn{a = 1/\phi} and \eqn{b = 1/(\phi\mu)}, which is what makes the variance
#' \eqn{\phi\mu^2}.
#'
#' @param theta A list with \code{mu} and \code{phi}.
#'
#' @return A list with \code{shape} and \code{rate}.
#'
#' @seealso \code{\link{gamma1_distrib}}
#'
#' @keywords internal
gamma1_shape_rate <- function(theta) {
  list(shape = 1 / theta[[2]], rate = 1 / (theta[[2]] * theta[[1]]))
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gamma Density in Mean and Dispersion
#' @name distrib_pdf.Gamma1Distrib
#' @description
#' The gamma density at shape \eqn{a = 1/\phi} and rate
#' \eqn{b = 1/(\phi\mu)}.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_pdf, Gamma1Distrib) <- function(distrib, y, theta, log = FALSE) {
  sr <- gamma1_shape_rate(theta)
  stats::dgamma(y, shape = sr$shape, rate = sr$rate, log = log)
}

#' @title Gamma Distribution Function in Mean and Dispersion
#' @name distrib_cdf.Gamma1Distrib
#' @description The gamma distribution function at the implied shape and rate.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns log-probabilities.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_cdf, Gamma1Distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log.p = FALSE, ...) {
  sr <- gamma1_shape_rate(theta)
  stats::pgamma(q, shape = sr$shape, rate = sr$rate,
                lower.tail = lower.tail, log.p = log.p)
}

#' @title Gamma Quantile Function in Mean and Dispersion
#' @name distrib_quantile.Gamma1Distrib
#' @description The gamma quantile function at the implied shape and rate.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as a log-probability.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_quantile, Gamma1Distrib) <- function(distrib, p, theta,
                                                         lower.tail = TRUE,
                                                         log.p = FALSE, ...) {
  sr <- gamma1_shape_rate(theta)
  stats::qgamma(p, shape = sr$shape, rate = sr$rate,
                lower.tail = lower.tail, log.p = log.p)
}

#' @title Gamma Random Generation in Mean and Dispersion
#' @name distrib_rng.Gamma1Distrib
#' @description Delegates to \code{\link[stats]{rgamma}}.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param n The number of draws.
#' @param theta A list with \code{mu} and \code{phi}.
#' @return A numeric vector.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_rng, Gamma1Distrib) <- function(distrib, n, theta) {
  sr <- gamma1_shape_rate(theta)
  stats::rgamma(n, shape = sr$shape, rate = sr$rate)
}

#' @title Gamma Analytical Gradient in Mean and Dispersion
#' @name distrib_gradient.Gamma1Distrib
#' @description
#' With \eqn{s = 1/\phi} and \eqn{z = y/\mu},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \dfrac{y-\mu}{\phi\mu^2}, \qquad
#'       \dfrac{\partial\ell}{\partial\phi}
#'         = -s^2\left\{\log s + 1 - \psi(s) + \log z - z\right\}}
#' the first being the score of a gamma generalized linear model.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_gradient, Gamma1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"), ...) {
  gamma1_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gamma Analytical Observed Hessian in Mean and Dispersion
#' @name distrib_hessian.Gamma1Distrib
#' @description
#' Closed form. The derivatives in \eqn{\phi} are those in \eqn{s = 1/\phi}
#' carried across by the one-variable chain rule, which is what keeps the
#' polygamma functions to one evaluation each.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_hessian, Gamma1Distrib) <- function(distrib, y, theta,
                                                        scale = c("parameter", "link"), ...) {
  gamma1_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gamma Analytical Expected Hessian in Mean and Dispersion
#' @name distrib_expected_hessian.Gamma1Distrib
#' @description
#' Closed form:
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\dfrac{1}{\phi\mu^2}, \qquad
#'       \mathbb{E}[\ell^{(\mu\phi)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\phi\phi)}] = s^4\left\{\dfrac{1}{s} - \psi'(s)\right\}}
#' The mean and the dispersion are orthogonal, which is what makes this the
#' natural parametrization for a generalized linear model. The expectation uses
#' \eqn{\mathbb{E}[\log(Y/\mu)] = \psi(s) - \log s}.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_expected_hessian, Gamma1Distrib) <- function(distrib, y, theta,
                                                                 scale = c("parameter", "link"),
                                                                 approx = c("bartlett", "integrate", "mc", "opg"),
                                                                 nsim = 10000, ...) {
  gamma1_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gamma Third-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv3.Gamma1Distrib
#' @description Closed form, observed or expected.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param expected Logical; if \code{TRUE}, returns the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_deriv3, Gamma1Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                       scale = c("parameter", "link"),
                                                       approx = c("integrate", "bartlett", "mc", "opg"),
                                                       nsim = 10000, ...) {
  if (expected) {
    gamma1_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  } else {
    gamma1_deriv3_cpp(y, theta[[1]], theta[[2]])
  }
}

#' @title Gamma Fourth-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv4.Gamma1Distrib
#' @description Closed form, observed or expected.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param expected Logical; if \code{TRUE}, returns the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_deriv4, Gamma1Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                       scale = c("parameter", "link"),
                                                       approx = c("integrate", "bartlett", "mc", "opg"),
                                                       nsim = 10000, ...) {
  if (expected) {
    gamma1_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  } else {
    gamma1_deriv4_cpp(y, theta[[1]], theta[[2]])
  }
}

#' @title Gamma Response Derivatives in Mean and Dispersion
#' @name distrib_grad_y.Gamma1Distrib
#' @description \eqn{\partial\ell/\partial y = (a-1)/y - b} at the implied
#'   shape and rate.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_grad_y, Gamma1Distrib) <- function(distrib, y, theta, ...) {
  sr <- gamma1_shape_rate(theta)
  (sr$shape - 1) / y - sr$rate
}

#' @title Gamma Second Response Derivative in Mean and Dispersion
#' @name distrib_hess_y.Gamma1Distrib
#' @description \eqn{\partial^2\ell/\partial y^2 = -(a-1)/y^2}.
#' @param distrib A \code{Gamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{phi}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gamma1_distrib}}
S7::method(distrib_hess_y, Gamma1Distrib) <- function(distrib, y, theta, ...) {
  sr <- gamma1_shape_rate(theta)
  -(sr$shape - 1) / y^2
}


#' Gamma Distribution in Mean and Dispersion
#'
#' @description
#' Creates a gamma distribution object parametrized by its mean and a
#' dispersion, with \eqn{\operatorname{Var}(Y) = \phi\mu^2}.
#'
#' @details
#' This is the parametrization a generalized linear model uses: the variance
#' function is \eqn{V(\mu) = \mu^2} and \eqn{\phi} is the dispersion that
#' multiplies it, so the mean and the dispersion are orthogonal and the score
#' in \eqn{\mu} is \eqn{(y-\mu)/(\phi\mu^2)}. The shape is \eqn{1/\phi} and the
#' rate \eqn{1/(\phi\mu)}.
#'
#' It is the same law as \code{\link{gamma2_distrib}}, which carries the mean
#' and the \emph{variance}: the two are related by
#' \eqn{\sigma^2 = \phi\mu^2}. They are separate families because the second
#' parameter is a different quantity in each, with its own interpretation,
#' standard error and interval.
#'
#' Derivatives are closed form to fourth order, observed and expected. Those in
#' \eqn{\phi} go through \eqn{s = 1/\phi} and the one-variable chain rule, so
#' each polygamma function is evaluated once.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{y^{1/\phi - 1}e^{-y/(\phi\mu)}}{(\phi\mu)^{1/\phi}\,\Gamma(1/\phi)}}
#' on \eqn{y \in (0, \infty)}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \phi\mu^{2}}
#'
#' @param link_mu Link function for \eqn{\mu}. Defaults to the log.
#' @param link_phi Link function for \eqn{\phi}. Defaults to the log.
#'
#' @return An S7 object of class \code{\link{Gamma1Distrib}}.
#'
#' @seealso \code{\link{gamma2_distrib}}
#'
#' @examples
#' d <- gamma1_distrib()
#' theta <- list(mu = 3, phi = 0.5)
#' distrib_pdf(d, c(1, 3, 5), theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' # the same law as gamma2 with sigma2 = phi * mu^2
#' distrib_pdf(gamma2_distrib(), 3, list(mu = 3, sigma2 = 0.5 * 9))
#'
#' @export
gamma1_distrib <- function(link_mu = log_link(), link_phi = log_link()) {
  Gamma1Distrib(
    distrib_name = "gamma1",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "phi"),
    params_interpretation = c(mu = "mean", phi = "dispersion"),
    n_params = 2,
    params_bounds = list(mu = c(0, Inf), phi = c(0, Inf)),
    link_params = list(mu = link_mu, phi = link_phi)
  )
}
