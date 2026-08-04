#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Weibull Distribution
#' @name WeibullDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the Weibull
#' distribution, parametrised by a scale and a shape.
#' @inheritParams distrib
#' @return An object of class \code{WeibullDistrib}.
#' @seealso \code{\link{weibull_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.WeibullDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.WeibullDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.WeibullDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.WeibullDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.WeibullDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.WeibullDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.WeibullDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.WeibullDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.WeibullDistrib]{distrib_rng()}},
#'   \code{\link[=kurtosis]{kurtosis()}},
#'   \code{\link[=mean.distrib]{mean()}},
#'   \code{\link[=skewness]{skewness()}},
#'   \code{\link[=variance]{variance()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
WeibullDistrib <- S7::new_class("WeibullDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' The Pieces a Weibull Evaluates From
#'
#' @description
#' Returns the standardised variable \eqn{z = y/\mu}, its logarithm and the
#' quantity \eqn{u = z^{\sigma}}, which every derivative of the Weibull
#' log-density is a polynomial in.
#'
#' @details
#' \eqn{u} is the substitution that makes the family tractable: under the model
#' \eqn{u \sim \mathrm{Exp}(1)} whatever the parameters are, so an expectation
#' of any polynomial in \eqn{u} and \eqn{\log u} is a derivative of the gamma
#' function at 2, which is what \code{\link{distrib_expected_hessian}} uses.
#'
#' @param y A numeric vector of observations.
#' @param mu The scale parameter.
#' @param sigma The shape parameter.
#'
#' @return A list with \code{z}, \code{lz} and \code{u}.
#'
#' @keywords internal
weibull_pieces <- function(y, mu, sigma) {
  z <- y / mu
  lz <- log(z)
  list(z = z, lz = lz, u = exp(sigma * lz))
}

#' @title Weibull Probability Density Function
#' @name distrib_pdf.WeibullDistrib
#' @description
#' Computes the probability density function for the Weibull distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{\sigma}{\mu}
#'   \left(\dfrac{y}{\mu}\right)^{\sigma - 1}
#'   \exp\left\{-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\}}
#' @param distrib A \code{WeibullDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_pdf, WeibullDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dweibull(x = y, shape = theta[[2]], scale = theta[[1]], log = log)
}

#' @title Weibull Cumulative Distribution Function
#' @name distrib_cdf.WeibullDistrib
#' @description
#' Computes the cumulative distribution function for the Weibull distribution:
#' \deqn{F(q; \mu, \sigma) = 1 - \exp\left\{-(q/\mu)^{\sigma}\right\}}
#' The survival function is exact on the log scale, which is what a censored
#' observation in the far tail needs.
#' @param distrib A \code{WeibullDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_cdf, WeibullDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pweibull(
    q = q, shape = theta[[2]], scale = theta[[1]],
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @title Weibull Quantile Function
#' @name distrib_quantile.WeibullDistrib
#' @description
#' Computes the quantile function for the Weibull distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu \left\{-\log(1 - p)\right\}^{1/\sigma}}
#' @param distrib A \code{WeibullDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_quantile, WeibullDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qweibull(
    p = p, shape = theta[[2]], scale = theta[[1]],
    lower.tail = lower.tail, log.p = log.p
  )
}

#' @title Weibull Random Number Generator
#' @name distrib_rng.WeibullDistrib
#' @description
#' Generates random numbers by inverse transform, which is exact here because
#' the quantile function is elementary.
#' @param distrib A \code{WeibullDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_rng, WeibullDistrib) <- function(distrib, n, theta) {
  stats::rweibull(n = n, shape = theta[[2]], scale = theta[[1]])
}

#' @title Weibull Analytical Gradient
#' @name distrib_gradient.WeibullDistrib
#' @description
#' Closed-form first derivatives of the Weibull log-density, written in
#' \eqn{u = (y/\mu)^{\sigma}} and \eqn{\log z = \log(y/\mu)}:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1),
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}
#'         + (1 - u)\log z}
#' @param distrib A \code{WeibullDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_gradient, WeibullDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  p <- weibull_pieces(y, mu, sigma)
  list(
    mu = sigma * (p$u - 1) / mu,
    sigma = 1 / sigma + (1 - p$u) * p$lz
  )
}

#' @title Weibull Analytical Observed Hessian
#' @name distrib_hessian.WeibullDistrib
#' @description
#' Closed-form second derivatives of the Weibull log-density:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2}
#'         = \dfrac{\sigma}{\mu^2}\left\{1 - (1 + \sigma) u\right\},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = -\dfrac{1}{\sigma^2} - u (\log z)^2,}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{1}{\mu}\left(u - 1 + \sigma u \log z\right).}
#' @param distrib A \code{WeibullDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_hessian, WeibullDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  p <- weibull_pieces(y, mu, sigma)
  list(
    mu_mu = sigma * (1 - (1 + sigma) * p$u) / mu^2,
    sigma_sigma = -1 / sigma^2 - p$u * p$lz^2,
    mu_sigma = (p$u - 1 + sigma * p$u * p$lz) / mu
  )
}

#' @title Weibull Analytical Expected Hessian
#' @name distrib_expected_hessian.WeibullDistrib
#' @description
#' Closed form. Under the model \eqn{u = (Y/\mu)^{\sigma}} is standard
#' exponential whatever the parameters are, so every expectation needed is a
#' derivative of \eqn{\Gamma} at 2: \eqn{E[u] = 1}, \eqn{E[u\log u] = 1 - \gamma}
#' and \eqn{E[u(\log u)^2] = (1-\gamma)^2 + \pi^2/6 - 1}. Hence
#' \deqn{E\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'         = -\dfrac{\sigma^2}{\mu^2},
#'       \qquad
#'       E\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right]
#'         = -\dfrac{1}{\sigma^2}\left\{(1-\gamma)^2 + \dfrac{\pi^2}{6}\right\},}
#' \deqn{E\left[\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}\right]
#'         = \dfrac{1 - \gamma}{\mu},}
#' with \eqn{\gamma} the Euler-Mascheroni constant. Because the closed form
#' exists, the \code{approx} argument is ignored.
#' @param distrib A \code{WeibullDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param approx Ignored; the expectation is exact.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_expected_hessian, WeibullDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  n <- length(y)
  eg <- -digamma(1)   # the Euler-Mascheroni constant
  list(
    mu_mu = rep(-sigma^2 / mu^2, length.out = n),
    sigma_sigma = rep(-((1 - eg)^2 + pi^2 / 6) / sigma^2, length.out = n),
    mu_sigma = rep((1 - eg) / mu, length.out = n)
  )
}

#' @title Weibull Response Derivative
#' @name distrib_grad_y.WeibullDistrib
#' @description
#' Closed form: \eqn{\partial \ell / \partial y = (\sigma - 1 - \sigma u)/y}.
#' @param distrib A \code{WeibullDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_grad_y, WeibullDistrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  p <- weibull_pieces(y, theta[[1]], sigma)
  (sigma - 1 - sigma * p$u) / y
}

#' @title Weibull Response Second Derivative
#' @name distrib_hess_y.WeibullDistrib
#' @description
#' Closed form:
#' \eqn{\partial^2 \ell / \partial y^2 = -(\sigma - 1)(1 + \sigma u)/y^2}.
#' @param distrib A \code{WeibullDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{weibull_distrib}}
S7::method(distrib_hess_y, WeibullDistrib) <- function(distrib, y, theta) {
  sigma <- theta[[2]]
  p <- weibull_pieces(y, theta[[1]], sigma)
  -(sigma - 1) * (1 + sigma * p$u) / y^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Weibull Distribution Object
#'
#' @description
#' Creates a distribution object for the Weibull distribution, parametrised by
#' a scale \eqn{\mu} and a shape \eqn{\sigma}, both positive.
#'
#' @param link_mu A link function object for the scale \eqn{\mu}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}}.
#' @param link_sigma A link function object for the shape \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}}.
#'
#' @details
#' \strong{Parametrisation.} \eqn{\mu} is the \strong{scale} and not the mean.
#' The mean is \eqn{\mu\,\Gamma(1 + 1/\sigma)}, which involves the shape, so a
#' mean parametrisation would make every derivative a derivative of the gamma
#' function and its inverse. The scale-shape form keeps the whole family
#' elementary, and \code{\link{mean.WeibullDistrib}} reports the mean.
#' This is the parametrisation of \code{WEI} in \pkg{gamlss}.
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma) = \dfrac{\sigma}{\mu}
#'   \left(\dfrac{y}{\mu}\right)^{\sigma - 1}
#'   \exp\left\{-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\}, \qquad y > 0}
#'
#' \strong{Cumulative distribution function:}
#' \deqn{F(q; \mu, \sigma) = 1 - \exp\left\{-(q/\mu)^{\sigma}\right\}}
#'
#' \strong{Quantile function:}
#' \deqn{Q(p; \mu, \sigma) = \mu\left\{-\log(1-p)\right\}^{1/\sigma}}
#'
#' \strong{Score}, with \eqn{u = (y/\mu)^{\sigma}} and \eqn{z = y/\mu}:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\sigma}{\mu}(u - 1),
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}
#'         + (1 - u)\log z}
#'
#' \strong{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2}
#'         = \dfrac{\sigma}{\mu^2}\left\{1 - (1 + \sigma) u\right\}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = -\dfrac{1}{\sigma^2} - u (\log z)^2, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{u - 1 + \sigma u \log z}{\mu}}
#'
#' \strong{Expected Hessian:} see
#' \code{\link{distrib_expected_hessian.WeibullDistrib}}. The substitution
#' \eqn{u \sim \mathrm{Exp}(1)} turns every expectation into a derivative of
#' \eqn{\Gamma} at 2.
#'
#' \strong{Moments.} With \eqn{g_k = \Gamma(1 + k/\sigma)}, the mean is
#' \eqn{\mu g_1} and the variance \eqn{\mu^2 (g_2 - g_1^2)}; the skewness and
#' the excess kurtosis follow from \eqn{g_3} and \eqn{g_4} and do not depend on
#' \eqn{\mu}.
#'
#' \strong{Special cases.} \eqn{\sigma = 1} is the exponential distribution with
#' mean \eqn{\mu}, and \eqn{\sigma = 2} the Rayleigh distribution. The hazard is
#' increasing for \eqn{\sigma > 1} and decreasing for \eqn{\sigma < 1}, which is
#' what the family is used for.
#'
#' \strong{Higher orders.} Third and fourth derivatives are not registered in
#' closed form and come from the numerical fallbacks
#' (\code{\link{numerical_deriv3}}, \code{\link{numerical_deriv4}}).
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' @return An S7 object of class \code{\link{WeibullDistrib}} (inheriting from
#'   \code{continuous_distrib}).
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994).
#' \emph{Continuous Univariate Distributions, Volume 1}, 2nd edition, chapter 21.
#' Wiley.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dweibull pweibull qweibull rweibull
#'
#' @examples
#' d <- weibull_distrib()
#' d@params
#'
#' theta <- list(mu = 2, sigma = 1.5)
#' distrib_pdf(d, c(0.5, 1, 2), theta)
#' distrib_gradient(d, c(0.5, 1, 2), theta)
#'
#' # the scale is not the mean
#' c(scale = theta$mu, mean = mean(d, theta))
#'
#' # shape 1 is the exponential distribution
#' max(abs(distrib_pdf(d, c(0.5, 1, 2), list(mu = 2, sigma = 1)) -
#'         stats::dexp(c(0.5, 1, 2), rate = 1 / 2)))
#'
#' @export
weibull_distrib <- function(link_mu = log_link(), link_sigma = log_link()) {
  WeibullDistrib(
    distrib_name = "weibull",
    dimension = "univariate",
    bounds = c(0, Inf),

    params = c("mu", "sigma"),
    params_interpretation = c(mu = "scale", sigma = "shape"),
    n_params = 2,

    params_bounds = list(
      mu = c(0, Inf),
      sigma = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma
    )
  )
}
