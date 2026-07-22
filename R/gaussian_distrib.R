#' @include distrib.R generics.R

#' @title S7 Class for Gaussian Distribution
#' @name GaussianDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Gaussian (Normal) distribution.
#' @inheritParams distrib
#' @seealso \code{\link{gaussian_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.GaussianDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.GaussianDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.GaussianDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.GaussianDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.GaussianDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.GaussianDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.GaussianDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.GaussianDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.GaussianDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.GaussianDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.GaussianDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
GaussianDistrib <- S7::new_class("GaussianDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Probability Density Function
#' @name distrib_pdf.GaussianDistrib
#' @description
#' Computes the probability density function for the Gaussian distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\}}
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_pdf, GaussianDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dnorm(
    x = y,
    mean = theta[[1]],
    sd = theta[[2]],
    log = log
  )
}

#' @title Gaussian Cumulative Distribution Function
#' @name distrib_cdf.GaussianDistrib
#' @description
#' Computes the cumulative distribution function for the Gaussian distribution:
#' \deqn{F(q; \mu, \sigma) = \Phi\left(\dfrac{q-\mu}{\sigma}\right)}
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_cdf, GaussianDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pnorm(
    q = q,
    mean = theta[[1]],
    sd = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gaussian Quantile Function
#' @name distrib_quantile.GaussianDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Gaussian distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \Phi^{-1}(p)}
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_quantile, GaussianDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qnorm(
    p = p,
    mean = theta[[1]],
    sd = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gaussian Random Number Generator
#' @name distrib_rng.GaussianDistrib
#' @description
#' Generates random numbers from the Gaussian distribution.
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_rng, GaussianDistrib) <- function(distrib, n, theta) {
  stats::rnorm(
    n = n,
    mean = theta[[1]],
    sd = theta[[2]]
  )
}

#' @title Gaussian Analytical Gradient
#' @name distrib_gradient.GaussianDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}}
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_gradient, GaussianDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  gaussian_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Analytical Observed Hessian
#' @name distrib_hessian.GaussianDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y - \mu)^2}{\sigma^4}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{2(y - \mu)}{\sigma^3}}
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_hessian, GaussianDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  gaussian_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Analytical Expected Hessian
#' @name distrib_expected_hessian.GaussianDistrib
#' @description
#' Computes the analytical expected Hessian of the Gaussian log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#' 
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_expected_hessian, GaussianDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  gaussian_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Analytical Third-Order Derivatives
#' @name distrib_deriv3.GaussianDistrib
#' @description
#' Computes the closed-form third-order partial derivatives of the Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\sigma} (observed, or expected when
#' \code{expected = TRUE}).
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_deriv3, GaussianDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    gaussian_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  } else {
    gaussian_deriv3_cpp(y, theta[[1]], theta[[2]])
  }
}

#' @title Gaussian Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.GaussianDistrib
#' @description
#' Computes the closed-form fourth-order partial derivatives of the Gaussian
#' log-density with respect to \eqn{\mu} and \eqn{\sigma} (observed, or expected when
#' \code{expected = TRUE}).
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_deriv4, GaussianDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    gaussian_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  } else {
    gaussian_deriv4_cpp(y, theta[[1]], theta[[2]])
  }
}

#' @title Gaussian Response Derivatives
#' @name distrib_grad_y.GaussianDistrib
#' @description
#' Closed-form derivatives of the Gaussian log-density with respect to the response:
#' \eqn{\partial \ell / \partial y = -(y-\mu)/\sigma^2} and
#' \eqn{\partial^2 \ell / \partial y^2 = -1/\sigma^2}.
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_grad_y, GaussianDistrib) <- function(distrib, y, theta) {
  -(y - theta[[1]]) / theta[[2]]^2
}

#' @title Gaussian Response Second Derivative
#' @name distrib_hess_y.GaussianDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -1/\sigma^2}.
#' @param distrib A \code{GaussianDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian_distrib}}
S7::method(distrib_hess_y, GaussianDistrib) <- function(distrib, y, theta) {
  rep(-1 / theta[[2]]^2, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Gaussian Distribution Object (Standard Deviation Parameterization)
#'
#' @description
#' Creates a distribution object for the Gaussian distribution parameterized by mean (\eqn{\mu}) and standard deviation (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The Gaussian (Normal) distribution is parameterized by its mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}.
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\,\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\}}
#'
#' \strong{Cumulative distribution function} (\eqn{\Phi} the standard normal CDF):
#' \deqn{F(q; \mu, \sigma) = \Phi\left(\dfrac{q-\mu}{\sigma}\right)}
#'
#' \strong{Quantile function:}
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma\,\Phi^{-1}(p)}
#'
#' \strong{Score} (gradient of the log-density):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}}
#'
#' \strong{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y-\mu)^2}{\sigma^4}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma} = -\dfrac{2(y-\mu)}{\sigma^3}}
#'
#' \strong{Expected Hessian} (negative Fisher information):
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0}
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\sigma^2}, skewness 0, excess kurtosis 0.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) and derivatives with respect to the response
#' (\code{\link{distrib_grad_y}}, \code{\link{distrib_hess_y}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.GaussianDistrib}} for the density function.
#'   \item \code{\link{distrib_cdf.GaussianDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.GaussianDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.GaussianDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.GaussianDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.GaussianDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.GaussianDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{GaussianDistrib} (inheriting from \code{continuous_distrib}) representing the Gaussian distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @export
gaussian_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  
  GaussianDistrib(
    distrib_name = "gaussian",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean", sigma = "standard deviation"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      sigma = link_sigma
    )
  )
  
}
