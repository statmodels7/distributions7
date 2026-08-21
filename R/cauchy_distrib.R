#' @include distrib.R generics.R
NULL

#' @title S7 Class for Cauchy Distribution
#' @name CauchyDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Cauchy distribution.
#' @inheritParams distrib
#' @return An object of class \code{CauchyDistrib}.
#' @seealso \code{\link{cauchy_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.CauchyDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.CauchyDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.CauchyDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.CauchyDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.CauchyDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.CauchyDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.CauchyDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.CauchyDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.CauchyDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.CauchyDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.CauchyDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
CauchyDistrib <- S7::new_class("CauchyDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Cauchy Probability Density Function
#' @name distrib_pdf.CauchyDistrib
#' @description
#' Computes the probability density function for the Cauchy distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right]}}
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_pdf, CauchyDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dcauchy(
    x = y,
    location = theta[[1]],
    scale = theta[[2]],
    log = log
  )
}

#' @title Cauchy Cumulative Distribution Function
#' @name distrib_cdf.CauchyDistrib
#' @description
#' Computes the cumulative distribution function for the Cauchy distribution:
#' \deqn{F(q; \mu, \sigma) = \dfrac{1}{2} + \dfrac{1}{\pi}\arctan\left(\dfrac{q-\mu}{\sigma}\right)}
#'
#' @param distrib A \code{CauchyDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_cdf, CauchyDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pcauchy(
    q = q,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Cauchy Quantile Function
#' @name distrib_quantile.CauchyDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Cauchy distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \tan\left(\pi\left(p - \dfrac{1}{2}\right)\right)}
#'
#' @param distrib A \code{CauchyDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_quantile, CauchyDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qcauchy(
    p = p,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Cauchy Random Number Generator
#' @name distrib_rng.CauchyDistrib
#' @description
#' Generates random numbers from the Cauchy distribution.
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_rng, CauchyDistrib) <- function(distrib, n, theta) {
  stats::rcauchy(
    n = n,
    location = theta[[1]],
    scale = theta[[2]]
  )
}

#' @title Cauchy Analytical Gradient
#' @name distrib_gradient.CauchyDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Cauchy log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{2(y-\mu)}{\sigma^2 + (y-\mu)^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y-\mu)^2 - \sigma^2}{\sigma(\sigma^2 + (y-\mu)^2)}}
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_gradient, CauchyDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  cauchy_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Analytical Observed Hessian
#' @name distrib_hessian.CauchyDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Cauchy log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(y-\mu)^2 - 2\sigma^2}{(\sigma^2 + (y-\mu)^2)^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 (y-\mu)^2 - (y-\mu)^4}{\sigma^2(\sigma^2 + (y-\mu)^2)^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{4\sigma (y-\mu)}{(\sigma^2 + (y-\mu)^2)^2}}
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_hessian, CauchyDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  cauchy_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Analytical Expected Hessian
#' @name distrib_expected_hessian.CauchyDistrib
#' @description
#' Computes the analytical expected Hessian of the Cauchy log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{2\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{2\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_expected_hessian, CauchyDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  cauchy_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Analytical Third-Order Derivatives
#' @name distrib_deriv3.CauchyDistrib
#' @description Closed-form third-order derivatives of the Cauchy log-density (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_deriv3, CauchyDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) cauchy_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else cauchy_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.CauchyDistrib
#' @description Closed-form fourth-order derivatives of the Cauchy log-density (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_deriv4, CauchyDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) cauchy_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else cauchy_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Cauchy Response Derivatives
#' @name distrib_grad_y.CauchyDistrib
#' @description
#' Closed-form derivatives of the Cauchy log-density with respect to the response.
#' Let \eqn{r = y - \mu} and \eqn{d = \sigma^2 + r^2}:
#' \eqn{\partial \ell / \partial y = -2r/d} and
#' \eqn{\partial^2 \ell / \partial y^2 = 2(r^2 - \sigma^2)/d^2}.
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_grad_y, CauchyDistrib) <- function(distrib, y, theta) {
  r <- y - theta[[1]]
  -2 * r / (theta[[2]]^2 + r^2)
}

#' @title Cauchy Response Second Derivative
#' @name distrib_hess_y.CauchyDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = 2(r^2 - \sigma^2)/(\sigma^2 + r^2)^2}, \eqn{r = y - \mu}.
#' @param distrib A \code{CauchyDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_hess_y, CauchyDistrib) <- function(distrib, y, theta) {
  r <- y - theta[[1]]
  s2 <- theta[[2]]^2
  2 * (r^2 - s2) / (s2 + r^2)^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Cauchy Distribution Object
#'
#' @description
#' Creates a distribution object for the Cauchy distribution, parameterized by location (\eqn{\mu})
#' and scale (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The Cauchy distribution is a heavy-tailed location-scale distribution with
#' location \eqn{\mu} and scale \eqn{\sigma}. Its moments (mean, variance, and
#' higher) are \strong{undefined}.
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right]}}
#'
#' \strong{Cumulative distribution function:}
#' \deqn{F(q; \mu, \sigma) = \dfrac{1}{2} + \dfrac{1}{\pi}\arctan\left(\dfrac{q-\mu}{\sigma}\right)}
#'
#' \strong{Quantile function:}
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \tan\left(\pi\left(p - \tfrac{1}{2}\right)\right)}
#'
#' \strong{Score} (with \eqn{d = \sigma^2 + (y-\mu)^2}):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{2(y-\mu)}{d}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y-\mu)^2 - \sigma^2}{\sigma d}}
#'
#' \strong{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(y-\mu)^2 - 2\sigma^2}{d^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 (y-\mu)^2 - (y-\mu)^4}{\sigma^2 d^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma} = -\dfrac{4\sigma(y-\mu)}{d^2}}
#'
#' \strong{Expected Hessian:}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] =
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{2\sigma^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0}
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) and response derivatives (\code{\link{distrib_grad_y}},
#' \code{\link{distrib_hess_y}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.CauchyDistrib}} for the density function.
#'   \item \code{\link{distrib_cdf.CauchyDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.CauchyDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.CauchyDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.CauchyDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.CauchyDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.CauchyDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{CauchyDistrib} (inheriting from \code{continuous_distrib}) representing the Cauchy distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dcauchy pcauchy qcauchy rcauchy
#' @examples
#' d <- cauchy_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' @export
cauchy_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  
  CauchyDistrib(
    distrib_name = "cauchy",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "location", sigma = "scale"),
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
