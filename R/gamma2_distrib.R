#' @include distrib.R generics.R
NULL

#' @title S7 Class for Gamma Distribution
#' @name Gamma2Distrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Gamma distribution 
#' under the mean-variance parameterization.
#' @inheritParams distrib
#' @return An object of class \code{Gamma2Distrib}.
#' @seealso \code{\link{gamma2_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.Gamma2Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.Gamma2Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Gamma2Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.Gamma2Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.Gamma2Distrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.Gamma2Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.Gamma2Distrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.Gamma2Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.Gamma2Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.Gamma2Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.Gamma2Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
Gamma2Distrib <- S7::new_class("Gamma2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gamma Probability Density Function
#' @name distrib_pdf.Gamma2Distrib
#' @description
#' Computes the probability density function for the Gamma distribution:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{\Gamma\left(\dfrac{\mu^2}{\sigma^2}\right)} \left(\dfrac{\mu}{\sigma^2}\right)^{\dfrac{\mu^2}{\sigma^2}} y^{\dfrac{\mu^2}{\sigma^2}-1} \exp\left(-\dfrac{\mu}{\sigma^2} y\right)}
#' 
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_pdf, Gamma2Distrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dgamma(
    x = y,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    log = log
  )
}

#' @title Gamma Cumulative Distribution Function
#' @name distrib_cdf.Gamma2Distrib
#' @description
#' Computes the cumulative distribution function for the Gamma distribution, using
#' the shape/rate reparameterization \eqn{\alpha = \mu^2/\sigma^2},
#' \eqn{\lambda = \mu/\sigma^2}:
#' \deqn{F(q; \mu, \sigma^2) = \dfrac{\gamma(\alpha, \lambda q)}{\Gamma(\alpha)}}
#' where \eqn{\gamma(\cdot, \cdot)} is the lower incomplete gamma function.
#'
#' @param distrib A \code{Gamma2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_cdf, Gamma2Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pgamma(
    q = q,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gamma Quantile Function
#' @name distrib_quantile.Gamma2Distrib
#' @description
#' Computes the quantile function for the Gamma distribution as the inverse of the
#' CDF, \eqn{Q(p; \mu, \sigma^2) = F^{-1}(p; \mu, \sigma^2)}. There is no elementary
#' closed form; it is obtained numerically (via \code{\link[stats]{qgamma}}) on the
#' shape/rate reparameterization \eqn{\alpha = \mu^2/\sigma^2}, \eqn{\lambda = \mu/\sigma^2}.
#'
#' @param distrib A \code{Gamma2Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_quantile, Gamma2Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qgamma(
    p = p,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gamma Random Number Generator
#' @name distrib_rng.Gamma2Distrib
#' @description
#' Generates random numbers from the Gamma distribution.
#' 
#' @param distrib A \code{Gamma2Distrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_rng, Gamma2Distrib) <- function(distrib, n, theta) {
  stats::rgamma(
    n = n,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]]
  )
}

#' @title Gamma Analytical Gradient
#' @name distrib_gradient.Gamma2Distrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Gamma log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{-2\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + 2\mu\log\left(\dfrac{\mu}{\sigma^2}\right) + \mu + 2\mu\log(y) - y}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} = -\dfrac{\mu\left[-\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + \mu + \mu\left(\log\left(\dfrac{\mu}{\sigma^2}\right) + \log(y)\right) - y\right]}{(\sigma^2)^2}}
#' 
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_gradient, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  gamma_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Analytical Observed Hessian
#' @name distrib_hessian.Gamma2Distrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Gamma log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_hessian, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  gamma_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Analytical Expected Hessian
#' @name distrib_expected_hessian.Gamma2Distrib
#' @description
#' Computes the analytical expected Hessian of the Gamma log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = \dfrac{3\sigma^2 - 4\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right)}{(\sigma^2)^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{\mu^2\left(\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - \sigma^2\right)}{(\sigma^2)^4}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2}\right] = \dfrac{2\mu\left(\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - \sigma^2\right)}{(\sigma^2)^3}}
#' 
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_expected_hessian, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  gamma_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Analytical Third-Order Derivatives
#' @name distrib_deriv3.Gamma2Distrib
#' @description Closed-form third-order derivatives of the Gamma log-density (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_deriv3, Gamma2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) gamma_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gamma_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.Gamma2Distrib
#' @description Closed-form fourth-order derivatives of the Gamma log-density (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_deriv4, Gamma2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) gamma_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else gamma_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Gamma Response Derivatives
#' @name distrib_grad_y.Gamma2Distrib
#' @description
#' Closed-form derivatives of the Gamma log-density with respect to the response,
#' using the shape/rate reparameterization \eqn{\alpha = \mu^2/\sigma^2},
#' \eqn{\lambda = \mu/\sigma^2}:
#' \eqn{\partial \ell / \partial y = (\alpha-1)/y - \lambda} and
#' \eqn{\partial^2 \ell / \partial y^2 = -(\alpha-1)/y^2}.
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A numeric vector.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_grad_y, Gamma2Distrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; s2 <- theta[[2]]
  (mu^2 / s2 - 1) / y - mu / s2
}

#' @title Gamma Response Second Derivative
#' @name distrib_hess_y.Gamma2Distrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -(\alpha-1)/y^2}, \eqn{\alpha = \mu^2/\sigma^2}.
#' @param distrib A \code{Gamma2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A numeric vector.
#' @seealso \code{\link{gamma2_distrib}}
S7::method(distrib_hess_y, Gamma2Distrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; s2 <- theta[[2]]
  -(mu^2 / s2 - 1) / y^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Gamma Distribution Object (Mean-Variance Parameterization)
#'
#' @description
#' Creates a distribution object for the Gamma distribution parameterized by mean (\eqn{\mu}) and variance (\eqn{\sigma^2}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#' @param link_sigma2 A link function object for the variance parameter \eqn{\sigma^2}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The Gamma distribution is given a mean/variance parameterization: \eqn{\mu} is
#' the mean and \eqn{\sigma^2} the variance. The standard shape \eqn{\alpha} and
#' rate \eqn{\lambda} are recovered as
#' \deqn{\alpha = \dfrac{\mu^2}{\sigma^2}, \qquad \lambda = \dfrac{\mu}{\sigma^2}}
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{\lambda^\alpha}{\Gamma(\alpha)}\, y^{\alpha-1} e^{-\lambda y}, \quad y > 0}
#'
#' \strong{Cumulative distribution function} (\eqn{\gamma} the lower incomplete gamma function):
#' \deqn{F(q; \mu, \sigma^2) = \dfrac{\gamma(\alpha, \lambda q)}{\Gamma(\alpha)}}
#'
#' \strong{Quantile function:} no closed form; the numerical inverse of the CDF.
#'
#' \strong{Score} (\eqn{\psi} the digamma function):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{-2\mu\psi(\alpha) + 2\mu\log\lambda + \mu + 2\mu\log y - y}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} = -\dfrac{\mu\left[-\mu\psi(\alpha) + \mu + \mu(\log\lambda + \log y) - y\right]}{(\sigma^2)^2}}
#'
#' \strong{Expected Hessian} (\eqn{\psi_1} the trigamma function):
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = \dfrac{3\sigma^2 - 4\mu^2\psi_1(\alpha)}{(\sigma^2)^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma^2}\right] = \dfrac{2\mu(\mu^2\psi_1(\alpha) - \sigma^2)}{(\sigma^2)^3}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{\mu^2(\mu^2\psi_1(\alpha) - \sigma^2)}{(\sigma^2)^4}}
#' The observed Hessian is available via \code{\link{distrib_hessian.Gamma2Distrib}}.
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\sigma^2},
#' skewness \eqn{2\sqrt{\sigma^2}/\mu}, excess kurtosis \eqn{6\sigma^2/\mu^2}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\sigma^2 \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) and response derivatives (\code{\link{distrib_grad_y}},
#' \code{\link{distrib_hess_y}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.Gamma2Distrib}} for the probability density function.
#'   \item \code{\link{distrib_cdf.Gamma2Distrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.Gamma2Distrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.Gamma2Distrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.Gamma2Distrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.Gamma2Distrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.Gamma2Distrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{Gamma2Distrib} (inheriting from \code{continuous_distrib}) representing the Gamma distribution.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dgamma pgamma qgamma rgamma
#' @examples
#' d <- gamma2_distrib()
#' d@params
#'
#' theta <- list(mu = 2, sigma2 = 1)
#' distrib_pdf(d, c(0.5, 1, 2), theta)
#' distrib_gradient(d, c(0.5, 1, 2), theta)
#'
#' @export
gamma2_distrib <- function(link_mu = log_link(), link_sigma2 = log_link()) {
  
  Gamma2Distrib(
    distrib_name = "gamma2",
    dimension = "univariate",
    bounds = c(0, Inf),
    
    params = c("mu", "sigma2"),
    params_interpretation = c(mu = "mean", sigma2 = "variance"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(0, Inf),
      sigma2 = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      sigma2 = link_sigma2
    )
  )
  
}
