#' @include distrib.R generics.R

#' @title S7 Class for Gamma Distribution
#' @name GammaDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Gamma distribution 
#' under the mean-variance parameterization.
#' @inheritParams distrib
#' @seealso \code{\link{gamma_distrib}}
GammaDistrib <- S7::new_class("GammaDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gamma Probability Density Function
#' @name distrib_pdf.GammaDistrib
#' @description
#' Computes the probability density function for the Gamma distribution:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{\Gamma\left(\dfrac{\mu^2}{\sigma^2}\right)} \left(\dfrac{\mu}{\sigma^2}\right)^{\dfrac{\mu^2}{\sigma^2}} y^{\dfrac{\mu^2}{\sigma^2}-1} \exp\left(-\dfrac{\mu}{\sigma^2} y\right)}
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_pdf, GammaDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dgamma(
    x = y,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    log = log
  )
}

#' @title Gamma Cumulative Distribution Function
#' @name distrib_cdf.GammaDistrib
#' @description
#' Computes the cumulative distribution function for the Gamma distribution.
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_cdf, GammaDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pgamma(
    q = q,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gamma Quantile Function
#' @name distrib_quantile.GammaDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Gamma distribution.
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_quantile, GammaDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qgamma(
    p = p,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Gamma Random Number Generator
#' @name distrib_rng.GammaDistrib
#' @description
#' Generates random numbers from the Gamma distribution.
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_rng, GammaDistrib) <- function(distrib, n, theta) {
  stats::rgamma(
    n = n,
    shape = theta[[1]]^2 / theta[[2]],
    rate = theta[[1]] / theta[[2]]
  )
}

#' @title Gamma Analytical Gradient
#' @name distrib_gradient.GammaDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Gamma log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{-2\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + 2\mu\log\left(\dfrac{\mu}{\sigma^2}\right) + \mu + 2\mu\log(y) - y}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} = -\dfrac{\mu\left[-\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + \mu + \mu\left(\log\left(\dfrac{\mu}{\sigma^2}\right) + \log(y)\right) - y\right]}{(\sigma^2)^2}}
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_gradient, GammaDistrib) <- function(distrib, y, theta) {
  gamma_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gamma Analytical Observed Hessian
#' @name distrib_hessian.GammaDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Gamma log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_hessian, GammaDistrib) <- function(distrib, y, theta) {
  gamma_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gamma Analytical Expected Hessian
#' @name distrib_expected_hessian.GammaDistrib
#' @description
#' Computes the analytical expected Hessian of the Gamma log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = \dfrac{3\sigma^2 - 4\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right)}{(\sigma^2)^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{\mu^2\left(\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - \sigma^2\right)}{(\sigma^2)^4}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2}\right] = \dfrac{2\mu\left(\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - \sigma^2\right)}{(\sigma^2)^3}}
#' 
#' @param distrib A \code{GammaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{gamma_distrib}}
S7::method(distrib_expected_hessian, GammaDistrib) <- function(distrib, y, theta) {
  gamma_expected_hessian_cpp(y, theta[[1]], theta[[2]])
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
#' The Gamma distribution is reparameterized from the standard shape \eqn{(\alpha)} and rate \eqn{(\lambda)} parameters using:
#' \deqn{\alpha = \dfrac{\mu^2}{\sigma^2}, \quad \lambda = \dfrac{\mu}{\sigma^2}}
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\sigma^2 \in (0, +\infty)}
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.GammaDistrib}} for the probability density function.
#'   \item \code{\link{distrib_cdf.GammaDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.GammaDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.GammaDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.GammaDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.GammaDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.GammaDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{GammaDistrib} (inheriting from \code{continuous_distrib}) representing the Gamma distribution.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dgamma pgamma qgamma rgamma
#' @export
gamma_distrib <- function(link_mu = log_link(), link_sigma2 = log_link()) {
  
  GammaDistrib(
    distrib_name = "gamma",
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