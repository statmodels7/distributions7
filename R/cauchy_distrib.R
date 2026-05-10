#' @include distrib.R generics.R

#' @title S7 Class for Cauchy Distribution
#' @name CauchyDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Cauchy distribution.
#' @inheritParams distrib
#' @seealso \code{\link{cauchy_distrib}}
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
S7::method(distrib_pdf, CauchyDistrib) <- function(distrib, y, theta, log = FALSE) {
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
#' Computes the cumulative distribution function for the Cauchy distribution.
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
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
#' Computes the quantile function (inverse CDF) for the Cauchy distribution.
#' 
#' @param distrib A \code{CauchyDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
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
S7::method(distrib_gradient, CauchyDistrib) <- function(distrib, y, theta) {
  cauchy_gradient_cpp(y, theta[[1]], theta[[2]])
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
S7::method(distrib_hessian, CauchyDistrib) <- function(distrib, y, theta) {
  cauchy_hessian_cpp(y, theta[[1]], theta[[2]])
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
S7::method(distrib_expected_hessian, CauchyDistrib) <- function(distrib, y, theta) {
  cauchy_expected_hessian_cpp(y, theta[[1]], theta[[2]])
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
#' The Cauchy distribution is a continuous probability distribution with heavy tails and no
#' defined mean or variance.
#'
#' \strong{Density function:}
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right]}}
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
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