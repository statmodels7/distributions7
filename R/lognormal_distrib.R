#' @include distrib.R generics.R

#' @title S7 Class for Lognormal Distribution
#' @name LognormalDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Lognormal distribution.
#' @inheritParams distrib
#' @seealso \code{\link{lognormal_distrib}}
LognormalDistrib <- S7::new_class("LognormalDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Lognormal Probability Density Function
#' @name distrib_pdf.LognormalDistrib
#' @description
#' Computes the probability density function for the Lognormal distribution:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\{-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\}}
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_pdf, LognormalDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dlnorm(
    x = y,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    log = log
  )
}

#' @title Lognormal Cumulative Distribution Function
#' @name distrib_cdf.LognormalDistrib
#' @description
#' Computes the cumulative distribution function for the Lognormal distribution.
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_cdf, LognormalDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::plnorm(
    q = q,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Lognormal Quantile Function
#' @name distrib_quantile.LognormalDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Lognormal distribution.
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_quantile, LognormalDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qlnorm(
    p = p,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Lognormal Random Number Generator
#' @name distrib_rng.LognormalDistrib
#' @description
#' Generates random numbers from the Lognormal distribution.
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_rng, LognormalDistrib) <- function(distrib, n, theta) {
  stats::rlnorm(
    n = n,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]])
  )
}

#' @title Lognormal Analytical Gradient
#' @name distrib_gradient.LognormalDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Lognormal log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_gradient, LognormalDistrib) <- function(distrib, y, theta) {
  lognormal_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Analytical Observed Hessian
#' @name distrib_hessian.LognormalDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Lognormal log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_hessian, LognormalDistrib) <- function(distrib, y, theta) {
  lognormal_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Analytical Expected Hessian
#' @name distrib_expected_hessian.LognormalDistrib
#' @description
#' Computes the analytical expected Hessian of the Lognormal log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{LognormalDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{lognormal_distrib}}
S7::method(distrib_expected_hessian, LognormalDistrib) <- function(distrib, y, theta) {
  lognormal_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

# --- CONSTRUCTOR WRAPPER ---

#' Lognormal Distribution Object (Log-Scale Parameterization)
#'
#' @description
#' Creates a distribution object for the Lognormal distribution parameterized by the mean (\eqn{\mu}) and the variance (\eqn{\sigma^2}) of the log-transformed variable.
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{identity_link}}.
#' @param link_sigma2 A link function object for the variance parameter \eqn{\sigma^2}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma^2 \in (0, +\infty)}
#' }
#'
#' @return An S7 object of class \code{LognormalDistrib} (inheriting from \code{continuous_distrib}) representing the Lognormal distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dlnorm plnorm qlnorm rlnorm
#' @export
lognormal_distrib <- function(link_mu = identity_link(), link_sigma2 = log_link()) {
  
  LognormalDistrib(
    distrib_name = "lognormal",
    dimension = "univariate",
    bounds = c(0, Inf),
    
    params = c("mu", "sigma2"),
    params_interpretation = c(mu = "mean (log scale)", sigma2 = "variance (log scale)"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma2 = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      sigma2 = link_sigma2
    )
  )
  
}