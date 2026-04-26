#' @include distrib.R generics.R

#' @title S7 Class for Logistic Distribution
#' @name LogisticDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Logistic distribution.
#' @inheritParams distrib
#' @seealso \code{\link{logistic_distrib}}
LogisticDistrib <- S7::new_class("LogisticDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Logistic Probability Density Function
#' @name distrib_pdf.LogisticDistrib
#' @description
#' Computes the probability density function for the Logistic distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left[1 + \exp\left(-\dfrac{y-\mu}{\sigma}\right)\right]^2}}
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_pdf, LogisticDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dlogis(
    x = y,
    location = theta[[1]],
    scale = theta[[2]],
    log = log
  )
}

#' @title Logistic Cumulative Distribution Function
#' @name distrib_cdf.LogisticDistrib
#' @description
#' Computes the cumulative distribution function for the Logistic distribution.
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_cdf, LogisticDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::plogis(
    q = q,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Logistic Quantile Function
#' @name distrib_quantile.LogisticDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Logistic distribution.
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_quantile, LogisticDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qlogis(
    p = p,
    location = theta[[1]],
    scale = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Logistic Random Number Generator
#' @name distrib_rng.LogisticDistrib
#' @description
#' Generates random numbers from the Logistic distribution.
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_rng, LogisticDistrib) <- function(distrib, n, theta) {
  stats::rlogis(
    n = n,
    location = theta[[1]],
    scale = theta[[2]]
  )
}

#' @title Logistic Analytical Gradient
#' @name distrib_gradient.LogisticDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Logistic log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right)}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma} \left[ 1 - \dfrac{y-\mu}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_gradient, LogisticDistrib) <- function(distrib, y, theta) {
  logistic_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Logistic Analytical Observed Hessian
#' @name distrib_hessian.LogisticDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Logistic log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{2\sigma^2} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right)}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1}{\sigma^2} \left[ 1 - \dfrac{(y-\mu)^2}{2\sigma^2} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) - \dfrac{2(y-\mu)}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{1}{\sigma^2} \left[ \tanh\left(\dfrac{y-\mu}{2\sigma}\right) + \dfrac{y-\mu}{2\sigma} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_hessian, LogisticDistrib) <- function(distrib, y, theta) {
  logistic_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Logistic Analytical Expected Hessian
#' @name distrib_expected_hessian.LogisticDistrib
#' @description
#' Computes the analytical expected Hessian of the Logistic log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{3\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{3+\pi^2}{9\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#' 
#' @param distrib A \code{LogisticDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_expected_hessian, LogisticDistrib) <- function(distrib, y, theta) {
  logistic_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

# --- CONSTRUCTOR WRAPPER ---

#' Logistic Distribution Object
#'
#' @description
#' Creates a distribution object for the Logistic distribution parameterized by location (\eqn{\mu})
#' and scale (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.LogisticDistrib}} for the probability density function.
#'   \item \code{\link{distrib_cdf.LogisticDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.LogisticDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.LogisticDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.LogisticDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.LogisticDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.LogisticDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{LogisticDistrib} (inheriting from \code{continuous_distrib}) representing the Logistic distribution.
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dlogis plogis qlogis rlogis
#' @export
logistic_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  LogisticDistrib(
    distrib_name = "logistic", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "mean", sigma = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma)
  )
}