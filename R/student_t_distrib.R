#' @include distrib.R generics.R

#' @title S7 Class for Student's t Distribution
#' @name StudentTDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Student's t distribution.
#' @inheritParams distrib
#' @seealso \code{\link{student_t_distrib}}
StudentTDistrib <- S7::new_class("StudentTDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Student's t Probability Density Function
#' @name distrib_pdf.StudentTDistrib
#' @description
#' Computes the probability density function for the Student's t distribution:
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\left(\dfrac{\nu}{2}\right)} \left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}}
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_pdf, StudentTDistrib) <- function(distrib, y, theta, log = FALSE) {
  val <- stats::dt(
    x = (y - theta[[1]]) / theta[[2]],
    df = theta[[3]],
    log = log
  )
  if (log) {
    val - log(theta[[2]])
  } else {
    val / theta[[2]]
  }
}

#' @title Student's t Cumulative Distribution Function
#' @name distrib_cdf.StudentTDistrib
#' @description
#' Computes the cumulative distribution function for the Student's t distribution.
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_cdf, StudentTDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pt(
    q = (q - theta[[1]]) / theta[[2]],
    df = theta[[3]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Student's t Quantile Function
#' @name distrib_quantile.StudentTDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Student's t distribution.
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_quantile, StudentTDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  theta[[1]] + theta[[2]] * stats::qt(
    p = p,
    df = theta[[3]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Student's t Random Number Generator
#' @name distrib_rng.StudentTDistrib
#' @description
#' Generates random numbers from the Student's t distribution.
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_rng, StudentTDistrib) <- function(distrib, n, theta) {
  theta[[1]] + theta[[2]] * stats::rt(
    n = n,
    df = theta[[3]]
  )
}

#' @title Student's t Analytical Gradient
#' @name distrib_gradient.StudentTDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Student's t log-density 
#' with respect to the parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\nu}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{(\nu+1)(y-\mu)}{\nu\sigma^2 + (y-\mu)^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{\nu\left[(y-\mu)^2 - \sigma^2\right]}{\sigma\left[\nu\sigma^2 + (y-\mu)^2\right]}}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left[ -\dfrac{1}{\nu} - \psi\left(\dfrac{\nu}{2}\right) + \psi\left(\dfrac{\nu+1}{2}\right) + \dfrac{(\nu+1)(y-\mu)^2}{\nu\left[\nu\sigma^2 + (y-\mu)^2\right]} - \log\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right) \right]}
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_gradient, StudentTDistrib) <- function(distrib, y, theta) {
  student_t_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Student's t Analytical Observed Hessian
#' @name distrib_hessian.StudentTDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Student's t log-density 
#' with respect to the parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\nu}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{(\nu+1)\left[(y-\mu)^2 - \nu\sigma^2\right]}{\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\nu\left[\nu\sigma^4 - (3\nu+1)\sigma^2(y-\mu)^2 - (y-\mu)^4\right]}{\sigma^2\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4}\left[ -\psi_1\left(\dfrac{\nu}{2}\right) + \psi_1\left(\dfrac{\nu+1}{2}\right) + \dfrac{2\left(\nu\sigma^4 + (y-\mu)^4\right)}{\nu\left[\nu\sigma^2 + (y-\mu)^2\right]^2} \right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{2\nu(\nu+1)\sigma(y-\mu)}{\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \nu} = \dfrac{(y-\mu)\left[(y-\mu)^2 - \sigma^2\right]}{\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu} = \dfrac{(y-\mu)^2\left[(y-\mu)^2 - \sigma^2\right]}{\sigma\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_hessian, StudentTDistrib) <- function(distrib, y, theta) {
  student_t_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Student's t Analytical Expected Hessian
#' @name distrib_expected_hessian.StudentTDistrib
#' @description
#' Computes the analytical expected Hessian of the Student's t log-density 
#' with respect to the parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\nu}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2\nu}{\sigma^2(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \nu^2}\right] = \dfrac{1}{4}\left[\psi_1\left(\dfrac{\nu+1}{2}\right) - \psi_1\left(\dfrac{\nu}{2}\right)\right] + \dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu}\right] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}}
#' 
#' The parameter \eqn{\mu} is orthogonal to \eqn{\sigma} and \eqn{\nu} (mixed expected derivatives are 0).
#' 
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma}, and \code{nu}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_expected_hessian, StudentTDistrib) <- function(distrib, y, theta) {
  student_t_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

# --- CONSTRUCTOR WRAPPER ---

#' Student's t Distribution Object (Location-Scale Parameterization)
#'
#' @description
#' Creates a distribution object for the Student's t distribution parameterized by location (\eqn{\mu}), scale (\eqn{\sigma}), and degrees of freedom (\eqn{\nu}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#' @param link_nu A link function object for the degrees of freedom parameter \eqn{\nu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The probability density function is given by:
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\left(\dfrac{\nu}{2}\right)} \left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}}
#' for \eqn{y \in (-\infty, +\infty)}.
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\nu \in (0, +\infty)}
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.StudentTDistrib}} for the probability density function.
#'   \item \code{\link{distrib_cdf.StudentTDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.StudentTDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.StudentTDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.StudentTDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.StudentTDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.StudentTDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{StudentTDistrib} (inheriting from \code{continuous_distrib}) representing the Student's t distribution.
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dt pt qt rt
#' @export
student_t_distrib <- function(link_mu = identity_link(), link_sigma = log_link(), link_nu = log_link()) {
  
  StudentTDistrib(
    distrib_name = "student t", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma", "nu"),
    params_interpretation = c(mu = "location", sigma = "scale", nu = "shape"),
    n_params = 3,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf), nu = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma, nu = link_nu)
  )
  
}