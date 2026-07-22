#' @include distrib.R generics.R

#' @title S7 Class for Student's t Distribution
#' @name StudentTDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Student's t distribution.
#' @inheritParams distrib
#' @seealso \code{\link{student_t_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.StudentTDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.StudentTDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.StudentTDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.StudentTDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.StudentTDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.StudentTDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.StudentTDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.StudentTDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.StudentTDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.StudentTDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.StudentTDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
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
#' Computes the cumulative distribution function for the (location-scale) Student's t
#' distribution:
#' \deqn{F(q; \mu, \sigma, \nu) = T_\nu\!\left(\dfrac{q-\mu}{\sigma}\right)}
#' where \eqn{T_\nu} is the CDF of the standard Student's t with \eqn{\nu} degrees of freedom.
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
#' Computes the quantile function (inverse CDF) for the (location-scale) Student's t
#' distribution:
#' \deqn{Q(p; \mu, \sigma, \nu) = \mu + \sigma\, T_\nu^{-1}(p)}
#' where \eqn{T_\nu^{-1}} is the standard Student's t quantile function with \eqn{\nu} degrees of freedom.
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
S7::method(distrib_gradient, StudentTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
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
S7::method(distrib_hessian, StudentTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
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
S7::method(distrib_expected_hessian, StudentTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  student_t_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Student's t Analytical Third-Order Derivatives
#' @name distrib_deriv3.StudentTDistrib
#' @description
#' Closed-form observed third-order derivatives of the Student's t log-density. The
#' expected third derivatives have no closed form, so \code{expected = TRUE} falls
#' back to the numerical \code{\link{expectation}} of the observed derivatives.
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param expected Logical; if \code{TRUE}, returns the (numerically integrated) expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_deriv3, StudentTDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    student_t_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

#' @title Student's t Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.StudentTDistrib
#' @description
#' Closed-form observed fourth-order derivatives of the Student's t log-density. The
#' expected fourth derivatives have no closed form, so \code{expected = TRUE} falls
#' back to the numerical \code{\link{expectation}} of the observed derivatives.
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param expected Logical; if \code{TRUE}, returns the (numerically integrated) expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_deriv4, StudentTDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    student_t_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

#' @title Student's t Response Derivatives
#' @name distrib_grad_y.StudentTDistrib
#' @description
#' Closed-form derivatives of the Student's t log-density with respect to the
#' response. Let \eqn{r = y - \mu} and \eqn{d = \nu\sigma^2 + r^2}:
#' \eqn{\partial \ell / \partial y = -(\nu+1)r/d} and
#' \eqn{\partial^2 \ell / \partial y^2 = (\nu+1)(r^2 - \nu\sigma^2)/d^2}.
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A numeric vector.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_grad_y, StudentTDistrib) <- function(distrib, y, theta) {
  r <- y - theta[[1]]
  nu <- theta[[3]]
  -(nu + 1) * r / (nu * theta[[2]]^2 + r^2)
}

#' @title Student's t Response Second Derivative
#' @name distrib_hess_y.StudentTDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = (\nu+1)(r^2 - \nu\sigma^2)/(\nu\sigma^2 + r^2)^2}, \eqn{r = y - \mu}.
#' @param distrib A \code{StudentTDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A numeric vector.
#' @seealso \code{\link{student_t_distrib}}
S7::method(distrib_hess_y, StudentTDistrib) <- function(distrib, y, theta) {
  r <- y - theta[[1]]
  nu <- theta[[3]]
  vs2 <- nu * theta[[2]]^2
  (nu + 1) * (r^2 - vs2) / (vs2 + r^2)^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Student's t Distribution Object (Location-Scale Parameterization)
#'
#' @description
#' Creates a distribution object for the Student's t distribution parameterized by location (\eqn{\mu}), scale (\eqn{\sigma}), and degrees of freedom (\eqn{\nu}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#' @param link_nu A link function object for the degrees of freedom parameter \eqn{\nu}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The (location-scale) Student's t distribution has location \eqn{\mu}, scale
#' \eqn{\sigma} and degrees of freedom \eqn{\nu}. Write \eqn{d = \nu\sigma^2 + (y-\mu)^2}.
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\left(\dfrac{\nu}{2}\right)} \left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}}
#'
#' \strong{Cumulative distribution function} (\eqn{T_\nu} the standard t CDF):
#' \deqn{F(q; \mu, \sigma, \nu) = T_\nu\!\left(\dfrac{q-\mu}{\sigma}\right)}
#'
#' \strong{Quantile function:}
#' \deqn{Q(p; \mu, \sigma, \nu) = \mu + \sigma\,T_\nu^{-1}(p)}
#'
#' \strong{Score} (\eqn{\psi} the digamma function):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{(\nu+1)(y-\mu)}{d}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{\nu\left[(y-\mu)^2 - \sigma^2\right]}{\sigma d}}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left[ -\dfrac{1}{\nu} - \psi\left(\dfrac{\nu}{2}\right) + \psi\left(\dfrac{\nu+1}{2}\right) + \dfrac{(\nu+1)(y-\mu)^2}{\nu d} - \log\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right) \right]}
#'
#' \strong{Expected Hessian} (\eqn{\psi_1} the trigamma function; \eqn{\mu} is orthogonal to \eqn{\sigma, \nu}):
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2\nu}{\sigma^2(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \nu^2}\right] = \dfrac{1}{4}\left[\psi_1\left(\dfrac{\nu+1}{2}\right) - \psi_1\left(\dfrac{\nu}{2}\right)\right] + \dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma\,\partial \nu}\right] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}}
#' The observed Hessian is available via \code{\link{distrib_hessian.StudentTDistrib}}.
#'
#' \strong{Moments} (defined for \eqn{\nu} large enough): mean \eqn{\mu} (\eqn{\nu>1}),
#' variance \eqn{\sigma^2\nu/(\nu-2)} (\eqn{\nu>2}), skewness 0 (\eqn{\nu>3}),
#' excess kurtosis \eqn{6/(\nu-4)} (\eqn{\nu>4}).
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\nu \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order observed derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}; the expected ones use the numerical fallback) and
#' response derivatives (\code{\link{distrib_grad_y}}, \code{\link{distrib_hess_y}})
#' are also available.
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
#' @importFrom linkfunctions7 identity_link log_link
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
