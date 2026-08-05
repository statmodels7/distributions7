#' @include distrib.R generics.R
NULL

#' @title S7 Class for Lognormal Distribution
#' @name Lognormal1Distrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Lognormal distribution.
#' @inheritParams distrib
#' @return An object of class \code{Lognormal1Distrib}.
#' @seealso \code{\link{lognormal1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.Lognormal1Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.Lognormal1Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Lognormal1Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.Lognormal1Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.Lognormal1Distrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.Lognormal1Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.Lognormal1Distrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.Lognormal1Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.Lognormal1Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.Lognormal1Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.Lognormal1Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
Lognormal1Distrib <- S7::new_class("Lognormal1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Lognormal Probability Density Function
#' @name distrib_pdf.Lognormal1Distrib
#' @description
#' Computes the probability density function for the Lognormal distribution:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\{-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\}}
#' 
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_pdf, Lognormal1Distrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dlnorm(
    x = y,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    log = log
  )
}

#' @title Lognormal Cumulative Distribution Function
#' @name distrib_cdf.Lognormal1Distrib
#' @description
#' Computes the cumulative distribution function for the Lognormal distribution
#' (with \eqn{\sigma = \sqrt{\sigma^2}} on the log scale):
#' \deqn{F(q; \mu, \sigma^2) = \Phi\left(\dfrac{\log q - \mu}{\sqrt{\sigma^2}}\right), \quad q > 0}
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_cdf, Lognormal1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::plnorm(
    q = q,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Lognormal Quantile Function
#' @name distrib_quantile.Lognormal1Distrib
#' @description
#' Computes the quantile function (inverse CDF) for the Lognormal distribution:
#' \deqn{Q(p; \mu, \sigma^2) = \exp\left(\mu + \sqrt{\sigma^2}\,\Phi^{-1}(p)\right)}
#' where \eqn{\Phi^{-1}} is the standard normal quantile function.
#'
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_quantile, Lognormal1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qlnorm(
    p = p,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]]),
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Lognormal Random Number Generator
#' @name distrib_rng.Lognormal1Distrib
#' @description
#' Generates random numbers from the Lognormal distribution.
#' 
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_rng, Lognormal1Distrib) <- function(distrib, n, theta) {
  stats::rlnorm(
    n = n,
    meanlog = theta[[1]],
    sdlog = sqrt(theta[[2]])
  )
}

#' @title Lognormal Analytical Gradient
#' @name distrib_gradient.Lognormal1Distrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Lognormal log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_gradient, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Analytical Observed Hessian
#' @name distrib_hessian.Lognormal1Distrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Lognormal log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_hessian, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  lognormal_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Analytical Expected Hessian
#' @name distrib_expected_hessian.Lognormal1Distrib
#' @description
#' Computes the analytical expected Hessian of the Lognormal log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\sigma^2}.
#' 
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_expected_hessian, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  lognormal_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Analytical Third-Order Derivatives
#' @name distrib_deriv3.Lognormal1Distrib
#' @description Closed-form third-order derivatives of the Lognormal log-density (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_deriv3, Lognormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) lognormal_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  else lognormal_deriv3_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.Lognormal1Distrib
#' @description Closed-form fourth-order derivatives of the Lognormal log-density (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_deriv4, Lognormal1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) lognormal_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  else lognormal_deriv4_cpp(y, theta[[1]], theta[[2]])
}

#' @title Lognormal Response Derivatives
#' @name distrib_grad_y.Lognormal1Distrib
#' @description
#' Closed-form derivatives of the Lognormal log-density with respect to the
#' response. With \eqn{r = \log y - \mu}:
#' \eqn{\partial \ell / \partial y = -\dfrac{1}{y}\left(1 + \dfrac{r}{\sigma^2}\right)} and
#' \eqn{\partial^2 \ell / \partial y^2 = \dfrac{1}{y^2}\left(1 + \dfrac{r-1}{\sigma^2}\right)}.
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A numeric vector.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_grad_y, Lognormal1Distrib) <- function(distrib, y, theta) {
  r <- log(y) - theta[[1]]
  s2 <- theta[[2]]
  -(1 + r / s2) / y
}

#' @title Lognormal Response Second Derivative
#' @name distrib_hess_y.Lognormal1Distrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = (1 + (\log y - \mu - 1)/\sigma^2)/y^2}.
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma2}.
#' @return A numeric vector.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_hess_y, Lognormal1Distrib) <- function(distrib, y, theta) {
  r <- log(y) - theta[[1]]
  s2 <- theta[[2]]
  (1 + (r - 1) / s2) / y^2
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
#' The Lognormal distribution describes a variable whose logarithm is Gaussian with
#' mean \eqn{\mu} and variance \eqn{\sigma^2} (both on the log scale). Write
#' \eqn{r = \log y - \mu}.
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\{-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\}, \quad y > 0}
#'
#' \strong{Cumulative distribution function} (\eqn{\Phi} the standard normal CDF):
#' \deqn{F(q; \mu, \sigma^2) = \Phi\left(\dfrac{\log q - \mu}{\sqrt{\sigma^2}}\right)}
#'
#' \strong{Quantile function:}
#' \deqn{Q(p; \mu, \sigma^2) = \exp\left(\mu + \sqrt{\sigma^2}\,\Phi^{-1}(p)\right)}
#'
#' \strong{Score:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma^2} = \dfrac{r^2 - \sigma^2}{2\sigma^4}}
#'
#' \strong{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial (\sigma^2)^2} = \dfrac{1}{2\sigma^4} - \dfrac{r^2}{\sigma^6}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma^2} = -\dfrac{r}{\sigma^4}}
#'
#' \strong{Expected Hessian:}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{1}{2\sigma^4}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma^2}\right] = 0}
#'
#' \strong{Moments:} mean \eqn{e^{\mu + \sigma^2/2}}, variance
#' \eqn{(e^{\sigma^2}-1)e^{2\mu+\sigma^2}}, skewness \eqn{(e^{\sigma^2}+2)\sqrt{e^{\sigma^2}-1}},
#' excess kurtosis \eqn{e^{4\sigma^2}+2e^{3\sigma^2}+3e^{2\sigma^2}-6}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma^2 \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) and response derivatives (\code{\link{distrib_grad_y}},
#' \code{\link{distrib_hess_y}}) are also available.
#'
#' @return An S7 object of class \code{Lognormal1Distrib} (inheriting from \code{continuous_distrib}) representing the Lognormal distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dlnorm plnorm qlnorm rlnorm
#' @examples
#' d <- lognormal1_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma2 = 1)
#' distrib_pdf(d, c(0.5, 1, 2), theta)
#' distrib_gradient(d, c(0.5, 1, 2), theta)
#'
#' @export
lognormal1_distrib <- function(link_mu = identity_link(), link_sigma2 = log_link()) {
  
  Lognormal1Distrib(
    distrib_name = "lognormal1",
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
