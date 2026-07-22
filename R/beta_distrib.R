#' @include distrib.R generics.R

#' @title S7 Class for Beta Distribution
#' @name BetaDistrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Beta distribution 
#' under the mean-precision parameterization.
#' @inheritParams distrib
#' @seealso \code{\link{beta_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.BetaDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.BetaDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.BetaDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.BetaDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.BetaDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.BetaDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.BetaDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.BetaDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.BetaDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.BetaDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.BetaDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
BetaDistrib <- S7::new_class("BetaDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta Probability Density Function
#' @name distrib_pdf.BetaDistrib
#' @description
#' Computes the probability density function for the Beta distribution:
#' \deqn{f(y; \mu, \phi) = \dfrac{\Gamma(\phi)}{\Gamma(\mu\phi)\Gamma((1-\mu)\phi)} y^{\mu\phi-1} (1-y)^{(1-\mu)\phi-1}}
#' 
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_pdf, BetaDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dbeta(
    x = y,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    log = log
  )
}

#' @title Beta Cumulative Distribution Function
#' @name distrib_cdf.BetaDistrib
#' @description
#' Computes the cumulative distribution function for the Beta distribution, using the
#' mean/precision reparameterization \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi}:
#' \deqn{F(q; \mu, \phi) = I_q(\alpha, \beta)}
#' where \eqn{I_q(\cdot, \cdot)} is the regularized incomplete beta function.
#'
#' @param distrib A \code{BetaDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_cdf, BetaDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbeta(
    q = q,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Beta Quantile Function
#' @name distrib_quantile.BetaDistrib
#' @description
#' Computes the quantile function for the Beta distribution as the inverse of the CDF,
#' \eqn{Q(p; \mu, \phi) = F^{-1}(p; \mu, \phi)}. There is no elementary closed form; it
#' is obtained numerically (via \code{\link[stats]{qbeta}}) on the mean/precision
#' reparameterization \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi}.
#'
#' @param distrib A \code{BetaDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_quantile, BetaDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbeta(
    p = p,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Beta Random Number Generator
#' @name distrib_rng.BetaDistrib
#' @description
#' Generates random numbers from the Beta distribution.
#' 
#' @param distrib A \code{BetaDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_rng, BetaDistrib) <- function(distrib, n, theta) {
  stats::rbeta(
    n = n,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]]
  )
}

#' @title Beta Analytical Gradient
#' @name distrib_gradient.BetaDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Beta log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \phi \left[ \log\left(\dfrac{y}{1-y}\right) - \psi(\mu\phi) + \psi((1-\mu)\phi) \right]}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\mu\phi) - (1-\mu)\psi((1-\mu)\phi) + \mu \log(y) + (1-\mu) \log(1-y)}
#' 
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_gradient, BetaDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  beta_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Beta Analytical Observed Hessian
#' @name distrib_hessian.BetaDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Beta log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_hessian, BetaDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  beta_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Beta Analytical Expected Hessian
#' @name distrib_expected_hessian.BetaDistrib
#' @description
#' Computes the analytical expected Hessian of the Beta log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\phi^2 \left[ \psi_1(\mu\phi) + \psi_1((1-\mu)\phi) \right]}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \phi^2}\right] = \psi_1(\phi) - \mu^2\psi_1(\mu\phi) - (1-\mu)^2\psi_1((1-\mu)\phi)}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \phi}\right] = -\phi \left[ \mu\psi_1(\mu\phi) - (1-\mu)\psi_1((1-\mu)\phi) \right]}
#' 
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_expected_hessian, BetaDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  beta_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Beta Analytical Third-Order Derivatives
#' @name distrib_deriv3.BetaDistrib
#' @description
#' Closed-form third-order derivatives of the Beta log-density. Because the
#' \eqn{y}-terms of the log-density are linear in the parameters, these derivatives
#' do not depend on \eqn{y}; the observed and expected values therefore coincide.
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param expected Logical; ignored (observed and expected coincide).
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_deriv3, BetaDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  beta_deriv3_cpp(y, theta[[1]], theta[[2]])
}

#' @title Beta Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.BetaDistrib
#' @description
#' Closed-form fourth-order derivatives of the Beta log-density. As for the third
#' order, they do not depend on \eqn{y}, so observed and expected coincide.
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param expected Logical; ignored (observed and expected coincide).
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_deriv4, BetaDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  beta_deriv4_cpp(y, theta[[1]], theta[[2]])
}

#' @title Beta Response Derivatives
#' @name distrib_grad_y.BetaDistrib
#' @description
#' Closed-form derivatives of the Beta log-density with respect to the response,
#' with \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi}:
#' \eqn{\partial \ell / \partial y = (\alpha-1)/y - (\beta-1)/(1-y)} and
#' \eqn{\partial^2 \ell / \partial y^2 = -(\alpha-1)/y^2 - (\beta-1)/(1-y)^2}.
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A numeric vector.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_grad_y, BetaDistrib) <- function(distrib, y, theta) {
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  (a - 1) / y - (b - 1) / (1 - y)
}

#' @title Beta Response Second Derivative
#' @name distrib_hess_y.BetaDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -(\alpha-1)/y^2 - (\beta-1)/(1-y)^2}.
#' @param distrib A \code{BetaDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A numeric vector.
#' @seealso \code{\link{beta_distrib}}
S7::method(distrib_hess_y, BetaDistrib) <- function(distrib, y, theta) {
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  -(a - 1) / y^2 - (b - 1) / (1 - y)^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Beta Distribution Object (Mean-Precision Parameterization)
#'
#' @description
#' Creates a distribution object for the Beta distribution parameterized by mean (\eqn{\mu}) and precision (\eqn{\phi}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{logit_link}} to ensure the parameter stays within (0, 1).
#' @param link_phi A link function object for the precision parameter \eqn{\phi}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The Beta distribution is given a mean/precision parameterization: \eqn{\mu} is
#' the mean and \eqn{\phi} a precision parameter. The standard shapes are
#' \deqn{\alpha = \mu\phi, \qquad \beta = (1-\mu)\phi}
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \phi) = \dfrac{\Gamma(\phi)}{\Gamma(\alpha)\Gamma(\beta)}\, y^{\alpha-1} (1-y)^{\beta-1}, \quad 0 < y < 1}
#'
#' \strong{Cumulative distribution function} (\eqn{I_q} the regularized incomplete beta function):
#' \deqn{F(q; \mu, \phi) = I_q(\alpha, \beta)}
#'
#' \strong{Quantile function:} no closed form; the numerical inverse of the CDF.
#'
#' \strong{Score} (\eqn{\psi} the digamma function):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \phi\left[\log\left(\dfrac{y}{1-y}\right) - \psi(\alpha) + \psi(\beta)\right]}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\alpha) - (1-\mu)\psi(\beta) + \mu\log y + (1-\mu)\log(1-y)}
#'
#' \strong{Expected Hessian} (\eqn{\psi_1} the trigamma function; the observed and
#' expected Hessians coincide, as they do not depend on \eqn{y}):
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\phi^2\left[\psi_1(\alpha) + \psi_1(\beta)\right], \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \phi} = -\phi\left[\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \phi^2} = \psi_1(\phi) - \mu^2\psi_1(\alpha) - (1-\mu)^2\psi_1(\beta)}
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\mu(1-\mu)/(\phi+1)},
#' skewness \eqn{\dfrac{2(1-2\mu)\sqrt{\phi+1}}{(\phi+2)\sqrt{\mu(1-\mu)}}}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#'   \item \eqn{\phi \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) and response derivatives (\code{\link{distrib_grad_y}},
#' \code{\link{distrib_hess_y}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.BetaDistrib}} for the probability density function.
#'   \item \code{\link{distrib_cdf.BetaDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.BetaDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.BetaDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.BetaDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.BetaDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.BetaDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{BetaDistrib} (inheriting from \code{continuous_distrib}) representing the Beta distribution.
#'
#' @importFrom linkfunctions7 logit_link log_link
#' @importFrom stats dbeta pbeta qbeta rbeta
#' @export
beta_distrib <- function(link_mu = logit_link(), link_phi = log_link()) {
  
  BetaDistrib(
    distrib_name = "beta",
    dimension = "univariate",
    bounds = c(0, 1),
    
    params = c("mu", "phi"),
    params_interpretation = c(mu = "mean", phi = "precision"),
    n_params = 2,
    
    params_bounds = list(
      mu = c(0, 1),
      phi = c(0, Inf)
    ),
    
    link_params = list(
      mu = link_mu,
      phi = link_phi
    )
  )
  
}
