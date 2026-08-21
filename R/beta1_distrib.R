#' @include distrib.R generics.R
NULL

#' @title S7 Class for Beta Distribution
#' @name Beta1Distrib
#' 
#' @description A subclass of \code{continuous_distrib} representing the Beta distribution 
#' under the mean-precision parameterization.
#' @inheritParams distrib
#' @return An object of class \code{Beta1Distrib}.
#' @seealso \code{\link{beta1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.Beta1Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.Beta1Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Beta1Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.Beta1Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.Beta1Distrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.Beta1Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.Beta1Distrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.Beta1Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.Beta1Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.Beta1Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.Beta1Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
Beta1Distrib <- S7::new_class("Beta1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta Probability Density Function
#' @name distrib_pdf.Beta1Distrib
#' @description
#' Computes the probability density function for the Beta distribution:
#' \deqn{f(y; \mu, \phi) = \dfrac{\Gamma(\phi)}{\Gamma(\mu\phi)\Gamma((1-\mu)\phi)} y^{\mu\phi-1} (1-y)^{(1-\mu)\phi-1}}
#' 
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_pdf, Beta1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbeta(
    x = y,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    log = log
  )
}

#' @title Beta Cumulative Distribution Function
#' @name distrib_cdf.Beta1Distrib
#' @description
#' Computes the cumulative distribution function for the Beta distribution, using the
#' mean/precision reparameterization \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi}:
#' \deqn{F(q; \mu, \phi) = I_q(\alpha, \beta)}
#' where \eqn{I_q(\cdot, \cdot)} is the regularized incomplete beta function.
#'
#' @param distrib A \code{Beta1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_cdf, Beta1Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbeta(
    q = q,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Beta Quantile Function
#' @name distrib_quantile.Beta1Distrib
#' @description
#' Computes the quantile function for the Beta distribution as the inverse of the CDF,
#' \eqn{Q(p; \mu, \phi) = F^{-1}(p; \mu, \phi)}. There is no elementary closed form; it
#' is obtained numerically (via \code{\link[stats]{qbeta}}) on the mean/precision
#' reparameterization \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi}.
#'
#' @param distrib A \code{Beta1Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_quantile, Beta1Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbeta(
    p = p,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Beta Random Number Generator
#' @name distrib_rng.Beta1Distrib
#' @description
#' Generates random numbers from the Beta distribution.
#' 
#' @param distrib A \code{Beta1Distrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_rng, Beta1Distrib) <- function(distrib, n, theta) {
  stats::rbeta(
    n = n,
    shape1 = theta[[1]] * theta[[2]],
    shape2 = (1 - theta[[1]]) * theta[[2]]
  )
}

#' @title Beta Analytical Gradient
#' @name distrib_gradient.Beta1Distrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Beta log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \phi \left[ \log\left(\dfrac{y}{1-y}\right) - \psi(\mu\phi) + \psi((1-\mu)\phi) \right]}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\mu\phi) - (1-\mu)\psi((1-\mu)\phi) + \mu \log(y) + (1-\mu) \log(1-y)}
#' 
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_gradient, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  beta_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Analytical Observed Hessian
#' @name distrib_hessian.Beta1Distrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Beta log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_hessian, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  beta_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Analytical Expected Hessian
#' @name distrib_expected_hessian.Beta1Distrib
#' @description
#' Computes the analytical expected Hessian of the Beta log-density 
#' with respect to the parameters \eqn{\mu} and \eqn{\phi}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\phi^2 \left[ \psi_1(\mu\phi) + \psi_1((1-\mu)\phi) \right]}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \phi^2}\right] = \psi_1(\phi) - \mu^2\psi_1(\mu\phi) - (1-\mu)^2\psi_1((1-\mu)\phi)}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \phi}\right] = -\phi \left[ \mu\psi_1(\mu\phi) - (1-\mu)\psi_1((1-\mu)\phi) \right]}
#' 
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_expected_hessian, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  beta_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Analytical Third-Order Derivatives
#' @name distrib_deriv3.Beta1Distrib
#' @description
#' Closed-form third-order derivatives of the Beta log-density. Because the
#' \eqn{y}-terms of the log-density are linear in the parameters, these derivatives
#' do not depend on \eqn{y}; the observed and expected values therefore coincide.
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param expected Logical; ignored (observed and expected coincide).
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_deriv3, Beta1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  beta_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.Beta1Distrib
#' @description
#' Closed-form fourth-order derivatives of the Beta log-density. As for the third
#' order, they do not depend on \eqn{y}, so observed and expected coincide.
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @param expected Logical; ignored (observed and expected coincide).
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_deriv4, Beta1Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  beta_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Beta Response Derivatives
#' @name distrib_grad_y.Beta1Distrib
#' @description
#' Closed-form derivatives of the Beta log-density with respect to the response,
#' with \eqn{\alpha = \mu\phi}, \eqn{\beta = (1-\mu)\phi}:
#' \eqn{\partial \ell / \partial y = (\alpha-1)/y - (\beta-1)/(1-y)} and
#' \eqn{\partial^2 \ell / \partial y^2 = -(\alpha-1)/y^2 - (\beta-1)/(1-y)^2}.
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A numeric vector.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_grad_y, Beta1Distrib) <- function(distrib, y, theta) {
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  (a - 1) / y - (b - 1) / (1 - y)
}

#' @title Beta Response Second Derivative
#' @name distrib_hess_y.Beta1Distrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -(\alpha-1)/y^2 - (\beta-1)/(1-y)^2}.
#' @param distrib A \code{Beta1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{phi}.
#' @return A numeric vector.
#' @seealso \code{\link{beta1_distrib}}
S7::method(distrib_hess_y, Beta1Distrib) <- function(distrib, y, theta) {
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
#'   \item \code{\link{distrib_pdf.Beta1Distrib}} for the probability density function.
#'   \item \code{\link{distrib_cdf.Beta1Distrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.Beta1Distrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.Beta1Distrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.Beta1Distrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.Beta1Distrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.Beta1Distrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{Beta1Distrib} (inheriting from \code{continuous_distrib}) representing the Beta distribution.
#'
#' @importFrom linkfunctions7 logit_link log_link
#' @importFrom stats dbeta pbeta qbeta rbeta
#' @examples
#' d <- beta1_distrib()
#' d@params
#'
#' theta <- list(mu = 0.4, phi = 5)
#' distrib_pdf(d, c(0.2, 0.5, 0.8), theta)
#' distrib_gradient(d, c(0.2, 0.5, 0.8), theta)
#'
#' @export
beta1_distrib <- function(link_mu = logit_link(), link_phi = log_link()) {
  
  Beta1Distrib(
    distrib_name = "beta1",
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
