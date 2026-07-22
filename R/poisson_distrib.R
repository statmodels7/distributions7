#' @include distrib.R generics.R

#' @title S7 Class for Poisson Distribution
#' @name PoissonDistrib
#' 
#' @description A subclass of \code{discrete_distrib} representing the Poisson distribution.
#' @inheritParams distrib
#' @seealso \code{\link{poisson_distrib}}
PoissonDistrib <- S7::new_class("PoissonDistrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Poisson Probability Mass Function
#' @name distrib_pdf.PoissonDistrib
#' @description
#' Computes the probability mass function for the Poisson distribution:
#' \deqn{P(Y=y; \mu) = \dfrac{\mu^y e^{-\mu}}{y!}}
#' 
#' @param distrib A \code{PoissonDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_pdf, PoissonDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dpois(
    x = y,
    lambda = theta[[1]],
    log = log
  )
}

#' @title Poisson Cumulative Distribution Function
#' @name distrib_cdf.PoissonDistrib
#' @description
#' Computes the cumulative distribution function for the Poisson distribution:
#' \deqn{F(q; \mu) = \sum_{k=0}^{\lfloor q \rfloor} \dfrac{\mu^k e^{-\mu}}{k!}}
#'
#' @param distrib A \code{PoissonDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_cdf, PoissonDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::ppois(
    q = q,
    lambda = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Poisson Quantile Function
#' @name distrib_quantile.PoissonDistrib
#' @description
#' Computes the quantile function for the Poisson distribution, defined as the
#' generalized inverse of the CDF:
#' \deqn{Q(p; \mu) = \min\left\{y \in \mathbb{N}_0 : F(y; \mu) \ge p\right\}}
#'
#' @param distrib A \code{PoissonDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_quantile, PoissonDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qpois(
    p = p,
    lambda = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Poisson Random Number Generator
#' @name distrib_rng.PoissonDistrib
#' @description
#' Generates random numbers from the Poisson distribution.
#' 
#' @param distrib A \code{PoissonDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameter \code{mu}.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_rng, PoissonDistrib) <- function(distrib, n, theta) {
  stats::rpois(
    n = n,
    lambda = theta[[1]]
  )
}

#' @title Poisson Analytical Gradient
#' @name distrib_gradient.PoissonDistrib
#' @description
#' Computes the analytical gradient (first derivative) of the Poisson log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu}}
#' 
#' @param distrib A \code{PoissonDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of first derivatives.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_gradient, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  poisson_gradient_cpp(y, theta[[1]])
}

#' @title Poisson Analytical Observed Hessian
#' @name distrib_hessian.PoissonDistrib
#' @description
#' Computes the analytical observed Hessian (second derivative) of the Poisson log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}}
#' 
#' @param distrib A \code{PoissonDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of second derivatives.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_hessian, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  poisson_hessian_cpp(y, theta[[1]])
}

#' @title Poisson Analytical Expected Hessian
#' @name distrib_expected_hessian.PoissonDistrib
#' @description
#' Computes the analytical expected Hessian of the Poisson log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu}}
#' 
#' @param distrib A \code{PoissonDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of expected second derivatives.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_expected_hessian, PoissonDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  poisson_expected_hessian_cpp(y, theta[[1]])
}

#' @title Poisson Analytical Third-Order Derivatives
#' @name distrib_deriv3.PoissonDistrib
#' @description Closed-form third-order derivative of the Poisson log-mass (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{PoissonDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivative.
#' @return A named list with the \code{mu_mu_mu} component.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_deriv3, PoissonDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) poisson_deriv3_expected_cpp(y, theta[[1]])
  else poisson_deriv3_cpp(y, theta[[1]])
}

#' @title Poisson Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.PoissonDistrib
#' @description Closed-form fourth-order derivative of the Poisson log-mass (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{PoissonDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivative.
#' @return A named list with the \code{mu_mu_mu_mu} component.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_deriv4, PoissonDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) poisson_deriv4_expected_cpp(y, theta[[1]])
  else poisson_deriv4_cpp(y, theta[[1]])
}

# --- CONSTRUCTOR WRAPPER ---

#' Poisson Distribution Object
#'
#' @description
#' Creates a distribution object for the Poisson distribution parameterized by the mean parameter \eqn{\mu}.
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The Poisson distribution models counts \eqn{y \in \{0, 1, 2, \dots\}} with mean
#' (and variance) \eqn{\mu}.
#'
#' \strong{Probability mass function:}
#' \deqn{P(Y=y; \mu) = \dfrac{\mu^y e^{-\mu}}{y!}}
#'
#' \strong{Cumulative distribution function:}
#' \deqn{F(q; \mu) = \sum_{k=0}^{\lfloor q \rfloor} \dfrac{\mu^k e^{-\mu}}{k!}}
#'
#' \strong{Quantile function:} the generalized inverse
#' \eqn{Q(p) = \min\{y \in \mathbb{N}_0 : F(y) \ge p\}}.
#'
#' \strong{Score, observed and expected Hessian:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu}}
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\mu}, skewness \eqn{1/\sqrt{\mu}},
#' excess kurtosis \eqn{1/\mu}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.PoissonDistrib}} for the probability mass function.
#'   \item \code{\link{distrib_cdf.PoissonDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.PoissonDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.PoissonDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.PoissonDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.PoissonDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.PoissonDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{PoissonDistrib} (inheriting from \code{discrete_distrib}) representing the Poisson distribution.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dpois ppois qpois rpois
#' @export
poisson_distrib <- function(link_mu = log_link()) {
  PoissonDistrib(
    distrib_name = "poisson", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)), link_params = list(mu = link_mu)
  )
}