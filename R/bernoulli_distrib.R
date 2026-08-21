#' @include distrib.R generics.R
NULL

#' @title S7 Class for Bernoulli Distribution
#' @name BernoulliDistrib
#' 
#' @description A subclass of \code{discrete_distrib} representing the Bernoulli distribution.
#' @inheritParams distrib
#' @return An object of class \code{BernoulliDistrib}.
#' @seealso \code{\link{bernoulli_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.BernoulliDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.BernoulliDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.BernoulliDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.BernoulliDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.BernoulliDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.BernoulliDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.BernoulliDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.BernoulliDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.BernoulliDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{discrete_distrib}}.
BernoulliDistrib <- S7::new_class("BernoulliDistrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Bernoulli Probability Mass Function
#' @name distrib_pdf.BernoulliDistrib
#' @description
#' Computes the probability mass function for the Bernoulli distribution:
#' \deqn{P(Y=y; \mu) = \mu^y (1-\mu)^{1-y}}
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param y A numeric vector of observations (\code{0} or \code{1}).
#' @param theta A list containing the parameter \code{mu}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_pdf, BernoulliDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbinom(
    x = y,
    size = 1,
    prob = theta[[1]],
    log = log
  )
}

#' @title Bernoulli Cumulative Distribution Function
#' @name distrib_cdf.BernoulliDistrib
#' @description
#' Computes the cumulative distribution function for the Bernoulli distribution:
#' \deqn{F(q; \mu) = \begin{cases} 0 & q < 0 \\ 1-\mu & 0 \le q < 1 \\ 1 & q \ge 1 \end{cases}}
#'
#' @param distrib A \code{BernoulliDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_cdf, BernoulliDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbinom(
    q = q,
    size = 1,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Bernoulli Quantile Function
#' @name distrib_quantile.BernoulliDistrib
#' @description
#' Computes the quantile function for the Bernoulli distribution, the generalized
#' inverse of the CDF:
#' \deqn{Q(p; \mu) = \begin{cases} 0 & p \le 1-\mu \\ 1 & p > 1-\mu \end{cases}}
#'
#' @param distrib A \code{BernoulliDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_quantile, BernoulliDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbinom(
    p = p,
    size = 1,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Bernoulli Random Number Generator
#' @name distrib_rng.BernoulliDistrib
#' @description
#' Generates random numbers from the Bernoulli distribution.
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameter \code{mu}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_rng, BernoulliDistrib) <- function(distrib, n, theta) {
  stats::rbinom(
    n = n,
    size = 1,
    prob = theta[[1]]
  )
}

#' @title Bernoulli Analytical Gradient
#' @name distrib_gradient.BernoulliDistrib
#' @description
#' Computes the analytical gradient (first derivative) of the Bernoulli log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1-\mu)}}
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of first derivatives.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_gradient, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  bernoulli_gradient_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Analytical Observed Hessian
#' @name distrib_hessian.BernoulliDistrib
#' @description
#' Computes the analytical observed Hessian (second derivative) of the Bernoulli log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{1-y}{(1-\mu)^2}}
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of second derivatives.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_hessian, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  bernoulli_hessian_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Analytical Expected Hessian
#' @name distrib_expected_hessian.BernoulliDistrib
#' @description
#' Computes the analytical expected Hessian of the Bernoulli log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu(1-\mu)}}
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of expected second derivatives.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_expected_hessian, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  bernoulli_expected_hessian_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Analytical Third-Order Derivatives
#' @name distrib_deriv3.BernoulliDistrib
#' @description Closed-form third-order derivative of the Bernoulli log-mass (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{BernoulliDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivative.
#' @return A named list with the \code{mu_mu_mu} component.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_deriv3, BernoulliDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) bernoulli_deriv3_expected_cpp(y, theta[[1]], threads)
  else bernoulli_deriv3_cpp(y, theta[[1]], threads)
}

#' @title Bernoulli Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.BernoulliDistrib
#' @description Closed-form fourth-order derivative of the Bernoulli log-mass (observed, or expected when \code{expected = TRUE}).
#' @param distrib A \code{BernoulliDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivative.
#' @return A named list with the \code{mu_mu_mu_mu} component.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_deriv4, BernoulliDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ..., threads = 1L) {
  if (expected) bernoulli_deriv4_expected_cpp(y, theta[[1]], threads)
  else bernoulli_deriv4_cpp(y, theta[[1]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Bernoulli Distribution Object
#'
#' @description
#' Creates a distribution object for the Bernoulli distribution parameterized by the probability of success \eqn{\mu}.
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu} (probability).
#'   Defaults to \code{\link[linkfunctions7]{logit_link}} to ensure the parameter stays within (0, 1).
#'
#' @details
#' The Bernoulli distribution models a binary outcome \eqn{y \in \{0, 1\}} with
#' success probability \eqn{\mu}.
#'
#' \strong{Probability mass function:}
#' \deqn{P(Y=y; \mu) = \mu^y (1-\mu)^{1-y}, \quad y \in \{0, 1\}}
#'
#' \strong{Cumulative distribution function:}
#' \deqn{F(q; \mu) = \begin{cases} 0 & q < 0 \\ 1-\mu & 0 \le q < 1 \\ 1 & q \ge 1 \end{cases}}
#'
#' \strong{Quantile function} (generalized inverse):
#' \deqn{Q(p; \mu) = \begin{cases} 0 & p \le 1-\mu \\ 1 & p > 1-\mu \end{cases}}
#'
#' \strong{Score, observed and expected Hessian:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu(1-\mu)}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{1-y}{(1-\mu)^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu(1-\mu)}}
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\mu(1-\mu)},
#' skewness \eqn{(1-2\mu)/\sqrt{\mu(1-\mu)}}, excess kurtosis
#' \eqn{(1 - 6\mu(1-\mu))/(\mu(1-\mu))}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.BernoulliDistrib}} for the probability mass function.
#'   \item \code{\link{distrib_cdf.BernoulliDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.BernoulliDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.BernoulliDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.BernoulliDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.BernoulliDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.BernoulliDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{BernoulliDistrib} (inheriting from \code{discrete_distrib}) representing the Bernoulli distribution.
#'
#' @importFrom linkfunctions7 logit_link
#' @importFrom stats dbinom pbinom qbinom rbinom
#' @examples
#' d <- bernoulli_distrib()
#' d@params
#'
#' theta <- list(mu = 0.3)
#' distrib_pdf(d, c(0, 1), theta)
#' distrib_gradient(d, c(0, 1), theta)
#'
#' @export
bernoulli_distrib <- function(link_mu = logit_link()) {
  
  BernoulliDistrib(
    distrib_name = "bernoulli",
    dimension = "univariate",
    bounds = c(0, 1),
    
    params = c("mu"),
    params_interpretation = c(mu = "probability"),
    n_params = 1,
    
    params_bounds = list(
      mu = c(0, 1)
    ),
    
    link_params = list(
      mu = link_mu
    )
  )
  
}
