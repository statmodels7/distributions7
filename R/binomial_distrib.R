#' @include distrib.R generics.R

#' @title S7 Class for Binomial Distribution
#' @name BinomialDistrib
#' 
#' @description A subclass of \code{discrete_distrib} representing the Binomial distribution.
#' @inheritParams distrib
#' @param size Integer or Numeric Vector. The number of trials \eqn{n}.
#'   Can be a single scalar (default 1) or a vector of the same length as the observations \eqn{y}.
#' @seealso \code{\link{binomial_distrib}}
BinomialDistrib <- S7::new_class("BinomialDistrib", 
  parent = discrete_distrib,
  properties = list(
    size = S7::class_numeric
  )
)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Binomial Probability Mass Function
#' @name distrib_pdf.BinomialDistrib
#' @description
#' Computes the probability mass function for the Binomial distribution:
#' \deqn{P(Y=y; \mu, n) = \dbinom{n}{y} \mu^y (1-\mu)^{n-y}}
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_pdf, BinomialDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dbinom(
    x = y,
    size = distrib@size,
    prob = theta[[1]],
    log = log
  )
}

#' @title Binomial Cumulative Distribution Function
#' @name distrib_cdf.BinomialDistrib
#' @description
#' Computes the cumulative distribution function for the Binomial distribution.
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_cdf, BinomialDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pbinom(
    q = q,
    size = distrib@size,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Binomial Quantile Function
#' @name distrib_quantile.BinomialDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Binomial distribution.
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_quantile, BinomialDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qbinom(
    p = p,
    size = distrib@size,
    prob = theta[[1]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Binomial Random Number Generator
#' @name distrib_rng.BinomialDistrib
#' @description
#' Generates random numbers from the Binomial distribution.
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameter \code{mu}.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_rng, BinomialDistrib) <- function(distrib, n, theta) {
  stats::rbinom(
    n = n,
    size = distrib@size,
    prob = theta[[1]]
  )
}

#' @title Binomial Analytical Gradient
#' @name distrib_gradient.BinomialDistrib
#' @description
#' Computes the analytical gradient (first derivative) of the Binomial log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - n\mu}{\mu(1-\mu)}}
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of first derivatives.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_gradient, BinomialDistrib) <- function(distrib, y, theta) {
  binomial_gradient_cpp(y, theta[[1]], distrib@size)
}

#' @title Binomial Analytical Observed Hessian
#' @name distrib_hessian.BinomialDistrib
#' @description
#' Computes the analytical observed Hessian (second derivative) of the Binomial log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{n-y}{(1-\mu)^2}}
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of second derivatives.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_hessian, BinomialDistrib) <- function(distrib, y, theta) {
  binomial_hessian_cpp(y, theta[[1]], distrib@size)
}

#' @title Binomial Analytical Expected Hessian
#' @name distrib_expected_hessian.BinomialDistrib
#' @description
#' Computes the analytical expected Hessian of the Binomial log-probability 
#' with respect to the parameter \eqn{\mu}.
#' 
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{n}{\mu(1-\mu)}}
#' 
#' @param distrib A \code{BinomialDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @return A list containing the vector of expected second derivatives.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_expected_hessian, BinomialDistrib) <- function(distrib, y, theta) {
  binomial_expected_hessian_cpp(y, theta[[1]], distrib@size)
}

# --- CONSTRUCTOR WRAPPER ---

#' Binomial Distribution Object
#'
#' @description
#' Creates a distribution object for the Binomial distribution parameterized by the probability of success \eqn{\mu} and a number of trials \eqn{n} (size).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu} (probability).
#'   Defaults to \code{\link[linkfunctions]{logit_link}} to ensure the parameter stays within (0, 1).
#' @param size Integer or Numeric Vector. The number of trials \eqn{n}.
#'   Can be a single scalar (default 1) or a vector of the same length as the observations \eqn{y}.
#'
#' @details
#' The Binomial distribution has the following probability mass function (PMF):
#' \deqn{P(Y=y; \mu, n) = \dbinom{n}{y} \mu^y (1-\mu)^{n-y}}
#' for \eqn{y \in \{0, 1, \dots, n\}}.
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#'   \item \eqn{n \in \mathbb{Z}^+} (fixed in constructor, can vary per observation)
#' }
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.BinomialDistrib}} for the probability mass function.
#'   \item \code{\link{distrib_cdf.BinomialDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.BinomialDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.BinomialDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.BinomialDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.BinomialDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.BinomialDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{BinomialDistrib} (inheriting from \code{discrete_distrib}) representing the Binomial distribution.
#'
#' @importFrom linkfunctions logit_link
#' @importFrom stats dbinom pbinom qbinom rbinom
#' @export
binomial_distrib <- function(link_mu = logit_link(), size = 1) {
  
  BinomialDistrib(
    distrib_name = "binomial",
    dimension = "univariate",
    bounds = c(0, max(size)),
    size = size,
    
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