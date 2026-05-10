#' @include distrib.R generics.R

#' @title S7 Class for Bernoulli Distribution
#' @name BernoulliDistrib
#' 
#' @description A subclass of \code{discrete_distrib} representing the Bernoulli distribution.
#' @inheritParams distrib
#' @seealso \code{\link{bernoulli_distrib}}
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
S7::method(distrib_pdf, BernoulliDistrib) <- function(distrib, y, theta, log = FALSE) {
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
#' Computes the cumulative distribution function for the Bernoulli distribution.
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
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
#' Computes the quantile function (inverse CDF) for the Bernoulli distribution.
#' 
#' @param distrib A \code{BernoulliDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
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
S7::method(distrib_gradient, BernoulliDistrib) <- function(distrib, y, theta) {
  bernoulli_gradient_cpp(y, theta[[1]])
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
S7::method(distrib_hessian, BernoulliDistrib) <- function(distrib, y, theta) {
  bernoulli_hessian_cpp(y, theta[[1]])
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
S7::method(distrib_expected_hessian, BernoulliDistrib) <- function(distrib, y, theta) {
  bernoulli_expected_hessian_cpp(y, theta[[1]])
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
#' The Bernoulli distribution has the following probability mass function (PMF):
#' \deqn{P(Y=y; \mu) = \mu^y (1-\mu)^{1-y}}
#' for \eqn{y \in \{0, 1\}}.
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#' }
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