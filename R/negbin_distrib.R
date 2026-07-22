#' @include distrib.R generics.R

#' @title S7 Class for Negative Binomial Distribution (NB2)
#' @name NegBinDistrib
#'
#' @description A subclass of \code{discrete_distrib} representing the Negative Binomial distribution (NB2 parameterization).
#' @inheritParams distrib
#' @seealso \code{\link{negbin_distrib}}
NegBinDistrib <- S7::new_class("NegBinDistrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Negative Binomial Probability Mass Function
#' @name distrib_pdf.NegBinDistrib
#' @description
#' Computes the probability mass function for the Negative Binomial distribution (NB2):
#' \deqn{P(Y=y; \mu, \theta) = \dfrac{\Gamma(y+\theta)}{y!\,\Gamma(\theta)} \left(\dfrac{\theta}{\theta+\mu}\right)^\theta \left(\dfrac{\mu}{\theta+\mu}\right)^y}
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_pdf, NegBinDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dnbinom(
    x = y,
    mu = theta[[1]],
    size = theta[[2]],
    log = log
  )
}

#' @title Negative Binomial Cumulative Distribution Function
#' @name distrib_cdf.NegBinDistrib
#' @description
#' Computes the cumulative distribution function for the Negative Binomial distribution:
#' \deqn{F(q; \mu, \theta) = \sum_{k=0}^{\lfloor q \rfloor} P(Y=k; \mu, \theta)}
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_cdf, NegBinDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::pnbinom(
    q = q,
    mu = theta[[1]],
    size = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Negative Binomial Quantile Function
#' @name distrib_quantile.NegBinDistrib
#' @description
#' Computes the quantile function for the Negative Binomial distribution, the
#' generalized inverse of the CDF:
#' \deqn{Q(p; \mu, \theta) = \min\left\{y \in \mathbb{N}_0 : F(y; \mu, \theta) \ge p\right\}}
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_quantile, NegBinDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qnbinom(
    p = p,
    mu = theta[[1]],
    size = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Negative Binomial Random Number Generator
#' @name distrib_rng.NegBinDistrib
#' @description
#' Generates random numbers from the Negative Binomial distribution.
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_rng, NegBinDistrib) <- function(distrib, n, theta) {
  stats::rnbinom(
    n = n,
    mu = theta[[1]],
    size = theta[[2]]
  )
}

#' @title Negative Binomial Analytical Gradient
#' @name distrib_gradient.NegBinDistrib
#' @description
#' Computes the analytical gradient (first derivatives) of the Negative Binomial log-probability
#' with respect to the parameters \eqn{\mu} and \eqn{\theta}.
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\theta}{\theta+\mu}\left(\dfrac{y}{\mu} - 1\right)}
#' \deqn{\dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) - \psi(\theta) + \log\left(\dfrac{\theta}{\theta+\mu}\right) + \dfrac{\mu - y}{\theta+\mu}}
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_gradient, NegBinDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  negbin_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Negative Binomial Analytical Observed Hessian
#' @name distrib_hessian.NegBinDistrib
#' @description
#' Computes the analytical observed Hessian (second derivatives) of the Negative Binomial log-probability
#' with respect to the parameters \eqn{\mu} and \eqn{\theta}.
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{y+\theta}{(\theta+\mu)^2} - \dfrac{y}{\mu^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \theta^2} = \psi_1(y+\theta) - \psi_1(\theta) + \dfrac{\mu}{\theta(\theta+\mu)} + \dfrac{y-\mu}{(\theta+\mu)^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \theta} = \dfrac{y-\mu}{(\theta+\mu)^2}}
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_hessian, NegBinDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  negbin_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Negative Binomial Analytical Expected Hessian
#' @name distrib_expected_hessian.NegBinDistrib
#' @description
#' Computes the analytical expected Hessian of the Negative Binomial log-probability
#' with respect to the parameters \eqn{\mu} and \eqn{\theta}.
#'
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\theta}{\mu(\theta+\mu)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \theta^2}\right] = \mathbb{E}[\psi_1(Y+\theta)] - \psi_1(\theta) + \dfrac{\mu}{\theta(\theta+\mu)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \theta}\right] = 0}
#'
#' The term \eqn{\mathbb{E}[\psi_1(Y+\theta)]} has no closed form and is evaluated by
#' summing over the support up to a far-tail quantile (with a tail-mass correction).
#'
#' @param distrib A \code{NegBinDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_expected_hessian, NegBinDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  negbin_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Negative Binomial Analytical Third-Order Derivatives
#' @name distrib_deriv3.NegBinDistrib
#' @description
#' Closed-form third-order derivatives of the Negative Binomial log-mass. The
#' expected pure-\eqn{\theta} derivative involves \eqn{\mathbb{E}[\psi_2(Y+\theta)]},
#' evaluated by summation over the support (with a far-tail correction).
#' @param distrib A \code{NegBinDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_deriv3, NegBinDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) negbin_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  else negbin_deriv3_cpp(y, theta[[1]], theta[[2]])
}

#' @title Negative Binomial Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.NegBinDistrib
#' @description
#' Closed-form fourth-order derivatives of the Negative Binomial log-mass. The
#' expected pure-\eqn{\theta} derivative involves \eqn{\mathbb{E}[\psi_3(Y+\theta)]},
#' evaluated by summation over the support (with a far-tail correction).
#' @param distrib A \code{NegBinDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{theta}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{negbin_distrib}}
S7::method(distrib_deriv4, NegBinDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) negbin_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  else negbin_deriv4_cpp(y, theta[[1]], theta[[2]])
}

# --- CONSTRUCTOR WRAPPER ---

#' Negative Binomial Distribution Object (NB2)
#'
#' @description
#' Creates a distribution object for the Negative Binomial distribution (NB2)
#' parameterized by mean (\eqn{\mu}) and dispersion (\eqn{\theta}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}}.
#' @param link_theta A link function object for the dispersion parameter \eqn{\theta}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}}.
#'
#' @details
#' The Negative Binomial distribution (NB2) is given a mean/dispersion
#' parameterization, with mean \eqn{\mu} and dispersion \eqn{\theta} (smaller
#' \eqn{\theta} means more overdispersion). Write \eqn{s = \theta + \mu}.
#'
#' \strong{Probability mass function:}
#' \deqn{P(Y=y; \mu, \theta) = \dfrac{\Gamma(y+\theta)}{y!\,\Gamma(\theta)} \left(\dfrac{\theta}{s}\right)^\theta \left(\dfrac{\mu}{s}\right)^y}
#'
#' \strong{Cumulative distribution function:}
#' \deqn{F(q; \mu, \theta) = \sum_{k=0}^{\lfloor q \rfloor} P(Y=k; \mu, \theta)}
#'
#' \strong{Quantile function:} the generalized inverse
#' \eqn{Q(p) = \min\{y \in \mathbb{N}_0 : F(y) \ge p\}}.
#'
#' \strong{Score} (\eqn{\psi} the digamma function):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\theta}{s}\left(\dfrac{y}{\mu} - 1\right), \qquad
#'       \dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) - \psi(\theta) + \log\left(\dfrac{\theta}{s}\right) + \dfrac{\mu - y}{s}}
#'
#' \strong{Observed Hessian} (\eqn{\psi_1} the trigamma function):
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{y+\theta}{s^2} - \dfrac{y}{\mu^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \theta} = \dfrac{y-\mu}{s^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \theta^2} = \psi_1(y+\theta) - \psi_1(\theta) + \dfrac{\mu}{\theta s} + \dfrac{y-\mu}{s^2}}
#'
#' \strong{Expected Hessian:}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\theta}{\mu s}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \theta^2}\right] = \mathbb{E}[\psi_1(Y+\theta)] - \psi_1(\theta) + \dfrac{\mu}{\theta s}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \theta}\right] = 0}
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{\mu + \mu^2/\theta},
#' skewness \eqn{(\theta + 2\mu)/\sqrt{\mu\theta(\theta+\mu)}}, excess kurtosis
#' \eqn{6/\theta + \theta/(\mu(\theta+\mu))}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\theta \in (0, +\infty)}
#' }
#'
#' Analytical third- and fourth-order derivatives (\code{\link{distrib_deriv3}},
#' \code{\link{distrib_deriv4}}) are also available.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.NegBinDistrib}} for the probability mass function.
#'   \item \code{\link{distrib_cdf.NegBinDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.NegBinDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.NegBinDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.NegBinDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.NegBinDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.NegBinDistrib}} for the analytical expected Hessian.
#' }
#'
#' @return An S7 object of class \code{NegBinDistrib} (inheriting from \code{discrete_distrib}) representing the Negative Binomial distribution.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dnbinom pnbinom qnbinom rnbinom
#' @export
negbin_distrib <- function(link_mu = log_link(), link_theta = log_link()) {

  NegBinDistrib(
    distrib_name = "negative binomial",
    dimension = "univariate",
    bounds = c(0, Inf),

    params = c("mu", "theta"),
    params_interpretation = c(mu = "mean", theta = "dispersion"),
    n_params = 2,

    params_bounds = list(
      mu = c(0, Inf),
      theta = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      theta = link_theta
    )
  )

}
