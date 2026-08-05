#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Beta-Binomial Distribution
#' @name BetaBinom1Distrib
#'
#' @description A subclass of \code{discrete_distrib} representing the
#'   beta-binomial distribution, written in its mean proportion and a
#'   dispersion parameter.
#' @inheritParams distrib
#' @param size The number of trials \eqn{n}, a constant of the distribution
#'   rather than a parameter, as for \code{\link{BinomialDistrib}}.
#' @return An object of class \code{BetaBinom1Distrib}.
#' @seealso \code{\link{betabinom1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.BetaBinom1Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.BetaBinom1Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.BetaBinom1Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.BetaBinom1Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.BetaBinom1Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.BetaBinom1Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.BetaBinom1Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{discrete_distrib}}; third and
#' fourth derivatives come from the numerical fallback, which for a family on
#' a finite support is differencing an exact mass function.
BetaBinom1Distrib <- S7::new_class("BetaBinom1Distrib",
  parent = discrete_distrib,
  properties = list(size = S7::class_numeric)
)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta-Binomial Probability Mass Function
#' @name distrib_pdf.BetaBinom1Distrib
#' @description
#' \deqn{P(Y = y) = \binom{n}{y}
#'       \dfrac{B(y + \alpha,\; n - y + \beta)}{B(\alpha, \beta)}, \qquad
#'       \alpha = \dfrac{\mu}{\sigma}, \quad \beta = \dfrac{1-\mu}{\sigma}}
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_pdf, BetaBinom1Distrib) <- function(distrib, y, theta, log = FALSE) {
  out <- betabinom_logpmf_cpp(y, theta[[1]], theta[[2]], distrib@size)
  if (log) out else exp(out)
}

#' @title Beta-Binomial Cumulative Distribution Function
#' @name distrib_cdf.BetaBinom1Distrib
#' @description The cumulative sum of the mass function over the finite
#'   support, which is exact rather than approximated.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_cdf, BetaBinom1Distrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE) {
  n <- distrib@size
  supp <- 0:n
  mass <- exp(betabinom_logpmf_cpp(supp, theta[[1]], theta[[2]], n))
  cum <- cumsum(mass)
  k <- floor(q)
  p <- ifelse(k < 0, 0, ifelse(k >= n, 1, cum[pmin(pmax(k, 0), n) + 1L]))
  p <- pmin(pmax(p, 0), 1)
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

#' @title Beta-Binomial Quantile Function
#' @name distrib_quantile.BetaBinom1Distrib
#' @description The generalized inverse
#'   \eqn{Q(p) = \min\{y : F(y) \ge p\}}, obtained from the exact cumulative
#'   sum.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), \code{p} is \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_quantile, BetaBinom1Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  n <- distrib@size
  cum <- cumsum(exp(betabinom_logpmf_cpp(0:n, theta[[1]], theta[[2]], n)))
  vapply(p, function(pp) {
    if (is.na(pp)) return(NA_real_)
    if (pp <= 0) return(0)
    if (pp >= 1) return(n)
    (0:n)[which(cum >= pp - 1e-12)[1L]]
  }, numeric(1))
}

#' @title Beta-Binomial Random Generation
#' @name distrib_rng.BetaBinom1Distrib
#' @description Draws a probability from the Beta and then a Binomial with it,
#'   which is the hierarchy the family is defined by.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_rng, BetaBinom1Distrib) <- function(distrib, n, theta) {
  a <- theta[[1]] / theta[[2]]
  b <- (1 - theta[[1]]) / theta[[2]]
  stats::rbinom(n, size = distrib@size, prob = stats::rbeta(n, a, b))
}

#' @title Beta-Binomial Analytical Gradient
#' @name distrib_gradient.BetaBinom1Distrib
#' @description
#' The chain rule from the shapes, where every derivative is a difference of
#' digammas:
#' \deqn{\dfrac{\partial \ell}{\partial \alpha}
#'       = \psi(y+\alpha) - \psi(\alpha) - \psi(n+S) + \psi(S), \quad S =
#'       \alpha + \beta}
#' and likewise in \eqn{\beta}, with \eqn{\alpha = \mu/\sigma} and
#' \eqn{\beta = (1-\mu)/\sigma}.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu} and \code{sigma} components.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_gradient, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  betabinom_gradient_cpp(y, theta[[1]], theta[[2]], distrib@size)
}

#' @title Beta-Binomial Analytical Observed Hessian
#' @name distrib_hessian.BetaBinom1Distrib
#' @description
#' The two-variable chain rule from the shapes, whose second derivatives are
#' differences of trigammas; the mixed shape component carries only the
#' \eqn{S = \alpha + \beta} part, the two shapes entering the mass function
#' separately otherwise.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_hessian, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  betabinom_hessian_cpp(y, theta[[1]], theta[[2]], distrib@size)
}

#' @title Beta-Binomial Analytical Expected Hessian
#' @name distrib_expected_hessian.BetaBinom1Distrib
#' @description
#' The observed Hessian averaged against the mass over \eqn{\{0, \dots, n\}}.
#' The support being finite, the expectation is an exact sum rather than a
#' quadrature, and \code{approx} is therefore ignored.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is an exact sum.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso \code{\link{betabinom1_distrib}}
S7::method(distrib_expected_hessian, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  betabinom_expected_hessian_cpp(y, theta[[1]], theta[[2]], distrib@size)
}

# --- CONSTRUCTOR WRAPPER ---

#' Beta-Binomial Distribution Object
#'
#' @description
#' Creates a distribution object for the beta-binomial distribution,
#' parametrised by the mean proportion \eqn{\mu} and a dispersion parameter
#' \eqn{\sigma}.
#'
#' @param size The number of trials \eqn{n}. A constant of the distribution
#'   rather than a parameter, as for \code{\link{binomial_distrib}}.
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   \code{\link[linkfunctions7]{logit_link}}, the mean being a proportion.
#' @param link_sigma A link function object for \eqn{\sigma}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The family is the binomial with its success probability drawn from a Beta,
#' which is what makes it the natural model for a proportion whose trials are
#' not independent. Writing \eqn{\alpha = \mu/\sigma} and
#' \eqn{\beta = (1-\mu)/\sigma}, the mean is \eqn{n\mu} and
#' \deqn{\operatorname{Var}(Y) = n\mu(1-\mu)
#'       \left(1 + (n-1)\dfrac{\sigma}{1+\sigma}\right),}
#' so the family is overdispersed relative to the binomial at every
#' \eqn{\sigma > 0} and approaches it as \eqn{\sigma \to 0}. The
#' intraclass correlation is \eqn{\sigma/(1+\sigma)}.
#'
#' \strong{Probability mass function:}
#' \deqn{P(Y = y) = \binom{n}{y}
#'       \dfrac{B(y+\alpha,\; n-y+\beta)}{B(\alpha, \beta)}}
#'
#' \strong{Score and information.} The parameters enter only through the two
#' shapes, where every derivative is a difference of polygammas, so the score
#' is that difference carried through the chain rule of
#' \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)}. The expected
#' information is an \strong{exact sum} over the finite support rather than a
#' quadrature, which is what a bounded count buys.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' The family is \strong{not} reachable from anything already in the package:
#' it is neither a binomial with a parameter held fixed nor a wrapper over
#' one, the mixing being over the success probability rather than over the
#' outcome.
#'
#' Third and fourth derivatives come from the numerical fallback. On a finite
#' support that fallback differences an exact mass function, so its accuracy
#' is the usual \code{1e-8} rather than the poorer figure a quadrature would
#' give.
#'
#' @return An S7 object of class \code{BetaBinom1Distrib}.
#'
#' @seealso \code{\link{binomial_distrib}}, \code{\link{beta1_distrib}},
#'   \code{\link{negbin2_distrib}}
#'
#' @importFrom linkfunctions7 logit_link log_link
#' @importFrom stats rbeta rbinom
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' d@params
#'
#' theta <- list(mu = 0.3, sigma = 0.5)
#' distrib_pdf(d, 0:10, theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' # overdispersed relative to the binomial with the same mean
#' variance(binomial_distrib(size = 10), list(mu = 0.3))
#'
#' @export
betabinom1_distrib <- function(size, link_mu = logit_link(),
                              link_sigma = log_link()) {
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
      size < 1 || size != round(size)) {
    stop("'size' must be a single positive integer.", call. = FALSE)
  }
  BetaBinom1Distrib(
    size = as.numeric(size),
    distrib_name = sprintf("beta-binomial [size=%g]", size),
    dimension = "univariate", bounds = c(0, size),
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean proportion", sigma = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, 1), sigma = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma)
  )
}
