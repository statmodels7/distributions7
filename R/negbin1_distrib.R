#' @include distrib.R generics.R
NULL

#' @title S7 Class for the NB1 Negative Binomial
#' @name NegBin1Distrib
#'
#' @description A subclass of `discrete_distrib` representing the
#'   negative binomial whose variance is **linear** in the mean.
#' @inheritParams distrib
#' @return An object of class `NegBin1Distrib`.
#' @seealso [negbin1_distrib()], [negbin2_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.NegBin1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.NegBin1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.NegBin1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.NegBin1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.NegBin1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.NegBin1Distrib],
#'   [`distrib_rng()`][distrib_rng.NegBin1Distrib]
#'
#' Everything else is inherited from [discrete_distrib()].
NegBin1Distrib <- S7::new_class("NegBin1Distrib", parent = discrete_distrib)

#' The Size Behind an NB1 Mean
#'
#' @description
#' The number of successes \eqn{r = \mu/\theta} that the base R negative
#' binomial functions take.
#'
#' @details
#' Requiring the variance to be \eqn{\mu(1+\theta)} fixes the success
#' probability at \eqn{1/(1+\theta)}, and the mean then determines the size.
#' It is this that puts \eqn{\mu} inside the gamma functions and makes the
#' family different from the quadratic form rather than a reparametrization of
#' it.
#'
#' @param mu The mean, a positive numeric vector.
#' @param theta The dispersion, a positive numeric vector.
#'
#' @return A numeric vector.
#'
#' @seealso [negbin1_distrib()], [nb1_prob()]
#'
#' @keywords internal
nb1_size <- function(mu, theta) mu / theta

#' The Success Probability Behind an NB1 Dispersion
#'
#' @description
#' \eqn{1/(1+\theta)}, the value that makes the variance-to-mean ratio
#' \eqn{1+\theta} at every mean.
#'
#' @param theta The dispersion, a positive numeric vector.
#'
#' @return A numeric vector in \eqn{(0, 1)}.
#'
#' @seealso [negbin1_distrib()], [nb1_size()]
#'
#' @keywords internal
nb1_prob <- function(theta) 1 / (1 + theta)

# --- S7 METHODS IMPLEMENTATION ---

#' @title NB1 Probability Mass Function
#' @name distrib_pdf.NegBin1Distrib
#' @description
#' The negative binomial mass at size \eqn{\mu/\theta} and success probability
#' \eqn{1/(1+\theta)}, which is what makes the variance \eqn{\mu(1+\theta)}.
#' @param distrib A `NegBin1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `theta`.
#' @param log Logical; if `TRUE`, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso [negbin1_distrib()]
S7::method(distrib_pdf, NegBin1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  out <- negbin1_logpmf_cpp(y, theta[[1]], theta[[2]])
  if (log) out else exp(out)
}

#' @title NB1 Cumulative Distribution Function
#' @name distrib_cdf.NegBin1Distrib
#' @description [stats::pnbinom()] at size \eqn{\mu/\theta} and
#'   probability \eqn{1/(1+\theta)}.
#' @param distrib A `NegBin1Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing `mu` and `theta`.
#' @param lower.tail Logical; if `TRUE` (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [negbin1_distrib()]
S7::method(distrib_cdf, NegBin1Distrib) <- function(distrib, q, theta,
                                                     lower.tail = TRUE,
                                                     log.p = FALSE) {
  stats::pnbinom(q, size = nb1_size(theta[[1]], theta[[2]]),
                 prob = nb1_prob(theta[[2]]),
                 lower.tail = lower.tail, log.p = log.p)
}

#' @title NB1 Quantile Function
#' @name distrib_quantile.NegBin1Distrib
#' @description [stats::qnbinom()] at the same size and probability.
#' @param distrib A `NegBin1Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing `mu` and `theta`.
#' @param lower.tail Logical; if `TRUE` (default), `p` is \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, `p` is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso [negbin1_distrib()]
S7::method(distrib_quantile, NegBin1Distrib) <- function(distrib, p, theta,
                                                          lower.tail = TRUE,
                                                          log.p = FALSE) {
  stats::qnbinom(p, size = nb1_size(theta[[1]], theta[[2]]),
                 prob = nb1_prob(theta[[2]]),
                 lower.tail = lower.tail, log.p = log.p)
}

#' @title NB1 Random Generation
#' @name distrib_rng.NegBin1Distrib
#' @description [stats::rnbinom()] at the same size and probability.
#' @param distrib A `NegBin1Distrib` object.
#' @param n The number of draws.
#' @param theta A list containing `mu` and `theta`.
#' @return A numeric vector of length `n`.
#' @seealso [negbin1_distrib()]
S7::method(distrib_rng, NegBin1Distrib) <- function(distrib, n, theta) {
  stats::rnbinom(n, size = nb1_size(theta[[1]], theta[[2]]),
                 prob = nb1_prob(theta[[2]]))
}

#' @title NB1 Analytical Gradient
#' @name distrib_gradient.NegBin1Distrib
#' @description
#' The chain rule through \eqn{r = \mu/\theta}. With
#' \eqn{P = \psi(y+r) - \psi(r) - \log(1+\theta)},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{P}{\theta}}
#' and the derivative in \eqn{\theta} adds the terms in which \eqn{\theta}
#' appears outside \eqn{r}.
#' @param distrib A `NegBin1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `theta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list with the `mu` and `theta` components.
#' @seealso [negbin1_distrib()]
S7::method(distrib_gradient, NegBin1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  negbin1_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title NB1 Analytical Observed Hessian
#' @name distrib_hessian.NegBin1Distrib
#' @description
#' The two-variable chain rule through \eqn{r = \mu/\theta}, whose second
#' derivative in \eqn{r} is \eqn{\psi'(y+r) - \psi'(r)} and whose mixed term is
#' \eqn{-1/(1+\theta)}.
#' @param distrib A `NegBin1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `theta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of second-derivative components.
#' @seealso [negbin1_distrib()]
S7::method(distrib_hessian, NegBin1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  negbin1_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title NB1 Analytical Expected Hessian
#' @name distrib_expected_hessian.NegBin1Distrib
#' @description
#' Every term carrying \eqn{P} drops out, its expectation vanishing by the
#' first Bartlett identity, and what remains needs only
#' \eqn{\mathbb{E}[\psi'(Y+r)]}. That has no closed form and is summed against
#' the exact mass to a far-tail quantile, as the NB2 kernel does; `approx`
#' is therefore ignored.
#' @param distrib A `NegBin1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `theta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @param threads How many threads the kernel may use; below the measured
#'   internal threshold it stays sequential whatever the count says.
#' @return A named list of expected second-derivative components.
#' @seealso [negbin1_distrib()]
S7::method(distrib_expected_hessian, NegBin1Distrib) <- function(distrib, y, theta,
                                                                  scale = c("parameter", "link"),
                                                                  approx = c("bartlett", "integrate", "mc", "opg"),
                                                                  nsim = 10000, ...,
                                       threads = 1L) {
  negbin1_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' NB1 Negative Binomial Distribution Object
#'
#' @description
#' Creates a distribution object for the negative binomial whose variance is
#' **linear** in the mean, \eqn{\operatorname{Var}(Y) = \mu(1+\theta)}.
#'
#' @section The distribution:
#' \deqn{P(Y=y) = \frac{\Gamma(y + \mu/\theta)}{\Gamma(\mu/\theta)\,y!}\left(\frac{1}{1+\theta}\right)^{\mu/\theta}\left(\frac{\theta}{1+\theta}\right)^{y}}
#' on \eqn{y \in \{0, 1, \dots\}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \mu(1+\theta)}
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   [linkfunctions7::log_link()].
#' @param link_theta A link function object for \eqn{\theta}. Defaults to
#'   [linkfunctions7::log_link()].
#'
#' @details
#' Two negative binomials are in common use and they are **different
#' families**, not two parametrizations of one. Here the variance is
#' \eqn{\mu(1+\theta)}, growing in proportion to the mean, so the dispersion
#' relative to a Poisson is the same at every mean;
#' [negbin2_distrib()] has \eqn{\mu + \mu^2/\theta}, growing
#' quadratically. Fitting one is not fitting the other, and a likelihood ratio
#' between them is not a test of nested models.
#'
#' The difference is visible in where the mean sits. The size is
#' \eqn{r = \mu/\theta} and the success probability \eqn{1/(1+\theta)}, so
#' \eqn{\mu} appears **inside** the gamma functions, while in the
#' quadratic form it stays outside them and \eqn{\theta} is the size.
#'
#' **Probability mass function:** the negative binomial mass at size
#' \eqn{\mu/\theta} and probability \eqn{1/(1+\theta)}.
#'
#' **Score.** Writing
#' \eqn{P = \psi(y+r) - \psi(r) - \log(1+\theta)},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{P}{\theta}, \qquad
#'       \dfrac{\partial \ell}{\partial \theta} = -\dfrac{\mu}{\theta^2}P
#'       - \dfrac{r}{1+\theta} + \dfrac{y}{\theta} - \dfrac{y}{1+\theta}}
#' and the Hessian is the same chain rule at second order.
#'
#' **Expected information.** Every term carrying \eqn{P} drops out, its
#' expectation vanishing by the first Bartlett identity, and only
#' \eqn{\mathbb{E}[\psi'(Y+r)]} remains. That has no closed form and is summed
#' against the exact mass to a far-tail quantile.
#'
#' **Parameter domains:**
#'
#' - \eqn{\mu \in (0, +\infty)}
#' - \eqn{\theta \in (0, +\infty)}
#'
#' As \eqn{\theta \to 0} the family approaches the Poisson, as the quadratic
#' form does when \eqn{\theta \to \infty}.
#'
#' @return An S7 object of class `NegBin1Distrib`.
#'
#' @seealso [negbin2_distrib()] for the quadratic variance,
#'   [poisson_distrib()], [geometric_distrib()]
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats pnbinom qnbinom rnbinom
#' @examples
#' d <- negbin1_distrib()
#' d@params
#'
#' theta <- list(mu = 4, theta = 4)
#' distrib_pdf(d, 0:6, theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' # the two negative binomials are different families: at the same (mu, theta)
#' # this one has variance mu(1+theta) = 20 and the other mu + mu^2/theta = 8
#' variance(negbin2_distrib(), list(mu = 4, theta = 4))
#'
#' @references
#' Cameron, A. C. and Trivedi, P. K. (1986). Econometric models based
#' on count data: comparisons and applications of some estimators and
#' tests. *Journal of Applied Econometrics* 1, 29-53.
#'
#' @export
negbin1_distrib <- function(link_mu = log_link(), link_theta = log_link()) {
  NegBin1Distrib(
    distrib_name = "negbin1", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu", "theta"),
    params_interpretation = c(mu = "mean", theta = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), theta = c(0, Inf)),
    link_params = list(mu = link_mu, theta = link_theta)
  )
}
