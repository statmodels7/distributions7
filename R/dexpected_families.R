#' @include dexpected_hessian.R bernoulli_distrib.R binomial_distrib.R
#' @include exponential_distrib.R geometric_distrib.R chisq_distrib.R
#' @include cauchy_distrib.R logistic_distrib.R gumbel_distrib.R
#' @include gaussian2_distrib.R gaussian3_distrib.R lognormal1_distrib.R
#' @include invgauss1_distrib.R invgauss2_distrib.R gamma2_distrib.R
#' @include beta2_distrib.R weibull1_distrib.R gpd_distrib.R
NULL

#' @title Derivatives of the Expected Information, Elementary Families
#' @name distrib_dexpected_hessian.elementary
#' @aliases distrib_dexpected_hessian.BernoulliDistrib
#'   distrib_d2expected_hessian.BernoulliDistrib
#'   distrib_dexpected_hessian.BinomialDistrib
#'   distrib_d2expected_hessian.BinomialDistrib
#'   distrib_dexpected_hessian.ExponentialDistrib
#'   distrib_d2expected_hessian.ExponentialDistrib
#'   distrib_dexpected_hessian.GeometricDistrib
#'   distrib_d2expected_hessian.GeometricDistrib
#'   distrib_dexpected_hessian.ChisqDistrib
#'   distrib_d2expected_hessian.ChisqDistrib
#'   distrib_dexpected_hessian.CauchyDistrib
#'   distrib_d2expected_hessian.CauchyDistrib
#'   distrib_dexpected_hessian.LogisticDistrib
#'   distrib_d2expected_hessian.LogisticDistrib
#'   distrib_dexpected_hessian.GumbelDistrib
#'   distrib_d2expected_hessian.GumbelDistrib
#'   distrib_dexpected_hessian.Gaussian2Distrib
#'   distrib_d2expected_hessian.Gaussian2Distrib
#'   distrib_dexpected_hessian.Gaussian3Distrib
#'   distrib_d2expected_hessian.Gaussian3Distrib
#'   distrib_dexpected_hessian.Lognormal1Distrib
#'   distrib_d2expected_hessian.Lognormal1Distrib
#'   distrib_dexpected_hessian.InvGauss1Distrib
#'   distrib_d2expected_hessian.InvGauss1Distrib
#'   distrib_dexpected_hessian.InvGauss2Distrib
#'   distrib_d2expected_hessian.InvGauss2Distrib
#'   distrib_dexpected_hessian.Gamma2Distrib
#'   distrib_d2expected_hessian.Gamma2Distrib
#'   distrib_dexpected_hessian.Beta2Distrib
#'   distrib_d2expected_hessian.Beta2Distrib
#'   distrib_dexpected_hessian.Weibull1Distrib
#'   distrib_d2expected_hessian.Weibull1Distrib
#'   distrib_dexpected_hessian.GPDDistrib
#'   distrib_d2expected_hessian.GPDDistrib
#' @description
#' The first and second derivatives of the expected information in the
#' parameters, from compiled kernels, for seventeen families whose expected
#' information is an elementary function of the parameters. Each component is
#' an ordinary derivative of the family's written-out
#' \eqn{\mathbb{E}[\ell_{ab}]}; on the link scale the result is carried across
#' by [dexpected_link()].
#'
#' @details
#' The expected informations differentiated, on the parameter scale:
#' \itemize{
#'   \item Bernoulli and binomial of size \eqn{n}:
#'     \eqn{-n/(\mu(1-\mu))}, with \eqn{n = 1} for the Bernoulli;
#'   \item exponential: \eqn{-1/\mu^2}; geometric: \eqn{-1/(\mu(1+\mu))};
#'   \item chi-squared: \eqn{-\psi'(\mu/2)/4};
#'   \item Cauchy, logistic and Gumbel: every component a constant over
#'     \eqn{\sigma^2}, so only \eqn{\sigma} moves them;
#'   \item gaussian by its variance, and lognormal in \eqn{(\mu, \sigma^2)}:
#'     \eqn{-1/\sigma^2} and \eqn{-1/(2\sigma^4)}; gaussian by its precision:
#'     \eqn{-\tau} and \eqn{-1/(2\tau^2)};
#'   \item inverse gaussian: \eqn{-1/(\phi\mu^3)} and \eqn{-1/(2\phi^2)} by the
#'     dispersion, \eqn{-\lambda/\mu^3} and \eqn{-1/(2\lambda^2)} by the shape;
#'   \item gamma by its mean and variance: with \eqn{a = \mu^2/\sigma^2} and
#'     \eqn{g(a) = a\,\psi'(a) - 1},
#'     \eqn{\mathbb{E}[\ell_{\mu\mu}] = -(1+4g)/\sigma^2},
#'     \eqn{\mathbb{E}[\ell_{\mu\sigma^2}] = 2\mu g/\sigma^4} and
#'     \eqn{\mathbb{E}[\ell_{\sigma^2\sigma^2}] = -a g/\sigma^4}; \eqn{g} and its
#'     derivatives are formed from the remainders of the polygamma functions'
#'     asymptotic series, so nothing cancels at a large shape, where
#'     \eqn{g \sim 1/(2a)};
#'   \item beta by its two shapes: \eqn{\psi'(\alpha+\beta) - \psi'(\alpha)},
#'     \eqn{\psi'(\alpha+\beta)} and \eqn{\psi'(\alpha+\beta) - \psi'(\beta)};
#'   \item Weibull by its scale and shape: \eqn{-\sigma^2/\mu^2},
#'     \eqn{(1-\gamma)/\mu} and \eqn{-C/\sigma^2}, with \eqn{\gamma} Euler's
#'     constant and \eqn{C = (1-\gamma)^2 + \pi^2/6};
#'   \item generalized Pareto by its scale and shape: with \eqn{d = 1 + 2\xi},
#'     \eqn{-1/(d\sigma^2)}, \eqn{-1/(d\sigma(1+\xi))} and
#'     \eqn{-2/(d(1+\xi))}, and `NA` for \eqn{\xi \le -1/2}, where the
#'     information does not exist.
#' }
#'
#' @param distrib A distribution object of one of the classes above.
#' @param y A numeric vector of observations, read for its length.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads A single positive integer, the kernel's thread count.
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#' @seealso [distrib_dexpected_hessian()], [distrib_d2expected_hessian()]
#' @keywords internal
NULL

#' Register a Family's Derivatives of the Expected Information
#'
#' @description
#' Registers methods of [distrib_dexpected_hessian()] and
#' [distrib_d2expected_hessian()] for one class, both reading one kernel on
#' the parameter scale through [dexpected_analytic()].
#'
#' @param cls An S7 class.
#' @param kernel A function of `(distrib, y, theta, order, threads)`
#'   returning the parameter-scale components keyed as [dexpected_names()]
#'   at order 1 and [d2expected_names()] at order 2.
#'
#' @return `NULL`, invisibly; called for the registration.
#'
#' @keywords internal
register_dexpected <- function(cls, kernel) {
  S7::method(distrib_dexpected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
      threads = 1L) {
    dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                       function(k) kernel(distrib, y, theta, k, threads))
  }
  S7::method(distrib_d2expected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
      threads = 1L) {
    dexpected_analytic(distrib, y, theta, match.arg(scale), 2L, threads,
                       function(k) kernel(distrib, y, theta, k, threads))
  }
  invisible(NULL)
}

register_dexpected(BernoulliDistrib, function(d, y, th, k, t)
  bernoulli_dexpected_cpp(y, th[[1]], k, t))
register_dexpected(BinomialDistrib, function(d, y, th, k, t)
  binomial_dexpected_cpp(y, th[[1]], d@size, k, t))
register_dexpected(ExponentialDistrib, function(d, y, th, k, t)
  exponential_dexpected_cpp(y, th[[1]], k, t))
register_dexpected(GeometricDistrib, function(d, y, th, k, t)
  geometric_dexpected_cpp(y, th[[1]], k, t))
register_dexpected(ChisqDistrib, function(d, y, th, k, t)
  chisq_dexpected_cpp(y, th[[1]], k, t))
register_dexpected(CauchyDistrib, function(d, y, th, k, t)
  cauchy_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(LogisticDistrib, function(d, y, th, k, t)
  logistic_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(GumbelDistrib, function(d, y, th, k, t)
  gumbel_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(Gaussian2Distrib, function(d, y, th, k, t)
  gaussian2_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(Gaussian3Distrib, function(d, y, th, k, t)
  gaussian3_dexpected_cpp(y, th[[1]], th[[2]], k, t))
# the lognormal's (mu, sigma2) information is the gaussian's by its variance
register_dexpected(Lognormal1Distrib, function(d, y, th, k, t)
  gaussian2_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(InvGauss1Distrib, function(d, y, th, k, t)
  invgauss1_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(InvGauss2Distrib, function(d, y, th, k, t)
  invgauss2_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(Gamma2Distrib, function(d, y, th, k, t)
  gamma2_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(Beta2Distrib, function(d, y, th, k, t)
  beta2_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(Weibull1Distrib, function(d, y, th, k, t)
  weibull1_dexpected_cpp(y, th[[1]], th[[2]], k, t))
register_dexpected(GPDDistrib, function(d, y, th, k, t)
  gpd_dexpected_cpp(y, th[[1]], th[[2]], k, t))
