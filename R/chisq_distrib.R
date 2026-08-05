#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Chi-Squared Distribution
#' @name ChisqDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the
#'   chi-squared distribution, whose degrees of freedom are its mean.
#' @inheritParams distrib
#' @return An object of class \code{ChisqDistrib}.
#' @seealso \code{\link{chisq_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.ChisqDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.ChisqDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.ChisqDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.ChisqDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.ChisqDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_grad_y.ChisqDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_hessian.ChisqDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_hess_y.ChisqDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_pdf.ChisqDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.ChisqDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.ChisqDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
ChisqDistrib <- S7::new_class("ChisqDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Chi-Squared Density
#' @name distrib_pdf.ChisqDistrib
#' @description
#' \deqn{f(y; \mu) = \dfrac{y^{\mu/2 - 1} e^{-y/2}}{2^{\mu/2}\Gamma(\mu/2)},
#'       \qquad y > 0}
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_pdf, ChisqDistrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dchisq(y, df = theta[[1]], log = log)
}

#' @title Chi-Squared Distribution Function
#' @name distrib_cdf.ChisqDistrib
#' @description The regularised lower incomplete gamma function
#'   \eqn{P(\mu/2, q/2)}.
#' @param distrib A \code{ChisqDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_cdf, ChisqDistrib) <- function(distrib, q, theta,
                                                  lower.tail = TRUE,
                                                  log.p = FALSE) {
  stats::pchisq(q, df = theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Chi-Squared Quantile Function
#' @name distrib_quantile.ChisqDistrib
#' @description The inverse of \code{\link[stats]{pchisq}}.
#' @param distrib A \code{ChisqDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameter \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), \code{p} is \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_quantile, ChisqDistrib) <- function(distrib, p, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE) {
  stats::qchisq(p, df = theta[[1]], lower.tail = lower.tail, log.p = log.p)
}

#' @title Chi-Squared Random Generation
#' @name distrib_rng.ChisqDistrib
#' @description Draws through \code{\link[stats]{rchisq}}.
#' @param distrib A \code{ChisqDistrib} object.
#' @param n The number of draws.
#' @param theta A list containing the parameter \code{mu}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_rng, ChisqDistrib) <- function(distrib, n, theta) {
  stats::rchisq(n, df = theta[[1]])
}

#' @title Chi-Squared Analytical Gradient
#' @name distrib_gradient.ChisqDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \mu} =
#'       \dfrac{\log y - \log 2 - \psi(\mu/2)}{2}}
#' the only order that involves the response.
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu} component.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_gradient, ChisqDistrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"), ...) {
  chisq_gradient_cpp(y, theta[[1]])
}

#' @title Chi-Squared Analytical Observed Hessian
#' @name distrib_hessian.ChisqDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\psi'(\mu/2)}{4}}
#' which does not involve the response, so it coincides with the expected
#' information.
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu} component.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_hessian, ChisqDistrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ...) {
  chisq_hessian_cpp(y, theta[[1]])
}

#' @title Chi-Squared Analytical Expected Hessian
#' @name distrib_expected_hessian.ChisqDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right]
#'       = -\dfrac{\psi'(\mu/2)}{4}}
#' identical to the observed Hessian on the parameter scale, which the second
#' derivative not depending on the response makes exact rather than
#' approximate. The two differ on the link scale, where the chain rule adds a
#' term proportional to the score; see \code{\link{chisq_distrib}}.
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu} component.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_expected_hessian, ChisqDistrib) <- function(distrib, y, theta,
                                                               scale = c("parameter", "link"),
                                                               approx = c("bartlett", "integrate", "mc", "opg"),
                                                               nsim = 10000, ...) {
  chisq_hessian_cpp(y, theta[[1]])
}

#' @title Chi-Squared Analytical Third-Order Derivative
#' @name distrib_deriv3.ChisqDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu)} = -\dfrac{\psi''(\mu/2)}{8}}
#' observed and expected alike.
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; the two coincide, so it changes nothing.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu_mu} component.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_deriv3, ChisqDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ...) {
  chisq_deriv3_cpp(y, theta[[1]])
}

#' @title Chi-Squared Analytical Fourth-Order Derivative
#' @name distrib_deriv4.ChisqDistrib
#' @description
#' \deqn{\ell^{(\mu\mu\mu\mu)} = -\dfrac{\psi'''(\mu/2)}{16}}
#' observed and expected alike.
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param expected Logical; the two coincide, so it changes nothing.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list with the \code{mu_mu_mu_mu} component.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_deriv4, ChisqDistrib) <- function(distrib, y, theta, expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ...) {
  chisq_deriv4_cpp(y, theta[[1]])
}

#' @title Chi-Squared Response Gradient
#' @name distrib_grad_y.ChisqDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\mu/2 - 1}{y} - \dfrac{1}{2}}
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_grad_y, ChisqDistrib) <- function(distrib, y, theta, ...) {
  (theta[[1]] / 2 - 1) / y - 0.5
}

#' @title Chi-Squared Response Hessian
#' @name distrib_hess_y.ChisqDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{\mu/2 - 1}{y^2}}
#' @param distrib A \code{ChisqDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameter \code{mu}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{chisq_distrib}}
S7::method(distrib_hess_y, ChisqDistrib) <- function(distrib, y, theta, ...) {
  -(theta[[1]] / 2 - 1) / (y * y)
}

# --- CONSTRUCTOR WRAPPER ---

#' Chi-Squared Distribution Object
#'
#' @description
#' Creates a distribution object for the chi-squared distribution, parametrised
#' by its mean \eqn{\mu}, which is the degrees of freedom.
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The degrees of freedom are treated as a continuous positive parameter, which
#' is what makes the family estimable; the mean is \eqn{\mu} and the variance
#' \eqn{2\mu}.
#'
#' \strong{Density:}
#' \deqn{f(y; \mu) = \dfrac{y^{\mu/2 - 1} e^{-y/2}}{2^{\mu/2}\Gamma(\mu/2)}}
#'
#' \strong{Score and information:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu}
#'         = \dfrac{\log y - \log 2 - \psi(\mu/2)}{2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\psi'(\mu/2)}{4}}
#'
#' The family is a one-parameter exponential family in \eqn{\log y}, so from
#' the second order on the derivatives do not involve the response at all,
#' \deqn{\ell^{(k)} = -\dfrac{\psi^{(k-2)}(\mu/2)}{2^{k}}, \qquad k \ge 2.}
#' On the parameter scale the observed information is therefore exactly the
#' expected information, and the same holds at third and fourth order: there is
#' nothing to average. \eqn{\mathbb{E}[\log y] = \psi(\mu/2) + \log 2} is what
#' makes the score have mean zero.
#'
#' That coincidence does not carry to the scale a fit optimises on. The
#' second-order chain rule of \code{\link{distrib_hessian}} adds a term
#' \eqn{h''(\eta)\,\partial\ell/\partial\mu} to the link-scale Hessian, and the
#' expected version drops it because the score has mean zero, while a sample
#' does not. Fisher scoring and Newton's method therefore take different steps
#' here, and agree at the optimum, where the summed score vanishes.
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{2\mu}, skewness
#' \eqn{2\sqrt{2/\mu}}, excess kurtosis \eqn{12/\mu}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#' }
#'
#' The family is a Gamma with shape \eqn{\mu/2} and scale 2, but it is
#' \strong{not} a Gamma with a fixed parameter: this package writes the Gamma
#' in \eqn{(\mu, \sigma^2)}, and a scale of 2 is the relation
#' \eqn{\sigma^2 = 2\mu} between two parameters rather than a value one of them
#' can be held at.
#'
#' @return An S7 object of class \code{ChisqDistrib}.
#'
#' @seealso \code{\link{gamma2_distrib}}, \code{\link{exponential_distrib}}
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dchisq pchisq qchisq rchisq
#' @examples
#' d <- chisq_distrib()
#' d@params
#'
#' theta <- list(mu = 4)
#' distrib_pdf(d, c(1, 4, 9), theta)
#' distrib_gradient(d, c(1, 4, 9), theta)
#'
#' # the observed and expected information coincide exactly
#' distrib_hessian(d, c(1, 4, 9), theta)
#' distrib_expected_hessian(d, c(1, 4, 9), theta)
#'
#' @export
chisq_distrib <- function(link_mu = log_link()) {
  ChisqDistrib(
    distrib_name = "chisq", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu"), params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = link_mu)
  )
}
