#' @include distrib.R generics.R pig1_distrib.R
NULL

#' @title S7 Class for the Poisson-Inverse Gaussian in Its Orthogonal Parametrization
#' @name Pig2Distrib
#'
#' @description A subclass of \code{discrete_distrib} representing the
#'   Poisson-inverse Gaussian distribution in the parametrization whose two
#'   parameters are orthogonal, gamlss's \code{PIG2}.
#' @inheritParams distrib
#' @return An object of class \code{Pig2Distrib}.
#' @seealso \code{\link{pig2_distrib}}, \code{\link{pig1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_pdf.Pig2Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_gradient.Pig2Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.Pig2Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_deriv3.Pig2Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Pig2Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_rng.Pig2Distrib]{distrib_rng()}}
#'
#' Everything else, the distribution function and the quantile included, is
#' inherited from \code{\link{discrete_distrib}}.
Pig2Distrib <- S7::new_class("Pig2Distrib", parent = discrete_distrib)

#' The Dispersion a Poisson-Inverse Gaussian Alpha Implies
#'
#' @description
#' Converts the orthogonal parametrization's \eqn{\alpha} into the
#' dispersion \eqn{\sigma} of \code{\link{pig1_distrib}}:
#' \eqn{\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2}, the positive
#' root of \eqn{\alpha^2\sigma^2 - 2\mu\sigma - 1 = 0}, with no
#' cancellation anywhere in the domain.
#'
#' @param mu The mean, a positive numeric vector.
#' @param alpha The orthogonal parameter, a positive numeric vector.
#' @return A numeric vector of dispersions.
#' @seealso \code{\link{pig2_distrib}}
#' @keywords internal
pig2_sigma <- function(mu, alpha) (mu + sqrt(mu^2 + alpha^2)) / alpha^2

# --- S7 METHODS IMPLEMENTATION ---

#' @title Orthogonal Poisson-Inverse Gaussian Probability Mass Function
#' @name distrib_pdf.Pig2Distrib
#' @description
#' The same law as \code{\link[=distrib_pdf.Pig1Distrib]{pig1's}} at
#' \eqn{\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2}: the parameter
#' \eqn{\alpha} of this parametrization is exactly the argument
#' \eqn{\sqrt{1 + 2\sigma\mu}/\sigma} the Bessel function is evaluated at.
#' @param distrib A \code{Pig2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{pig2_distrib}}
S7::method(distrib_pdf, Pig2Distrib) <- function(distrib, y, theta, log = FALSE) {
  out <- pig_hd_block(y, theta, c(l = "l"), pig2_hd_cpp)$l
  out[is.nan(out)] <- -Inf
  if (log) out else exp(out)
}

#' @title Orthogonal Poisson-Inverse Gaussian Analytical Gradient
#' @name distrib_gradient.Pig2Distrib
#' @description The exact score in \eqn{(\mu, \alpha)}, from the compiled
#' fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu} and \code{alpha} components.
#' @seealso \code{\link{pig2_distrib}}
S7::method(distrib_gradient, Pig2Distrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ...) {
  pig_hd_block(y, theta, c(mu = "d10", alpha = "d01"), pig2_hd_cpp)
}

#' @title Orthogonal Poisson-Inverse Gaussian Analytical Observed Hessian
#' @name distrib_hessian.Pig2Distrib
#' @description The exact second derivatives in \eqn{(\mu, \alpha)}, from
#' the compiled fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{pig2_distrib}}
S7::method(distrib_hessian, Pig2Distrib) <- function(distrib, y, theta,
                                                     scale = c("parameter", "link"), ...) {
  pig_hd_block(y, theta,
               c(mu_mu = "d20", alpha_alpha = "d02", mu_alpha = "d11"),
               pig2_hd_cpp)
}

#' @title Orthogonal Poisson-Inverse Gaussian Analytical Third Derivatives
#' @name distrib_deriv3.Pig2Distrib
#' @description The exact third derivatives in \eqn{(\mu, \alpha)}, from the
#' compiled fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @param expected Logical; the expected version goes through the generic's
#'   strategies.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx,nsim Passed on when \code{expected} is \code{TRUE}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{pig2_distrib}}
S7::method(distrib_deriv3, Pig2Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig_hd_block(y, theta,
               c(mu_mu_mu = "d30", mu_mu_alpha = "d21",
                 mu_alpha_alpha = "d12", alpha_alpha_alpha = "d03"),
               pig2_hd_cpp)
}

#' @title Orthogonal Poisson-Inverse Gaussian Analytical Fourth Derivatives
#' @name distrib_deriv4.Pig2Distrib
#' @description The exact fourth derivatives in \eqn{(\mu, \alpha)}, from
#' the compiled fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @param expected Logical; the expected version goes through the generic's
#'   strategies.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx,nsim Passed on when \code{expected} is \code{TRUE}.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{pig2_distrib}}
S7::method(distrib_deriv4, Pig2Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig_hd_block(y, theta,
               c(mu_mu_mu_mu = "d40", mu_mu_mu_alpha = "d31",
                 mu_mu_alpha_alpha = "d22", mu_alpha_alpha_alpha = "d13",
                 alpha_alpha_alpha_alpha = "d04"),
               pig2_hd_cpp)
}

#' @title Orthogonal Poisson-Inverse Gaussian Random Generation
#' @name distrib_rng.Pig2Distrib
#' @description The exact mixture sampler of
#' \code{\link[=distrib_rng.Pig1Distrib]{pig1}}, at the dispersion
#' \code{\link{pig2_sigma}} implies.
#' @param distrib A \code{Pig2Distrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{mu} and \code{alpha}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{pig2_distrib}}
S7::method(distrib_rng, Pig2Distrib) <- function(distrib, n, theta) {
  mu <- rep_len(theta[[1]], n)
  sg <- pig2_sigma(mu, rep_len(theta[[2]], n))
  lam <- statmod::rinvgauss(n, mean = mu, shape = mu / sg)
  stats::rpois(n, lam)
}

# --- CONSTRUCTOR WRAPPER ---

#' Poisson-Inverse Gaussian Distribution in Its Orthogonal Parametrization
#'
#' @description
#' Creates a Poisson-inverse Gaussian distribution in the parametrization
#' whose parameters are orthogonal -- the expected information is diagonal
#' -- which is gamlss's \code{PIG2}. The mean stays \eqn{\mu}; the second
#' parameter \eqn{\alpha} is exactly the argument the Bessel function of the
#' mass function is evaluated at, related to the dispersion of
#' \code{\link{pig1_distrib}} by \eqn{\alpha = \sqrt{1 + 2\sigma\mu}/\sigma}
#' (gamlss states the inverse of the same map,
#' \eqn{\alpha = 1/(\sqrt{\mu^2 + \sigma_2^2} - \mu)} with its own
#' \eqn{\sigma_2}, which coincides with this \eqn{\alpha}).
#'
#' @details
#' Orthogonality makes the maximum likelihood estimates of \eqn{\mu} and
#' \eqn{\alpha} asymptotically independent, so Fisher scoring steps in one
#' parameter do not disturb the other; the cost is that \eqn{\alpha} has no
#' moment reading of its own. Derivatives to fourth order are exact,
#' computed by the same compiled kernel as \code{\link{pig1_distrib}} with
#' \eqn{\alpha} a seed variable of the jet, so the Bessel argument needs no
#' chain rule at all.
#'
#' @section The distribution:
#' \deqn{P(Y=y) = \sqrt{\frac{2\alpha}{\pi}}\,\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\,y!}\,K_{y-1/2}(\alpha), \qquad \sigma = \frac{1}{\sqrt{\mu^{2}+\alpha^{2}} - \mu}}
#' on \eqn{y \in \{0, 1, \dots\}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \mu + \sigma\mu^{2}}
#'
#' @param link_mu The link for \eqn{\mu}; defaults to \code{log_link()}.
#' @param link_alpha The link for \eqn{\alpha}; defaults to \code{log_link()}.
#'
#' @return A \code{Pig2Distrib} object.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. \emph{Applied Statistics} 54(3), 507--554.
#'
#' Heller, G. Z., Couturier, D.-L., and Heritier, S. R. (2019). Beyond mean
#' modelling: bias due to misspecification of dispersion in Poisson-inverse
#' Gaussian regression. \emph{Biometrical Journal} 61(2), 333--342.
#'
#' @examples
#' d <- pig2_distrib()
#' theta <- list(mu = 3, alpha = 1.2)
#' distrib_pdf(d, 0:5, theta)
#' mean(d, theta)
#'
#' @seealso \code{\link{pig1_distrib}}, \code{\link{negbin2_distrib}}
#' @export
pig2_distrib <- function(link_mu = log_link(), link_alpha = log_link()) {
  Pig2Distrib(
    distrib_name = "poisson-inverse gaussian (orthogonal)",
    dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "alpha"),
    params_interpretation = c(mu = "mean", alpha = "bessel argument"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), alpha = c(0, Inf)),
    link_params = list(mu = link_mu, alpha = link_alpha)
  )
}
