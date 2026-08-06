#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Poisson-Inverse Gaussian Distribution
#' @name Pig1Distrib
#'
#' @description A subclass of \code{discrete_distrib} representing the
#'   Poisson-inverse Gaussian distribution on \eqn{\{0, 1, 2, \dots\}} in its
#'   mean-dispersion parametrization, gamlss's \code{PIG}.
#' @inheritParams distrib
#' @return An object of class \code{Pig1Distrib}.
#' @seealso \code{\link{pig1_distrib}}, \code{\link{pig2_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_pdf.Pig1Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_gradient.Pig1Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.Pig1Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_deriv3.Pig1Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Pig1Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_rng.Pig1Distrib]{distrib_rng()}}
#'
#' Everything else, the distribution function and the quantile included, is
#' inherited from \code{\link{discrete_distrib}}.
Pig1Distrib <- S7::new_class("Pig1Distrib", parent = discrete_distrib)

#' Log-Likelihood Derivatives of the Poisson-Inverse Gaussian
#'
#' @description
#' Evaluates the log-likelihood and its fourteen partial derivatives to
#' fourth order at once, through the compiled kernel, and returns the block
#' of components an order asks for.
#'
#' @details
#' With \eqn{c = 1 + 2\sigma\mu} and \eqn{\alpha = \sqrt{c}/\sigma}, the
#' half-integer order collapses the Bessel function to a finite sum and the
#' log-likelihood to
#' \deqn{\ell(y) = y\log\mu - \tfrac{y}{2}\log c + \tfrac{1}{\sigma}
#'   + \psi(\alpha) - \log y!,\qquad
#'   \psi(\alpha) = -\alpha + \log S_y(\alpha),}
#' where \eqn{S_y} sums \eqn{y} positive terms on the log scale. The kernel
#' carries a bivariate jet truncated at total order four through this
#' expression, so every partial is exact and no chain rule is transcribed.
#'
#' @param y A numeric vector of observations.
#' @param theta The aligned parameter list.
#' @param cols The kernel columns wanted, by name.
#' @param kernel The compiled kernel, \code{pig1_hd_cpp} or
#'   \code{pig2_hd_cpp}.
#' @return A named list of derivative component vectors.
#' @keywords internal
pig_hd_block <- function(y, theta, cols, kernel) {
  n <- max(length(y), length(theta[[1]]), length(theta[[2]]))
  y <- rep_len(y, n)
  p1 <- rep_len(theta[[1]], n)
  p2 <- rep_len(theta[[2]], n)
  ok <- is.finite(y) & y >= 0 & y == floor(y)
  m <- matrix(NaN, n, 15)
  if (any(ok)) m[ok, ] <- kernel(y[ok], p1[ok], p2[ok])
  colnames(m) <- c("l", "d10", "d01", "d20", "d11", "d02",
                   "d30", "d21", "d12", "d03",
                   "d40", "d31", "d22", "d13", "d04")
  stats::setNames(lapply(cols, function(cc) m[, cc]), names(cols))
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Poisson-Inverse Gaussian Probability Mass Function
#' @name distrib_pdf.Pig1Distrib
#' @description
#' With \eqn{c = 1 + 2\sigma\mu} and \eqn{\alpha = \sqrt{c}/\sigma},
#' \deqn{P(Y = y) = \sqrt{\dfrac{2\alpha}{\pi}}\,
#'   \dfrac{\mu^y e^{1/\sigma}}{(\alpha\sigma)^y\, y!}\, K_{y-1/2}(\alpha),}
#' evaluated through the finite half-integer Bessel sum, in which the
#' prefactors cancel down to
#' \eqn{\ell(y) = y\log\mu - (y/2)\log c + (1-\sqrt{c})/\sigma
#'   + \log S_y(\alpha) - \log y!}. A non-integer or negative \eqn{y} has
#' probability zero.
#' @param distrib A \code{Pig1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probability values.
#' @seealso \code{\link{pig1_distrib}}
S7::method(distrib_pdf, Pig1Distrib) <- function(distrib, y, theta, log = FALSE) {
  out <- pig_hd_block(y, theta, c(l = "l"), pig1_hd_cpp)$l
  out[is.nan(out)] <- -Inf
  if (log) out else exp(out)
}

#' @title Poisson-Inverse Gaussian Analytical Gradient
#' @name distrib_gradient.Pig1Distrib
#' @description The exact score in \eqn{(\mu, \sigma)}, from the compiled
#' fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu} and \code{sigma} components.
#' @seealso \code{\link{pig1_distrib}}
S7::method(distrib_gradient, Pig1Distrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ...) {
  pig_hd_block(y, theta, c(mu = "d10", sigma = "d01"), pig1_hd_cpp)
}

#' @title Poisson-Inverse Gaussian Analytical Observed Hessian
#' @name distrib_hessian.Pig1Distrib
#' @description The exact second derivatives in \eqn{(\mu, \sigma)}, from
#' the compiled fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{pig1_distrib}}
S7::method(distrib_hessian, Pig1Distrib) <- function(distrib, y, theta,
                                                     scale = c("parameter", "link"), ...) {
  pig_hd_block(y, theta,
               c(mu_mu = "d20", sigma_sigma = "d02", mu_sigma = "d11"),
               pig1_hd_cpp)
}

#' @title Poisson-Inverse Gaussian Analytical Third Derivatives
#' @name distrib_deriv3.Pig1Distrib
#' @description The exact third derivatives in \eqn{(\mu, \sigma)}, from the
#' compiled fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param expected Logical; the expected version goes through the generic's
#'   strategies.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx,nsim Passed on when \code{expected} is \code{TRUE}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{pig1_distrib}}
S7::method(distrib_deriv3, Pig1Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig_hd_block(y, theta,
               c(mu_mu_mu = "d30", mu_mu_sigma = "d21",
                 mu_sigma_sigma = "d12", sigma_sigma_sigma = "d03"),
               pig1_hd_cpp)
}

#' @title Poisson-Inverse Gaussian Analytical Fourth Derivatives
#' @name distrib_deriv4.Pig1Distrib
#' @description The exact fourth derivatives in \eqn{(\mu, \sigma)}, from
#' the compiled fourth-order kernel described in \code{\link{pig_hd_block}}.
#' @param distrib A \code{Pig1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param expected Logical; the expected version goes through the generic's
#'   strategies.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx,nsim Passed on when \code{expected} is \code{TRUE}.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{pig1_distrib}}
S7::method(distrib_deriv4, Pig1Distrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  pig_hd_block(y, theta,
               c(mu_mu_mu_mu = "d40", mu_mu_mu_sigma = "d31",
                 mu_mu_sigma_sigma = "d22", mu_sigma_sigma_sigma = "d13",
                 sigma_sigma_sigma_sigma = "d04"),
               pig1_hd_cpp)
}

#' @title Poisson-Inverse Gaussian Random Generation
#' @name distrib_rng.Pig1Distrib
#' @description
#' Exact mixture sampling: \eqn{\lambda} is drawn from the inverse Gaussian
#' with mean \eqn{\mu} and shape \eqn{\mu/\sigma}, whose variance
#' \eqn{\sigma\mu^2} is exactly what the mixing construction requires, and
#' \eqn{Y \mid \lambda} from the Poisson.
#' @param distrib A \code{Pig1Distrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{pig1_distrib}}
S7::method(distrib_rng, Pig1Distrib) <- function(distrib, n, theta) {
  lam <- statmod::rinvgauss(n, mean = rep_len(theta[[1]], n),
                            shape = rep_len(theta[[1]], n) / rep_len(theta[[2]], n))
  stats::rpois(n, lam)
}

# --- CONSTRUCTOR WRAPPER ---

#' Poisson-Inverse Gaussian Distribution Object
#'
#' @description
#' Creates a Poisson-inverse Gaussian distribution in its mean-dispersion
#' parametrization: mean \eqn{\mu} and variance \eqn{\mu + \sigma\mu^2},
#' gamlss's \code{PIG} (Rigby and Stasinopoulos, 2005). The family is the
#' Poisson mixed over an inverse Gaussian rate, an overdispersed count model
#' with a heavier tail than the negative binomial at the same variance.
#'
#' @details
#' The mass function carries the modified Bessel function
#' \eqn{K_{y-1/2}}, which at half-integer order is a finite sum; the
#' log-likelihood and its derivatives to fourth order are exact, computed by
#' a compiled kernel that carries a bivariate jet through the closed
#' expression. The expected information has no closed form and goes through
#' the summation strategies of
#' \code{\link{expected_derivative_methods}}. For the parametrization in
#' which \eqn{\mu} and the second parameter are orthogonal, see
#' \code{\link{pig2_distrib}}.
#'
#' @param link_mu The link for \eqn{\mu}; defaults to \code{log_link()}.
#' @param link_sigma The link for \eqn{\sigma}; defaults to \code{log_link()}.
#'
#' @return A \code{Pig1Distrib} object.
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. \emph{Applied Statistics} 54(3), 507--554.
#'
#' Dean, C., Lawless, J. F., and Willmot, G. E. (1989). A mixed
#' Poisson-inverse-Gaussian regression model. \emph{Canadian Journal of
#' Statistics} 17(2), 171--181.
#'
#' @examples
#' d <- pig1_distrib()
#' theta <- list(mu = 3, sigma = 0.8)
#' distrib_pdf(d, 0:5, theta)
#' mean(d, theta)
#' variance(d, theta)
#'
#' @export
pig1_distrib <- function(link_mu = log_link(), link_sigma = log_link()) {
  Pig1Distrib(
    distrib_name = "poisson-inverse gaussian", dimension = "univariate",
    bounds = c(0, Inf),
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean", sigma = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), sigma = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma)
  )
}
