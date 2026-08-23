#' @include distrib.R generics.R vonmises1_distrib.R moments.R
NULL

# The von Mises in its mean resultant length. The concentration and the
# resultant length are related by rho = A(kappa) = I_1(kappa)/I_0(kappa), a
# strictly increasing bijection from (0, Inf) onto (0, 1) whose inverse has no
# closed form. numericals7::bessel_i_ratio_inverse() obtains it by root
# finding and differentiates it by the inverse function rule, with A' to
# A'''' from the Bessel recurrences rather than four more evaluations.

#' @title S7 Class for the von Mises Distribution in Its Resultant Length
#' @name VonMises2Distrib
#'
#' @description A subclass of \code{continuous_distrib} for the von Mises
#'   written in its mean direction and its mean resultant length.
#' @inheritParams distrib
#' @return An object of class \code{VonMises2Distrib}.
#' @seealso \code{\link{vonmises2_distrib}}, \code{\link{vonmises1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_expected_hessian.VonMises2Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.VonMises2Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.VonMises2Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.VonMises2Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_rng.VonMises2Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
VonMises2Distrib <- S7::new_class("VonMises2Distrib", parent = continuous_distrib)

#' The Pieces a von Mises Derivative in rho Needs
#'
#' @description
#' The concentration, its four derivatives in \eqn{\rho}, and the derivatives
#' of \eqn{A} at that concentration, computed once per call.
#'
#' @param theta A list with \code{mu} and \code{rho}.
#'
#' @return A list with \code{kappa}, \code{kd} and \code{ad}.
#'
#' @seealso \code{\link{vonmises2_distrib}}
#'
#' @keywords internal
vm2_parts <- function(theta) {
  kd <- numericals7::bessel_i_ratio_inverse(theta[[2]])
  list(kappa = kd$kappa, kd = kd, ad = numericals7::bessel_i_ratio_derivs(kd$kappa))
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title von Mises Density in the Resultant Length
#' @name distrib_pdf.VonMises2Distrib
#' @description
#' The von Mises density at the concentration \eqn{\kappa = A^{-1}(\rho)}.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list with \code{mu} and \code{rho}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector.
#' @seealso \code{\link{vonmises2_distrib}}
S7::method(distrib_pdf, VonMises2Distrib) <- function(distrib, y, theta,
                                                     log = FALSE, ...,
                                                     threads = 1L) {
  distrib_pdf(vonmises1_distrib(), y,
              list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa), log = log,
              threads = threads)
}

#' @title von Mises Random Generation in the Resultant Length
#' @name distrib_rng.VonMises2Distrib
#' @description Delegates to the concentration parametrization.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param n The number of draws.
#' @param theta A list with \code{mu} and \code{rho}.
#' @return A numeric vector of angles.
#' @seealso \code{\link{vonmises2_distrib}}
S7::method(distrib_rng, VonMises2Distrib) <- function(distrib, n, theta) {
  distrib_rng(vonmises1_distrib(), n,
              list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa))
}

#' @title von Mises Analytical Gradient in the Resultant Length
#' @name distrib_gradient.VonMises2Distrib
#' @description
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
#'       \dfrac{\partial\ell}{\partial\rho}
#'         = \left\{\cos(y-\mu) - A(\kappa)\right\}\kappa'(\rho)}
#' The map touches only the second parameter, so the chain rule is the
#' one-variable one and no cancellation is involved.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list with \code{mu} and \code{rho}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{vonmises2_distrib}}
S7::method(distrib_gradient, VonMises2Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  list(mu = p$kappa * sin(d),
       rho = (cos(d) - p$ad$A) * p$kd$d1)
}

#' @title von Mises Analytical Observed Hessian in the Resultant Length
#' @name distrib_hessian.VonMises2Distrib
#' @description
#' The concentration parametrization's second derivatives carried through the
#' one-variable chain rule,
#' \deqn{\ell^{(\rho\rho)} = \ell^{(\kappa\kappa)}(\kappa')^2
#'                          + \ell^{(\kappa)}\kappa'',}
#' with \eqn{\ell^{(\kappa\kappa)} = -A'(\kappa)} and
#' \eqn{\ell^{(\mu\kappa)} = \sin(y-\mu)}.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list with \code{mu} and \code{rho}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{vonmises2_distrib}}
S7::method(distrib_hessian, VonMises2Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...) {
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  list(
    mu_mu = -p$kappa * cos(d),
    rho_rho = -p$ad$d1 * p$kd$d1^2 + (cos(d) - p$ad$A) * p$kd$d2,
    mu_rho = sin(d) * p$kd$d1
  )
}

#' @title von Mises Analytical Expected Hessian in the Resultant Length
#' @name distrib_expected_hessian.VonMises2Distrib
#' @description
#' Closed form. Since \eqn{\mathbb{E}[\cos(Y-\mu)] = A} and
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0}, the term in the second derivative of the
#' map drops and
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\kappa A, \qquad
#'       \mathbb{E}[\ell^{(\mu\rho)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\rho\rho)}] = -\dfrac{1}{A'(\kappa)}.}
#' The last is the inverse of the information in \eqn{\kappa}, which is what a
#' one-to-one reparametrization of a single parameter must give, and the two
#' parameters stay orthogonal.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list with \code{mu} and \code{rho}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{vonmises2_distrib}}
S7::method(distrib_expected_hessian, VonMises2Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  p <- vm2_parts(theta)
  n <- length(y)
  list(
    mu_mu = rep(-p$kappa * p$ad$A, length.out = n),
    rho_rho = rep(-p$ad$d1 * p$kd$d1^2, length.out = n),
    mu_rho = rep(0, length.out = n)
  )
}

#' @title Moments of the von Mises in the Resultant Length
#' @name mean.VonMises2Distrib
#' @description
#' The ordinary moments of \eqn{Y} as a number on \eqn{[-\pi, \pi)}, obtained
#' numerically as they are for \code{\link{vonmises1_distrib}}: \eqn{\mu} is
#' the mean \emph{direction} and \eqn{\rho} the mean resultant length, and
#' neither is an ordinary moment.
#' @param x A \code{VonMises2Distrib} object.
#' @param theta A list with \code{mu} and \code{rho}.
#' @param ... Passed on.
#' @return A numeric vector.
#' @seealso \code{\link{vonmises2_distrib}}
#' @keywords internal
S7::method(mean, VonMises2Distrib) <- function(x, theta, ...) {
  mean(vonmises1_distrib(), list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa), ...)
}


#' von Mises Distribution in the Mean Resultant Length
#'
#' @description
#' Creates a von Mises distribution object parametrized by its mean direction
#' and its \strong{mean resultant length} \eqn{\rho = A(\kappa)}, which lives
#' in \eqn{(0, 1)}.
#'
#' @details
#' The concentration \eqn{\kappa} of \code{\link{vonmises1_distrib}} is
#' unbounded and hard to read; the resultant length is bounded, is the quantity
#' circular statistics reports, and is one minus the circular variance. The two
#' are related by \eqn{\rho = I_1(\kappa)/I_0(\kappa)}, a strictly increasing
#' bijection.
#'
#' That bijection has no closed inverse, which is why this is a family of its
#' own rather than a \code{\link{reparametrize}} of the other: \eqn{\kappa} is
#' obtained by root finding on \eqn{\log\kappa}, and its four derivatives come
#' from the inverse function rule applied to \eqn{A' \dots A''''}, which the
#' Bessel recurrences give from the same two evaluations \eqn{A} already needs.
#'
#' The map touches the second parameter only, so the chain rule is the
#' one-variable one and the derivatives are exact. The expected information is
#' closed form and the two parameters are orthogonal, as they are in the
#' concentration parametrization.
#'
#' \strong{The moments are not the parameters.} \code{\link{mean}} returns the
#' ordinary expectation of \eqn{Y} on \eqn{[-\pi, \pi)}, which differs from
#' \eqn{\mu} whenever \eqn{\mu \ne 0}; see \code{\link{vonmises1_distrib}}.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{e^{\kappa\cos(y-\mu)}}{2\pi I_{0}(\kappa)}, \qquad \kappa = A^{-1}(\rho)}
#' on \eqn{y \in [-\pi, \pi)}.
#'
#' \deqn{\mathbb{E}[\cos(Y-\mu)] = \rho}
#'
#' @param link_mu Link function for the mean direction. Defaults to a link
#'   bounded to \eqn{(-\pi, \pi)}.
#' @param link_rho Link function for the resultant length. Defaults to the
#'   logit, the natural link onto \eqn{(0, 1)}.
#'
#' @return An S7 object of class \code{\link{VonMises2Distrib}}.
#'
#' @seealso \code{\link{vonmises1_distrib}}
#'
#' @examples
#' d <- vonmises2_distrib()
#' theta <- list(mu = 0.5, rho = 0.7)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#'
#' # rho is bounded, which is what makes it readable
#' d@params_bounds$rho
#'
#' @export
vonmises2_distrib <- function(link_mu = bounded_link(lwr = -pi, upr = pi),
                              link_rho = logit_link()) {
  VonMises2Distrib(
    distrib_name = "von mises2",
    dimension = "univariate",
    bounds = c(-pi, pi),
    params = c("mu", "rho"),
    params_interpretation = c(mu = "mean direction",
                              rho = "mean resultant length"),
    n_params = 2,
    params_bounds = list(mu = c(-pi, pi), rho = c(0, 1)),
    link_params = list(mu = link_mu, rho = link_rho)
  )
}


#' @title von Mises Third and Fourth Derivatives in the Resultant Length
#' @name distrib_deriv3.VonMises2Distrib
#' @description
#' Closed form at both orders. Every derivative that carries at least one
#' \eqn{\mu} is \eqn{D_a \kappa^{(b)}}, where \eqn{D_a} is the \eqn{a}-th
#' \eqn{\mu}-derivative of \eqn{\cos(y-\mu)} and \eqn{\kappa^{(b)}} the
#' \eqn{b}-th derivative of \eqn{A^{-1}}: the concentration parametrization's
#' \eqn{\mu}-derivatives are linear in \eqn{\kappa}, so the composition
#' collapses to a single term. The pure \eqn{\rho} components carry the full
#' one-variable chain rule on \eqn{\log I_0}.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{rho}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @param approx The approximation used when \code{expected} is \code{TRUE}.
#' @param nsim Monte Carlo draws when \code{approx = "mc"}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{vonmises2_distrib}}
S7::method(distrib_deriv3, VonMises2Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         expected = FALSE,
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  A <- p$ad
  k <- p$kd
  # d^3/drho^3 of log I_0(kappa(rho)), Faa di Bruno written out
  phi3 <- A$d2 * k$d1^3 + 3 * A$d1 * k$d1 * k$d2 + A$A * k$d3
  list(mu_mu_mu = -sin(d) * p$kappa,
       mu_mu_rho = -cos(d) * k$d1,
       mu_rho_rho = sin(d) * k$d2,
       rho_rho_rho = cos(d) * k$d3 - phi3)
}

#' @rdname distrib_deriv3.VonMises2Distrib
#' @name distrib_deriv4.VonMises2Distrib
#' @return A named list of fourth-derivative components.
S7::method(distrib_deriv4, VonMises2Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         expected = FALSE,
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  A <- p$ad
  k <- p$kd
  phi4 <- A$d3 * k$d1^4 + 6 * A$d2 * k$d1^2 * k$d2 +
    A$d1 * (3 * k$d2^2 + 4 * k$d1 * k$d3) + A$A * k$d4
  list(mu_mu_mu_mu = cos(d) * p$kappa,
       mu_mu_mu_rho = -sin(d) * k$d1,
       mu_mu_rho_rho = -cos(d) * k$d2,
       mu_rho_rho_rho = sin(d) * k$d3,
       rho_rho_rho_rho = cos(d) * k$d4 - phi4)
}


#' @title von Mises Distribution Function in the Resultant Length
#' @name distrib_cdf.VonMises2Distrib
#' @description
#' The concentration parametrization's series, read at the concentration this
#' one's resultant length implies.
#' @details
#' The map touches the second parameter only and the response not at all, so
#' the distribution function is the other family's at
#' \eqn{\kappa = A^{-1}(\rho)}. What it replaces is the base class's
#' quadrature, one per observation; see \code{\link{vm_cdf}} for the series
#' and for how many terms it takes.
#' @param distrib A \code{VonMises2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with \code{mu} and \code{rho}.
#' @param ... Unused.
#' @return The distribution function at \code{q}.
#' @seealso \code{\link{vm_cdf}}, \code{\link{vonmises2_distrib}}
#' @keywords internal
S7::method(distrib_cdf, VonMises2Distrib) <- function(distrib, q, theta,
                                                      ...) {
  vm_cdf(q, theta[[1]], vm2_parts(theta)$kappa)
}
