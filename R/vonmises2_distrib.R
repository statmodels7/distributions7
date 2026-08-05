#' @include distrib.R generics.R vonmises1_distrib.R moments.R
NULL

# The von Mises in its mean resultant length. The concentration and the
# resultant length are related by rho = A(kappa) = I_1(kappa)/I_0(kappa), a
# strictly increasing bijection from (0, Inf) onto (0, 1) whose inverse has no
# closed form. It is obtained by root finding and differentiated by the inverse
# function rule, which needs A' to A'''' -- and those come from the Bessel
# recurrences rather than from four more Bessel evaluations.

#' Higher Derivatives of the Mean Resultant Length
#'
#' @description
#' \eqn{A''}, \eqn{A'''} and \eqn{A''''} at \eqn{\kappa}, obtained by
#' differentiating \eqn{A' = 1 - A/\kappa - A^2} repeatedly.
#'
#' @details
#' Each order is written in the orders below it, so the whole table costs the
#' two Bessel functions \code{\link{vm_A}} already evaluates and nothing more:
#' \deqn{A'' = -\dfrac{A'}{\kappa} + \dfrac{A}{\kappa^2} - 2AA'}
#' and so on. The alternative, evaluating a Bessel function of higher order for
#' each derivative, costs more and is less accurate at large \eqn{\kappa},
#' where the functions themselves overflow and only their ratio does not.
#'
#' @param kappa The concentration.
#'
#' @return A named list with \code{A} and its four derivatives.
#'
#' @seealso \code{\link{vonmises2_distrib}}, \code{\link{vm_A}}
#'
#' @keywords internal
vm_A_derivs <- function(kappa) {
  k <- kappa
  A <- vm_A(k)
  d1 <- 1 - A / k - A * A
  d2 <- -d1 / k + A / k^2 - 2 * A * d1
  d3 <- -d2 / k + 2 * d1 / k^2 - 2 * A / k^3 - 2 * d1^2 - 2 * A * d2
  d4 <- -d3 / k + 3 * d2 / k^2 - 6 * d1 / k^3 + 6 * A / k^4 -
    6 * d1 * d2 - 2 * A * d3
  list(A = A, d1 = d1, d2 = d2, d3 = d3, d4 = d4)
}

#' The Concentration a Mean Resultant Length Implies
#'
#' @description
#' \eqn{\kappa = A^{-1}(\rho)}, by root finding, together with the four
#' derivatives of the inverse.
#'
#' @details
#' \eqn{A} has no elementary inverse, so \eqn{\kappa} is found by bisection on
#' \eqn{\log\kappa}, where the function is well conditioned over the whole
#' range. The derivatives then come from the inverse function rule, which needs
#' no further root finding:
#' \deqn{\kappa' = \dfrac{1}{A'}, \qquad
#'       \kappa'' = -\dfrac{A''}{(A')^3}, \qquad
#'       \kappa''' = \dfrac{3(A'')^2 - A'A'''}{(A')^5},}
#' and the fourth in the same pattern. \eqn{A'} is the variance of
#' \eqn{\cos(Y-\mu)} and therefore strictly positive, so none of these divides
#' by zero in the interior.
#'
#' @param rho The mean resultant length, in \eqn{(0, 1)}.
#'
#' @return A named list with \code{kappa} and its four derivatives in
#'   \code{rho}.
#'
#' @seealso \code{\link{vonmises2_distrib}}
#'
#' @keywords internal
vm_kappa_of_rho <- function(rho) {
  one <- function(r) {
    if (!is.finite(r) || r <= 0 || r >= 1) return(NA_real_)
    # A starting value from the usual approximation, then a bracket widened
    # around it. A blind bracket is not an option: the scaled Bessel functions
    # underflow past about kappa = 1e15 and their ratio comes back NaN, so the
    # search has to stay where A can be evaluated at all.
    g <- if (r < 0.53) {
      2 * r + r^3 + 5 * r^5 / 6
    } else if (r < 0.85) {
      -0.4 + 1.39 * r + 0.43 / (1 - r)
    } else {
      1 / (r^3 - 4 * r^2 + 3 * r)
    }
    g <- min(max(g, 1e-8), 1e12)
    f <- function(k) vm_A(k) - r
    lo <- g / 2
    hi <- g * 2
    it <- 0L
    while (f(lo) > 0 && lo > 1e-10 && it < 60L) {
      lo <- lo / 2
      it <- it + 1L
    }
    it <- 0L
    while (f(hi) < 0 && hi < 1e13 && it < 60L) {
      hi <- hi * 2
      it <- it + 1L
    }
    stats::uniroot(f, c(lo, hi), tol = .Machine$double.eps^0.75)$root
  }
  k <- vapply(rho, one, numeric(1))
  a <- vm_A_derivs(k)
  p1 <- a$d1
  list(
    kappa = k,
    d1 = 1 / p1,
    d2 = -a$d2 / p1^3,
    d3 = (3 * a$d2^2 - p1 * a$d3) / p1^5,
    d4 = (-15 * a$d2^3 + 10 * p1 * a$d2 * a$d3 - p1^2 * a$d4) / p1^7
  )
}

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
  kd <- vm_kappa_of_rho(theta[[2]])
  list(kappa = kd$kappa, kd = kd, ad = vm_A_derivs(kd$kappa))
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
S7::method(distrib_pdf, VonMises2Distrib) <- function(distrib, y, theta, log = FALSE) {
  distrib_pdf(vonmises1_distrib(), y,
              list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa), log = log)
}

#' @title von Mises Random Generation in the Resultant Length
#' @name distrib_rng.VonMises2Distrib
#' @description Delegates to the concentration parametrisation.
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
#' The concentration parametrisation's second derivatives carried through the
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
#' one-to-one reparametrisation of a single parameter must give, and the two
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
#' Creates a von Mises distribution object parametrised by its mean direction
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
#' concentration parametrisation.
#'
#' \strong{The moments are not the parameters.} \code{\link{mean}} returns the
#' ordinary expectation of \eqn{Y} on \eqn{[-\pi, \pi)}, which differs from
#' \eqn{\mu} whenever \eqn{\mu \ne 0}; see \code{\link{vonmises1_distrib}}.
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
