#' @include distrib.R generics.R
NULL

#' @title S7 Class for the von Mises Distribution
#' @name VonMisesDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the von
#'   Mises distribution on the circle.
#' @inheritParams distrib
#' @return An object of class \code{VonMisesDistrib}.
#' @seealso \code{\link{vonmises_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_expected_hessian.VonMisesDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.VonMisesDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_grad_y.VonMisesDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_hessian.VonMisesDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_hess_y.VonMisesDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_pdf.VonMisesDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_rng.VonMisesDistrib]{distrib_rng()}}
#'
#' The distribution function and the quantile come from
#' \code{\link{continuous_distrib}}, by quadrature and root finding over the
#' bounded support.
VonMisesDistrib <- S7::new_class("VonMisesDistrib", parent = continuous_distrib)

#' The Mean Resultant Length of a von Mises
#'
#' @description
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}, the expected cosine of the
#' deviation from the location.
#'
#' @details
#' Both Bessel functions are taken exponentially scaled, so the factor
#' \eqn{e^{\kappa}} they share cancels in the ratio and the result stays
#' finite for a concentration of any size, where the unscaled functions
#' overflow past about \eqn{\kappa = 700}.
#'
#' @param kappa The concentration, a positive numeric vector.
#'
#' @return A numeric vector in \eqn{(0, 1)}.
#'
#' @seealso \code{\link{vonmises_distrib}}, \code{\link{vm_dA}}
#'
#' @keywords internal
vm_A <- function(kappa) {
  besselI(kappa, 1, expon.scaled = TRUE) /
    besselI(kappa, 0, expon.scaled = TRUE)
}

#' The Derivative of the Mean Resultant Length
#'
#' @description
#' \eqn{A'(\kappa) = 1 - A(\kappa)/\kappa - A(\kappa)^2}, which is also the
#' variance of \eqn{\cos(Y - \mu)} and therefore positive.
#'
#' @details
#' The identity follows from \eqn{I_0' = I_1} and
#' \eqn{I_1' = I_0 - I_1/\kappa}, so no further Bessel evaluation is needed.
#'
#' @param kappa The concentration, a positive numeric vector.
#' @param A The value of \code{\link{vm_A}} at \code{kappa}, passed in when it
#'   has already been computed.
#'
#' @return A numeric vector.
#'
#' @seealso \code{\link{vonmises_distrib}}, \code{\link{vm_A}}
#'
#' @keywords internal
vm_dA <- function(kappa, A = vm_A(kappa)) 1 - A / kappa - A * A

# --- S7 METHODS IMPLEMENTATION ---

#' @title von Mises Density
#' @name distrib_pdf.VonMisesDistrib
#' @description
#' \deqn{f(y) = \dfrac{e^{\kappa \cos(y - \mu)}}{2\pi I_0(\kappa)},
#'       \qquad y \in [-\pi, \pi)}
#' @param distrib A \code{VonMisesDistrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_pdf, VonMisesDistrib) <- function(distrib, y, theta, log = FALSE) {
  k <- theta[[2]]
  # log I_0 through the scaled Bessel: the exponent is added back rather than
  # left to overflow, which it does past about kappa = 700.
  log_i0 <- log(besselI(k, 0, expon.scaled = TRUE)) + k
  out <- k * cos(y - theta[[1]]) - log(2 * pi) - log_i0
  out[y < -pi | y >= pi] <- -Inf
  if (log) out else exp(out)
}

#' @title von Mises Random Generation
#' @name distrib_rng.VonMisesDistrib
#' @description
#' The rejection algorithm of Best and Fisher (1979), which draws from a
#' wrapped Cauchy envelope and accepts with a probability that does not
#' involve the Bessel function at all.
#' @param distrib A \code{VonMisesDistrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @return A numeric vector of length \code{n}, in \eqn{[-\pi, \pi)}.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_rng, VonMisesDistrib) <- function(distrib, n, theta) {
  mu <- theta[[1]]
  k <- theta[[2]]
  a <- 1 + sqrt(1 + 4 * k * k)
  b <- (a - sqrt(2 * a)) / (2 * k)
  r <- (1 + b * b) / (2 * b)

  out <- numeric(0)
  while (length(out) < n) {
    m <- 2 * (n - length(out)) + 16
    u1 <- stats::runif(m); u2 <- stats::runif(m); u3 <- stats::runif(m)
    z <- cos(pi * u1)
    f <- (1 + r * z) / (r + z)
    c0 <- k * (r - f)
    ok <- (c0 * (2 - c0) - u2 > 0) | (log(c0 / u2) + 1 - c0 >= 0)
    ang <- sign(u3 - 0.5) * acos(pmin(pmax(f, -1), 1))
    out <- c(out, ang[ok])
  }
  # wrapped back into the declared support
  ((out[seq_len(n)] + mu + pi) %% (2 * pi)) - pi
}

#' @title von Mises Analytical Gradient
#' @name distrib_gradient.VonMisesDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \kappa \sin(y - \mu), \qquad
#'       \dfrac{\partial \ell}{\partial \kappa} = \cos(y - \mu) - A(\kappa)}
#' with \eqn{A = I_1/I_0}, the derivative of \eqn{\log I_0}.
#' @param distrib A \code{VonMisesDistrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{mu} and \code{kappa} components.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_gradient, VonMisesDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  d <- y - theta[[1]]
  list(mu = theta[[2]] * sin(d), kappa = cos(d) - vm_A(theta[[2]]))
}

#' @title von Mises Analytical Observed Hessian
#' @name distrib_hessian.VonMisesDistrib
#' @description
#' \deqn{\ell^{(\mu\mu)} = -\kappa\cos(y-\mu), \qquad
#'       \ell^{(\mu\kappa)} = \sin(y-\mu), \qquad
#'       \ell^{(\kappa\kappa)} = -A'(\kappa)}
#' the last free of the data, \eqn{\log I_0} being the only place
#' \eqn{\kappa} appears alone.
#' @param distrib A \code{VonMisesDistrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_hessian, VonMisesDistrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...) {
  d <- y - theta[[1]]
  k <- theta[[2]]
  list(mu_mu = -k * cos(d), mu_kappa = sin(d),
       kappa_kappa = rep_len(-vm_dA(k), length(d)))
}

#' @title von Mises Analytical Expected Hessian
#' @name distrib_expected_hessian.VonMisesDistrib
#' @description
#' \eqn{\mathbb{E}[\cos(Y-\mu)] = A(\kappa)} and
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0} by symmetry, so
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\kappa A(\kappa), \qquad
#'       \mathbb{E}[\ell^{(\mu\kappa)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\kappa\kappa)}] = -A'(\kappa)}
#' The location and the concentration are therefore orthogonal.
#' @param distrib A \code{VonMisesDistrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_expected_hessian, VonMisesDistrib) <- function(distrib, y, theta,
                                                                   scale = c("parameter", "link"),
                                                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                                                   nsim = 10000, ...) {
  k <- theta[[2]]
  A <- vm_A(k)
  n <- length(y)
  list(mu_mu = rep_len(-k * A, n), mu_kappa = rep_len(0, n),
       kappa_kappa = rep_len(-vm_dA(k, A), n))
}

#' @title von Mises Response Gradient
#' @name distrib_grad_y.VonMisesDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\kappa \sin(y - \mu)}
#' @param distrib A \code{VonMisesDistrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_grad_y, VonMisesDistrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * sin(y - theta[[1]])
}

#' @title von Mises Response Hessian
#' @name distrib_hess_y.VonMisesDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\kappa \cos(y - \mu)}
#' @param distrib A \code{VonMisesDistrib} object.
#' @param y A numeric vector of angles.
#' @param theta A list containing \code{mu} and \code{kappa}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{vonmises_distrib}}
S7::method(distrib_hess_y, VonMisesDistrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * cos(y - theta[[1]])
}

# --- CONSTRUCTOR WRAPPER ---

#' von Mises Distribution Object
#'
#' @description
#' Creates a distribution object for the von Mises distribution, the natural
#' family for an angle, parametrised by a mean direction \eqn{\mu} and a
#' concentration \eqn{\kappa}.
#'
#' @param link_mu A link function object for \eqn{\mu}. Defaults to
#'   \code{\link[linkfunctions7]{bounded_link}} on \eqn{(-\pi, \pi)}.
#' @param link_kappa A link function object for \eqn{\kappa}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The observation is an angle and the support is a circle, written here as
#' \eqn{[-\pi, \pi)}. This is the first family in the package whose support has
#' that topology: the two ends of the interval are the same point, so a density
#' need not vanish at either.
#'
#' \strong{Density:}
#' \deqn{f(y) = \dfrac{e^{\kappa\cos(y-\mu)}}{2\pi I_0(\kappa)}}
#' The normalising constant is a modified Bessel function, and it is evaluated
#' exponentially scaled with the exponent added back, so a concentration past
#' \eqn{\kappa = 700} does not overflow.
#'
#' \strong{Score, observed and expected Hessian.} Writing
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} for the derivative of
#' \eqn{\log I_0},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
#'       \dfrac{\partial\ell}{\partial\kappa} = \cos(y-\mu) - A(\kappa),}
#' and every second derivative is closed form, with
#' \eqn{A'(\kappa) = 1 - A(\kappa)/\kappa - A(\kappa)^2} obtained from the
#' Bessel recurrences rather than from a further evaluation. Since
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0} by symmetry, \strong{the two parameters
#' are orthogonal}: the expected information is diagonal, and Fisher scoring
#' updates the direction and the concentration independently.
#'
#' \eqn{A(\kappa)} is the mean resultant length, so \eqn{A'(\kappa)} is the
#' variance of \eqn{\cos(Y-\mu)} and is positive, which is what makes the
#' information positive definite.
#'
#' \strong{The mean direction is carried on a bounded chart}, the default link
#' mapping the free scale onto \eqn{(-\pi, \pi)}. That keeps the parameter
#' identified, at the cost that a fit cannot walk across the boundary: data
#' concentrated near \eqn{\pm\pi} are better rotated before fitting than
#' handed to the optimiser as they are. Leaving \eqn{\mu} unbounded instead
#' would make the likelihood periodic and every maximum one of infinitely
#' many.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\pi, \pi)}
#'   \item \eqn{\kappa \in (0, +\infty)}
#' }
#'
#' The distribution function has no elementary form and comes from the base
#' class by quadrature over the bounded support, with the quantile by root
#' finding on it.
#'
#' @return An S7 object of class \code{VonMisesDistrib}.
#'
#' @references
#' Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von Mises
#' distribution. \emph{Applied Statistics} 28, 152-157.
#'
#' Mardia, K. V. and Jupp, P. E. (2000). \emph{Directional Statistics}. Wiley.
#'
#' @seealso \code{\link{gaussian_distrib}} for the analogous family on the line
#'
#' @importFrom linkfunctions7 bounded_link log_link
#' @importFrom stats runif
#' @examples
#' d <- vonmises_distrib()
#' d@params
#'
#' theta <- list(mu = 0.5, kappa = 2)
#' distrib_pdf(d, c(-1, 0, 0.5, 2), theta)
#'
#' # the expected information is diagonal: direction and concentration are
#' # orthogonal, the sine having mean zero
#' distrib_expected_hessian(d, 0, theta)
#'
#' @export
vonmises_distrib <- function(link_mu = bounded_link(lwr = -pi, upr = pi),
                             link_kappa = log_link()) {
  VonMisesDistrib(
    distrib_name = "von mises", dimension = "univariate",
    bounds = c(-pi, pi),
    params = c("mu", "kappa"),
    params_interpretation = c(mu = "mean direction", kappa = "concentration"),
    n_params = 2,
    params_bounds = list(mu = c(-pi, pi), kappa = c(0, Inf)),
    link_params = list(mu = link_mu, kappa = link_kappa)
  )
}
