#' @include gpd_distrib.R reparam_maps.R
NULL

# The generalized Pareto's third and fourth derivatives, written out.
#
# Splitting the log-density as
#
#   l = -log(sigma) - L - W,   L = log(t),  W = L/xi,  t = 1 + xi y/sigma,
#
# separates the part that is analytic at xi = 0 from the part that is not
# written that way. L is a logarithm of a function whose partials are
# elementary and vanish beyond first order in xi, so fdb2()'s written-out
# template covers it. W = L/xi has a removable singularity: the Leibniz form
#
#   d^{a,b} W = sum_j C(b,j) d^{a,j} L * d^{b-j}(1/xi),
#
# carries terms of size xi^-(b+1) that cancel against each other, so below a
# measured threshold it is replaced by the series the cancellation leaves,
#
#   W = sum_{k >= 0} (-1)^k xi^k z^(k+1) / (k+1),   z = y/sigma,
#
# differentiated term by term. The threshold is on xi*z and not on xi: at
# order b the Leibniz terms are of size z xi^-b against an answer of size
# z^(b+1), so the relative cancellation is (xi z)^-b and the accuracy is
# eps (xi z)^-b, which at order four reaches 0.3 by xi z = 1.7e-4 -- measured,
# and matching that prediction is what fixed the switch. The series needs
# |xi z| < 1 to converge and is cut at 0.2, where forty terms leave 1e-28 and
# the Leibniz form still carries 1e-13. The two agree in the overlap.

#' Derivative Components of the Generalized Pareto
#'
#' @description
#' The components of \eqn{\partial^{a+b}\ell/\partial\sigma^a\partial\xi^b} at
#' any order from one to four.
#'
#' @details
#' The shape direction goes through the series of \eqn{\log(1+\xi z)/\xi}
#' wherever \eqn{\lvert\xi z\rvert} is below \code{cut}, which is the
#' region where the Leibniz form's terms of size \eqn{\xi^{-(b+1)}} cancel
#' against each other; elsewhere the Leibniz form is used directly. The
#' threshold is exposed so that a test can force either route where both are
#' accurate and compare them.
#'
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param order The derivative order, 1 to 4.
#' @param cut The value of \eqn{\lvert\xi z\rvert} below which the series
#'   is used.
#'
#' @return A named list of component vectors, keyed as
#'   \code{\link{deriv_names}}.
#'
#' @seealso \code{\link{gpd_distrib}}, \code{\link{fdb2}}
#' @keywords internal
gpd_components <- function(y, theta, order, cut = 0.2) {
  sg <- theta[[1]]
  xi <- theta[[2]]
  n <- max(length(y), lengths(theta[1:2]))
  y <- rep_len(y, n)
  sg <- rep_len(sg, n)
  xi <- rep_len(xi, n)
  z <- y / sg
  t <- 1 + xi * z
  one <- rep_len(1, n)

  # L = log(t); t is affine in xi and a power series in 1/sigma, so every
  # partial of t with two or more xi vanishes. x stands for sigma, z for xi.
  Ld <- fdb2(
    list(one / t, -one / t^2, 2 / t^3, -6 / t^4),
    list(x = -xi * y / sg^2, z = y / sg,
         xx = 2 * xi * y / sg^3, xz = -y / sg^2,
         xxx = -6 * xi * y / sg^4, xxz = 2 * y / sg^3,
         xxxx = 24 * xi * y / sg^5, xxxz = -6 * y / sg^4)
  )
  Lget <- function(a, b) {
    if (a == 0L && b == 0L) return(base::log(t))
    v <- Ld[[paste0(strrep("x", a), strrep("z", b))]]
    if (is.null(v)) rep_len(0, n) else v * one
  }

  # W = L/xi, by Leibniz away from zero and by its own series near it
  small <- abs(xi * z) < cut
  rising <- function(m, a) {
    out <- rep_len(1, length(m))
    if (a >= 1L) for (i in 0:(a - 1L)) out <- out * (m + i)
    out
  }
  Wget <- function(a, b) {
    out <- numeric(n)
    if (any(!small)) {
      acc <- numeric(n)
      for (j in 0:b) {
        m <- b - j
        acc <- acc + choose(b, j) * Lget(a, j) * (-1)^m * factorial(m) / xi^(m + 1L)
      }
      out[!small] <- acc[!small]
    }
    if (any(small)) {
      acc <- numeric(n)
      for (k in b:40L) {
        acc <- acc + (-1)^k * (factorial(k) / factorial(k - b)) * xi^(k - b) *
          (-1)^a * rising(k + 1, a) * z^(k + 1) / ((k + 1) * sg^a)
      }
      out[small] <- acc[small]
    }
    out
  }

  comp <- function(a, b) {
    out <- -Lget(a, b) - Wget(a, b)
    if (b == 0L) {
      out <- out - if (a == 0L) base::log(sg) else
        (-1)^(a - 1L) * factorial(a - 1L) / sg^a
    }
    out
  }

  nms <- deriv_names(c("sigma", "xi"), order)
  stats::setNames(lapply(nms, function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    comp(sum(parts == "sigma"), sum(parts == "xi"))
  }), nms)
}


#' @title Generalized Pareto Third and Fourth Derivatives
#' @name distrib_deriv3.GPDDistrib
#' @description
#' Closed form at both orders, from \code{\link{gpd_components}}: the
#' log-density splits into \eqn{-\log\sigma}, \eqn{-\log t} and
#' \eqn{-\log(t)/\xi}, and the last is taken from its series where the
#' Leibniz form's terms cancel.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx The approximation used when \code{expected} is \code{TRUE}.
#' @param nsim Monte Carlo draws when \code{approx = "mc"}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_deriv3, GPDDistrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  gpd_components(y, theta, 3L)
}

#' @rdname distrib_deriv3.GPDDistrib
#' @name distrib_deriv4.GPDDistrib
#' @return A named list of fourth-derivative components.
S7::method(distrib_deriv4, GPDDistrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  gpd_components(y, theta, 4L)
}
