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
#' wherever \eqn{\lvert\xi z\rvert} is below `cut`, which is the
#' region where the Leibniz form's terms of size \eqn{\xi^{-(b+1)}} cancel
#' against each other; elsewhere the Leibniz form is used directly. The
#' threshold is exposed so that a test can force either route where both are
#' accurate and compare them.
#'
#' @param y A numeric vector of observations.
#' @param theta A list containing `sigma` and `xi`.
#' @param order The derivative order, 1 to 4.
#' @param cut The value of \eqn{\lvert\xi z\rvert} below which the series
#'   is used.
#' @param threads How many threads the series kernel may use.
#'
#' @return A named list of component vectors, keyed as
#'   [deriv_names()].
#'
#' @seealso [gpd_distrib()], [fdb2()]
#' @keywords internal
gpd_components <- function(y, theta, order, cut = 0.2, threads = 1L) {
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
  # Each branch runs on ITS OWN elements. Both used to be evaluated over the
  # whole vector and subset afterwards, so a sample straddling the cut paid
  # for both in full -- the arithmetic is elementwise, so subsetting first
  # returns the same numbers and does a fraction of the work.
  big_i <- which(!small)
  small_i <- which(small)
  Wget <- function(a, b) {
    out <- numeric(n)
    if (length(big_i)) {
      xi_i <- xi[big_i]
      acc <- numeric(length(big_i))
      for (j in 0:b) {
        m <- b - j
        acc <- acc + choose(b, j) * Lget(a, j)[big_i] * (-1)^m * factorial(m) /
          xi_i^(m + 1L)
      }
      out[big_i] <- acc
    }
    if (length(small_i)) {
      z_i <- z[small_i]
      u_i <- xi[small_i] * z_i
      # The two powers the loop used to raise PER ELEMENT are one power:
      # with u = xi z, xi^(k-b) z^(k+1) = u^(k-b) z^(b+1), so z^(b+1)/sigma^a
      # leaves the loop entirely and what is left is a POLYNOMIAL IN u whose
      # coefficients are scalar in k. Evaluating it is a scalar recursion of
      # forty-one steps per element, which is what gpd_poly_cpp() does by
      # Horner -- from the highest power down, so a decaying series is summed
      # smallest-first where the loop this replaces summed largest-first.
      # |u| < cut on this branch, so the powers decay and cannot overflow.
      ck <- vapply(b:40L, function(k) {
        (-1)^k * (factorial(k) / factorial(k - b)) * rising(k + 1, a) / (k + 1)
      }, 0)
      pref <- (-1)^a * z_i^(b + 1L) / sg[small_i]^a
      out[small_i] <- gpd_poly_cpp(u_i, ck, threads) * pref
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
#' Closed form at both orders, from [gpd_components()]: the
#' log-density splits into \eqn{-\log\sigma}, \eqn{-\log t} and
#' \eqn{-\log(t)/\xi}, and the last is taken from its series where the
#' Leibniz form's terms cancel.
#' @param distrib A `GPDDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `sigma` and `xi`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx The approximation used when `expected` is `TRUE`.
#' @param nsim Monte Carlo draws when `approx = "mc"`.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [gpd_distrib()]
S7::method(distrib_deriv3, GPDDistrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...,
                                                    threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  gpd_components(y, theta, 3L, threads = threads)
}

#' @rdname distrib_deriv3.GPDDistrib
#' @name distrib_deriv4.GPDDistrib
#' @return A named list of fourth-derivative components.
S7::method(distrib_deriv4, GPDDistrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...,
                                                    threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  gpd_components(y, theta, 4L, threads = threads)
}
