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
#' Returns the components of
#' \eqn{\partial^{a+b}\ell/\partial\sigma^a\partial\xi^b} at any order from one
#' to four, by splitting the log-density into a part that is analytic at
#' \eqn{\xi = 0} and a part that carries a removable singularity there.
#'
#' @details
#' # The split
#'
#' With \eqn{z = y/\sigma} and \eqn{t = 1 + \xi z},
#' \deqn{\ell = -\log\sigma - L - W, \qquad L = \log t, \qquad W = L/\xi.}
#' \eqn{t} is affine in \eqn{\xi}, so every partial of it carrying two or more
#' \eqn{\xi} vanishes and the written-out template of [fdb2()] covers \eqn{L}
#' outright.
#'
#' # Why W has two routes
#'
#' \eqn{W = L/\xi} has a removable singularity at \eqn{\xi = 0}, where it tends
#' to \eqn{z}. Differentiating it by Leibniz,
#' \deqn{\partial^{a,b} W = \sum_{j} \binom{b}{j}\, \partial^{a,j}L \cdot
#'       \partial^{b-j}(1/\xi),}
#' produces terms of size \eqn{\xi^{-(b+1)}} that cancel against each other, so
#' near zero the form loses every digit. Below a threshold it is replaced by the
#' series the cancellation leaves,
#' \deqn{W = \sum_{k \ge 0} \frac{(-1)^k \xi^k z^{k+1}}{k+1},}
#' differentiated term by term.
#'
#' The threshold is on \eqn{\lvert\xi z\rvert} and **not** on \eqn{\xi}. At
#' order \eqn{b} the Leibniz terms are of size \eqn{z\xi^{-b}} against an answer
#' of size \eqn{z^{b+1}}, so the relative cancellation is
#' \eqn{(\xi z)^{-b}} and the accuracy is \eqn{\varepsilon(\xi z)^{-b}}, which
#' at order four reaches 0.3 by \eqn{\xi z = 1.7\times10^{-4}}. That prediction
#' was measured, and matching it is what fixed the switch. The series needs
#' \eqn{\lvert\xi z\rvert < 1} to converge and is cut at 0.2, where forty terms
#' leave \eqn{10^{-28}} while the Leibniz form still carries \eqn{10^{-13}}, so
#' the two overlap.
#'
#' # What the series branch evaluates
#'
#' The two elementwise powers the sum appears to need are algebraically one:
#' with \eqn{u = \xi z}, \eqn{\xi^{k-b}z^{k+1} = u^{k-b}z^{b+1}}, so
#' \eqn{z^{b+1}/\sigma^a} leaves the sum and what remains is a **polynomial in**
#' \eqn{u} whose coefficients are scalar in \eqn{k}. `gpd_poly_cpp()` evaluates
#' it by Horner from the highest power down, which sums a decaying series
#' smallest first. Each branch also runs on its own elements only, so a sample
#' straddling the cut does not pay for both in full.
#'
#' @param y A numeric vector of observations on the support.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`. `sigma` must be strictly
#'   positive; `xi` may be of either sign, including zero, which is the
#'   exponential limit. Shorter components are recycled.
#' @param order The derivative order, an integer from 1 to 4.
#' @param cut The value of \eqn{\lvert\xi z\rvert} below which the series is
#'   used, defaulting to `0.2`. It is exposed so that a test can force either
#'   route in the region where both are accurate and compare them; a caller has
#'   no reason to change it.
#' @param threads A single positive integer, how many threads the polynomial
#'   kernel may use. Defaults to `1L`.
#'
#' @return A named list of component vectors, one per distinct multi-index of
#'   the given order and keyed as [deriv_names()] keys them: two at order 1,
#'   three at order 2, four at order 3 and five at order 4. Each has the
#'   recycled length of the inputs.
#'
#' @seealso [distrib_deriv3.GPDDistrib()] and [distrib_deriv4.GPDDistrib()],
#'   which call this; [fdb2()] for the two-variable composition template; and
#'   [gpd_distrib()] for the family.
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


#' @title Generalized Pareto Third-Order Derivatives
#' @name distrib_deriv3.GPDDistrib
#' @description
#' Computes the four distinct third derivatives of the generalized Pareto
#' log-density in \eqn{\sigma} and \eqn{\xi}, **in closed form**, through
#' [gpd_components()]. The log-density splits as
#' \eqn{-\log\sigma - \log t - \log(t)/\xi} with \eqn{t = 1 + \xi y/\sigma}, and
#' the last piece is taken from its own series wherever the Leibniz form's
#' terms of size \eqn{\xi^{-(b+1)}} cancel against each other.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead: the
#' expected third derivatives have no closed form. That is the one place on
#' this page where `approx` and `nsim` are read.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector of observations on the support. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`. `sigma` must be strictly
#'   positive; `xi` may be of either sign, including zero.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the polynomial
#'   kernel of the series branch may use. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `sigma_sigma_sigma`,
#'   `sigma_sigma_xi`, `sigma_xi_xi` and `xi_xi_xi`, each of length
#'   `max(length(y), lengths(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\sigma > 0} the
#' scale, \eqn{\xi} the shape, \eqn{z = y/\sigma} and \eqn{t = 1 + \xi z}.
#'
#' @seealso [distrib_hessian.GPDDistrib()] for the order below,
#'   [distrib_deriv4.GPDDistrib()] for the order above,
#'   [gpd_components()] for the two routes and the measured threshold, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' y <- c(0.5, 2, 8)
#' th <- list(sigma = 1.5, xi = 0.3)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-scale component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(sigma = 1.5 + eps, xi = 0.3))$sigma_sigma
#' dn <- distrib_hessian(d, y, list(sigma = 1.5 - eps, xi = 0.3))$sigma_sigma
#' all.equal((up - dn) / (2 * eps), d3$sigma_sigma_sigma, tolerance = 1e-6)
#'
#' # At a shape near zero the family is the exponential, and the scale
#' # components converge onto that family's at rate O(xi).
#' de <- exponential_distrib()
#' vapply(c(1e-6, 1e-9), function(x)
#'   max(abs(distrib_deriv3(d, y, list(sigma = 1.5, xi = x))$sigma_sigma_sigma -
#'           distrib_deriv3(de, y, list(mu = 1.5))$mu_mu_mu)), numeric(1))
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

#' @title Generalized Pareto Fourth-Order Derivatives
#' @name distrib_deriv4.GPDDistrib
#' @description
#' Computes the five distinct fourth derivatives of the generalized Pareto
#' log-density in \eqn{\sigma} and \eqn{\xi}, **in closed form**, by the
#' construction [distrib_deriv3.GPDDistrib()] describes carried one order
#' further. This is the order at which the choice between the two routes bites
#' hardest: the Leibniz form's relative cancellation is \eqn{(\xi z)^{-b}}, so
#' at \eqn{b = 4} it has lost a third of the answer by
#' \eqn{\xi z = 1.7\times10^{-4}} and everything below that.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead: the
#' expected fourth derivatives have no closed form.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param y A numeric vector of observations on the support. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `sigma` and `xi`, each a numeric
#'   vector of length 1 or of the length of `y`. `sigma` must be strictly
#'   positive; `xi` may be of either sign, including zero.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the polynomial
#'   kernel of the series branch may use. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors, `sigma_sigma_sigma_sigma`,
#'   `sigma_sigma_sigma_xi`, `sigma_sigma_xi_xi`, `sigma_xi_xi_xi` and
#'   `xi_xi_xi_xi`, each of length `max(length(y), lengths(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\sigma > 0} the
#' scale, \eqn{\xi} the shape, \eqn{z = y/\sigma} and \eqn{t = 1 + \xi z}.
#'
#' @seealso [distrib_deriv3.GPDDistrib()] for the order below and the
#'   construction, [gpd_components()] for the two routes and the measured
#'   threshold, and [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- gpd_distrib()
#' y <- c(0.5, 2, 8)
#' th <- list(sigma = 1.5, xi = 0.3)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The two routes agree wherever both are accurate: forcing the cut either
#' # way at xi * z = 0.05, 0.15 and 0.30.
#' gc <- distributions7:::gpd_components
#' gap <- function(xz) {
#'   xi <- xz / (2 / 1.5)
#'   a <- gc(2, list(sigma = 1.5, xi = xi), 4L, cut = 1e-9)
#'   b <- gc(2, list(sigma = 1.5, xi = xi), 4L, cut = 10)
#'   max(abs(unlist(a) - unlist(b)) / pmax(abs(unlist(a)), 1e-10))
#' }
#' vapply(c(0.05, 0.15, 0.30), gap, numeric(1))
#'
#' # And below the cut the Leibniz form loses the answer entirely, which is
#' # the reason the series branch exists.
#' leib <- function(x) gc(2, list(sigma = 1.5, xi = x), 4L, cut = 1e-12)$xi_xi_xi_xi
#' ser <- function(x) gc(2, list(sigma = 1.5, xi = x), 4L, cut = 10)$xi_xi_xi_xi
#' rbind(xi = c(1e-2, 1e-4, 1e-6),
#'       leibniz = vapply(c(1e-2, 1e-4, 1e-6), leib, numeric(1)),
#'       series = vapply(c(1e-2, 1e-4, 1e-6), ser, numeric(1)))
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
