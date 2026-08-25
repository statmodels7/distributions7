#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R numerical_functions.R
NULL

# Higher-order (3rd/4th) derivatives of the log-density.
#
# The generics distrib_deriv3 / distrib_deriv4 dispatch to closed-form C++ kernels
# for the distributions that have them; every other distribution (and every
# wrapper / transformed / user-defined one) is served by the finite-difference
# fallbacks below, which differentiate the Hessian returned by distrib_hessian
# (analytical when available, itself a finite-difference fallback otherwise).

#' Numerical Third-Order Derivatives of the Log-Density
#'
#' @description
#' Computes the unique third-order partial derivatives of the log-density by central
#' finite differences of [distrib_hessian()]. This powers the default
#' [distrib_deriv3()] method for distributions without a closed-form
#' implementation, and is the reference used to validate the analytical kernels.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or `length(y)`).
#' @param h_rel Numeric. Relative finite-difference step. Defaults to
#'   `.Machine$double.eps^(1/3)`.
#'
#' @return A named list of third-derivative component vectors, keyed as in
#'   [`deriv_names(distrib@params, 3)`][deriv_names].
#'
#' @details
#' Each component \eqn{\partial^3 \ell / \partial\theta_i\partial\theta_j\partial\theta_k}
#' (with \eqn{i \le j \le k}) is obtained by differentiating the Hessian entry
#' \eqn{(i, j)} along \eqn{\theta_k}. Steps are scaled by `max(1, |theta|)` and
#' shrunk near parameter-domain boundaries.
#'
#' @seealso [numerical_deriv4()], [distrib_deriv3()]
#' @examples
#' numerical_deriv3(gaussian1_distrib(), 0, list(mu = 0, sigma = 1))
#'
#' @export
numerical_deriv3 <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3)) {
  params <- distrib@params
  bounds <- distrib@params_bounds
  nms <- deriv_names(params, 3)
  # Taken from the same enumeration that produced the names, never recovered by
  # splitting them: see deriv_indices().
  idx_of <- deriv_indices(params, 3)

  out <- vector("list", length(nms))
  names(out) <- nms

  for (t in seq_along(nms)) {
    nm <- nms[t]
    idx <- idx_of[[t]]
    i <- idx[1]; j <- idx[2]; k <- idx[3]
    hk <- fd_steps(theta[[k]], bounds[[params[k]]], h_rel)
    hcomp <- paste(params[c(i, j)], collapse = "_")

    tp <- tm <- theta
    tp[[k]] <- theta[[k]] + hk
    tm[[k]] <- theta[[k]] - hk

    out[[nm]] <- (distrib_hessian(distrib, y, tp)[[hcomp]] -
      distrib_hessian(distrib, y, tm)[[hcomp]]) / (2 * hk)
  }

  out
}

#' Numerical Fourth-Order Derivatives of the Log-Density
#'
#' @description
#' Computes the unique fourth-order partial derivatives of the log-density by second
#' central differences of [distrib_hessian()]. This powers the default
#' [distrib_deriv4()] method for distributions without a closed-form
#' implementation.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or `length(y)`).
#' @param h_rel Numeric. Relative finite-difference step. Defaults to
#'   `.Machine$double.eps^(1/4)`.
#'
#' @return A named list of fourth-derivative component vectors, keyed as in
#'   [`deriv_names(distrib@params, 4)`][deriv_names].
#'
#' @details
#' Each component
#' \eqn{\partial^4 \ell / \partial\theta_i\partial\theta_j\partial\theta_k\partial\theta_l}
#' (with \eqn{i \le j \le k \le l}) is obtained as the second derivative of the
#' Hessian entry \eqn{(i, j)} along \eqn{(\theta_k, \theta_l)}: a three-point stencil
#' when \eqn{k = l}, a four-point cross stencil otherwise.
#'
#' @seealso [numerical_deriv3()], [distrib_deriv4()]
#' @examples
#' numerical_deriv4(gaussian1_distrib(), 0, list(mu = 0, sigma = 1))
#'
#' @export
numerical_deriv4 <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 4)) {
  params <- distrib@params
  bounds <- distrib@params_bounds
  nms <- deriv_names(params, 4)
  idx_of <- deriv_indices(params, 4)

  out <- vector("list", length(nms))
  names(out) <- nms

  H <- function(th) distrib_hessian(distrib, y, th)
  # The center point of the three-point stencil does not move, so it is computed
  # once rather than once per component with k == l.
  H0 <- H(theta)

  for (t in seq_along(nms)) {
    nm <- nms[t]
    idx <- idx_of[[t]]
    i <- idx[1]; j <- idx[2]; k <- idx[3]; l <- idx[4]
    hcomp <- paste(params[c(i, j)], collapse = "_")
    hk <- fd_steps(theta[[k]], bounds[[params[k]]], h_rel)
    hl <- fd_steps(theta[[l]], bounds[[params[l]]], h_rel)

    if (k == l) {
      tp <- tm <- theta
      tp[[k]] <- theta[[k]] + hk
      tm[[k]] <- theta[[k]] - hk
      out[[nm]] <- (H(tp)[[hcomp]] - 2 * H0[[hcomp]] + H(tm)[[hcomp]]) / (hk^2)
    } else {
      shift <- function(a, b) {
        th <- theta
        th[[k]] <- theta[[k]] + a * hk
        th[[l]] <- theta[[l]] + b * hl
        th
      }
      out[[nm]] <- (H(shift(1, 1))[[hcomp]] - H(shift(1, -1))[[hcomp]] -
        H(shift(-1, 1))[[hcomp]] + H(shift(-1, -1))[[hcomp]]) / (4 * hk * hl)
    }
  }

  out
}

# --- DEFAULT (FALLBACK) METHODS ---

#' @title Default Third-Order Derivatives for `distrib` Objects
#' @name distrib_deriv3.distrib
#'
#' @description
#' The fallback for a family that registers no third-order method. Observed
#' derivatives come from [numerical_deriv3()], one central difference of
#' [distrib_hessian()] along each parameter; expected ones from
#' [expected_derivative()] at the strategy `approx` names.
#'
#' **No family shipped in this package reaches this method for its observed
#' derivatives.** All 46 write the third order out, 24 of them in compiled
#' kernels. It exists for a family defined outside the package, which gets four
#' orders from a density alone, and it is the reference the analytical kernels
#' are validated against.
#'
#' @param distrib An object inheriting from `distrib` that registers no method
#'   of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param expected Logical of length 1. `FALSE`, the default, differences the
#'   Hessian; `TRUE` takes the expectation of the observed derivatives through
#'   [expected_derivative()].
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param approx Which strategy takes the expectation, read only when
#'   `expected` is `TRUE`: `"integrate"` (the default at this order),
#'   `"bartlett"`, `"opg"` or `"mc"`. See [expected_derivative_methods()].
#' @param nsim Monte Carlo sample size, a single positive number, read only
#'   under `approx = "mc"`. Defaults to 10000.
#' @param ... Unused.
#'
#' @return A named list of third-derivative component vectors, each of length
#'   `length(y)`, keyed lexicographically as
#'   [`deriv_names(distrib@params, 3)`][deriv_names] gives them.
#'
#' @seealso [numerical_deriv3()], which does the differencing;
#'   [distrib_deriv4.distrib()] for the order above;
#'   [expected_derivative_methods()] for what `approx` selects.
#' @keywords internal
S7::method(distrib_deriv3, distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    numerical_deriv3(distrib, y, theta)
  }
}

#' @title Default Fourth-Order Derivatives for `distrib` Objects
#' @name distrib_deriv4.distrib
#'
#' @description
#' The fallback for a family that registers no fourth-order method. Observed
#' derivatives come from [numerical_deriv4()], a **second** difference of
#' [distrib_hessian()] rather than a difference of the third order, so the
#' package's rule against nesting one difference inside another holds here;
#' expected ones from [expected_derivative()] at the strategy `approx` names.
#'
#' As at the order below, no family shipped in this package reaches it. Being a
#' second difference it is the least accurate route the package offers, and its
#' step is chosen accordingly.
#'
#' @param distrib An object inheriting from `distrib` that registers no method
#'   of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param expected Logical of length 1. `FALSE`, the default, differences the
#'   Hessian twice; `TRUE` takes the expectation of the observed derivatives.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param approx Which strategy takes the expectation, read only when
#'   `expected` is `TRUE`: `"integrate"` (the default at this order),
#'   `"bartlett"`, `"opg"` or `"mc"`. See [expected_derivative_methods()].
#' @param nsim Monte Carlo sample size, a single positive number, read only
#'   under `approx = "mc"`. Defaults to 10000.
#' @param ... Unused.
#'
#' @return A named list of fourth-derivative component vectors, each of length
#'   `length(y)`, keyed lexicographically as
#'   [`deriv_names(distrib@params, 4)`][deriv_names] gives them.
#'
#' @seealso [numerical_deriv4()], which does the differencing;
#'   [distrib_deriv3.distrib()] for the order below;
#'   [expected_derivative_methods()] for what `approx` selects.
#' @keywords internal
S7::method(distrib_deriv4, distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    numerical_deriv4(distrib, y, theta)
  }
}
