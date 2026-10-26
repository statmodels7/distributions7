#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R numerical_functions.R
NULL

# Higher-order (3rd/4th) derivatives of the log-density.
#
# The generics distrib_deriv3 / distrib_deriv4 dispatch to closed-form C++ kernels
# for the distributions that have them; every other distribution is served by
# the finite-difference fallbacks below. A family with its own Hessian has that
# Hessian differenced; a family without one has its score or its log-density
# differenced by one tensor-product stencil (tensor_derivatives()), because
# differencing a Hessian that is itself a difference nests one difference in
# another.

#' Numerical Third-Order Derivatives of the Log-Density
#'
#' @description
#' Computes the unique third-order partial derivatives of the log-density by
#' central finite differences of the highest order that the family implements
#' itself. This powers the default [distrib_deriv3()] method for distributions
#' without a closed-form implementation, and is the reference used to validate
#' the analytical kernels.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or `length(y)`).
#' @param h_rel Numeric; the relative finite-difference step. Defaults to
#'   `.Machine$double.eps^(1/3)` for a family with its own Hessian, and to the
#'   step that [tensor_derivatives()] chooses for the order of its difference
#'   otherwise; a value given explicitly is used on either route.
#' @param skip Character vector of component names, or `NULL`, the default.
#'   A named component is left `NULL` in the result rather than computed, for
#'   a caller that supplies it in closed form. The names and their order are
#'   unchanged, so a downstream caller can read the result without checking
#'   which components were skipped.
#'
#' @return A named list of third-derivative component vectors, keyed as in
#'   [`deriv_names(distrib@params, 3)`][deriv_names].
#'
#' @details
#' The route depends on what the family implements, as [analytic_order()]
#' reports it.
#'
#' For a family with its own Hessian, each component
#' \eqn{\partial^3 \ell / \partial\theta_i\partial\theta_j\partial\theta_k}
#' (with \eqn{i \le j \le k}) is obtained by differentiating the Hessian entry
#' \eqn{(i, j)} along \eqn{\theta_k}, one central difference of an analytical
#' quantity. Steps are scaled by `max(1, |theta|)` and shrunk near
#' parameter-domain boundaries.
#'
#' For a family without its own Hessian, that Hessian is itself a difference,
#' and differencing it would nest one difference in another: on a diagonal
#' component the same parameter would be differenced twice in succession. Each
#' component is then one tensor-product central stencil of order two on the
#' score, or of order three on the log-density when the family has no score
#' either; see [tensor_derivatives()].
#'
#' @seealso [numerical_deriv4()], [distrib_deriv3()]
#' @examples
#' numerical_deriv3(gaussian1_distrib(), 0, list(mu = 0, sigma = 1))
#'
#' @export
numerical_deriv3 <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3), skip = NULL) {
  base <- analytic_order(distrib, 2L)
  if (base < 2L) {
    return(tensor_derivatives(distrib, y, theta, order = 3L, base = base,
                              h_rel = if (missing(h_rel)) NULL else h_rel,
                              skip = skip))
  }
  params <- distrib@params
  bounds <- distrib@params_bounds
  nms <- deriv_names(params, 3)
  # Taken from the same enumeration that produced the names, never recovered by
  # splitting them: see deriv_indices().
  idx_of <- deriv_indices(params, 3)

  out <- vector("list", length(nms))
  names(out) <- nms

  # The step of each parameter is chosen once, on the whole Hessian differenced
  # along it, and then used by every component that differences along that
  # parameter: the choice costs four evaluations per parameter rather than per
  # component, and it does not depend on which component happens to come first.
  hsteps <- vector("list", length(params))
  step_for <- function(k) {
    if (is.null(hsteps[[k]])) {
      quotient <- function(hk) {
        tp <- tm <- theta
        tp[[k]] <- theta[[k]] + hk
        tm[[k]] <- theta[[k]] - hk
        mapply(function(a, b) (a - b) / (2 * hk),
               distrib_hessian(distrib, y, tp),
               distrib_hessian(distrib, y, tm), SIMPLIFY = FALSE)
      }
      hsteps[[k]] <<- fd_stable_step(quotient, theta[[k]],
                                     bounds[[params[k]]], h_rel,
                                     need_value = FALSE)$h
    }
    hsteps[[k]]
  }

  for (t in seq_along(nms)) {
    nm <- nms[t]
    if (nm %in% skip) next
    idx <- idx_of[[t]]
    i <- idx[1]; j <- idx[2]; k <- idx[3]
    hk <- step_for(k)
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
#' Computes the unique fourth-order partial derivatives of the log-density by
#' central finite differences of the highest order that the family implements
#' itself, up to the Hessian. This powers the default [distrib_deriv4()]
#' method for distributions without a closed-form implementation.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or `length(y)`).
#' @param h_rel Numeric; the relative finite-difference step. Defaults to
#'   `.Machine$double.eps^(1/4)` for a family with its own Hessian, and to the
#'   step that [tensor_derivatives()] chooses for the order of its difference
#'   otherwise; a value given explicitly is used on either route.
#' @param skip Character vector of component names, or `NULL`, the default.
#'   A named component is left `NULL` in the result rather than computed, for
#'   a caller that supplies it in closed form. The names and their order are
#'   unchanged, so a downstream caller can read the result without checking
#'   which components were skipped.
#'
#' @return A named list of fourth-derivative component vectors, keyed as in
#'   [`deriv_names(distrib@params, 4)`][deriv_names].
#'
#' @details
#' For a family with its own Hessian, each component
#' \eqn{\partial^4 \ell / \partial\theta_i\partial\theta_j\partial\theta_k\partial\theta_l}
#' (with \eqn{i \le j \le k \le l}) is obtained as the second derivative of the
#' Hessian entry \eqn{(i, j)} along \eqn{(\theta_k, \theta_l)}: a three-point stencil
#' when \eqn{k = l}, a four-point cross stencil otherwise.
#'
#' For a family without its own Hessian, each component is one tensor-product
#' central stencil of order three on the score, or of order four on the
#' log-density when the family has no score either; see
#' [tensor_derivatives()] and the same paragraph of [numerical_deriv3()].
#'
#' @seealso [numerical_deriv3()], [distrib_deriv4()]
#' @examples
#' numerical_deriv4(gaussian1_distrib(), 0, list(mu = 0, sigma = 1))
#'
#' @export
numerical_deriv4 <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 4), skip = NULL) {
  base <- analytic_order(distrib, 2L)
  if (base < 2L) {
    return(tensor_derivatives(distrib, y, theta, order = 4L, base = base,
                              h_rel = if (missing(h_rel)) NULL else h_rel,
                              skip = skip))
  }
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

  # As in numerical_deriv3(), the step of each parameter is chosen once, here on
  # the second difference of the whole Hessian along that parameter, and reused
  # by every component that steps along it.
  hsteps <- vector("list", length(params))
  step_for <- function(k) {
    if (is.null(hsteps[[k]])) {
      quotient <- function(hk) {
        tp <- tm <- theta
        tp[[k]] <- theta[[k]] + hk
        tm[[k]] <- theta[[k]] - hk
        Hp <- H(tp); Hm <- H(tm)
        stats::setNames(lapply(names(H0), function(nm)
          (Hp[[nm]] - 2 * H0[[nm]] + Hm[[nm]]) / (hk^2)), names(H0))
      }
      hsteps[[k]] <<- fd_stable_step(quotient, theta[[k]],
                                     bounds[[params[k]]], h_rel,
                                     need_value = FALSE, order = 2L)$h
    }
    hsteps[[k]]
  }

  for (t in seq_along(nms)) {
    nm <- nms[t]
    if (nm %in% skip) next
    idx <- idx_of[[t]]
    i <- idx[1]; j <- idx[2]; k <- idx[3]; l <- idx[4]
    hcomp <- paste(params[c(i, j)], collapse = "_")
    hk <- step_for(k)
    hl <- step_for(l)

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
#' Provides the fallback for a family that registers no third-order method.
#' Observed derivatives come from [numerical_deriv3()]: one central difference
#' of the family's own [distrib_hessian()] along each parameter, or, for a
#' family without its own Hessian, one tensor-product stencil on its score or
#' its log-density. Expected ones come from [expected_derivative()] at the
#' approximation selected by `approx`.
#'
#' Every family shipped in this package writes the observed third order out,
#' so this method serves a family defined outside the package, and five orders
#' are available for such a family from its density without any other method.
#' It is also the reference against which the
#' analytical kernels are validated.
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
#' @param approx An object inheriting from [ExpectApprox], [expect_integrate()]
#'   by default, that gives the approximation of the expectation. It is read
#'   only when `expected` is `TRUE`. The Monte Carlo sample size is the `nsim`
#'   of [expect_mc()]. See [expected_derivative_methods].
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
S7::method(distrib_deriv3, distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = expect_integrate(), ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = approx)
  } else {
    numerical_deriv3(distrib, y, theta)
  }
}

#' @title Default Fourth-Order Derivatives for `distrib` Objects
#' @name distrib_deriv4.distrib
#'
#' @description
#' Provides the fallback for a family that registers no fourth-order method.
#' Observed
#' derivatives come from [numerical_deriv4()]: a second difference of the
#' family's own [distrib_hessian()] rather than a difference of the third
#' order, or, for a family without its own Hessian, one tensor-product stencil
#' on its score or its log-density. On both routes each component is a single
#' stencil, so the package's rule against nesting one difference inside
#' another holds. Expected ones come from [expected_derivative()] at the
#' approximation that `approx` names.
#'
#' As at the order below, the families shipped in this package write the
#' fourth order out themselves. A second difference of the Hessian is less
#' accurate than the first difference of the third order, and its step is
#' chosen accordingly.
#'
#' @param distrib An object inheriting from `distrib` that registers no method
#'   of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param expected Logical of length 1. `FALSE`, the default, differences the
#'   Hessian twice; `TRUE` takes the expectation of the observed derivatives.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param approx An object inheriting from [ExpectApprox], [expect_integrate()]
#'   by default, that gives the approximation of the expectation. It is read
#'   only when `expected` is `TRUE`. The Monte Carlo sample size is the `nsim`
#'   of [expect_mc()]. See [expected_derivative_methods].
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
S7::method(distrib_deriv4, distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = expect_integrate(), ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = approx)
  } else {
    numerical_deriv4(distrib, y, theta)
  }
}


#' Fifth-Order Derivatives by One Tensor Stencil
#'
#' @description
#' Computes the fifth-order derivative components, where each component is a
#' single tensor-product central difference of the highest analytic order
#' below it, through [tensor_derivatives()].
#'
#' @details
#' The quantity differenced is the fourth order where it is analytic, as
#' [has_exact_deriv4()] reads it, and each component is then one central
#' difference of a fourth-order component in one parameter: the component of
#' a non-decreasing five-tuple is the four-tuple left when one index is
#' dropped, and the whole order costs one stencil per parameter. Where the
#' fourth order is not analytic, a family that registers no wrapper reads
#' the highest order that it implements itself, as [analytic_order()] gives
#' it, and a wrapper of a parent without an analytic fourth order differences
#' its log-density, because its lower orders are partition sums over the
#' parent's fallbacks. A difference of a difference is never formed.
#'
#' The step of a difference of order \eqn{r} is
#' \eqn{\varepsilon^{1/(r+a)}\max(1, |\theta_j|)}, with \eqn{a} the
#' accuracy, reduced so that every node stays inside the parameter's domain.
#' On the link scale the perturbation is applied to \eqn{\eta} and the
#' differenced quantity is read on the link scale as well, so the result is
#' the derivative of a link-scale quantity. This route requires neither the
#' fifth derivative of a link nor an order-5 entry in [bell_partial()], both
#' of which the chain rule route uses.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale `"parameter"` or `"link"`; the scale on which the
#'   differenced quantity is read and on which the perturbation is applied.
#' @param accuracy The order of accuracy of the central rules, passed to
#'   [numericals7::fd_tensor()]. The default 4 uses four evaluations of the
#'   fourth order per parameter where it is analytic; 2 uses two.
#'
#' @return A named list of fifth-derivative component vectors, each of length
#'   `length(y)`, keyed lexicographically as
#'   [`deriv_names(distrib@params, 5)`][deriv_names] gives them.
#'
#' @examples
#' numerical_deriv5(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @seealso [distrib_deriv5()], the generic; [numerical_deriv4()], the order
#'   below; [tensor_derivatives()], which does the differencing.
#' @export
numerical_deriv5 <- function(distrib, y, theta, scale = c("parameter", "link"),
                             accuracy = 4L) {
  scale <- match.arg(scale)
  tensor_derivatives(distrib, y, theta, order = 5L, base = deriv5_base(distrib),
                     scale = scale, accuracy = accuracy)
}


#' The Order from Which a Fifth Derivative Is Differenced
#'
#' @description
#' Returns 4 where [has_exact_deriv4()] is `TRUE`; 0 for a wrapper whose
#' parent has no analytic fourth order, where a wrapper is a class that
#' overrides [has_exact_deriv4()] to read its parent; and otherwise the
#' highest order that the family registers itself, as [analytic_order()] gives
#' it up to 3.
#'
#' @param x A distribution object.
#'
#' @return An integer from 0 to 4.
#'
#' @seealso [numerical_deriv5()], which passes it to [tensor_derivatives()].
#' @keywords internal
deriv5_base <- function(x) {
  if (has_exact_deriv4(x)) return(4L)
  m <- S7::method(has_exact_deriv4, S7::S7_class(x))
  if (!is_class(attr(m, "signature")[[1]], distrib)) return(0L)
  analytic_order(x, 3L)
}


#' @title Default Fifth-Order Derivatives for `distrib` Objects
#' @name distrib_deriv5.distrib
#'
#' @description
#' Defines the route that a family takes at the fifth order when it registers
#' no method of its own: one tensor-product difference of its highest analytic
#' order, through [numerical_deriv5()].
#'
#' @details
#' This is registered on the base class. The families that write the fifth
#' order out themselves (see [distrib_deriv5()]) and the wrappers override it.
#' Among the families that reach it, the quantity that is differenced varies:
#' [has_exact_deriv4()] reports whether that quantity is analytic, and
#' [check_distrib()] emits no order-5 row where it is not.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale Passed through to [numerical_deriv5()], which reads the fourth
#'   order on that scale and perturbs on it.
#' @param ... Passed to [numerical_deriv5()], which accepts `accuracy`.
#'
#' @return A named list of fifth-derivative component vectors, each of length
#'   `length(y)`, keyed lexicographically as
#'   [`deriv_names(distrib@params, 5)`][deriv_names] gives them.
#'
#' @seealso [numerical_deriv5()], which does the differencing;
#'   [distrib_deriv4.distrib()] for the order below.
#' @keywords internal
S7::method(distrib_deriv5, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  numerical_deriv5(distrib, y, theta, scale = scale, ...)
}
