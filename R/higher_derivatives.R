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
#' @param skip Character vector of component names, or `NULL`, the default.
#'   A named component is left `NULL` in the result rather than computed, for
#'   a caller that supplies it in closed form. The names and their order are
#'   unchanged, so nothing downstream has to know which were skipped.
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
numerical_deriv3 <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3), skip = NULL) {
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
#' @param skip Character vector of component names, or `NULL`, the default.
#'   A named component is left `NULL` in the result rather than computed, for
#'   a caller that supplies it in closed form. The names and their order are
#'   unchanged, so nothing downstream has to know which were skipped.
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
numerical_deriv4 <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 4), skip = NULL) {
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
                                     need_value = FALSE)$h
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


#' Fifth-Order Derivatives by One Central Difference
#'
#' @description
#' The fifth-order derivative components, obtained by differentiating the
#' fourth-order ones once in each parameter.
#'
#' @details
#' One central difference of the analytic fourth order, and never a difference
#' of a difference. That is the rule the rest of this package's numerical
#' surface obeys and it is worth restating here, because the sibling
#' [numerical_deriv4()] does **not** obey it: it differences the analytic
#' Hessian twice, which it can afford because a second difference of an exact
#' quantity is still only two orders removed from one. A fifth order built the
#' same way -- differencing the analytic third three times -- would not be.
#'
#' The component of a non-decreasing five-tuple is the four-tuple that remains
#' when its last index is dropped, differentiated in the parameter that index
#' names. The tuple is non-decreasing, so dropping the last index leaves a
#' valid order-4 key, and mixed partials commute, so which index is dropped
#' does not matter. Components are grouped by that last index, so the whole
#' order costs `2p` evaluations of the fourth rather than two per component.
#'
#' # The stencil, and the step
#'
#' The nodes, the weights and the step are all \pkg{numericals7}'s, so the
#' rule lives in one place and raising `accuracy` is an argument rather than
#' new arithmetic. [numericals7::fd_derivative()] is not called directly for
#' the reason [numDeriv_grad()] already records: its `f` maps a vector of
#' points to the values at those points, while this one reads a whole named
#' list at each node and keeps every component that shares the differentiated
#' index. A node of zero weight is skipped, so the default central rule costs
#' two evaluations of the fourth order per parameter and not three.
#'
#' [numericals7::fd_step()] at order 1 is \eqn{\epsilon^{1/3}} scaled by the
#' magnitude of the evaluation point and shrunk to keep the stencil inside the
#' parameter's domain. It is the right rule here because the quantity being
#' differenced is evaluated to machine precision. The same rule would be a
#' thousand times too small in statmodels7's outer stencil, where what limits
#' the step is the reproducibility of refitting the mode rather than the
#' rounding of an arithmetic expression. The two cases look alike and take
#' opposite answers.
#'
#' On the link scale the perturbation is applied to \eqn{\eta} and the
#' fourth order is asked for on the link scale as well, so the result is the
#' derivative of an analytic link-scale quantity. Nothing here needs the fifth
#' derivative of a link, nor an order-5 entry in [bell_partial()], both of
#' which the chain rule route would have required.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale `"parameter"` or `"link"`; the scale the fourth order is read
#'   on and the scale the perturbation is applied on.
#' @param accuracy The order of accuracy of the central rule, passed to
#'   [numericals7::fd_offsets()] and [numericals7::fd_step()]. The default 2
#'   is the three-point rule; 4 costs four evaluations of the fourth order per
#'   parameter instead of two.
#'
#' @return A named list of fifth-derivative component vectors, each of length
#'   `length(y)`, keyed lexicographically as
#'   [`deriv_names(distrib@params, 5)`][deriv_names] gives them.
#'
#' @examples
#' numerical_deriv5(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @seealso [distrib_deriv5()], the generic; [numerical_deriv4()], the order
#'   below and the one shape not to copy.
#' @export
numerical_deriv5 <- function(distrib, y, theta, scale = c("parameter", "link"),
                             accuracy = 2L) {
  scale <- match.arg(scale)
  params <- distrib@params
  p <- length(params)
  nms <- deriv_names(params, 5)
  idx_of <- deriv_indices(params, 5)

  # The stencil is numericals7's, nodes and weights together. A zero weight is
  # dropped rather than evaluated: at the default accuracy the centre node
  # carries one, so the rule costs two calls per parameter.
  s <- numericals7::fd_offsets(1L, accuracy = accuracy)$central
  w <- numericals7::fd_weights(s, 1L)
  nodes <- which(w != 0)

  # The point the difference is taken at, and the map back to a theta the
  # fourth-order method can be called with. On the link scale the two differ;
  # on the parameter scale the map is the identity.
  if (scale == "link") {
    links <- distrib@link_params
    x0 <- lapply(params, function(nm)
      linkfunctions7::linkfun(links[[nm]], theta[[nm]]))
    bnds <- lapply(params, function(nm)
      linkfunctions7::eta_bounds(links[[nm]]))
    to_theta <- function(x) {
      th <- theta
      for (i in seq_len(p)) {
        th[[params[i]]] <- linkfunctions7::linkinv(links[[params[i]]], x[[i]])
      }
      th
    }
  } else {
    x0 <- lapply(params, function(nm) theta[[nm]])
    bnds <- lapply(params, function(nm) distrib@params_bounds[[nm]])
    to_theta <- function(x) {
      th <- theta
      for (i in seq_len(p)) th[[params[i]]] <- x[[i]]
      th
    }
  }

  D4 <- function(th) distrib_deriv4(distrib, y, th, scale = scale)

  out <- vector("list", length(nms))
  names(out) <- nms

  # The component of a non-decreasing five-tuple is the four-tuple left when
  # its last index is dropped, differentiated in the parameter that index
  # names. Grouping by that index is what makes the whole order cost one
  # stencil per parameter rather than one per component.
  last <- vapply(idx_of, function(r) r[5L], integer(1))
  key4 <- vapply(idx_of, function(r) paste(params[r[1:4]], collapse = "_"),
                 character(1))

  for (m in seq_len(p)) {
    who <- which(last == m)
    if (!length(who)) next

    h <- numericals7::fd_step(x0[[m]], 1L, accuracy = accuracy,
                              bounds = bnds[[m]])

    acc <- vector("list", length(who))
    for (j in nodes) {
      x <- x0
      x[[m]] <- x0[[m]] + s[j] * h
      D <- D4(to_theta(x))
      for (t in seq_along(who)) {
        term <- w[j] * D[[key4[who[t]]]]
        acc[[t]] <- if (is.null(acc[[t]])) term else acc[[t]] + term
      }
    }

    for (t in seq_along(who)) out[[nms[who[t]]]] <- acc[[t]] / h
  }

  out
}


#' @title Default Fifth-Order Derivatives for `distrib` Objects
#' @name distrib_deriv5.distrib
#'
#' @description
#' The route every family takes at the fifth order: one central difference of
#' the fourth, through [numerical_deriv5()].
#'
#' @details
#' This is registered on the base class and nothing overrides it, so it is not
#' a fallback in the usual sense -- no family writes the fifth order out yet.
#' What differs between families is the quantity being differenced.
#' [has_exact_deriv4()] answers whether that quantity is analytic, and
#' [check_distrib()] reports the order-5 row as unchecked where it is not.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale Passed through to [numerical_deriv5()], which reads the fourth
#'   order on that scale and perturbs on it.
#' @param ... Unused.
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
