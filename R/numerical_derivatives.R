#' @include distrib.R generics.R utility_functions.R numerical_functions.R
NULL

#' Finite-Difference Steps That Respect a Parameter's Domain
#'
#' @description
#' Builds the step \eqn{h} for a central difference in one parameter: scaled by
#' the parameter's magnitude and kept strictly inside the parameter's
#' mathematical domain, either by cutting it to under half the distance to the
#' nearest finite bound or by scaling it on that distance.
#'
#' @details
#' The domain clamp makes it possible to offer a finite-difference fallback at
#' all. Parameter domains here are open: a scale parameter is positive, not
#' non-negative. A step chosen only from the magnitude therefore takes a small
#' \eqn{\sigma} straight through zero, and the log-density comes back `NaN`
#' for reasons that look like a bug in the density. Clamping to 49\% of the
#' distance to each finite boundary keeps both evaluation points inside.
#'
#' `to_bound = "scale"` takes \eqn{h = h_{rel}\min(\max(1,|\theta|), d)}, with
#' \eqn{d} the distance to the nearest finite bound, so the step shrinks in
#' proportion as the parameter approaches the bound rather than sitting at
#' \eqn{0.49\,d}. Where the log-density is singular in that parameter at the
#' bound the clamped step's error stops falling once the clamp binds, and where
#' the log-density instead saturates the scaled step is the worse of the two,
#' because its shorter step lets the rounding dominate. Each of the two is wrong
#' somewhere, and [fd_stable_step()] chooses between them.
#'
#' The two give the same step wherever \eqn{d \ge \max(1, |\theta|)}, which is
#' every parameter of a family with no finite bound, and also a positive
#' parameter at or above one.
#'
#' A parameter already on or outside its boundary cannot be rescued this way, and
#' is reported rather than differentiated.
#'
#' The clamp applies to the stencil's OUTERMOST node and not to the step, so it
#' divides by the reach, the largest offset (in units of \eqn{h}) at which the
#' rule evaluates. That number belongs to the stencil and is read from
#' [numericals7::fd_offsets()] rather than assumed: a central difference and a
#' three-point second difference both reach one step out, which is every use
#' that this package makes of the clamped branch today, while the five-point
#' rules of accuracy four reach two. Without the division a step cut to
#' \eqn{0.49 d} puts the outermost node of a reach-two rule at \eqn{0.98 d},
#' which is inside the domain but lies farther out than the \eqn{0.49 d} that
#' the clamp is meant to enforce.
#'
#' When the reach is accounted for, the clamped branch at the default `rel_step` is
#' [numericals7::fd_step()] with `bounds`, and a test asserts the two are
#' `identical()` so they cannot drift. It is written here rather than
#' delegated because `rel_step` is an argument that the callers vary:
#' [fd_stable_step()] halves it and `check_distrib()` reads the quotient at two
#' steps, neither of which `fd_step()` can express.
#'
#' The scaled branch needs no such division. Its step is at most
#' \eqn{h_{rel} d}, so the outermost node sits at \eqn{d(1 - r\,h_{rel})} for
#' a reach \eqn{r}, inside the domain for any stencil at any root of machine
#' epsilon.
#'
#' @param theta_j A numeric vector, the values of one parameter.
#' @param bounds_j A length-2 numeric vector giving that parameter's domain, or
#'   `NULL`.
#' @param rel_step The relative step size, typically a root of machine epsilon
#'   chosen for the stencil in use.
#' @param to_bound `"clamp"`, the default, or `"scale"`, as above.
#' @param order,accuracy The stencil for which the step is built; they fix its
#'   reach through [numericals7::fd_offsets()]. The defaults are the central
#'   difference, whose reach is one.
#'
#' @return A numeric vector of steps, the same length as `theta_j`.
#'
#' @seealso [fd_stable_step()], which chooses between the two, and
#'   [numerical_gradient()] and [numerical_hessian()], which take the chosen
#'   step. [fd_steps_y()] is the response counterpart.
#'   [numericals7::fd_step()] is the same rule at the default `rel_step`.
#' @keywords internal
fd_steps <- function(theta_j, bounds_j, rel_step, to_bound = c("clamp", "scale"),
                     order = 1L, accuracy = 2L) {
  to_bound <- match.arg(to_bound)
  if (identical(to_bound, "scale")) {
    s <- pmax(1, abs(theta_j))
    if (!is.null(bounds_j)) {
      if (is.finite(bounds_j[1])) s <- pmin(s, theta_j - bounds_j[1])
      if (is.finite(bounds_j[2])) s <- pmin(s, bounds_j[2] - theta_j)
    }
    h <- rel_step * s
  } else {
    reach <- numericals7::fd_offsets(order, accuracy = accuracy)$half_width
    h <- rel_step * pmax(1, abs(theta_j))
    if (!is.null(bounds_j)) {
      if (is.finite(bounds_j[1]))
        h <- pmin(h, 0.49 * (theta_j - bounds_j[1]) / reach)
      if (is.finite(bounds_j[2]))
        h <- pmin(h, 0.49 * (bounds_j[2] - theta_j) / reach)
    }
  }
  if (any(!is.finite(h) | h <= 0)) {
    stop("Cannot build finite-difference steps: some parameter values lie on or outside their domain boundary.", call. = FALSE)
  }
  h
}

#' @title The Self-Consistency of a Difference Quotient Under Step Halving
#'
#' @description
#' Returns \eqn{\max|q(h/2) - q(h)| / \max|q(h)|}, the reading on which
#' [fd_stable_step()] compares its two candidate steps. It is `Inf` where
#' the half-step quotient is not finite and `NA` where the full-step one is
#' identically zero, and in both cases the reading is not usable.
#'
#' @param x_half The quotient at half the step, a numeric vector.
#' @param x_full The quotient at the full step, the same length.
#'
#' @return A single number, possibly `Inf` or `NA`.
#'
#' @seealso [fd_stable_step()].
#' @keywords internal
fd_self_consistency <- function(x_half, x_full) {
  if (!all(is.finite(x_half))) return(Inf)
  m <- max(abs(x_full))
  if (!is.finite(m) || m == 0) return(NA_real_)
  max(abs(x_half - x_full)) / m
}

#' @title Selecting the Finite-Difference Step in One Parameter
#'
#' @description
#' Evaluates a difference quotient in one parameter at both steps of
#' [fd_steps()] and returns the one that agrees better with itself at half its
#' own step, together with the quotient there. The numerical gradient and
#' Hessian, the third and fourth orders that difference a Hessian, and the
#' mixed response derivatives take their step from it;
#' [tensor_derivatives()] and [numerical_deriv5()] use rules of their own.
#'
#' @details
#' For each of the two steps the quotient \eqn{q(h)} is also taken at
#' \eqn{h/2}, and the step kept is the one whose
#' [fd_self_consistency()] reading is strictly the smaller. Ties and missing
#' readings keep the clamped step, and an infinite reading on the clamped side
#' gives way to a finite one on the scaled side, so the answer is the clamped
#' step unless the scaled step is measurably better.
#'
#' The choice is made once for the whole vector rather than observation by
#' observation, which is the opposite of what [fd_stable_quotient()] does in
#' the response direction. There each
#' observation carries its own \eqn{y} and therefore its own step, so a
#' per-observation choice is a choice between two quantities that genuinely
#' differ; here the parameter is usually one number for the whole sample, the
#' two candidate quotients estimate the same thing, and the variation of the
#' self-consistency reading across observations is the rounding rather than a
#' signal.
#'
#' It costs four quotients where one would do, except where the two steps
#' coincide at every entry, which is every call on a parameter with no finite
#' bound and every positive parameter at or above one: there the clamped step
#' is returned at once and each quotient is evaluated only once.
#'
#' Each order selects its own step, with its own `rel_step` and its own stencil;
#' a choice made at first order is not reused at higher orders.
#'
#' @param quotient A function of one step, returning the difference quotient at
#'   that step: a numeric vector, or a list of them where the caller
#'   differences several components at once, in which case the reading is taken
#'   over all of them together.
#' @param theta_j A numeric vector, the values of one parameter.
#' @param bounds_j That parameter's domain, a length-2 numeric vector, or
#'   `NULL`.
#' @param rel_step The relative step size.
#' @param need_value Whether the caller will use the quotient at the chosen
#'   step. A caller that needs only the step sets it to `FALSE`, and the
#'   function then skips every evaluation where the two candidates coincide: the
#'   higher orders read only `h`, and their quotient is a difference of whole
#'   Hessians.
#' @param order,accuracy The stencil that the quotient implements; it is passed
#'   to [fd_steps()] so that the clamp applies to the outermost node. The
#'   defaults are the central difference; a caller whose quotient is a second
#'   difference passes `order = 2`.
#' @return A list of two elements: `h`, the step chosen, and `value`, the
#'   quotient at it, or `NULL` when `need_value` is `FALSE` and no quotient had
#'   to be taken.
#'
#' @seealso [fd_steps()] for the two candidates, [fd_self_consistency()] for
#'   the reading they are compared on, and [fd_stable_quotient()] for the
#'   response direction.
#' @keywords internal
fd_stable_step <- function(quotient, theta_j, bounds_j, rel_step, need_value = TRUE,
                           order = 1L, accuracy = 2L) {
  hA <- fd_steps(theta_j, bounds_j, rel_step, "clamp", order, accuracy)
  hB <- fd_steps(theta_j, bounds_j, rel_step, "scale", order, accuracy)
  if (isTRUE(all(hA == hB))) {
    return(list(h = hA, value = if (need_value) quotient(hA) else NULL))
  }
  flat <- function(x) unlist(x, use.names = FALSE)
  a1 <- quotient(hA)
  b1 <- quotient(hB)
  kA <- fd_self_consistency(flat(quotient(hA / 2)), flat(a1))
  kB <- fd_self_consistency(flat(quotient(hB / 2)), flat(b1))
  if (isTRUE(kB < kA)) list(h = hB, value = b1) else list(h = hA, value = a1)
}

#' Numerical Gradient of the Log-Density
#'
#' @description
#' Computes the gradient of the log-density with respect to each parameter by
#' central finite differences of `distrib_pdf(..., log = TRUE)`. This powers
#' the default [distrib_gradient()] method for distributions that do not
#' implement an analytical gradient: any `distrib` subclass that defines only
#' `distrib_pdf` receives its score function from this default at no further
#' cost to its author.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each element has length 1 or
#'   `length(y)`).
#' @param rel_step Numeric, the relative step size. Defaults to
#'   `.Machine$double.eps^(1/3)` (optimal for central differences).
#'
#' @return A named list (one element per parameter) of gradient vectors.
#'
#' @details
#' Each component is one central difference of the log-density in its own
#' parameter,
#'
#' \deqn{l^{(i)} = \frac{\partial}{\partial \theta_i} \log f(y; \theta)
#'   \method \frac{\log f(y; \theta + h_i e_i)
#'     - \log f(y; \theta - h_i e_i)}{2 h_i},}
#'
#' so the cost is two density evaluations per parameter when the two candidate
#' steps of [fd_steps()] coincide, and more when [fd_stable_step()] compares
#' them. Truncation is of
#' order \eqn{h^{2}} and rounding of order \eqn{\varepsilon / h}, which the
#' default \eqn{h \propto \varepsilon^{1/3}} balances.
#'
#' Steps are scaled by `max(1, |theta|)` and kept inside the boundaries of
#' `distrib@params_bounds`, by [fd_stable_step()], which reads both of
#' [fd_steps()]'s candidates and keeps the one that agrees better with itself
#' at half its own step. Accuracy is roughly `eps^(2/3)` (about 10 significant
#' digits at best): sufficient for optimization, but slower and less precise
#' than an analytical implementation.
#'
#' @seealso [numerical_hessian()], [distrib_gradient()]
#' @examples
#' numerical_gradient(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_gradient <- function(distrib, y, theta, rel_step = .Machine$double.eps^(1 / 3)) {
  params <- distrib@params
  out <- vector("list", length(params))
  names(out) <- params

  for (j in seq_along(params)) {
    quotient <- function(h) {
      tp <- tm <- theta
      tp[[j]] <- theta[[j]] + h
      tm[[j]] <- theta[[j]] - h
      (distrib_pdf(distrib, y, tp, log = TRUE) -
        distrib_pdf(distrib, y, tm, log = TRUE)) / (2 * h)
    }
    out[[j]] <- fd_stable_step(quotient, theta[[j]],
                               distrib@params_bounds[[params[j]]], rel_step)$value
  }

  out
}

#' Numerical Hessian of the Log-Density
#'
#' @description
#' Computes the observed Hessian of the log-density with respect to the parameters
#' by central finite differences of `distrib_pdf(..., log = TRUE)`. The
#' default [distrib_hessian()] method calls it for a family that has neither an
#' analytical Hessian nor an analytical score; for a family with an analytical
#' score, the method differences that score instead.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each element has length 1 or
#'   `length(y)`).
#' @param rel_step Numeric, the relative step size. Defaults to
#'   `.Machine$double.eps^(1/4)` (optimal for second differences).
#'
#' @return A named list of Hessian component vectors, in [hess_names()] order
#'   (diagonal elements first, then the upper-triangular mixed derivatives).
#'
#' @details
#' Diagonal components use the three-point stencil
#' \eqn{(\ell(\theta+h) - 2\ell(\theta) + \ell(\theta-h))/h^2}; mixed components use
#' the four-point cross stencil. Steps are scaled and clamped as in
#' [numerical_gradient()]. Accuracy is roughly `sqrt(eps)`.
#'
#' @seealso [numerical_gradient()], [distrib_hessian()]
#' @examples
#' numerical_hessian(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_hessian <- function(distrib, y, theta, rel_step = .Machine$double.eps^(1 / 4)) {
  params <- distrib@params
  p <- length(params)
  lp <- function(th) distrib_pdf(distrib, y, th, log = TRUE)

  out <- list()
  lp0 <- lp(theta)

  # The step of each parameter is chosen on that parameter's own diagonal
  # component and then used for every component it enters, the mixed ones
  # included, so the choice is paid once per parameter rather than once per
  # component.
  chosen <- lapply(seq_len(p), function(j) {
    quotient <- function(hj) {
      tp <- tm <- theta
      tp[[j]] <- theta[[j]] + hj
      tm[[j]] <- theta[[j]] - hj
      (lp(tp) - 2 * lp0 + lp(tm)) / (hj^2)
    }
    fd_stable_step(quotient, theta[[j]], distrib@params_bounds[[params[j]]],
                   rel_step, order = 2L)
  })
  h <- lapply(chosen, `[[`, "h")

  shift <- function(j, a, k = NULL, b = 0) {
    th <- theta
    th[[j]] <- th[[j]] + a * h[[j]]
    if (!is.null(k)) th[[k]] <- th[[k]] + b * h[[k]]
    th
  }

  # Diagonal: three-point stencil, already evaluated at the chosen step
  for (j in seq_len(p)) {
    out[[paste0(params[j], "_", params[j])]] <- chosen[[j]]$value
  }

  # Mixed: four-point cross stencil
  if (p > 1) {
    for (j in 1:(p - 1)) {
      for (k in (j + 1):p) {
        out[[paste0(params[j], "_", params[k])]] <-
          (lp(shift(j, 1, k, 1)) - lp(shift(j, 1, k, -1)) -
            lp(shift(j, -1, k, 1)) + lp(shift(j, -1, k, -1))) / (4 * h[[j]] * h[[k]])
      }
    }
  }

  out[hess_names(params)]
}


#' Derivatives by One Stencil on the Highest Analytical Order
#'
#' @description
#' Computes the derivative components of order `order` of the log-density
#' from the highest order that the family implements itself: the log-density
#' (`base = 0`), the score (`base = 1`), the Hessian (`base = 2`), the third
#' order (`base = 3`) or the fourth (`base = 4`). Each component is a single
#' tensor-product central difference of order `order - base`, never a
#' difference of a difference.
#'
#' @details
#' The computation is [numericals7::fd_tensor()] over the coordinates
#' \eqn{\theta_1, \dots, \theta_p}, where the log-density is the function and
#' the derivatives of order `base` are declared as the known quantity. A
#' component is a multiset of parameters: \eqn{(\mu, \mu, \phi)} is
#' \eqn{\partial^3\ell / \partial\mu^2\partial\phi}. `base` of its indices
#' select the analytical component that is differenced, and these indices are
#' taken from the parameters that occur most often. The remaining
#' multiplicities are the orders of the one-dimensional central stencils whose
#' tensor product is applied to it. Each parameter is differentiated only once,
#' however many times it occurs.
#'
#' The step in \eqn{\theta_j} is
#' \eqn{h_j = \varepsilon^{1/(r+a)}\max(1, |\theta_j|)}, with
#' \eqn{r = } `order - base` the order of the difference and \eqn{a} the
#' accuracy, reduced so that every node stays inside the parameter's bounds.
#' On the link scale the coordinates are the linear predictors
#' \eqn{\eta_j}, the bounds are those of the links, and the known quantity is
#' read on the link scale, so the result is the derivative of a link-scale
#' quantity.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param order The order of the derivatives, a whole number above `base`.
#' @param base The highest order that the family implements itself, 0 to 4.
#' @param rel_step `NULL`, the default, for the step above, or a relative step
#'   that replaces \eqn{\varepsilon^{1/(r+a)}}.
#' @param skip Character vector of component names left `NULL`, or `NULL`.
#' @param scale `"parameter"` or `"link"`, the scale of the coordinates and
#'   of the result.
#' @param accuracy The order of accuracy of the one-dimensional stencils,
#'   passed to [numericals7::fd_tensor()].
#'
#' @return A named list of component vectors, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names] gives them.
#'
#' @seealso [analytic_order()] for how `base` is found; [numerical_deriv3()],
#'   [numerical_deriv4()] and [numerical_deriv5()], which call this.
#' @keywords internal
tensor_derivatives <- function(distrib, y, theta, order, base, rel_step = NULL,
                               skip = NULL, scale = c("parameter", "link"),
                               accuracy = 4L) {
  scale <- match.arg(scale)
  order <- as.integer(order)
  base <- as.integer(base)
  if (order - base < 1L) stop("'order' must exceed 'base'.", call. = FALSE)
  params <- distrib@params
  p <- length(params)
  nms <- deriv_names(params, order)
  out <- vector("list", length(nms))
  names(out) <- nms
  want <- which(!nms %in% skip)
  if (!length(want)) return(out)

  cs <- tensor_coords(distrib, theta, scale)
  comp <- t(vapply(deriv_indices(params, order)[want], tabulate, integer(p),
                   nbins = p))
  if (p == 1L) comp <- matrix(comp, ncol = 1L)
  f <- function(x) distrib_pdf(distrib, y, cs$to_theta(x), log = TRUE)
  known <- list()
  if (base > 0L) {
    gen <- list(distrib_gradient, distrib_hessian, distrib_deriv3,
                distrib_deriv4)[[base]]
    keys <- if (base == 2L) hess_names(params) else deriv_names(params, base)
    known <- list(numericals7::fd_known(function(x) {
      unname(gen(distrib, y, cs$to_theta(x), scale = scale)[keys])
    }, order = base))
  }
  got <- numericals7::fd_tensor(f, cs$x, components = comp, known = known,
                                bounds = cs$bounds, accuracy = accuracy,
                                rel_step = rel_step)
  out[want] <- unname(got)
  out
}


#' The Coordinates of a Tensor Stencil over the Parameters
#'
#' @description
#' Returns the coordinates on which [tensor_derivatives()],
#' [mixed_tensor_derivatives()] and [numerical_deriv5()] difference: the
#' parameters themselves on the parameter scale, or the linear predictors
#' \eqn{\eta = g(\theta)} on the link scale, with the bounds of each and the
#' map back to a `theta` that the family's methods accept.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale `"parameter"` or `"link"`.
#'
#' @return A list with `x`, the named list of coordinates; `bounds`, a list of
#'   `c(lower, upper)` per coordinate; and `to_theta`, a function from a list
#'   of coordinates to `theta`.
#' @keywords internal
tensor_coords <- function(distrib, theta, scale = c("parameter", "link")) {
  scale <- match.arg(scale)
  params <- distrib@params
  if (scale == "link") {
    links <- distrib@link_params
    x <- lapply(params, function(nm)
      linkfunctions7::linkfun(links[[nm]], theta[[nm]]))
    bounds <- lapply(params, function(nm)
      linkfunctions7::eta_bounds(links[[nm]]))
    to_theta <- function(x) {
      th <- theta
      for (nm in params) {
        th[[nm]] <- linkfunctions7::linkinv(links[[nm]], x[[nm]])
      }
      th
    }
  } else {
    x <- lapply(params, function(nm) theta[[nm]])
    bounds <- lapply(params, function(nm) distrib@params_bounds[[nm]])
    to_theta <- function(x) {
      th <- theta
      for (nm in params) th[[nm]] <- x[[nm]]
      th
    }
  }
  names(x) <- params
  names(bounds) <- params
  list(x = x, bounds = bounds, to_theta = to_theta)
}


#' The Highest Derivative Order That a Family Implements Itself
#'
#' @description
#' Returns the highest order, up to `upto`, at which the family registers its
#' own method: 0 for a family that has only a density, 1 with a score, 2 with a
#' Hessian, 3 with third derivatives, 4 with fourth derivatives. A method
#' inherited from one of the package's base classes is a numerical fallback
#' and does not count; one inherited from a family class does.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param upto The highest order to be examined, 1 to 4.
#'
#' @return An integer from 0 to `upto`.
#'
#' @seealso [tensor_derivatives()], which differences that order.
#' @keywords internal
analytic_order <- function(distrib, upto) {
  gens <- list(distrib_gradient, distrib_hessian, distrib_deriv3,
               distrib_deriv4)
  m <- 0L
  for (o in seq_len(upto)) if (owns_method(distrib, gens[[o]])) m <- o
  m
}

# --- DEFAULT (FALLBACK) METHODS ---
# Registered on the base `distrib` class: any subclass that implements only
# distrib_pdf automatically gets a score function, an observed Hessian and an
# expected Hessian. Subclasses with analytical implementations override these
# through normal S7 dispatch.

#' @title Default Numerical Gradient for `distrib` Objects
#' @name distrib_gradient.distrib
#'
#' @description
#' This method is the fallback for a family that implements no analytical
#' score. It computes the gradient of `distrib_pdf(..., log = TRUE)` by one
#' central difference per parameter, through [numerical_gradient()]. For this
#' reason [distrib_pdf()]
#' is the only compulsory method of the package: a family that defines only the
#' density has a score, an information, five orders of derivative and a fit
#' available.
#'
#' @details
#' # The stencil, the step and the cost
#' Each component is \eqn{[\ell(\theta_i + h) - \ell(\theta_i - h)]/(2h)}, so
#' one gradient costs \eqn{2p} evaluations of the log-density when the two
#' candidate steps coincide. The step is
#' \eqn{h = \varepsilon^{1/3}\max(1, |\theta_i|) \method 6.06\times10^{-6}}
#' at a parameter of order one, which balances the \eqn{O(h^2)} truncation of a
#' central difference against a rounding term growing as \eqn{1/h}. Parameter
#' domains here are open and a step through zero returns `NaN` from the density
#' for reasons that look like a defect in the family, so near a finite boundary
#' the step is kept inside: [fd_steps()] offers it cut to 49\% of the distance
#' or scaled on that distance, and [fd_stable_step()] keeps whichever agrees
#' better with itself at half its own step, which removes most of the error of
#' the clamped step close to a boundary.
#'
#' @param distrib An object inheriting from `distrib` that registers no method
#'   of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`; each vector has length `length(y)`.
#'
#' @seealso [numerical_gradient()], which does the differencing;
#'   [fd_steps()] for the boundary rule;
#'   [distrib_hessian.distrib()] for the order above;
#'   [distrib_gradient()] for the generic.
#' @keywords internal
S7::method(distrib_gradient, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  numerical_gradient(distrib, y, theta)
}

#' @title Default Numerical Hessian for `distrib` Objects
#' @name distrib_hessian.distrib
#'
#' @description
#' This method is the fallback for a family that implements no analytical
#' Hessian. For a family that has only a density it takes second differences of
#' `distrib_pdf(..., log = TRUE)` through [numerical_hessian()]: a diagonal
#' component takes the three-point stencil
#' \eqn{[\ell(\theta_i+h) - 2\ell(\theta_i) + \ell(\theta_i-h)]/h^2} and an
#' off-diagonal one the four-point mixed stencil, so both are a single
#' difference of the log-density. For a family with its own score it takes
#' one central difference of that score instead, through
#' [tensor_derivatives()] with `base = 1`: the component \eqn{(i, j)} is the
#' score in \eqn{\theta_i} differenced along \eqn{\theta_j}, a first
#' difference of an analytical quantity, with the step
#' \eqn{\varepsilon^{1/3}\max(1, |\theta_j|)}.
#'
#' @details
#' # The step and the cost
#' For a family that has only a density, the step is \eqn{h =
#' \varepsilon^{1/4}\max(1, |\theta_i|) \method 1.22\times10^{-4}}, twenty times
#' the gradient's: a second difference divides by \eqn{h^2}, so rounding grows
#' as \eqn{1/h^2} and the optimum moves out. [fd_steps()] applies the same
#' boundary clamp. One Hessian costs \eqn{2p} evaluations for the diagonal,
#' \eqn{4} per distinct pair and one at the center, which is 9 in all for a
#' two-parameter family, and more when [fd_stable_step()] compares two candidate
#' steps. A second difference is about two digits less accurate than the first
#' difference of the gradient.
#'
#' @param distrib An object inheriting from `distrib` that registers no method
#'   of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param ... Unused.
#'
#' @return A named list of Hessian component vectors, each of which has length
#'   `length(y)`, keyed by [hess_names()], which puts the diagonal first and is
#'   not the lexicographic keying that [deriv_names()] uses above order 2.
#'
#' @seealso [numerical_hessian()] and [tensor_derivatives()], which do the
#'   differencing; [fd_steps()] for the boundary rule;
#'   [distrib_gradient.distrib()] for the order below;
#'   [distrib_expected_hessian()] for the expectation of this.
#' @keywords internal
S7::method(distrib_hessian, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  if (analytic_order(distrib, 1L) == 1L) {
    out <- tensor_derivatives(distrib, y, theta, order = 2L, base = 1L)
    return(out[hess_names(distrib@params)])
  }
  numerical_hessian(distrib, y, theta)
}

#' @title Default Expected Hessian for `distrib` Objects
#' @name distrib_expected_hessian.distrib
#' @description
#' This method is the fallback for a family that does not write its expected
#' information out. It rests on the second Bartlett identity,
#' \eqn{\mathbb{E}[\ell^{(ij)}] = -\mathbb{E}[\ell^{(i)}\ell^{(j)}]}, which
#' holds for a regular model and, unlike \eqn{\mathbb{E}[\ell^{(ij)}]} read
#' directly, survives a log-likelihood that is not differentiable in a
#' parameter (the location of a Laplace, where the observed Hessian is
#' degenerate while the score variance is still the information).
#'
#' @details
#' The argument `method` selects how the right-hand side is obtained, and the
#' methods differ in cost. With `method = NULL` the default is [expect_opg()],
#' which reads \eqn{-\ell^{(i)}\ell^{(j)}} at each observation and takes no
#' expectation, so it costs one call to [distrib_gradient()];
#' `expect_integrate(identity = "score")` evaluates the expectation itself,
#' which is a sum over the support for a discrete family and a quadrature for a
#' continuous one, and is orders of magnitude dearer. See [expected_by_opg()]
#' for the properties that the default loses and those that it retains.
#'
#' The score is taken from [distrib_gradient()], so it uses the analytical
#' gradient where the family has one and finite differences otherwise.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`, the scale on which the components
#'   are reported. The transformation is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param method `NULL` (the default), which selects [expect_opg()], or an
#'   object inheriting from [ExpectMethod] that [expected_derivative()]
#'   applies. See [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the
#'   generic's.
#' @return A named list of expected Hessian component vectors.
#' @seealso [expected_by_opg()] and [expected_by_bartlett()] for the two
#'   readings of the identity, and [expected_hessian_exact()] for the
#'   predicate that determines whether a family reaches this method at all.
#' @keywords internal
S7::method(distrib_expected_hessian, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), method = NULL, ...) {
  expected_derivative(distrib, y, theta, order = 2L,
                      method = method)
}
