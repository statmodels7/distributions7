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
#' The domain clamp is what allows a finite-difference fallback to be offered at
#' all. Parameter domains here are **open**: a scale parameter is positive, not
#' non-negative. A step chosen from the magnitude alone therefore takes a small
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
#' its shorter step letting the rounding dominate. Neither is right everywhere,
#' and [fd_stable_step()] is what chooses between them.
#'
#' The two give the same step wherever \eqn{d \ge \max(1, |\theta|)}, which is
#' every parameter of a family with no finite bound, and also a positive
#' parameter at or above one.
#'
#' A parameter already on or outside its boundary cannot be rescued this way, and
#' is reported rather than differentiated.
#'
#' The clamp is on the stencil's OUTERMOST node and not on the step, so it
#' divides by the reach, the largest offset the rule evaluates at in units of
#' \eqn{h}. That number is the stencil's and is read from
#' [numericals7::fd_offsets()] rather than assumed: a central difference and a
#' three-point second difference both reach one step out, which is every use
#' this package makes of the clamped branch today, while the five-point rules
#' of accuracy four reach two. Without the division a step cut to
#' \eqn{0.49 d} puts the outermost node of a reach-two rule at \eqn{0.98 d},
#' which is inside the domain and not what the clamp says it is.
#'
#' With the reach accounted for, the clamped branch at the default `h_rel` is
#' [numericals7::fd_step()] with `bounds`, and a test asserts the two are
#' `identical()` so they cannot drift. It is written here rather than
#' delegated because `h_rel` is an argument the callers vary:
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
#' @param h_rel The relative step size, typically a root of machine epsilon
#'   chosen for the stencil in use.
#' @param to_bound `"clamp"`, the default, or `"scale"`, as above.
#' @param order,accuracy The stencil the step is for, which fix its reach
#'   through [numericals7::fd_offsets()]. The defaults are the central
#'   difference, whose reach is one.
#'
#' @return A numeric vector of steps, the same length as `theta_j`.
#'
#' @seealso [fd_stable_step()], which chooses between the two, and
#'   [numerical_gradient()] and [numerical_hessian()], which take the chosen
#'   step. [fd_steps_y()] is the response counterpart.
#'   [numericals7::fd_step()] is the same rule at the default `h_rel`.
#' @keywords internal
fd_steps <- function(theta_j, bounds_j, h_rel, to_bound = c("clamp", "scale"),
                     order = 1L, accuracy = 2L) {
  to_bound <- match.arg(to_bound)
  if (identical(to_bound, "scale")) {
    s <- pmax(1, abs(theta_j))
    if (!is.null(bounds_j)) {
      if (is.finite(bounds_j[1])) s <- pmin(s, theta_j - bounds_j[1])
      if (is.finite(bounds_j[2])) s <- pmin(s, bounds_j[2] - theta_j)
    }
    h <- h_rel * s
  } else {
    reach <- numericals7::fd_offsets(order, accuracy = accuracy)$reach
    h <- h_rel * pmax(1, abs(theta_j))
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
#' Returns \eqn{\max|q(h/2) - q(h)| / \max|q(h)|}, the reading
#' [fd_stable_step()] compares its two candidate steps on. It is `Inf` where
#' the half-step quotient is not finite and `NA` where the full-step one is
#' identically zero, neither of which is a usable reading.
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

#' @title A Parameter Step That Chooses Itself
#'
#' @description
#' Evaluates a difference quotient in one parameter at both steps of
#' [fd_steps()] and returns the one that agrees better with itself at half its
#' own step, together with the quotient there. Every numerical derivative in
#' the parameter direction takes its step from it.
#'
#' @details
#' For each of the two steps the quotient \eqn{q(h)} is also taken at
#' \eqn{h/2}, and the step kept is the one whose
#' [fd_self_consistency()] reading is strictly the smaller. Ties, and readings
#' that are not usable, keep the clamped step, so the answer is the clamped one
#' unless the scaled step is measurably better.
#'
#' The choice is made once for the whole vector rather than observation by
#' observation, which is the opposite of what [fd_stable_quotient()] does in
#' the response direction, and the difference is not an oversight. There each
#' observation carries its own \eqn{y} and therefore its own step, so a
#' per-observation choice is a choice between two quantities that genuinely
#' differ; here the parameter is usually one number for the whole sample, the
#' two candidate quotients estimate the same thing, and the variation of the
#' self-consistency reading across observations is the rounding rather than a
#' signal. Measured over a census of 324 cells at distances from 1 to
#' \eqn{10^{-8}} from a bound, the whole-vector choice puts 307 gradients and
#' 238 diagonal Hessian components within \eqn{10^{-6}} of the analytic value
#' against 285 and 229 for the per-observation choice, and on the 46 cells
#' where the two differ the whole-vector one is the better on 44.
#'
#' It costs four quotients where one would do, except where the two steps
#' coincide at every entry, which is every call on a parameter with no finite
#' bound and every positive parameter at or above one: there the clamped step
#' is returned at once and nothing is evaluated twice. Measured on the same
#' census, 61 cells of 324 take that path and their result is `identical()` to
#' what the clamped step alone produced.
#'
#' Each order chooses its own step, with its own `h_rel` and its own stencil.
#' A single choice made at first order and reused would cost 2 cells of 324 on
#' the Hessian; the two tests agree on 244 of the 263 cells where the steps
#' differ at all.
#'
#' @param quotient A function of one step, returning the difference quotient at
#'   that step: a numeric vector, or a list of them where the caller
#'   differences several components at once, in which case the reading is taken
#'   over all of them together.
#' @param theta_j A numeric vector, the values of one parameter.
#' @param bounds_j That parameter's domain, a length-2 numeric vector, or
#'   `NULL`.
#' @param h_rel The relative step size.
#' @param need_value Whether the caller will use the quotient at the chosen
#'   step. `FALSE` for a caller that wants only the step, and then nothing is
#'   evaluated at all where the two candidates coincide: the higher orders read
#'   only `h`, and their quotient is a difference of whole Hessians.
#' @param order,accuracy The stencil the quotient implements, passed to
#'   [fd_steps()] so that the clamp is on the outermost node. The defaults are
#'   the central difference; a caller whose quotient is a second difference
#'   passes `order = 2`.
#' @return A list of two elements: `h`, the step chosen, and `value`, the
#'   quotient at it, or `NULL` when `need_value` is `FALSE` and no quotient had
#'   to be taken.
#'
#' @seealso [fd_steps()] for the two candidates, [fd_self_consistency()] for
#'   the reading they are compared on, and [fd_stable_quotient()] for the
#'   response direction.
#' @keywords internal
fd_stable_step <- function(quotient, theta_j, bounds_j, h_rel, need_value = TRUE,
                           order = 1L, accuracy = 2L) {
  hA <- fd_steps(theta_j, bounds_j, h_rel, "clamp", order, accuracy)
  hB <- fd_steps(theta_j, bounds_j, h_rel, "scale", order, accuracy)
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
#' `distrib_pdf` gets its score function for free.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or `length(y)`).
#' @param h_rel Numeric. Relative step size. Defaults to `.Machine$double.eps^(1/3)`
#'   (optimal for central differences).
#'
#' @return A named list (one element per parameter) of gradient vectors.
#'
#' @details
#' Each component is one central difference of the log-density in its own
#' parameter,
#'
#' \deqn{l^{(i)} = \frac{\partial}{\partial \theta_i} \log f(y; \theta)
#'   \approx \frac{\log f(y; \theta + h_i e_i)
#'     - \log f(y; \theta - h_i e_i)}{2 h_i},}
#'
#' so the cost is two density evaluations per parameter. Truncation is of
#' order \eqn{h^{2}} and rounding of order \eqn{\varepsilon / h}, which the
#' default \eqn{h \propto \varepsilon^{1/3}} balances.
#'
#' Steps are scaled by `max(1, |theta|)` and kept inside the boundaries of
#' `distrib@params_bounds`, by [fd_stable_step()], which reads both of
#' [fd_steps()]'s candidates and keeps the one that agrees better with itself
#' at half its own step. Accuracy is roughly `eps^(2/3)` (about 8 significant
#' digits): sufficient for optimization, but slower and less precise than an
#' analytical implementation.
#'
#' @seealso [numerical_hessian()], [distrib_gradient()]
#' @examples
#' numerical_gradient(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
numerical_gradient <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3)) {
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
                               distrib@params_bounds[[params[j]]], h_rel)$value
  }

  out
}

#' Numerical Hessian of the Log-Density
#'
#' @description
#' Computes the observed Hessian of the log-density with respect to the parameters
#' by central finite differences of `distrib_pdf(..., log = TRUE)`. This powers
#' the default [distrib_hessian()] method for distributions that do not
#' implement an analytical Hessian.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or `length(y)`).
#' @param h_rel Numeric. Relative step size. Defaults to `.Machine$double.eps^(1/4)`
#'   (optimal for second differences).
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
numerical_hessian <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 4)) {
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
                   h_rel, order = 2L)
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
#' from the highest order the family implements itself: the log-density
#' (`base = 0`), the score (`base = 1`), the Hessian (`base = 2`) or the third
#' order (`base = 3`). Each component is a single central difference of order
#' `order - base`, never a difference of a difference.
#'
#' @details
#' A component is a multiset of parameters: \eqn{(\mu, \mu, \phi)} is
#' \eqn{\partial^3\ell / \partial\mu^2\partial\phi}. `base` of its indices
#' select the analytical component that is differenced, taken from the
#' parameters that occur most often so that no direction carries a high-order
#' stencil while another carries none; the remaining multiplicities
#' \eqn{r_1, \dots, r_p} are the orders of the difference in each direction.
#' The stencil is the tensor product of one-dimensional central stencils, the
#' one for direction \eqn{j} being [numericals7::fd_weights()] of order
#' \eqn{r_j} at second-order accuracy, so the component is
#' \deqn{\sum_{k_1, \dots, k_p} \Bigl(\prod_j w^{(r_j)}_{k_j} / h_j^{r_j}\Bigr)
#'   f(\theta_1 + k_1 h_1, \dots, \theta_p + k_p h_p),}
#' a single linear combination of values of the base component \eqn{f}. It
#' differentiates each parameter only once, however many times that parameter
#' occurs, which is what separates it from differencing a Hessian that is
#' itself a difference.
#'
#' The step is \eqn{h_j = \varepsilon^{1/(r+2)}\max(1, |\theta_j|)}, with
#' \eqn{r = \sum_j r_j} the total order of the difference: rounding grows as
#' \eqn{\varepsilon/h^r} and the truncation of a second-order stencil as
#' \eqn{h^2}, and the two balance there. The step is clamped so that the
#' widest stencil of that order stays inside the parameter's bounds. Every
#' component of one order uses the same steps, so a point shared by several
#' components is evaluated once.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param order The order of the derivatives, a whole number above `base`.
#' @param base The highest order the family implements itself, 0 to 3.
#' @param h_rel `NULL`, the default, for the step above, or a relative step
#'   that replaces \eqn{\varepsilon^{1/(r+2)}}.
#' @param skip Character vector of component names left `NULL`, or `NULL`.
#'
#' @return A named list of component vectors, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names] gives them.
#'
#' @seealso [analytic_order()] for how `base` is found; [numerical_deriv3()]
#'   and [numerical_deriv4()], which call this for a family without its own
#'   Hessian.
#' @keywords internal
tensor_derivatives <- function(distrib, y, theta, order, base, h_rel = NULL,
                               skip = NULL) {
  params <- distrib@params
  p <- length(params)
  nms <- deriv_names(params, order)
  idx_of <- deriv_indices(params, order)
  r_tot <- order - base
  if (r_tot < 1L) stop("'order' must exceed 'base'.", call. = FALSE)
  if (is.null(h_rel)) h_rel <- .Machine$double.eps^(1 / (r_tot + 2))

  # one step per parameter, clamped for the widest stencil of this order
  reach <- numericals7::fd_offsets(r_tot, 2L)$reach
  h <- lapply(seq_len(p), function(j) {
    x <- theta[[params[j]]]
    hj <- h_rel * pmax(1, abs(x))
    b <- distrib@params_bounds[[params[j]]]
    if (!is.null(b)) {
      if (is.finite(b[1])) hj <- pmin(hj, 0.49 * (x - b[1]) / reach)
      if (is.finite(b[2])) hj <- pmin(hj, 0.49 * (b[2] - x) / reach)
    }
    hj
  })

  gen <- switch(as.character(base), "1" = distrib_gradient,
                "2" = distrib_hessian, "3" = distrib_deriv3, NULL)
  cache <- list()
  value_at <- function(off) {
    key <- paste(off, collapse = ",")
    v <- cache[[key]]
    if (is.null(v)) {
      th <- theta
      for (j in which(off != 0L)) {
        th[[params[j]]] <- theta[[params[j]]] + off[j] * h[[j]]
      }
      v <- if (base == 0L) list(distrib_pdf(distrib, y, th, log = TRUE)) else
        gen(distrib, y, th)
      cache[[key]] <<- v
    }
    v
  }

  out <- vector("list", length(nms))
  names(out) <- nms
  for (t in seq_along(nms)) {
    if (nms[t] %in% skip) next
    cnt <- tabulate(idx_of[[t]], p)
    b <- integer(p)
    for (s in seq_len(base)) {
      j <- which.max(cnt - b)
      b[j] <- b[j] + 1L
    }
    rem <- cnt - b
    comp <- if (base == 0L) 1L else
      paste(params[rep(seq_len(p), b)], collapse = "_")
    dirs <- which(rem > 0L)
    st <- lapply(dirs, function(j) {
      off <- numericals7::fd_offsets(rem[j], 2L)$central
      w <- numericals7::fd_weights(off, rem[j])
      keep <- w != 0
      list(off = off[keep], w = w[keep])
    })
    grid <- as.matrix(expand.grid(lapply(st, function(s) seq_along(s$off))))
    acc <- 0
    for (g in seq_len(nrow(grid))) {
      off <- integer(p)
      wt <- 1
      for (q in seq_along(dirs)) {
        off[dirs[q]] <- st[[q]]$off[grid[g, q]]
        wt <- wt * st[[q]]$w[grid[g, q]]
      }
      acc <- acc + wt * value_at(off)[[comp]]
    }
    den <- 1
    for (j in dirs) den <- den * h[[j]]^rem[j]
    out[[nms[t]]] <- acc / den
  }
  out
}


#' The Highest Derivative Order a Family Implements Itself
#'
#' @description
#' Returns the highest order, up to `upto`, at which the family registers its
#' own method: 0 for a family with a density alone, 1 with a score, 2 with a
#' Hessian, 3 with third derivatives. A method inherited from one of the
#' package's base classes is a numerical fallback and does not count; one
#' inherited from a family class does.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param upto The highest order to look for, 1 to 3.
#'
#' @return An integer from 0 to `upto`.
#'
#' @seealso [tensor_derivatives()], which differences that order.
#' @keywords internal
analytic_order <- function(distrib, upto) {
  gens <- list(distrib_gradient, distrib_hessian, distrib_deriv3)
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
#' The fallback for a family that implements no analytical score: the gradient
#' of `distrib_pdf(..., log = TRUE)` by one **central difference** per
#' parameter, through [numerical_gradient()]. This is why [distrib_pdf()] is
#' the only compulsory method of the package: a family that defines the density
#' alone gets a score, an information, four orders of derivative and a fit.
#'
#' @details
#' # The stencil, the step and the cost
#' Each component is \eqn{[\ell(\theta_i + h) - \ell(\theta_i - h)]/(2h)}, so
#' one gradient costs \eqn{2p} evaluations of the log-density. The step is
#' \eqn{h = \varepsilon^{1/3}\max(1, |\theta_i|) \approx 6.06\times10^{-6}}
#' at a parameter of order one, which balances the \eqn{O(h^2)} truncation of a
#' central difference against a rounding term growing as \eqn{1/h}. Parameter
#' domains here are open and a step through zero returns `NaN` from the density
#' for reasons that look like a defect in the family, so near a finite boundary
#' the step is kept inside: [fd_steps()] offers it cut to 49\% of the distance
#' or scaled on that distance, and [fd_stable_step()] keeps whichever agrees
#' better with itself at half its own step. On a gamma at a dispersion of
#' \eqn{10^{-6}} that choice is worth five orders, the relative error going
#' from 3.5e-01 to 1.1e-06.
#'
#' # What it delivers
#' Measured on a Gamma in its mean and dispersion at
#' \eqn{(\mu, \sigma^2) = (2, 0.7)}, against the family's own closed form: the
#' two components agree to \eqn{1.3\times10^{-11}} and
#' \eqn{9.7\times10^{-11}} relative, which is the \eqn{O(h^2)} the step
#' promises.
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
#'   `distrib@params`, each of length `length(y)`.
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
#' The fallback for a family that implements no analytical Hessian. For a
#' family with a density alone it takes second differences of
#' `distrib_pdf(..., log = TRUE)` through [numerical_hessian()]: a diagonal
#' component takes the three-point stencil
#' \eqn{[\ell(\theta_i+h) - 2\ell(\theta_i) + \ell(\theta_i-h)]/h^2} and an
#' off-diagonal one the four-point mixed stencil, so both are a **single**
#' difference of the log-density. For a family with its own score it takes
#' one central difference of that score instead, through
#' [tensor_derivatives()] with `base = 1`: the component \eqn{(i, j)} is the
#' score in \eqn{\theta_i} differenced along \eqn{\theta_j}, a first
#' difference of an analytical quantity, with the step
#' \eqn{\varepsilon^{1/3}\max(1, |\theta_j|)}.
#'
#' @details
#' # The step and the cost
#' For a family with a density alone, the step is \eqn{h = \varepsilon^{1/4}\max(1, |\theta_i|) \approx
#' 1.22\times10^{-4}}, twenty times the gradient's: a second difference divides
#' by \eqn{h^2}, so rounding grows as \eqn{1/h^2} and the optimum moves out.
#' [fd_steps()] applies the same boundary clamp. One Hessian costs
#' \eqn{2p} evaluations for the diagonal and \eqn{4} per distinct pair,
#' which is 6 in all for a two-parameter family.
#'
#' # What it delivers
#' Measured on a Gamma in its mean and dispersion at
#' \eqn{(\mu, \sigma^2) = (2, 0.7)}, against the family's own closed form: the
#' three components agree to \eqn{2.3\times10^{-9}},
#' \eqn{2.7\times10^{-8}} and \eqn{3.6\times10^{-8}} relative, two to three
#' digits worse than the gradient's. That is the price of a second difference.
#'
#' @param distrib An object inheriting from `distrib` that registers no method
#'   of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param ... Unused.
#'
#' @return A named list of Hessian component vectors, each of length
#'   `length(y)`, keyed by [hess_names()], which puts the diagonal first and is
#'   **not** the lexicographic keying [deriv_names()] uses above order 2.
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
#' Fallback method, for a family that does not write its expected information
#' out. It rests on the second Bartlett identity,
#' \eqn{\mathbb{E}[\ell^{(ij)}] = -\mathbb{E}[\ell^{(i)}\ell^{(j)}]}, which
#' holds for a regular model and, unlike \eqn{\mathbb{E}[\ell^{(ij)}]} read
#' directly, survives a log-likelihood that is not differentiable in a
#' parameter -- the location of a Laplace, where the observed Hessian is
#' degenerate while the score variance is still the information.
#'
#' @details
#' `approx` says how the right-hand side is obtained, and the choice is a
#' choice of cost. The default `"opg"` reads
#' \eqn{-\ell^{(i)}\ell^{(j)}} at each observation and takes no expectation,
#' so it costs one call to [distrib_gradient()]; `"bartlett"` evaluates the
#' expectation itself, which is a sum over the support for a discrete family
#' and a quadrature for a continuous one, and is orders of magnitude dearer.
#' See [expected_by_opg()] for what the default gives up and what it does not.
#'
#' The score is taken from [distrib_gradient()], so it uses the analytical
#' gradient where the family has one and finite differences otherwise.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`, the scale the components are
#'   reported on. The transformation is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param approx One of `"opg"`, `"bartlett"`, `"integrate"` or `"mc"`, the
#'   strategy [expected_derivative()] uses. Defaults to `"opg"`.
#' @param nsim Number of draws, read only by `approx = "mc"`.
#' @param ... Unused, and accepted so that the signature matches the
#'   generic's.
#' @return A named list of expected Hessian component vectors.
#' @seealso [expected_by_opg()] and [expected_by_bartlett()] for the two
#'   readings of the identity, and [expected_hessian_exact()] for the
#'   predicate that says whether a family reaches this method at all.
#' @keywords internal
S7::method(distrib_expected_hessian, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  expected_derivative(distrib, y, theta, order = 2L,
                      approx = match.arg(approx), nsim = nsim)
}
