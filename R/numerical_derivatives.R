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
#' @param theta_j A numeric vector, the values of one parameter.
#' @param bounds_j A length-2 numeric vector giving that parameter's domain, or
#'   `NULL`.
#' @param h_rel The relative step size, typically a root of machine epsilon
#'   chosen for the stencil in use.
#' @param to_bound `"clamp"`, the default, or `"scale"`, as above.
#'
#' @return A numeric vector of steps, the same length as `theta_j`.
#'
#' @seealso [fd_stable_step()], which chooses between the two, and
#'   [numerical_gradient()] and [numerical_hessian()], which take the chosen
#'   step. [fd_steps_y()] is the response counterpart.
#' @keywords internal
fd_steps <- function(theta_j, bounds_j, h_rel, to_bound = c("clamp", "scale")) {
  to_bound <- match.arg(to_bound)
  if (identical(to_bound, "scale")) {
    s <- pmax(1, abs(theta_j))
    if (!is.null(bounds_j)) {
      if (is.finite(bounds_j[1])) s <- pmin(s, theta_j - bounds_j[1])
      if (is.finite(bounds_j[2])) s <- pmin(s, bounds_j[2] - theta_j)
    }
    h <- h_rel * s
  } else {
    h <- h_rel * pmax(1, abs(theta_j))
    if (!is.null(bounds_j)) {
      if (is.finite(bounds_j[1])) h <- pmin(h, 0.49 * (theta_j - bounds_j[1]))
      if (is.finite(bounds_j[2])) h <- pmin(h, 0.49 * (bounds_j[2] - theta_j))
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
#'
#' @return A list of two elements: `h`, the step chosen, and `value`, the
#'   quotient at it, or `NULL` when `need_value` is `FALSE` and no quotient had
#'   to be taken.
#'
#' @seealso [fd_steps()] for the two candidates, [fd_self_consistency()] for
#'   the reading they are compared on, and [fd_stable_quotient()] for the
#'   response direction.
#' @keywords internal
fd_stable_step <- function(quotient, theta_j, bounds_j, h_rel, need_value = TRUE) {
  hA <- fd_steps(theta_j, bounds_j, h_rel, "clamp")
  hB <- fd_steps(theta_j, bounds_j, h_rel, "scale")
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
    fd_stable_step(quotient, theta[[j]], distrib@params_bounds[[params[j]]], h_rel)
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
#' The fallback for a family that implements no analytical Hessian: second
#' differences of `distrib_pdf(..., log = TRUE)` through
#' [numerical_hessian()]. A diagonal component takes the three-point stencil
#' \eqn{[\ell(\theta_i+h) - 2\ell(\theta_i) + \ell(\theta_i-h)]/h^2} and an
#' off-diagonal one the four-point mixed stencil, so both are a **single**
#' difference of the log-density and neither is a difference of the gradient.
#'
#' @details
#' # The step and the cost
#' The step is \eqn{h = \varepsilon^{1/4}\max(1, |\theta_i|) \approx
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
#' @seealso [numerical_hessian()], which does the differencing;
#'   [fd_steps()] for the boundary rule;
#'   [distrib_gradient.distrib()] for the order below;
#'   [distrib_expected_hessian()] for the expectation of this.
#' @keywords internal
S7::method(distrib_hessian, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
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
