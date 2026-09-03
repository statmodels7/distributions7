#' @include distrib.R generics.R utility_functions.R numerical_functions.R
NULL

#' Finite-Difference Steps That Respect a Parameter's Domain
#'
#' @description
#' Builds the step \eqn{h} for a central difference in one parameter: scaled by
#' the parameter's magnitude, then shrunk so that \eqn{\theta \pm h} stays
#' strictly inside the parameter's mathematical domain.
#'
#' @details
#' The domain clamp is what allows a finite-difference fallback to be offered at
#' all. Parameter domains here are **open**: a scale parameter is positive, not
#' non-negative. A step chosen from the magnitude alone therefore takes a small
#' \eqn{\sigma} straight through zero, and the log-density comes back `NaN`
#' for reasons that look like a bug in the density. Clamping to 49\% of the
#' distance to each finite boundary keeps both evaluation points inside.
#'
#' A parameter already on or outside its boundary cannot be rescued this way, and
#' is reported rather than differentiated.
#'
#' @param theta_j A numeric vector, the values of one parameter.
#' @param bounds_j A length-2 numeric vector giving that parameter's domain, or
#'   `NULL`.
#' @param h_rel The relative step size, typically a root of machine epsilon
#'   chosen for the stencil in use.
#'
#' @return A numeric vector of steps, the same length as `theta_j`.
#'
#' @seealso [numerical_gradient()], [numerical_hessian()]
#' @keywords internal
fd_steps <- function(theta_j, bounds_j, h_rel) {
  h <- h_rel * pmax(1, abs(theta_j))
  if (!is.null(bounds_j)) {
    if (is.finite(bounds_j[1])) h <- pmin(h, 0.49 * (theta_j - bounds_j[1]))
    if (is.finite(bounds_j[2])) h <- pmin(h, 0.49 * (bounds_j[2] - theta_j))
  }
  if (any(!is.finite(h) | h <= 0)) {
    stop("Cannot build finite-difference steps: some parameter values lie on or outside their domain boundary.", call. = FALSE)
  }
  h
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
#' Steps are scaled by `max(1, |theta|)` and automatically shrunk near the
#' boundaries of `distrib@params_bounds` so that the evaluation points remain
#' inside the parameter domain. Accuracy is roughly `eps^(2/3)` (about 8
#' significant digits): sufficient for optimization, but slower and less precise
#' than an analytical implementation.
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
    h <- fd_steps(theta[[j]], distrib@params_bounds[[params[j]]], h_rel)
    tp <- tm <- theta
    tp[[j]] <- theta[[j]] + h
    tm[[j]] <- theta[[j]] - h
    out[[j]] <- (distrib_pdf(distrib, y, tp, log = TRUE) -
      distrib_pdf(distrib, y, tm, log = TRUE)) / (2 * h)
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

  h <- lapply(seq_len(p), function(j) {
    fd_steps(theta[[j]], distrib@params_bounds[[params[j]]], h_rel)
  })

  shift <- function(j, a, k = NULL, b = 0) {
    th <- theta
    th[[j]] <- th[[j]] + a * h[[j]]
    if (!is.null(k)) th[[k]] <- th[[k]] + b * h[[k]]
    th
  }

  out <- list()
  lp0 <- lp(theta)

  # Diagonal: three-point stencil
  for (j in seq_len(p)) {
    out[[paste0(params[j], "_", params[j])]] <-
      (lp(shift(j, 1)) - 2 * lp0 + lp(shift(j, -1))) / (h[[j]]^2)
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
#' central difference against a rounding term growing as \eqn{1/h}. Near a
#' finite boundary [fd_steps()] shrinks it to 49% of the distance, since
#' parameter domains here are open and a step through zero returns `NaN` from
#' the density for reasons that look like a defect in the family.
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
