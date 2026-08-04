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
#' all. Parameter domains here are \strong{open} -- a scale parameter is positive,
#' not non-negative -- so a step chosen from the magnitude alone will step a small
#' \eqn{\sigma} straight through zero, and the log-density comes back \code{NaN}
#' for reasons that look like a bug in the density. Clamping to 49\% of the
#' distance to each finite boundary keeps both evaluation points inside.
#'
#' A parameter already on or outside its boundary cannot be rescued this way, and
#' is reported rather than differentiated.
#'
#' @param theta_j A numeric vector, the values of one parameter.
#' @param bounds_j A length-2 numeric vector giving that parameter's domain, or
#'   \code{NULL}.
#' @param h_rel The relative step size, typically a root of machine epsilon
#'   chosen for the stencil in use.
#'
#' @return A numeric vector of steps, the same length as \code{theta_j}.
#'
#' @seealso \code{\link{numerical_gradient}}, \code{\link{numerical_hessian}}
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
#' central finite differences of \code{distrib_pdf(..., log = TRUE)}. This powers
#' the default \code{\link{distrib_gradient}} method for distributions that do not
#' implement an analytical gradient: any \code{distrib} subclass that defines only
#' \code{distrib_pdf} gets its score function for free.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or \code{length(y)}).
#' @param h_rel Numeric. Relative step size. Defaults to \code{.Machine$double.eps^(1/3)}
#'   (optimal for central differences).
#'
#' @return A named list (one element per parameter) of gradient vectors.
#'
#' @details
#' Steps are scaled by \code{max(1, |theta|)} and automatically shrunk near the
#' boundaries of \code{distrib@params_bounds} so that the evaluation points remain
#' inside the parameter domain. Accuracy is roughly \code{eps^(2/3)} (about 8
#' significant digits): sufficient for optimization, but slower and less precise
#' than an analytical implementation.
#'
#' @seealso \code{\link{numerical_hessian}}, \code{\link{distrib_gradient}}
#' @examples
#' numerical_gradient(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
#' by central finite differences of \code{distrib_pdf(..., log = TRUE)}. This powers
#' the default \code{\link{distrib_hessian}} method for distributions that do not
#' implement an analytical Hessian.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters (each of length 1 or \code{length(y)}).
#' @param h_rel Numeric. Relative step size. Defaults to \code{.Machine$double.eps^(1/4)}
#'   (optimal for second differences).
#'
#' @return A named list of Hessian component vectors, in \code{\link{hess_names}} order
#'   (diagonal elements first, then the upper-triangular mixed derivatives).
#'
#' @details
#' Diagonal components use the three-point stencil
#' \eqn{(\ell(\theta+h) - 2\ell(\theta) + \ell(\theta-h))/h^2}; mixed components use
#' the four-point cross stencil. Steps are scaled and clamped as in
#' \code{\link{numerical_gradient}}. Accuracy is roughly \code{sqrt(eps)}.
#'
#' @seealso \code{\link{numerical_gradient}}, \code{\link{distrib_hessian}}
#' @examples
#' numerical_hessian(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
#' @description
#' Fallback method: distributions that do not implement an analytical gradient get
#' one computed by central finite differences of the log-density
#' (see \code{\link{numerical_gradient}}).
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @return A named list of gradient vectors.
#' @keywords internal
S7::method(distrib_gradient, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  numerical_gradient(distrib, y, theta)
}

#' @title Default Numerical Hessian for `distrib` Objects
#' @name distrib_hessian.distrib
#' @description
#' Fallback method: distributions that do not implement an analytical Hessian get
#' one computed by finite differences of the log-density
#' (see \code{\link{numerical_hessian}}).
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @return A named list of Hessian component vectors in \code{\link{hess_names}} order.
#' @keywords internal
S7::method(distrib_hessian, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  numerical_hessian(distrib, y, theta)
}

#' @title Default Expected Hessian for `distrib` Objects
#' @name distrib_expected_hessian.distrib
#' @description
#' Fallback method: the expected Hessian is computed as the negative Fisher
#' information via the outer product of the score (gradient),
#' \eqn{-\mathbb{E}[\nabla\ell\,\nabla\ell^\top]}, obtained by numerical
#' \code{\link{expectation}}. This first-Bartlett form equals \eqn{\mathbb{E}[H]}
#' for regular (twice-differentiable) models but, unlike \eqn{\mathbb{E}[H]},
#' remains valid when the log-likelihood is non-differentiable in a parameter
#' (e.g. the location of a Laplace distribution), where the observed Hessian is
#' degenerate. The score is taken from \code{\link{distrib_gradient}}, so it uses
#' the analytical gradient when available and finite differences otherwise.
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @return A named list of expected Hessian component vectors.
#' @keywords internal
S7::method(distrib_expected_hessian, distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  expected_derivative(distrib, y, theta, order = 2L,
                      approx = match.arg(approx), nsim = nsim)
}
