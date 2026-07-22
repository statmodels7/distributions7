#' Align Parameters to the Distribution's Parameter Order
#'
#' @description
#' Coerces `theta` to a list and, when it is named, reorders it to match
#' `distrib@params`, so that methods can safely access parameters by position.
#' This prevents silently wrong results when a user supplies a named list in a
#' different order (e.g. `list(sigma = 2, mu = 0)` for a Gaussian).
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param theta A named list (or named numeric vector) of parameters, or an
#'   unnamed list/vector given in the order of `distrib@params`.
#'
#' @return A list whose first `distrib@n_params` elements correspond to
#'   `distrib@params`, in that order.
#'
#' @keywords internal
#' @noRd
align_theta <- function(distrib, theta) {
  if (!is.list(theta)) theta <- as.list(theta)

  # Every generic funnels through here, so this runs on every call and its fixed
  # cost is the package's fixed cost. Each S7 property read costs a couple of
  # microseconds, so the two that are needed are read once and passed on rather
  # than reached for again inside the bounds check.
  params <- distrib@params
  bounds <- distrib@params_bounds
  nms <- names(theta)

  # Strip names from the values themselves. A parameter that has been through a
  # link function comes back carrying its own name, which is meaningless on a
  # numeric value and leaks into the results: distrib_pdf() would return a
  # density labelled "nu". Only worth a pass over the list when there is
  # something to strip.
  if (any(vapply(theta, function(x) !is.null(names(x)), logical(1)))) {
    theta <- lapply(theta, unname)
  }

  # Unnamed: trust positional order, but require enough elements
  if (is.null(nms) || !any(nzchar(nms))) {
    if (length(theta) < distrib@n_params) {
      stop(sprintf(
        "'theta' has %d element(s) but %d parameter(s) are expected (%s).",
        length(theta), distrib@n_params, paste(params, collapse = ", ")
      ), call. = FALSE)
    }
    check_bounds_fast(params, bounds, theta, distrib@distrib_name)
    return(theta)
  }

  # Named: all parameters must be present, then reorder (extras are dropped).
  # The overwhelmingly common case is a list already in the right order, and
  # identical() settles that in a fraction of a microsecond; setdiff() costs
  # twenty even on two names, so it is kept for the branch that reports the error.
  if (!identical(nms, params)) {
    if (!all(params %in% nms)) {
      stop(sprintf(
        "Missing parameter(s) in 'theta': %s. Expected: %s.",
        paste(setdiff(params, nms), collapse = ", "), paste(params, collapse = ", ")
      ), call. = FALSE)
    }
    theta <- theta[params]
  }
  check_bounds_fast(params, bounds, theta, distrib@distrib_name)
  theta
}

#' Check Parameter Values Against Their Domains
#'
#' @description
#' Verifies that every supplied parameter value lies strictly inside the
#' distribution's \code{params_bounds}, and is finite. Domains are treated as
#' \strong{open} intervals: for instance a Gaussian requires \eqn{\sigma > 0} and a
#' Bernoulli requires \eqn{0 < \mu < 1}, since the log-likelihood and its
#' derivatives are not defined at the boundary.
#'
#' This is called automatically by every generic (through the internal
#' \code{align_theta()}), so passing an out-of-domain value raises an informative
#' error instead of silently producing \code{NaN}. It is exported so that it can
#' also be used directly, e.g. when writing an optimizer.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param theta A list of parameter values, ordered as \code{distrib@params}.
#'
#' @return Invisibly \code{NULL}. Raises an error listing every offending
#'   parameter, the offending value(s) and the expected domain.
#'
#' @examples
#' d <- gaussian_distrib()
#' check_theta_bounds(d, list(mu = 0, sigma = 1))
#' \dontrun{
#' check_theta_bounds(d, list(mu = 0, sigma = -1)) # error: sigma outside (0, Inf)
#' }
#'
#' @export
check_theta_bounds <- function(distrib, theta) {
  check_bounds_fast(distrib@params, distrib@params_bounds, theta, distrib@distrib_name)
}

# The same check, taking the properties it needs as arguments. align_theta calls
# it on every generic invocation, and reading an S7 property costs a couple of
# microseconds each time, which is material against a total budget of a few tens.
check_bounds_fast <- function(params, param_bounds, theta, distrib_name) {
  problems <- character()

  for (i in seq_along(params)) {
    p <- params[i]
    b <- param_bounds[[p]]
    if (is.null(b) || i > length(theta)) next
    v <- theta[[i]]

    if (!is.numeric(v)) {
      problems <- c(problems, sprintf("  '%s' must be numeric (got %s)", p, class(v)[1]))
      next
    }

    ok <- is.finite(v) & v > b[1] & v < b[2]
    if (!all(ok)) {
      offend <- unique(v[!ok])
      shown <- paste(format(utils::head(offend, 3), trim = TRUE), collapse = ", ")
      if (length(offend) > 3) shown <- paste0(shown, ", ...")
      problems <- c(problems, sprintf(
        "  '%s' = %s is outside its domain (%s, %s)",
        p, shown, format(b[1]), format(b[2])
      ))
    }
  }

  if (length(problems) > 0) {
    stop("Invalid parameter value(s) for the '", distrib_name,
         "' distribution:\n", paste(problems, collapse = "\n"), call. = FALSE)
  }

  invisible(NULL)
}


#' Check Consistency of Parameter Dimensions
#'
#' Validates that all elements in the provided parameter list have compatible lengths.
#' Each parameter must have a length of either 1 (scalar) or exactly equal to
#' `n`. This ensures safe vector recycling and dimensional consistency.
#'
#' @param theta A named list of vectors (parameters). Each element represents
#'   a parameter of a distribution (e.g., `mu`, `sigma`).
#' @param n (Optional) An integer specifying the required maximum length.
#'   If not provided, it defaults to the maximum length found among the elements
#'   of `theta`. Providing this argument allows validation against an external
#'   dimension (e.g., sample size `n`).
#'
#' @return Returns `NULL` invisibly if the check passes.
#'
#' @section Errors:
#' The function throws an error (`stop`) if it detects any parameter with a length
#' that is neither 1 nor `n`. The error message lists the specific parameters
#' causing the mismatch.
#'
#' @examples
#' # --- Case 1: Implicit max length ---
#' # Valid: all scalars
#' check_params_dim(list(mu = 1, sigma = 2))
#'
#' # Valid: mixing scalar and vector
#' check_params_dim(list(mu = 1:5, sigma = 1))
#'
#' # Invalid: incompatible lengths (2 vs 3)
#' \dontrun{
#' check_params_dim(list(mu = 1:2, sigma = 1:3))
#' }
#'
#' # --- Case 2: Explicit n ---
#' # Valid: vector matches n (5)
#' check_params_dim(list(mu = 1:5, sigma = 1), n = 5)
#'
#' # Invalid: vector length (3) does not match required n (5)
#' # This is useful to enforce consistency with a dataset size n = 5
#' \dontrun{
#' check_params_dim(list(mu = 1:3, sigma = 1), n = 5)
#' }
#'
#' @export
check_params_dim <- function(theta, n) {
  len_theta <- lengths(theta)

  if (missing(n)) {
    n <- max(len_theta)
  }

  # Check: length must be 1 OR exactly n
  mismatch_idx <- which(len_theta != 1 & len_theta != n)

  if (length(mismatch_idx) > 0) {
    bad_params <- names(theta)[mismatch_idx]
    if (is.null(bad_params)) bad_params <- paste0("[[", mismatch_idx, "]]")
    bad_lens <- len_theta[mismatch_idx]

    stop(
      "Parameter dimension mismatch. All parameters should have length 1 or ", n, ".\n",
      "  Offending: ",
      paste0(bad_params, " (length ", bad_lens, ")", collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Expand Parameters to Common Length
#'
#' Expands scalar parameters in a list to match the maximum length found (or a specified length),
#' ensuring all vectors are ready for element-wise operations.
#'
#' @param theta A named list of parameters.
#' @param n (Optional) The target length. If missing, defaults to `max(lengths(theta))`.
#'
#' @return A list where all elements have length `n`.
#' @export
expand_params <- function(theta, n) {
  lens <- lengths(theta)

  if (missing(n)) {
    n <- max(lens)
  }

  if (all(lens == n)) {
    return(theta)
  }

  check_params_dim(theta, n = n)

  idx_to_expand <- which(lens == 1)
  theta[idx_to_expand] <- lapply(theta[idx_to_expand], rep, times = n)

  theta
}


#' Generate Names for Hessian Matrix Components
#'
#' @description
#' Generates the names of the unique second-order partial derivatives (Hessian
#' components) for a vector of parameter names: first the diagonal elements
#' (\code{"mu_mu"}, ...), then the upper-triangular off-diagonal elements in
#' row-major order (\code{"mu_sigma"}, ...).
#'
#' @param params A character vector of parameter names (e.g., \code{c("mu", "sigma")}).
#'
#' @return A character vector of length \eqn{n + n(n-1)/2}.
#'
#' @examples
#' hess_names(c("mu", "sigma"))
#' # "mu_mu" "sigma_sigma" "mu_sigma"
#'
#' @export
hess_names <- function(params) {
  n_params <- length(params)
  diagonal <- paste0(params, "_", params)
  if (n_params < 2) {
    return(diagonal)
  }
  off_diagonal <- character(0.5 * n_params * (n_params - 1))
  k <- 1
  for (i in 1:(n_params - 1)) {
    for (j in (i + 1):n_params) {
      off_diagonal[k] <- paste0(params[i], "_", params[j])
      k <- k + 1
    }
  }
  c(diagonal, off_diagonal)
}


#' Generate Names for Higher-Order Derivative Components
#'
#' @description
#' Generates the names of the unique partial derivatives of a given \code{order}
#' with respect to a vector of parameters. Because mixed partial derivatives are
#' symmetric, only one representative per multi-index is listed, using
#' non-decreasing parameter order (e.g. \code{"mu_mu_sigma"} but not
#' \code{"mu_sigma_mu"}). For \code{order = 2} this coincides with the set of
#' \code{\link{hess_names}} (though possibly in a different order).
#'
#' @param params A character vector of parameter names (e.g., \code{c("mu", "sigma")}).
#' @param order A positive integer, the derivative order (e.g. \code{3} or \code{4}).
#'
#' @return A character vector of the \eqn{\binom{p + \text{order} - 1}{\text{order}}}
#'   unique component names, where \eqn{p} is the number of parameters.
#'
#' @examples
#' deriv_names(c("mu", "sigma"), 3)
#' # "mu_mu_mu" "mu_mu_sigma" "mu_sigma_sigma" "sigma_sigma_sigma"
#'
#' @export
deriv_names <- function(params, order) {
  p <- length(params)
  # Non-decreasing index tuples of length `order` over 1:p (combinations with
  # repetition), enumerated in lexicographic order.
  idx <- as.matrix(do.call(expand.grid, rep(list(seq_len(p)), order)))
  idx <- idx[, rev(seq_len(order)), drop = FALSE]           # lexicographic
  keep <- apply(idx, 1L, function(r) all(diff(r) >= 0))
  idx <- idx[keep, , drop = FALSE]
  apply(idx, 1L, function(r) paste(params[r], collapse = "_"))
}


#' Transpose and Simplify Parameter List Structure
#'
#' @description
#' Transposes a list structure (swapping "columns" and "rows") and simplifies the
#' inner elements into atomic vectors.
#'
#' Turns a list of \code{k} equal-length vectors into a list of \code{n} vectors of
#' length \code{k}, one per observation, keeping the names.
#'
#' @param theta A list to be transposed.
#' @return A \code{list} where each element has been transposed and simplified to an atomic vector.
#'
#' @export
transpose_params <- function(theta) {
  if (!length(theta)) return(list())

  # A true transpose: the names of the outer list become the names inside each
  # row, and the names inside the first element become the names of the result.
  # Applying the function twice therefore returns the original structure.
  #
  # Deliberately independent of the *inner* names when splitting: parameter
  # values that have travelled through a link function carry their own parameter
  # name, and keying the rows by those names collapses a multi-parameter theta
  # onto its first column.
  n <- length(theta[[1L]])
  inner <- names(theta)
  outer <- names(theta[[1L]])

  out <- lapply(seq_len(n), function(i) {
    stats::setNames(vapply(theta, function(x) x[[i]], numeric(1)), inner)
  })
  names(out) <- outer
  out
}