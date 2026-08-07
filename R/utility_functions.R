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
  # density labeled "nu". Only worth a pass over the list when there is
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
#' d <- gaussian1_distrib()
#' check_theta_bounds(d, list(mu = 0, sigma = 1))
#' \dontrun{
#' check_theta_bounds(d, list(mu = 0, sigma = -1)) # error: sigma outside (0, Inf)
#' }
#'
#' @export
check_theta_bounds <- function(distrib, theta) {
  check_bounds_fast(distrib@params, distrib@params_bounds, theta, distrib@distrib_name)
}

#' Check Parameter Domains, Taking the Properties as Arguments
#'
#' @description
#' The body behind \code{\link{check_theta_bounds}}, differing only in that it
#' receives the distribution's properties instead of reaching for them.
#'
#' @details
#' \code{align_theta()} calls this on every generic invocation, so its fixed cost
#' is the package's fixed cost. Reading an S7 property costs a couple of
#' microseconds, which is material against a budget of a few tens, and the caller
#' has already read the properties it needs.
#'
#' Every offending parameter is collected before anything is raised, so a call
#' with two bad values reports both rather than one at a time. At most three
#' distinct offending values are shown.
#'
#' @param params A character vector of parameter names.
#' @param param_bounds A named list of length-2 domain vectors.
#' @param theta A list of parameter values, ordered as \code{params}.
#' @param distrib_name The distribution's name, used in the message.
#'
#' @return Invisibly \code{NULL}; raises an error if any value is outside its
#'   open domain or not finite.
#'
#' @seealso \code{\link{check_theta_bounds}}
#' @keywords internal
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
#' @examples
#' expand_params(list(mu = 0, sigma = 1), n = 3)
#'
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


#' Is an S7 Class the Given Base Class?
#'
#' @description
#' Compares two S7 class objects, treating a class re-created from the same
#' definition as the same class.
#'
#' @details
#' Several places here ask "did this method come from the base class, or did the
#' subclass register its own?", and answer it with the documented S7 trick of
#' reading \code{attr(m, "signature")[[1]]}. The comparison that follows must not
#' be \code{identical()}: on S7 class objects that is object identity, so it
#' returns \code{FALSE} for a class re-created from the same definition — which
#' is what happens whenever the package's code is re-evaluated rather than
#' loaded, as it is under coverage instrumentation.
#'
#' The failure is silent and can be severe. In \pkg{linkfunctions7} the same
#' mistake made a numerical fallback look like an analytic method and sent the
#' fourth derivative of the log link wrong by a factor of 900, under coverage
#' only, while every other check stayed green. Identity is kept as a fast path
#' and name-with-package as the answer.
#'
#' @param cls The S7 class recorded on a method.
#' @param base The base class to compare against.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{has_analytic_quantile}}, \code{\link{has_exact_cdf_deriv}}
#' @keywords internal
is_class <- function(cls, base) {
  if (identical(cls, base)) return(TRUE)
  identical(attr(cls, "name"), attr(base, "name")) &&
    identical(attr(cls, "package"), attr(base, "package"))
}


#' Invert the Hessian Component Names
#'
#' @description
#' Maps each name produced by \code{\link{hess_names}} back to the pair of
#' parameters it differentiates with respect to.
#'
#' @details
#' The wrappers need to go from \code{"mu_sigma"} back to
#' \code{c("mu", "sigma")} in order to combine the parent's score with its
#' Hessian. Splitting the string on \code{"_"} is the obvious way and it is
#' wrong: a parameter whose own name contains an underscore
#' (\code{"log_scale"}) makes \code{"log_scale_log_scale"} split into four
#' pieces, and taking the first and the last silently yields
#' \code{c("log", "scale")}. Building the map from the parameter vector cannot
#' be fooled.
#'
#' \code{\link{deriv_indices}} is the same idea for orders above two.
#'
#' @param params A character vector of parameter names.
#'
#' @return A named list, parallel to \code{hess_names(params)}, each element a
#'   character pair.
#'
#' @seealso \code{\link{hess_names}}, \code{\link{deriv_indices}}
#' @keywords internal
hess_pairs <- function(params) {
  nms <- hess_names(params)
  n <- length(params)
  pairs <- c(
    lapply(params, function(p) c(p, p)),
    if (n >= 2) {
      unlist(lapply(seq_len(n - 1L), function(i) {
        lapply((i + 1L):n, function(j) params[c(i, j)])
      }), recursive = FALSE)
    }
  )
  stats::setNames(pairs, nms)
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
  vapply(deriv_indices(params, order),
         function(r) paste(params[r], collapse = "_"), character(1))
}


#' Index Tuples Behind the Higher-Order Derivative Names
#'
#' @description
#' The multi-indices \code{\link{deriv_names}} names, in exactly the same order:
#' non-decreasing tuples of length \code{order} over \code{seq_along(params)},
#' enumerated lexicographically.
#'
#' @details
#' This exists so that nothing has to recover an index tuple by splitting a
#' component name back apart. Splitting fails for a parameter whose own name
#' contains an underscore: \code{"mu_log_scale_log_scale"} splits into five
#' pieces, and matching those against \code{params} yields \code{NA}s. The
#' multivariate families carry such names by construction
#' (\code{sigma_log_L1}), so the parsing route is wrong for shipped families
#' as well as user-defined ones. Generating the indices and the names from the
#' same enumeration removes the possibility of disagreement.
#'
#' Note that this is \strong{not} interchangeable with
#' \code{deriv_index_list()} in \code{link_scale.R}: that one is ordered to match
#' \code{\link{hess_names}} at order 2, which puts the diagonal first, whereas
#' this one is lexicographic throughout to match \code{\link{deriv_names}}. At
#' order 2 use \code{\link{hess_pairs}}.
#'
#' @param params A character vector of parameter names.
#' @param order A positive integer, the derivative order.
#'
#' @return A list of integer vectors, each of length \code{order}, parallel to
#'   \code{deriv_names(params, order)}.
#'
#' @seealso \code{\link{deriv_names}}, \code{\link{hess_pairs}}
#' @keywords internal
deriv_indices <- function(params, order) {
  p <- length(params)
  # Non-decreasing index tuples of length `order` over 1:p (combinations with
  # repetition), enumerated in lexicographic order.
  idx <- as.matrix(do.call(expand.grid, rep(list(seq_len(p)), order)))
  idx <- idx[, rev(seq_len(order)), drop = FALSE]           # lexicographic
  keep <- apply(idx, 1L, function(r) all(diff(r) >= 0))
  idx <- idx[keep, , drop = FALSE]
  lapply(seq_len(nrow(idx)), function(k) as.integer(idx[k, ]))
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
#' @examples
#' transpose_params(list(mu = c(0, 1), sigma = c(1, 2)))
#'
#' @export
transpose_params <- function(theta) {
  if (!length(theta)) return(list())

  # A true transpose: the names of the outer list become the names inside each
  # row, and the names inside the first element become the names of the result.
  # Applying the function twice therefore returns the original structure.
  #
  # Deliberately independent of the *inner* names when splitting: parameter
  # values that have traveled through a link function carry their own parameter
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