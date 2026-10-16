#' Align Parameters to the Distribution's Parameter Order
#'
#' @description
#' Coerces `theta` to a list and, when it is named, reorders it to match
#' `distrib@params`, so that methods can safely access parameters by position.
#' This prevents silently wrong results when a user supplies a named list in a
#' different order (e.g. `list(sigma = 2, mu = 0)` for a Gaussian).
#'
#' @param distrib An object inheriting from class `"distrib"`.
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
#' distribution's `params_bounds`, and is finite. Domains are treated as
#' open intervals: for instance a Gaussian requires \eqn{\sigma > 0} and a
#' Bernoulli requires \eqn{0 < \mu < 1}, since the log-likelihood and its
#' derivatives are not defined at the boundary.
#'
#' This is called automatically by every generic (through the internal
#' `align_theta()`), so passing an out-of-domain value raises an informative
#' error instead of silently producing `NaN`. It is exported so that it can
#' also be used directly, e.g. when writing an optimizer.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param theta A list of parameter values, ordered as `distrib@params`.
#'
#' @return Invisibly `NULL`. Raises an error listing every offending
#'   parameter, the offending value(s) and the expected domain.
#'
#' @examples
#' d <- gaussian1_distrib()
#' check_theta_bounds(d, list(mu = 0, sigma = 1))
#' \dontrun{
#' check_theta_bounds(d, list(mu = 0, sigma = -1)) # error: sigma outside (0, Inf)
#' }
#'
#' @seealso [deriv_names()], [hess_names()], [expand_params()], [transpose_params()], [check_params_dim()], [param_smoothness()], [generate_random_theta()]
#' @export
check_theta_bounds <- function(distrib, theta) {
  check_bounds_fast(distrib@params, distrib@params_bounds, theta, distrib@distrib_name)
}

#' Check Parameter Domains, Taking the Properties as Arguments
#'
#' @description
#' The body behind [check_theta_bounds()], differing only in that it
#' receives the distribution's properties instead of reaching for them.
#'
#' @details
#' `align_theta()` calls this on every generic invocation, so its fixed cost
#' is the package's fixed cost. Each read of an S7 property has a cost of its
#' own, so the caller reads the properties it needs once and passes them on.
#'
#' Every offending parameter is collected before anything is raised, so a call
#' with two bad values reports both rather than one at a time. At most three
#' distinct offending values are shown.
#'
#' @param params A character vector of parameter names.
#' @param param_bounds A named list of length-2 domain vectors.
#' @param theta A list of parameter values, ordered as `params`.
#' @param distrib_name The distribution's name, used in the message.
#'
#' @return Invisibly `NULL`; raises an error if any value is outside its
#'   open domain or not finite.
#'
#' @seealso [check_theta_bounds()]
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
#' The function signals an error (`stop`) if any parameter has a length
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
#' @seealso [deriv_names()], [hess_names()], [expand_params()], [transpose_params()], [check_theta_bounds()], [param_smoothness()], [generate_random_theta()]
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
#' @seealso [deriv_names()], [hess_names()], [transpose_params()], [check_params_dim()], [check_theta_bounds()], [param_smoothness()], [generate_random_theta()]
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
#' (`"mu_mu"`, ...), then the upper-triangular off-diagonal elements in
#' row-major order (`"mu_sigma"`, ...).
#'
#' @param params A character vector of parameter names (e.g., `c("mu", "sigma")`).
#'
#' @return A character vector of length \eqn{n + n(n-1)/2}.
#'
#' @examples
#' hess_names(c("mu", "sigma"))
#' # "mu_mu" "sigma_sigma" "mu_sigma"
#'
#' @seealso [deriv_names()], [expand_params()], [transpose_params()], [check_params_dim()], [check_theta_bounds()], [param_smoothness()], [generate_random_theta()]
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


#' Whether an S7 Class Is the Given Base Class
#'
#' @description
#' Compares two S7 class objects, treating a class re-created from the same
#' definition as the same class.
#'
#' @details
#' Several functions here decide whether a method comes from the base class or
#' from a subclass that registered its own, by reading
#' `attr(m, "signature")[[1]]`, the documented S7 route. The comparison that
#' follows must not use `identical()`: on S7 class objects that tests object
#' identity, so it returns `FALSE` for a class re-created from the same
#' definition. A class is re-created whenever the code of the package is
#' re-evaluated instead of being loaded, as under coverage instrumentation.
#'
#' When `identical()` fails in this way the comparison is silent, and a
#' numerical fallback is taken for an analytic method. The function therefore
#' returns `TRUE` for identical objects as a fast path, and otherwise compares
#' the class name together with the package name.
#'
#' @param cls The S7 class recorded on a method.
#' @param base The base class to compare against.
#'
#' @return A single logical.
#'
#' @seealso [has_analytic_quantile()], [has_exact_cdf_deriv()]
#' @keywords internal
is_class <- function(cls, base) {
  if (identical(cls, base)) return(TRUE)
  identical(attr(cls, "name"), attr(base, "name")) &&
    identical(attr(cls, "package"), attr(base, "package"))
}


#' Invert the Hessian Component Names
#'
#' @description
#' Maps each name produced by [hess_names()] back to the pair of
#' parameters it differentiates with respect to.
#'
#' @details
#' The wrappers need to go from `"mu_sigma"` back to
#' `c("mu", "sigma")` in order to combine the parent's score with its
#' Hessian. Splitting the string on `"_"` is the obvious way and it is
#' wrong: a parameter whose own name contains an underscore
#' (`"log_scale"`) makes `"log_scale_log_scale"` split into four
#' pieces, and taking the first and the last silently yields
#' `c("log", "scale")`. The map is built from the parameter vector, so it does
#' not depend on the characters in the names.
#'
#' [deriv_indices()] is the same idea for orders above two.
#'
#' @param params A character vector of parameter names.
#'
#' @return A named list, parallel to `hess_names(params)`, each element a
#'   character pair.
#'
#' @seealso [hess_names()], [deriv_indices()]
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
#' Generates the names of the unique partial derivatives of a given `order`
#' with respect to a vector of parameters. Because mixed partial derivatives are
#' symmetric, only one representative per multi-index is listed, using
#' non-decreasing parameter order (e.g. `"mu_mu_sigma"` but not
#' `"mu_sigma_mu"`). For `order = 2` this coincides with the set of
#' [hess_names()] (though possibly in a different order).
#'
#' @param params A character vector of parameter names (e.g., `c("mu", "sigma")`).
#' @param order A positive integer, the derivative order (e.g. `3` or `4`).
#'
#' @return A character vector of the \eqn{\binom{p + \text{order} - 1}{\text{order}}}
#'   unique component names, where \eqn{p} is the number of parameters.
#'
#' @examples
#' deriv_names(c("mu", "sigma"), 3)
#' # "mu_mu_mu" "mu_mu_sigma" "mu_sigma_sigma" "sigma_sigma_sigma"
#'
#' @seealso [hess_names()], [expand_params()], [transpose_params()], [check_params_dim()], [check_theta_bounds()], [param_smoothness()], [generate_random_theta()]
#' @export
deriv_names <- function(params, order) {
  vapply(deriv_indices(params, order),
         function(r) paste(params[r], collapse = "_"), character(1))
}


#' Index Tuples Behind the Higher-Order Derivative Names
#'
#' @description
#' The multi-indices that [deriv_names()] names, in exactly the same order:
#' non-decreasing tuples of length `order` over `seq_along(params)`,
#' enumerated lexicographically.
#'
#' @details
#' This function exists so that no caller needs to recover an index tuple by
#' splitting a component name. Splitting fails for a parameter whose own name
#' contains an underscore: `"mu_log_scale_log_scale"` splits into five
#' pieces, and matching those against `params` yields `NA`s. The
#' multivariate families carry such names by construction
#' (`sigma_log_L1`), so the parsing route is wrong for shipped families
#' as well as user-defined ones. Generating the indices and the names from the
#' same enumeration removes the possibility of disagreement.
#'
#' Note that this is not interchangeable with
#' `deriv_index_list()` in `link_scale.R`: that one is ordered to match
#' [hess_names()] at order 2, which puts the diagonal first, whereas
#' this one is lexicographic throughout to match [deriv_names()]. At
#' order 2 use [hess_pairs()].
#'
#' @param params A character vector of parameter names.
#' @param order A positive integer, the derivative order.
#'
#' @return A list of integer vectors, each of length `order`, parallel to
#'   `deriv_names(params, order)`.
#'
#' @seealso [deriv_names()], [hess_pairs()]
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
#' Turns a list of `k` equal-length vectors into a list of `n` vectors of
#' length `k`, one per observation, keeping the names.
#'
#' @param theta A list to be transposed.
#' @return A `list` where each element has been transposed and simplified to an atomic vector.
#'
#' @examples
#' transpose_params(list(mu = c(0, 1), sigma = c(1, 2)))
#'
#' @seealso [deriv_names()], [hess_names()], [expand_params()], [check_params_dim()], [check_theta_bounds()], [param_smoothness()], [generate_random_theta()]
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

#' @title Tail Probabilities on the Log Scale
#' @name tail_helpers
#'
#' @description
#' Three internal helpers for the distribution and quantile functions of the
#' families with a closed form. `log1mexp_neg()` computes
#' \eqn{\log(1 - e^a)} for \eqn{a \le 0}, through `log(-expm1(a))` when
#' \eqn{a > -\log 2} and through `log1p(-exp(a))` otherwise, so that the
#' result keeps its relative accuracy at both ends (Maechler, 2012).
#' `cdf_from_log_tails()` returns the probability that `lower.tail` and
#' `log.p` select from the logarithms of the two tails, both computed by the
#' caller. `quantile_log_tails()` converts a probability given with
#' `lower.tail` and `log.p` into the logarithms of the lower and the upper
#' tail, without forming \eqn{1 - p} or exponentiating a logarithm; a
#' probability outside \eqn{[0, 1]} gives `NaN` for both.
#'
#' @param a A numeric vector of non-positive numbers.
#' @param log_lower,log_upper Numeric vectors, the logarithms of
#'   \eqn{P(Y \le q)} and \eqn{P(Y > q)}.
#' @param p A numeric vector of probabilities, or of their logarithms when
#'   `log.p = TRUE`.
#' @param lower.tail,log.p Logical of length 1, with the meaning they have in
#'   [distrib_cdf()] and [distrib_quantile()].
#'
#' @return `log1mexp_neg()` and `cdf_from_log_tails()` return a numeric vector
#'   of the length of their input; `quantile_log_tails()` returns a list with
#'   the numeric vectors `lower` and `upper`.
#'
#' @references Maechler, M. (2012). Accurately computing
#'   \eqn{\log(1 - \exp(-|a|))}. Vignette of the R package Rmpfr.
#'
#' @examples
#' # The helper keeps the digits that the direct formula loses (it gives -Inf).
#' distributions7:::log1mexp_neg(-1e-20)
#' log(1 - exp(-1e-20))
#'
#' # The upper tail of a probability given on the log scale.
#' distributions7:::quantile_log_tails(log(1e-30), lower.tail = FALSE,
#'                                     log.p = TRUE)
#'
#' @keywords internal
log1mexp_neg <- function(a) {
  # by indexing, so that a NaN argument gives NaN and not NA
  out <- log1p(-exp(a))
  big <- !is.na(a) & a > -log(2)
  out[big] <- log(-expm1(a[big]))
  out
}

#' The Logarithm of minus log(1 - e^p)
#'
#' @description
#' Returns \eqn{\log\{-\log(1 - e^{p})\}} for \eqn{p \le 0} without
#' forming \eqn{1 - e^{p}}: as \eqn{\log(-\log(-\mathrm{expm1}(p)))} where
#' \eqn{p > -\log 2}, and otherwise as
#' \eqn{p + \log\{-\mathrm{log1p}(-u)/u\}}, \eqn{u = e^{p}}, which is
#' \eqn{p} itself once \eqn{u} underflows. With \eqn{p = \log S} it is the
#' \eqn{\log(-\log F)} of the Gumbel quantile; with \eqn{p = \log F} it is
#' the logarithm of the Weibull cumulative hazard at the quantile.
#'
#' @param p A numeric vector of log probabilities.
#'
#' @return A numeric vector of the length of `p`.
#'
#' @examples
#' # at p = -800 the value is -800 to double precision
#' distributions7:::log_neg_log1mexp(c(-0.1, -50, -800))
#'
#' @keywords internal
log_neg_log1mexp <- function(p) {
  u <- exp(p)
  far <- !is.na(p) & p <= -log(2)
  out <- log(-log(-expm1(p)))
  out[far] <- ifelse(u[far] == 0, p[far], p[far] + log(-log1p(-u[far]) / u[far]))
  out
}

#' Lower-Tail Logarithms From the Logarithm of the Argument
#'
#' @description
#' `weibull_log_lower()` returns \eqn{\log F = \log(1 - e^{-H})} from
#' \eqn{\log H}, and `gamma_log_lower()` returns the logarithm of the
#' regularized lower incomplete gamma function \eqn{P(k, w)} from
#' \eqn{\log w}, so that a cumulative hazard or an argument that underflows
#' still gives a finite logarithm. Below \eqn{\log H = -30} the first is
#' \eqn{\log H - H/2}; below \eqn{\log w = -50} the second is
#' \eqn{k\log w - \log\Gamma(k + 1) - w\,k/(k + 1)}, the first terms of
#' the series of \eqn{P}; above, they are `log(-expm1(-H))` and
#' [stats::pgamma()] on the log scale.
#'
#' @param lh,lw Numeric vectors, the logarithm of the cumulative hazard or of
#'   the argument (`-Inf` at zero).
#' @param k The shape of the gamma, recycled.
#'
#' @return A numeric vector of log probabilities.
#'
#' @examples
#' # a hazard of e^-800 and an argument of e^-800, both below the double range
#' distributions7:::weibull_log_lower(c(-3, -800))
#' distributions7:::gamma_log_lower(c(-3, -800), 0.7)
#'
#' @keywords internal
weibull_log_lower <- function(lh) {
  ifelse(lh < -30, lh - exp(lh) / 2, log(-expm1(-exp(lh))))
}

#' @rdname weibull_log_lower
gamma_log_lower <- function(lw, k) {
  n <- max(length(lw), length(k))
  lw <- rep_len(lw, n)
  k <- rep_len(k, n)
  far <- !is.na(lw) & lw < -50
  out <- stats::pgamma(exp(lw), shape = k, rate = 1, log.p = TRUE)
  w <- exp(lw[far])
  out[far] <- k[far] * lw[far] - lgamma(k[far] + 1) - w * k[far] / (k[far] + 1)
  out
}

#' @rdname tail_helpers
cdf_from_log_tails <- function(log_lower, log_upper, lower.tail, log.p) {
  lp <- if (lower.tail) log_lower else log_upper
  if (log.p) lp else exp(lp)
}

#' @rdname tail_helpers
quantile_log_tails <- function(p, lower.tail, log.p) {
  lp <- if (log.p) p else log(p)
  # a probability above one: both tails NaN
  lp[!is.na(lp) & lp > 0] <- NaN
  lq <- log1mexp_neg(lp)
  if (lower.tail) list(lower = lp, upper = lq) else list(lower = lq, upper = lp)
}
