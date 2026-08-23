# The quadrature and summation engines live in numericals7; this file owns
# what is distributions7's business -- the split of the domain at quantiles,
# the handling of the support's endpoints, and the contract expectation()
# offers its integrand.

#' Batched Quadrature with Rejection
#'
#' @description
#' Calls [numericals7::quad_vec()] with its warning muffled and
#' returns the per-row integrals, `NA` where the accuracy was not
#' reached; the caller decides what a failed row means and names it.
#'
#' @param integrand The integrand, in `quad_vec`'s matrix contract.
#' @param lower,upper Numeric vectors of panel endpoints.
#' @return A numeric vector of integrals, `NA` for failed rows.
#' @keywords internal
quad_rows <- function(integrand, lower, upper) {
  withCallingHandlers(
    numericals7::quad_vec(integrand, lower, upper),
    warning = function(w) {
      if (grepl("quad_vec", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Summation over an Integer Support
#'
#' @description
#' Sums `term(k, i)` over the integers of `[from, to]` for each of
#' `n_rows` parameter rows at once. A finite support is summed directly
#' in one matrix evaluation; a support unbounded above goes through
#' [numericals7::series_vec()]; one unbounded below is reflected;
#' one unbounded on both sides is folded around zero. A row whose series
#' does not converge raises an error naming it.
#'
#' @param term A function `term(k, i)` of two equal-length integer
#'   vectors, returning the terms elementwise.
#' @param from,to The endpoints of the support, either possibly infinite.
#' @param n_rows The number of parameter rows.
#' @return A numeric vector of sums, one per row.
#' @keywords internal
discrete_support_sum <- function(term, from, to, n_rows) {
  if (is.finite(from) && is.finite(to)) {
    ks <- seq.int(from, to)
    K <- rep(ks, each = n_rows)
    I <- rep.int(seq_len(n_rows), length(ks))
    return(rowSums(matrix(term(K, I), n_rows, length(ks))))
  }
  if (is.finite(from)) {
    return(series_rows(term, from, n_rows))
  }
  if (is.finite(to)) {
    return(series_rows(function(k, i) term(-k, i), -to, n_rows))
  }
  center <- term(rep.int(0L, n_rows), seq_len(n_rows))
  center + series_rows(function(k, i) term(k, i) + term(-k, i), 1L, n_rows)
}

#' Batched Series Summation with Rejection
#'
#' @description
#' Calls [numericals7::series_vec()] with its warning muffled and
#' promotes a row that did not converge to an error naming it: a series that
#' does not converge is a failure of the request, not a number.
#'
#' @param term A function `term(k, i)` in `series_vec`'s contract.
#' @param from The first summation index.
#' @param n_rows The number of parameter rows.
#' @return A numeric vector of sums, one per row.
#' @keywords internal
series_rows <- function(term, from, n_rows) {
  s <- withCallingHandlers(
    numericals7::series_vec(term, n = n_rows, from = from),
    warning = function(w) {
      if (grepl("series_vec", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (anyNA(s)) {
    stop(sprintf(
      "The series did not converge for parameter combination(s) %s.",
      paste(which(is.na(s)), collapse = ", ")
    ), call. = FALSE)
  }
  s
}

#' Expected Value of a Function of a Random Variable
#'
#' @description
#' Computes \eqn{E[f(Y)]} under the distribution: by numerical integration for a
#' continuous distribution and by series summation for a discrete one, with the
#' methods described on their own pages.
#'
#' @param distrib An object inheriting from [distrib()].
#' @param f The function whose expectation is taken, with signature
#'   `f(y, theta, ...)`. It is evaluated elementwise: `y` arrives as
#'   a numeric vector and every component of `theta`, like every argument
#'   passed through `...`, as a vector of the same length, so `f`
#'   must be vectorized in all of them jointly -- which any expression built
#'   from arithmetic and the `distrib_*` generics already is.
#' @param theta A named list of parameters. Vectors are supported and are
#'   recycled against any vectors in `...`, so several parameter values can
#'   be handled in one call; all combinations share one batched evaluation.
#' @param ... Further arguments passed to `f`; their names must not clash
#'   with those of `theta`.
#'
#' @return A numeric vector of expected values, one per parameter combination.
#'
#' @examples
#' expectation(poisson_distrib(), function(y, theta, k = 1) y^k,
#'             theta = list(mu = 10), k = 2)
#'
#' @seealso [moment()], [variance()], [std_dev()], [skewness()], [kurtosis()]
#' @export
expectation <- S7::new_generic("expectation", "distrib", fun = function(distrib, f, theta, ...) {
  S7::S7_dispatch()
})

#' Aligned Parameter Columns for an Expectation
#'
#' @description
#' Shared preparation for the [expectation()] methods: checks that
#' the names in `...` do not collide with those of `theta`, then
#' expands every component to one aligned column per parameter combination.
#'
#' @param f_env_theta The named list of parameters.
#' @param dots The list of further arguments destined for `f`.
#' @return A list with the theta columns `th`, the dot columns
#'   `dots` and the number of combinations `n`.
#' @keywords internal
expectation_columns <- function(f_env_theta, dots) {
  if (any(names(dots) %in% names(f_env_theta))) {
    stop("Arguments in '...' cannot have the same names as parameters in 'theta'.")
  }
  all_params <- expand_params(c(f_env_theta, dots))
  n_theta <- length(f_env_theta)
  list(
    th = all_params[seq_len(n_theta)],
    dots = all_params[-seq_len(n_theta)],
    n = length(all_params[[1L]])
  )
}

#' @title Expectation of a Continuous Distribution
#' @name expectation.continuous_distrib
#' @description Evaluates \eqn{E[f(Y)] = \int f(y)\,p(y;\theta)\,dy} by the
#' batched adaptive quadrature of [numericals7::quad_vec()]: the
#' panels of every parameter combination are refined in one call, so a vector
#' `theta` costs matrix evaluations rather than one adaptive run per
#' value. The domain of each combination is split at its 0.1, 0.5 and 0.9
#' quantiles, which anchors the quadrature on the probability mass wherever it
#' sits. A combination the batched quadrature rejects -- an integrable
#' endpoint singularity too harsh for bisection -- is rescued by one scalar
#' [stats::integrate()] run, whose extrapolation reaches it; an
#' error naming the combination is raised only when both routes fail.
#' @param distrib A `continuous_distrib`.
#' @param f The function whose expectation is taken.
#' @param theta A named list of parameters.
#' @param ... Further arguments passed to `f`.
#' @return A numeric vector of expected values.
#' @keywords internal
S7::method(expectation, continuous_distrib) <- function(distrib, f, theta, ...) {
  cols <- expectation_columns(theta, list(...))
  b <- distrib@bounds

  # quantile knots for every combination in one elementwise call
  pr <- c(0.1, 0.5, 0.9)
  th_rep <- lapply(cols$th, function(v) rep(v, each = length(pr)))
  qs <- suppressWarnings(distrib_quantile(distrib, rep(pr, times = cols$n), th_rep))
  qm <- matrix(qs, nrow = length(pr))

  lower <- upper <- numeric(0)
  comb <- integer(0)
  for (j in seq_len(cols$n)) {
    kj <- qm[, j]
    kj <- unique(kj[is.finite(kj) & kj > b[1L] & kj < b[2L]])
    knots <- sort(c(b[1L], kj, b[2L]))
    lower <- c(lower, knots[-length(knots)])
    upper <- c(upper, knots[-1L])
    comb <- c(comb, rep.int(j, length(knots) - 1L))
  }

  integrand <- function(x, i) {
    xv <- as.numeric(x)
    idx <- rep(comb[i], times = ncol(x))
    th_e <- lapply(cols$th, function(v) v[idx])
    dt_e <- lapply(cols$dots, function(v) v[idx])
    val_f <- do.call(f, c(list(y = xv, theta = th_e), dt_e))
    val_f * distrib_pdf(distrib, xv, th_e, log = FALSE)
  }

  panels <- quad_rows(integrand, lower, upper)
  out <- as.numeric(rowsum(panels, comb, reorder = FALSE))

  # A combination the batched quadrature refuses gets one scalar rescue
  # through stats::integrate, whose extrapolation reaches integrable
  # endpoint singularities that bisection cannot: a gamma-weighted Hessian
  # component at shape below one behaves like y^(shape-2+k) at zero, where
  # quad_vec would need a bisection depth in the hundreds. The batch stays
  # the fast path for every combination that converges; the rescue costs
  # one adaptive run per refused combination, exactly the old engine's
  # price, and an error is raised only when both routes fail.
  if (anyNA(out)) {
    for (j in which(is.na(out))) {
      th_j <- lapply(cols$th, `[`, j)
      dt_j <- lapply(cols$dots, `[`, j)
      one <- function(y) {
        val_f <- do.call(f, c(list(y = y, theta = lapply(th_j, rep_len, length(y))),
                              lapply(dt_j, rep_len, length(y))))
        val_f * distrib_pdf(distrib, y, th_j, log = FALSE)
      }
      ks <- sort(unique(c(b, lower[comb == j], upper[comb == j])))
      ks <- ks[!is.na(ks)]
      out[j] <- tryCatch(
        sum(vapply(seq_len(length(ks) - 1L), function(k)
          stats::integrate(one, ks[k], ks[k + 1L])$value, numeric(1))),
        error = function(e) NA_real_
      )
    }
  }
  if (anyNA(out)) {
    stop(sprintf(
      "The quadrature did not reach the requested accuracy for parameter combination(s) %s.",
      paste(which(is.na(out)), collapse = ", ")
    ), call. = FALSE)
  }
  unname(out)
}

#' @title Expectation of a Discrete Distribution
#' @name expectation.discrete_distrib
#' @description Evaluates \eqn{E[f(Y)] = \sum_y f(y)\,P(Y = y;\theta)} by
#' summation over the support, every parameter combination in one batched
#' pass: a finite support is summed exactly in a single matrix evaluation,
#' an infinite one through [numericals7::series_vec()], whose rows
#' retire as they converge.
#' @param distrib A `discrete_distrib`.
#' @param f The function whose expectation is taken.
#' @param theta A named list of parameters.
#' @param ... Further arguments passed to `f`.
#' @return A numeric vector of expected values.
#' @keywords internal
S7::method(expectation, discrete_distrib) <- function(distrib, f, theta, ...) {
  cols <- expectation_columns(theta, list(...))
  b <- distrib@bounds

  term <- function(k, i) {
    th_e <- lapply(cols$th, function(v) v[i])
    dt_e <- lapply(cols$dots, function(v) v[i])
    val_f <- do.call(f, c(list(y = k, theta = th_e), dt_e))
    val_f * distrib_pdf(distrib, k, th_e, log = FALSE)
  }

  unname(discrete_support_sum(term, b[1L], b[2L], cols$n))
}
