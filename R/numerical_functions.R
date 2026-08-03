#' Numerical Summation of Discrete Series
#'
#' Calculates the sum of a function `f(x)` over a sequence of integers from `start` to `end`.
#' The function is designed to handle finite sums, one-sided infinite series, and
#' doubly infinite series by automatically adapting its summation strategy.
#'
#' @param f A function taking a vector of integers `x` and returning a vector of numeric values.
#'   **Must be vectorized**.
#' @param start Numeric. Starting value. Can be finite, `Inf`, or `-Inf`. Defaults to `0`.
#' @param end Numeric. Ending value. Can be finite, `Inf`, or `-Inf`. Defaults to `Inf`.
#' @param step Integer. Number of terms to calculate in a single vectorized batch. Defaults to `1000`.
#' @param tol Numeric. Tolerance threshold for convergence. Defaults to `1e-10`.
#' @param maxit Integer. Safety limit for the maximum number of batch iterations. Defaults to `1000000`.
#' @param reltol Logical. If `TRUE` (default), uses a hybrid relative tolerance.
#'
#' @details
#' **1. Summation Strategies:**
#' The function automatically detects the domain topology based on `start` and `end`:
#' * **Forward (Standard):** If `start <= end` (e.g., `1` to `Inf`).
#' * **Backward (Reflection):** If `start > end` (e.g., `-1` to `-Inf`), evaluates `f(-x)`.
#' * **Doubly Infinite (Folding):** If `start == -Inf` and `end == Inf`, folds around 0.
#'
#' **2. Speed and Convergence:**
#' It monitors convergence using the sum of absolute values in the current chunk, preventing premature stops on alternating series while maintaining high precision.
#' 
#' **3. Underflow & Divergence Detection:**
#' Includes heuristics to stop early if the sequence starts growing in absolute terms
#' (divergence), or skips up to 50 empty chunks to protect from premature stopping
#' when `f(x)` evaluates exactly to `0` at the start.
#'
#' @return A numeric scalar representing the calculated sum.
#'
#' @export
numerical_series <- function(f, start = 0, end = Inf, step = 10000, tol = 1e-10, maxit = 1000000L, reltol = TRUE) {
  
  # --- Setup Range and Direction ---
  if (is.infinite(start) && start < 0 && is.infinite(end) && end > 0) {
    s_init <- f(0)
    start_internal <- 1
    end_internal <- Inf
    f_internal <- function(x) f(x) + f(-x)
  } else if (end < start) {
    s_init <- 0
    start_internal <- -start
    end_internal <- -end
    f_internal <- function(x) f(-x)
  } else {
    s_init <- 0
    start_internal <- start
    end_internal <- end
    f_internal <- f
  }
  
  s <- 0.0
  it <- 0L
  climbing <- TRUE
  prev_max_abs <- Inf
  divergence_counter <- 0L
  flat_counter <- 0L
  
  while (climbing && it < maxit) {
    it <- it + 1L
    
    # Upper limit of the current block
    upper_limit <- min(start_internal + step - 1, end_internal)
    
    # Evaluate function in vectorized chunk (using ALTREP for the sequence)
    x <- start_internal:upper_limit
    vals <- f_internal(x)
    
    # Single pass to evaluate the chunk sum
    chunk_sum <- sum(vals)
    
    s <- s + chunk_sum
    
    # 1. Explicit divergence check
    if (is.infinite(s) || is.na(s)) {
      warning("The series reached Inf or NaN. Stopping.")
      return(s_init + s)
    }
    
    # 2. Divergence Early-Exit (check only the last term to avoid full array scans)
    last_val_abs <- abs(vals[length(vals)])
    if (last_val_abs > prev_max_abs && last_val_abs > tol) {
      divergence_counter <- divergence_counter + 1L
      if (divergence_counter >= 3L) {
        warning("The series seems to be divergent (terms are growing). Stopping early.")
        return(s_init + s)
      }
    } else {
      divergence_counter <- 0L
    }
    prev_max_abs <- last_val_abs
    
    # 3. Convergence or Underflow Protection
    scaled_tol <- if (reltol) tol * max(abs(s), 1.0) else tol
    
    if (abs(chunk_sum) < scaled_tol) {
      # Short-circuit logic: only perform full vector allocation and absolute sum 
      # if the chunk sum is small, preventing false positives from alternating series.
      chunk_abs_sum <- sum(abs(vals))
      if (chunk_abs_sum < scaled_tol) {
        flat_counter <- flat_counter + 1L
        if (abs(s) > tol || flat_counter >= 50L) {
          climbing <- FALSE
        }
      } else {
        flat_counter <- 0L
      }
    } else {
      flat_counter <- 0L
    }
    
    # 4. Advancement
    if (climbing) {
      start_internal <- upper_limit + 1
      if (start_internal > end_internal) {
        climbing <- FALSE
      }
    }
  }
  
  if (it >= maxit) {
    warning("Maximum number of iterations reached.")
  }
  
  s_init + s
}

#' Expected Value of a Function of a Random Variable
#'
#' @description
#' Computes \eqn{E[f(Y)]} under the distribution: by numerical integration for a
#' continuous distribution and by series summation for a discrete one, with the
#' methods described on their own pages.
#'
#' @param distrib An object inheriting from \code{\link{distrib}}.
#' @param f The function whose expectation is taken, with signature
#'   \code{f(y, theta, ...)}.
#' @param theta A named list of parameters. Vectors are supported and are
#'   recycled against any vectors in \code{...}, so several parameter values can
#'   be handled in one call.
#' @param ... Further arguments passed to \code{f}; their names must not clash
#'   with those of \code{theta}.
#'
#' @return A numeric vector of expected values, one per parameter combination.
#'
#' @importFrom stats integrate
#'
#' @examples
#' expectation(poisson_distrib(), function(y, theta, k = 1) y^k,
#'             theta = list(mu = 10), k = 2)
#'
#' @export
expectation <- S7::new_generic("expectation", "distrib", fun = function(distrib, f, theta, ...) {
  S7::S7_dispatch()
})

#' @title Expectation of a Continuous Distribution
#' @name expectation.continuous_distrib
#' @description Evaluates \eqn{E[f(Y)] = \int f(y)\,p(y;\theta)\,dy} by adaptive quadrature (\code{\link[stats]{integrate}}). The domain is split at the 0.1, 0.5 and 0.9 quantiles of the distribution and each panel is integrated separately, which anchors the quadrature nodes on the probability mass wherever it sits; the panels are then summed.
#' @param distrib A \code{continuous_distrib}.
#' @param f The function whose expectation is taken.
#' @param theta A named list of parameters.
#' @param ... Further arguments passed to \code{f}.
#' @return A numeric vector of expected values.
#' @keywords internal
#' @export
S7::method(expectation, continuous_distrib) <- function(distrib, f, theta, ...) {
  # Capture extra arguments and check for name collisions
  dots <- list(...)
  if (any(names(dots) %in% names(theta))) {
    stop("Arguments in '...' cannot have the same names as parameters in 'theta'.")
  }

  # Combine all parameters to handle vectorization
  all_params <- c(theta, dots)
  n_theta <- length(theta) 

  # Define the worker function for a single set of parameters
  compute_single <- function(params) {
    p_theta <- as.list(params[1:n_theta])
    p_dots <- if (length(params) > n_theta) as.list(params[-(1:n_theta)]) else list()

    integrand <- function(y) {
      val_f <- do.call(f, c(list(y = y, theta = p_theta), p_dots))
      val_p <- distrib_pdf(distrib, y, p_theta, log = FALSE)
      val_f * val_p
    }

    # `integrate` over (-Inf, Inf) transforms the domain around the origin and
    # can silently return 0 when the probability mass sits far from it.
    # Splitting the domain at the 0.1/0.5/0.9 quantiles anchors every panel on
    # the mass. Exactly three knots: more panels extend the range of extreme
    # shapes that integrate at all, but convert loud failures into silent wrong
    # answers and add each panel's quadrature error in the ordinary range.
    b <- distrib@bounds
    qs <- suppressWarnings(distrib_quantile(distrib, c(0.1, 0.5, 0.9), p_theta))
    qs <- unique(qs[is.finite(qs) & qs > b[1] & qs < b[2]])
    knots <- sort(c(b[1], qs, b[2]))

    sum(vapply(
      seq_len(length(knots) - 1L),
      function(k) stats::integrate(
        integrand,
        lower = knots[k], upper = knots[k + 1L]
      )$value,
      numeric(1)
    ))
  }
  unname(sapply(transpose_params(expand_params(all_params)), compute_single))
}

#' @title Expectation of a Discrete Distribution
#' @name expectation.discrete_distrib
#' @description Evaluates \eqn{E[f(Y)] = \sum_y f(y)\,P(Y = y;\theta)} by direct summation over the support, truncating the series once the accumulated tail contribution falls below tolerance.
#' @param distrib A \code{discrete_distrib}.
#' @param f The function whose expectation is taken.
#' @param theta A named list of parameters.
#' @param ... Further arguments passed to \code{f}.
#' @return A numeric vector of expected values.
#' @keywords internal
#' @export
S7::method(expectation, discrete_distrib) <- function(distrib, f, theta, ...) {
  # Capture extra arguments and check for name collisions
  dots <- list(...)
  if (any(names(dots) %in% names(theta))) {
    stop("Arguments in '...' cannot have the same names as parameters in 'theta'.")
  }

  # Combine all parameters to handle vectorization
  all_params <- c(theta, dots)
  n_theta <- length(theta) 

  # Define the worker function for a single set of parameters
  compute_single <- function(params) {
    p_theta <- as.list(params[1:n_theta])
    p_dots <- if (length(params) > n_theta) as.list(params[-(1:n_theta)]) else list()

    integrand <- function(y) {
      val_f <- do.call(f, c(list(y = y, theta = p_theta), p_dots))
      val_p <- distrib_pdf(distrib, y, p_theta, log = FALSE)
      val_f * val_p
    }

    numerical_series(integrand, start = distrib@bounds[1], end = distrib@bounds[2])
  }
  unname(sapply(transpose_params(expand_params(all_params)), compute_single))
}