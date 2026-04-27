#' Numerical Summation of Discrete Series (Rcpp Accelerated)
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
#' Powered by Rcpp. The summation is performed in vectorized blocks of size `step`.
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
numerical_series <- function(f, start = 0, end = Inf, step = 1000, tol = 1e-10, maxit = 1000000L, reltol = TRUE) {
  
  # --- Setup Range and Direction in R to minimize C++ context switches ---
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
  
  s_init + series_cpp(
    f = f_internal, start = as.numeric(start_internal), end = as.numeric(end_internal), 
    step = as.integer(step), tol = as.numeric(tol), 
    maxit = as.integer(maxit), reltol = as.logical(reltol)
  )
}

#' Calculate the Expected Value of a Function
#'
#' Computes the expected value of a given function \eqn{f(y)} with respect to a probability distribution defined by \code{distrib}.
#' It automatically handles continuous distributions (via numerical integration) and discrete distributions (via series summation).
#'
#' @param distrib An object of class \code{"distrib"}
#' @param f A function representing the transformation of the random variable \eqn{y}.
#'   **Signature:** It must accept arguments \code{y}, \code{theta}, and \code{...} (see Details).
#' @param theta A named list of parameters for the distribution (e.g., \code{list(mu=10, sigma=2)}).
#'   Vectors inside this list allow computing expectations for multiple distribution parametrizations at once.
#' @param ... Additional arguments passed directly to the function \code{f}.
#'   **Vectorization:** These arguments are fully vectorized. If vectors are provided, they are recycled
#'   against the parameters in \code{theta} according to standard R recycling rules.
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item \eqn{E[f(Y)] = \int_{lb}^{ub} f(y, \theta, \dots) \cdot p(y|\theta) \, dy} (Continuous)
#'   \item \eqn{E[f(Y)] = \sum_{y=lb}^{ub} f(y, \theta, \dots) \cdot P(y|\theta)} (Discrete)
#' }
#'
#' **Vectorization:**
#' The function iterates over the longest vector found among \code{theta} and \code{...}.
#' For example, if \code{theta$mu} has length 2 and you pass a vector of length 2 to \code{...},
#' the function computes the expectation for the paired values. If lengths differ, standard R recycling applies.
#'
#' **Requirements for `f`:**
#' The user-provided function \code{f} must be defined with the signature:
#' \code{f(y, theta, ...)}
#'
#' @return A numeric vector containing the expected values. The length corresponds to the
#'   maximum length among all vectors in \code{theta} and \code{...}.
#'
#' @importFrom stats integrate
#'
#' @examples
#' \dontrun{
#' distrib <- poisson_distrib()
#'
#' # Define f accepting y, theta, and extra parameter gamma
#' f_pow <- function(y, theta, gamma = 1) {
#'   y^gamma
#' }
#'
#' # --- Example 1: Basic usage ---
#' expectation(distrib, f_pow, theta = list(mu = 10), gamma = 2)
#' }
#'
#' @export
expectation <- S7::new_generic("expectation", "distrib", fun = function(distrib, f, theta, ...) {
  S7::S7_dispatch()
})

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

    stats::integrate(integrand, lower = distrib@bounds[1], upper = distrib@bounds[2])$value
  }
  unname(sapply(transpose_params(expand_params(all_params)), compute_single))
}

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