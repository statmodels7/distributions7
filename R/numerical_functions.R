#' Numerical Summation of Discrete Series (Rcpp Accelerated)
#'
#' Calculates the sum of a function `f(x)` over a sequence of integers from `start` to `end`.
#' The function is designed to handle finite sums, one-sided infinite series, and
#' doubly infinite series by automatically adapting its summation strategy.
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
#' @param f A function taking a vector of integers `x` and returning a vector of numeric values.
#'   **Must be vectorized**.
#' @param start Numeric. Starting value. Can be finite, `Inf`, or `-Inf`. Defaults to `0`.
#' @param end Numeric. Ending value. Can be finite, `Inf`, or `-Inf`. Defaults to `Inf`.
#' @param step Integer. Number of terms to calculate in a single vectorized batch. Defaults to `1000`.
#' @param tol Numeric. Tolerance threshold for convergence. Defaults to `1e-10`.
#' @param maxit Integer. Safety limit for the maximum number of batch iterations. Defaults to `1000000`.
#' @param reltol Logical. If `TRUE` (default), uses a hybrid relative tolerance.
#'
#' @return A numeric scalar representing the calculated sum.
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