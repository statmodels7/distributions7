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
    bad_lens <- len_theta[mismatch_idx]

    error_msg <- paste0(
      "Parameter dimension mismatch. All parameters should have length 1 or ", n, ".\n"
    )

    stop(error_msg, call. = FALSE)
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


#' Transpose and Simplify Parameter List Structure
#'
#' @description
#' Transposes a list structure (swapping "columns" and "rows") and simplifies the
#' inner elements into atomic vectors.
#'
#' @param theta A list to be transposed.
#' @return A \code{list} where each element has been transposed and simplified to an atomic vector.
#'
#' @importFrom purrr transpose
#' @export
transpose_params <- function(theta) {
  lapply(purrr::transpose(theta), unlist)
}