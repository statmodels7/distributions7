## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib distributions7, .registration = TRUE
## usethis namespace: end
NULL

#' Register the Package's S7 Methods on Load
#'
#' @description
#' Calls \code{S7::methods_register()}, which is what makes methods registered on
#' generics owned by other packages take effect.
#'
#' @details
#' It matters here for the S3 generics the package extends -- \code{print},
#' \code{plot}, \code{mean}, \code{simulate} -- whose owners are \pkg{base} and
#' \pkg{stats}. Standard R load hook; not called directly.
#'
#' @param ... Ignored; R calls the hook with the library path and package name.
#'
#' @return Called for its side effect; the return value is discarded by R.
#'
#' @keywords internal
#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}
