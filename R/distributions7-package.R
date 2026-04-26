## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib distributions7, .registration = TRUE
## usethis namespace: end
NULL

.onLoad <- function(...) {
  S7::methods_register()
}
