#' @include distrib.R binomial_distrib.R
NULL

#' The Scalar Route of a Distribution
#'
#' @description
#' Returns how the package's compiled scalar entry points address a
#' distribution: the name that the C function `d7_scalar_id()` recognizes,
#' and the constants of the distribution that follow its parameters in the
#' parameter vector those entry points read. A consumer that resolves the
#' entry points with `R_GetCCallable()` uses the name to obtain the
#' distribution's identifier once, and appends the constants to the
#' parameters of every observation.
#'
#' @details
#' The entry points are `d7_scalar_id(name)`, `d7_score_curv(id, k, y, th,
#' out)`, which writes the score and the \eqn{(k, k)} second derivative of
#' the log-density in parameter \eqn{k}, and `d7_info_dinfo(id, k, y, th,
#' out)`, which writes the \eqn{(k, k)} expected second derivative and its
#' derivative in the same parameter, all on the parameter scale, with `k`
#' counted from zero over `distrib@params` and `th` the parameters of one
#' observation followed by its constants. `d7_scalar_thread_safe(id)` returns
#' 1 when the entries of a distribution never reach the R API, so that they
#' may be called from a worker thread, and 0 otherwise.
#'
#' The default method returns the class name when the registry covers the
#' class and no constants; the binomial carries its size as a constant.
#'
#' @param distrib A distribution object inheriting from `distrib`.
#' @param ... Passed to methods.
#'
#' @return `NULL` when the registry does not cover the distribution;
#'   otherwise a list with `name`, a single character string, and
#'   `constants`, a named list of numeric vectors, each of length 1 or of the
#'   number of observations, in the order they follow the parameters.
#'
#' @examples
#' distrib_scalar_route(gaussian1_distrib())
#' distrib_scalar_route(binomial_distrib(size = 10))
#' is.null(distrib_scalar_route(mvgaussian1_distrib(n_dim = 2)))
#'
#' @export
distrib_scalar_route <- S7::new_generic(
  "distrib_scalar_route", "distrib",
  function(distrib, ...) S7::S7_dispatch()
)

S7::method(distrib_scalar_route, distrib) <- function(distrib, ...) {
  name <- attr(S7::S7_class(distrib), "name")
  if (!name %in% d7_scalar_classes_covered()) return(NULL)
  list(name = name, constants = list())
}

S7::method(distrib_scalar_route, BinomialDistrib) <- function(distrib, ...) {
  list(name = "BinomialDistrib", constants = list(size = distrib@size))
}
