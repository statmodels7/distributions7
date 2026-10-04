#' @include distrib.R binomial_distrib.R betabinom1_distrib.R
#'   betabinom2_distrib.R fixed.R
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
#' observation followed by its constants. `d7_logpdf(id, y, th)` returns the
#' log-density of one observation. `d7_scalar_thread_safe(id)` returns
#' 1 when the entries of a distribution never reach the R API, so that they
#' may be called from a worker thread, and 0 otherwise.
#'
#' The default method returns the class name when the registry covers the
#' class and no constants; the binomial and the two beta-binomials carry
#' their size as a constant.
#'
#' A wrapped family is named `"<wrapper class>[:aux]|<inner class>"`, `aux`
#' being what the wrapper needs besides its constants, and its constants are
#' the inner family's followed by the wrapper's. For [fixed()], `aux` is the
#' mask of the fixed parameters (bit `j - 1` for the `j`-th parameter of the
#' inner family) and the wrapper's constants are the fixed values in the
#' inner family's order. A wrapper of a wrapper, or a fixed value that varies
#' by observation, has no route.
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

S7::method(distrib_scalar_route, BetaBinom1Distrib) <- function(distrib, ...) {
  list(name = "BetaBinom1Distrib", constants = list(size = distrib@size))
}

S7::method(distrib_scalar_route, BetaBinom2Distrib) <- function(distrib, ...) {
  list(name = "BetaBinom2Distrib", constants = list(size = distrib@size))
}

#' The Scalar Route of a Fixed Family
#'
#' @description
#' Returns the route of a family built by [fixed()]: the name
#' `"<class>:<mask>|<inner name>"` and the inner family's constants followed
#' by the fixed values, as [distrib_scalar_route()] describes, or `NULL` when
#' the parent has no route of its own or a fixed value varies by
#' observation.
#'
#' @param distrib A fixed family.
#' @param ... Unused.
#'
#' @return A list with components `name` and `constants`, or `NULL`.
#'
#' @keywords internal
fixed_scalar_route <- function(distrib, ...) {
  inner <- distrib_scalar_route(distrib@parent_distrib)
  if (is.null(inner) || grepl("|", inner$name, fixed = TRUE)) return(NULL)
  fp <- distrib@fixed_params
  if (any(lengths(fp) != 1L)) return(NULL)
  P <- distrib@parent_distrib@params
  is_fixed <- P %in% names(fp)
  mask <- sum(2^(which(is_fixed) - 1L))
  cls <- attr(S7::S7_class(distrib), "name")
  list(name = paste0(cls, ":", mask, "|", inner$name),
       constants = c(inner$constants, unname(fp[P[is_fixed]])))
}
S7::method(distrib_scalar_route, FixedContinuousDistrib) <- fixed_scalar_route
S7::method(distrib_scalar_route, FixedDiscreteDistrib) <- fixed_scalar_route
