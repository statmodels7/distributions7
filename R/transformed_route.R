#' @include transformations.R scalar_route.R
NULL

#' The Scalar Code of a Ready-Made Transformer
#'
#' @description
#' Identifies a transformer built by one of the twelve ready-made
#' constructors and returns its code in the compiled registry and its
#' parameters. A candidate is rebuilt from the parameters found in the
#' transformer's closures, and the transformer is recognized when every one
#' of its functions has the candidate's code and, where it is a closure of
#' the constructor, the transformer's own environment; any other transformer
#' has no code.
#'
#' @param tr A [transformer()] object.
#'
#' @return `NULL`, or a list with `code`, a single integer, and `par`, a
#'   named list of the transformer's parameters.
#'
#' @keywords internal
transformer_scalar_code <- function(tr) {
  e <- environment(tr@trans_abs_jac)
  if (is.null(e)) return(NULL)
  kinds <- list(
    list(1L, character(), function() log_transform()),
    list(2L, character(), function() exp_transform()),
    list(3L, character(), function() inverse_transform()),
    list(4L, character(), function() sqrt_transform()),
    list(5L, "p", function(p) power_transform(p)),
    list(6L, character(), function() asinh_transform()),
    list(7L, "lambda", function(lambda) bc_transform(lambda)),
    list(8L, "lambda", function(lambda) yj_transform(lambda)),
    list(9L, c("loc", "scale"), function(loc, scale) affine_transform(loc, scale)),
    list(10L, character(), function() logit_transform()),
    list(11L, character(), function() expit_transform()),
    list(12L, "a", function(a) softplus_transform(a))
  )
  slots <- c("trans_fun", "trans_inv", "trans_abs_jac", "trans_inv_hessian",
             "grad_log_jac", "hess_log_jac", "bounds_fun", "valid_support")
  for (kd in kinds) {
    par <- lapply(kd[[2]], function(nm) get0(nm, envir = e, inherits = FALSE))
    names(par) <- kd[[2]]
    if (!all(vapply(par, function(v) is.numeric(v) && length(v) == 1L &&
                      is.finite(v), logical(1)))) next
    cand <- do.call(kd[[3]], par)
    ce <- environment(cand@trans_abs_jac)
    same <- identical(tr@decreasing, cand@decreasing)
    for (sl in slots) {
      if (!same) break
      f <- S7::prop(tr, sl)
      g <- S7::prop(cand, sl)
      same <- identical(deparse(f), deparse(g)) &&
        if (identical(environment(g), ce)) identical(environment(f), e)
        else identical(f, g)
    }
    if (same) return(list(code = kd[[1]], par = par))
  }
  NULL
}

#' The Scalar Route of a Transformed Family
#'
#' @description
#' Returns the route of a family built by [transformation()]: the name
#' `"TransformedDistrib:<code>|<inner name>"` and the inner family's
#' constants followed by the transformer's parameters, as
#' [distrib_scalar_route()] describes, or `NULL` when the parent has no route
#' of its own or the transformer is not one of the ready-made ones.
#'
#' @param distrib A transformed family.
#' @param ... Unused.
#'
#' @return A list with components `name` and `constants`, or `NULL`.
#'
#' @keywords internal
transformed_scalar_route <- function(distrib, ...) {
  inner <- distrib_scalar_route(distrib@parent_distrib)
  if (is.null(inner) || grepl("|", inner$name, fixed = TRUE)) return(NULL)
  tc <- transformer_scalar_code(distrib@transformer)
  if (is.null(tc)) return(NULL)
  list(name = paste0("TransformedDistrib:", tc$code, "|", inner$name),
       constants = c(inner$constants, tc$par))
}
S7::method(distrib_scalar_route, TransformedDistrib) <- transformed_scalar_route
