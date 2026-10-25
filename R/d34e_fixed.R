#' @include fixed.R dexpected_families.R d34expected_analytic.R
NULL

#' The Third and Fourth Derivatives of a Fixed Family's Expected Information
#'
#' @description
#' Returns the parent family's third (`k = 3`) or fourth (`k = 4`) derivatives
#' of the expected information at the full parameter vector, subset to the free
#' parameters, on the parameter scale. The function is registered through
#' [register_d34expected()], which carries the derivatives to the link scale
#' with the wrapper's links.
#'
#' @details
#' The expected information of the wrapper is the parent's submatrix over the
#' free parameters, evaluated at the spliced vector, and the fixed values do
#' not move. A derivative in free parameters of an entry in free parameters is
#' therefore the parent's component with every index free. The free set keeps
#' the parent's order, so the subset by [d3expected_names()] or
#' [d4expected_names()] of the free parameters pairs each component with the
#' same indices as in the parent. The parent's own route is used, closed form
#' where the parent has one and the stencil of
#' [numerical_dexpected_higher()] otherwise, so the accuracy of each component
#' is the parent's. A fully fixed family has no free parameter and returns an
#' empty list.
#'
#' @param d A fixed family, from [fixed()].
#' @param y The response, read for its length.
#' @param th The free parameters.
#' @param k The order, `3L` or `4L`.
#' @param t The thread count passed to the parent's methods.
#'
#' @return A named list keyed as [d3expected_names()] (`k = 3`) or
#'   [d4expected_names()] (`k = 4`) over the free parameters.
#'
#' @keywords internal
fixed_d34expected <- function(d, y, th, k, t) {
  if (d@n_params == 0L) return(list())
  parent <- d@parent_distrib
  full <- fixed_full_theta(d, th)
  if (k == 3L) {
    res <- distrib_d3expected_hessian(parent, y, full, scale = "parameter",
                                      threads = t)
    return(res[d3expected_names(d@params)])
  }
  res <- distrib_d4expected_hessian(parent, y, full, scale = "parameter",
                                    threads = t)
  res[d4expected_names(d@params)]
}
register_d34expected(FixedContinuousDistrib, fixed_d34expected)
register_d34expected(FixedDiscreteDistrib, fixed_d34expected)
