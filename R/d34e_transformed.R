#' @include transformations.R dexpected_families.R d34expected_analytic.R
NULL

#' The Third and Fourth Derivatives of a Transformed Family's Expected Information
#'
#' @description
#' Returns the parent family's third (`k = 3`) or fourth (`k = 4`) derivatives
#' of the expected information on the parameter scale. The function is
#' registered through [register_d34expected()], which carries the derivatives
#' to the link scale with the wrapper's links.
#'
#' @details
#' A transformation of the response leaves the expected information unchanged
#' as a function of the parameters: the log-density of \eqn{Y = g(X)} differs
#' from the parent's at \eqn{g^{-1}(y)} by \eqn{\log|d g^{-1}/dy|}, which does
#' not depend on the parameters, and \eqn{Y} and \eqn{X} have the same
#' distribution of the score. The expected Hessian
#' [distrib_expected_hessian()] of the wrapper is accordingly the parent's, and
#' so are its derivatives in the parameters at every order. The response is
#' carried to the parent's scale by the inverse transformation, as
#' [distrib_expected_hessian()] does, and is read only for its length.
#'
#' @param d A transformed family, from [transformation()].
#' @param y The response, read for its length.
#' @param th The parameters.
#' @param k The order, `3L` or `4L`.
#' @param t The thread count passed to the parent's methods.
#'
#' @return A named list keyed as [d3expected_names()] (`k = 3`) or
#'   [d4expected_names()] (`k = 4`).
#'
#' @keywords internal
transformed_d34expected <- function(d, y, th, k, t) {
  y <- d@transformer@trans_inv(y)
  if (k == 3L) {
    return(distrib_d3expected_hessian(d@parent_distrib, y, th,
                                      scale = "parameter", threads = t))
  }
  distrib_d4expected_hessian(d@parent_distrib, y, th, scale = "parameter",
                             threads = t)
}
register_d34expected(TransformedDistrib, transformed_d34expected)
