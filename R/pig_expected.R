#' @include pig1_distrib.R pig2_distrib.R dexpected_hessian.R expected_loc_scale.R
NULL

#' Poisson-Inverse Gaussian Expected Information
#'
#' @name distrib_expected_hessian.Pig2Distrib
#' @aliases distrib_dexpected_hessian.Pig2Distrib
#'   distrib_d2expected_hessian.Pig2Distrib
#'   distrib_expected_hessian.Pig1Distrib
#'   distrib_dexpected_hessian.Pig1Distrib
#'   distrib_d2expected_hessian.Pig1Distrib
#'
#' @description
#' The expected information of the Poisson-inverse Gaussian and its first and
#' second derivatives in the parameters, computed exactly by one pass over the
#' support for each observation.
#'
#' @details
#' For [pig2_distrib()] the expectations are sums over the support of the mass
#' times products of the observed derivatives, written through the second
#' Bartlett identity with the measure moving:
#' \deqn{\mathbb{E}[\ell_{ab}] = -\mathbb{E}[\ell_a\ell_b],}
#' \deqn{\partial_c\,\mathbb{E}[\ell_{ab}] = -\mathbb{E}[\ell_{ac}\ell_b
#'   + \ell_a\ell_{bc} + \ell_a\ell_b\ell_c],}
#' and the corresponding expression at the second order. The forms in
#' \eqn{\mathbb{E}[\ell_{ab}]} itself sum terms whose mean is of a higher order
#' in \eqn{1/\alpha} than the terms, and lose their digits in the Poisson
#' limit; the forms in the score do not.
#'
#' Each observed derivative reads \eqn{\log S_y(\alpha)} and the moments of the
#' finite sum \eqn{S_y}. Computed afresh at each \eqn{y} they cost \eqn{y}
#' terms apiece; here they come from the Bessel recurrence
#' \eqn{S_{y+1} = S_{y-1} + \{(2y-1)/\alpha\}S_y}, differentiated in
#' \eqn{\alpha} and written on positive quantities, so a whole pass costs a
#' number of terms linear in the length of the support. The pass stops past
#' the mean once the geometric bound on the remaining mass, the mass at
#' \eqn{y} times \eqn{2\sigma\mu}, is below \eqn{10^{-17}} of the mass
#' accumulated.
#'
#' For [pig1_distrib()] the recurrence runs in \eqn{w = 1/(2\alpha) =
#' \sigma/(2\sqrt{1 + 2\sigma\mu})} instead, \eqn{S_{y+1} = S_{y-1} +
#' 2(2y-1)\,w\,S_y}, because the row of that parametrization reads
#' \eqn{\log S_y} and its derivatives in \eqn{w}: all of them are finite as
#' \eqn{\sigma \to 0}, where those in \eqn{\alpha} carry
#' \eqn{\sigma^{-k}} and cancel.
#'
#' The `approx` and `nsim` arguments are accepted and ignored.
#'
#' @param distrib A `Pig1Distrib` or `Pig2Distrib` object.
#' @param y A numeric vector, read for its length.
#' @param theta A named list with `mu` and `sigma` (pig1) or `alpha` (pig2).
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads A single positive integer, the kernel's thread count.
#'
#' @return A named list keyed as [hess_names()], [dexpected_names()] or
#'   [d2expected_names()].
#'
#' @seealso [distrib_expected_hessian()]; `pig1_expected_cpp()` and
#'   `pig2_expected_cpp()` for the kernels.
#'
#' @examples
#' d <- pig2_distrib()
#' th <- list(mu = 3, alpha = 2)
#' distrib_expected_hessian(d, 1, th)
#'
#' # Against the outer product of the scores summed over the support.
#' y <- 0:400
#' f <- distrib_pdf(d, y, th)
#' g <- distrib_gradient(d, y, th)
#' c(kernel = distrib_expected_hessian(d, 1, th)$alpha_alpha,
#'   sum = -sum(f * g$alpha^2))
#' @keywords internal
S7::method(distrib_expected_hessian, Pig2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  pig2_expected_cpp(y, theta[[1]], theta[[2]], 0L, threads)
}

S7::method(distrib_dexpected_hessian, Pig2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                     function(k) pig2_expected_cpp(y, theta[[1]], theta[[2]],
                                                   k, threads))
}

S7::method(distrib_d2expected_hessian, Pig2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  dexpected_analytic(distrib, y, theta, match.arg(scale), 2L, threads,
                     function(k) pig2_expected_cpp(y, theta[[1]], theta[[2]],
                                                   k, threads))
}

S7::method(distrib_expected_hessian, Pig1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  pig1_expected_cpp(y, theta[[1]], theta[[2]], 0L, threads)
}

S7::method(distrib_dexpected_hessian, Pig1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                     function(k) pig1_expected_cpp(y, theta[[1]], theta[[2]],
                                                   k, threads))
}

S7::method(distrib_d2expected_hessian, Pig1Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  dexpected_analytic(distrib, y, theta, match.arg(scale), 2L, threads,
                     function(k) pig1_expected_cpp(y, theta[[1]], theta[[2]],
                                                   k, threads))
}

# The sum over the support is exact and costs far more than the observed
# information; see expected_hessian_costly().
S7::method(expected_hessian_costly, Pig1Distrib) <- function(x, ...) TRUE
S7::method(expected_hessian_costly, Pig2Distrib) <- function(x, ...) TRUE
