
#' Distribution Function and Its Derivatives from the Compiled Kernels
#'
#' @description
#' Method bodies for [distrib_cdf()], [distrib_grad_cdf()] and
#' [distrib_hess_cdf()] that read the compiled distribution functions of the
#' package's scalar registry, for the families whose distribution function
#' has no closed derivative in a shape parameter: the gamma (both
#' parametrizations), the chi-squared, the two generalized gammas, the two
#' betas, the two von Mises, the two Student t, the two pseudo-Huber and the
#' skew t.
#'
#' @details
#' A derivative in a shape parameter is the integral of the density's own
#' derivative,
#' \deqn{\partial_k F(q) = \int_{-\infty}^{q} f\, s_k \, dy, \qquad
#'   \partial_{kl} F(q) = \int_{-\infty}^{q} f\,(\ell_{kl} + s_k s_l)\, dy,}
#' with \eqn{s} and \eqn{\ell} the score and the Hessian of the log-density.
#' The integral is taken by a fixed rule (tanh-sinh panels doubling away from
#' the family's center, an exp-sinh tail beyond), in \eqn{\log y} on
#' \eqn{(0, \infty)} and in \eqn{\mathrm{logit}(y)} on \eqn{(0, 1)}, so that
#' a density that behaves like \eqn{y^{k-1}} at an end of its support is
#' integrated as an exponential tail. Above the family's center the integral
#' runs over \eqn{[q, \infty)} with the opposite sign. For the families that
#' are location-scale in \eqn{(\mu, \sigma)} only the shape parameters are
#' integrated: with \eqn{z = (q - \mu)/\sigma} the location and scale
#' components are \eqn{-f}, \eqn{-z f} and, at second order, \eqn{-f s_\mu},
#' \eqn{f(-z^2 s_\mu + 2z/\sigma)} and \eqn{f(-z s_\mu + 1/\sigma)}, and the
#' mixed components with a shape \eqn{k} are \eqn{-f s_k} and
#' \eqn{-z f s_k}, all at \eqn{q}. The distribution function itself is
#' `pt()`, `pgamma()`, `pchisq()` or `pbeta()` where the family has one, and
#' the same integral of the density otherwise. The tail and the log scale
#' are applied as for every family, by the conversion of the lower-tail
#' derivatives.
#'
#' @param distrib A distribution object of one of the families above.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameter values.
#' @param lower.tail Is the lower tail wanted? A single logical.
#' @param log,log.p Are derivatives of (or values of) the log probability
#'   wanted? A single logical.
#' @param ... Unused.
#'
#' @return For `compiled_cdf()`, a numeric vector; for the two derivative
#'   bodies, a named list of numeric vectors as [distrib_grad_cdf()] and
#'   [distrib_hess_cdf()] return.
#'
#' @keywords internal
compiled_cdf <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE,
                         ...) {
  a <- compiled_cdf_args(distrib, q, theta)
  out <- d7_cdf_cpp(a$cls, a$q, a$tm, lower.tail)
  if (log.p) log(out) else out
}

#' @rdname compiled_cdf
compiled_grad_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE,
                              ...) {
  a <- compiled_cdf_args(distrib, q, theta)
  G <- d7_cdf_grad_cpp(a$cls, a$q, a$tm)
  d1 <- stats::setNames(lapply(seq_len(ncol(G)), function(j) G[, j]),
                        distrib@params)
  cdf_tail_scale(distrib, distrib_cdf(distrib, a$q, theta), d1, NULL,
                 lower.tail, log)
}

#' @rdname compiled_cdf
compiled_hess_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE,
                              ...) {
  a <- compiled_cdf_args(distrib, q, theta)
  G <- d7_cdf_grad_cpp(a$cls, a$q, a$tm)
  H <- d7_cdf_hess_cpp(a$cls, a$q, a$tm)
  d1 <- stats::setNames(lapply(seq_len(ncol(G)), function(j) G[, j]),
                        distrib@params)
  d2 <- stats::setNames(lapply(seq_len(ncol(H)), function(j) H[, j]),
                        hess_names(distrib@params))
  cdf_tail_scale(distrib, distrib_cdf(distrib, a$q, theta), d1, d2,
                 lower.tail, log)
}

#' @rdname compiled_cdf
#' @return `compiled_cdf_args()` returns the class name, the recycled
#'   quantiles and the parameter matrix the kernels read.
compiled_cdf_args <- function(distrib, q, theta) {
  n <- max(length(q), lengths(theta))
  th <- lapply(align_theta(distrib, theta)[distrib@params], rep_len,
               length.out = n)
  list(cls = attr(S7::S7_class(distrib), "name"),
       q = rep_len(as.numeric(q), n),
       tm = matrix(unlist(th, use.names = FALSE), nrow = n))
}
