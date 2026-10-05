
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


#' Third and Fourth Derivatives of a Compiled Distribution Function
#'
#' @description
#' The derivatives of \eqn{F} of order three or four in the parameters, for
#' the families whose compiled distribution function integrates its shape
#' derivatives ([compiled_cdf()]): gamma1, gamma2, chisq, both generalized
#' gammas, both betas, both von Mises, both Student t, both pseudo-Huber and
#' the skew t.
#'
#' @details
#' With \eqn{f = e^{\ell}}, the derivative of \eqn{f} in a multi-index
#' \eqn{I} is \eqn{f B_I}, where \eqn{B_I} is the complete Bell polynomial
#' in the derivatives of \eqn{\ell} ([bell_f_ratio()]), so
#' \eqn{\partial_I F(q) = \int_{lo}^{q} f B_I}, taken as
#' \eqn{-\int_{q}^{hi} f B_I} on the side [cdf_rule_cpp()] chooses, by the
#' nodes and weights of the rule of the first two orders. The family's own
#' third and fourth derivatives of \eqn{\ell} supply \eqn{B_I}.
#'
#' A family location-scale in its first two parameters (both Student t,
#' both pseudo-Huber, the skew t) reads its components in the location and
#' the scale at \eqn{q}: with \eqn{\partial_\mu F = -f(q)} and
#' \eqn{\partial_\sigma F = -(q - \mu) f(q)/\sigma},
#' \eqn{\partial_{\mu J} F = -f(q) B_J} and, for a multi-index with
#' \eqn{m} scale indices and shape indices \eqn{S},
#' \eqn{\sum_{j=0}^{m-1} \binom{m-1}{j} c_j f(q) B_{\sigma^{m-1-j} S}}
#' with \eqn{c_j = -(q - \mu)(-1)^j j!/\sigma^{j+1}}.
#'
#' @param distrib A family named in the description.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The derivative order, 3 or 4.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} on the
#'   natural scale and the lower tail, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names].
#'
#' @seealso [numerical_cdf_deriv_k()], the stencil this replaces for these
#'   families; [continuous_cdf_deriv_k()], which chooses between them.
#'
#' @keywords internal
compiled_cdf_deriv_k <- function(distrib, q, theta, order) {
  params <- distrib@params
  a <- compiled_cdf_args(distrib, q, theta)
  n <- length(a$q)
  th <- stats::setNames(lapply(seq_along(params), function(j) a$tm[, j]),
                        params)
  rule <- cdf_rule_cpp(a$cls, a$q, a$tm)
  ix <- rule$idx
  thp <- lapply(th, `[`, ix)
  any_pts <- length(ix) > 0L
  fw <- if (any_pts) exp(distrib_pdf(distrib, rule$y, thp, log = TRUE)) * rule$w
  live <- if (any_pts) is.finite(fw) & fw != 0
  ell <- if (any_pts) parent_ell(distrib, rule$y, thp, order, params)
  sgn <- ifelse(rule$upper, -1, 1)
  ls <- a$cls %in% compiled_ls_classes
  if (ls) {
    fq <- distrib_pdf(distrib, a$q, th)
    ellq <- parent_ell(distrib, a$q, th, order - 1L, params)
    bq <- function(J) if (!length(J)) 1 else bell_f_ratio(params[J], ellq)
  }
  out <- lapply(deriv_indices(params, order), function(I) {
    if (ls && any(I <= 2L)) {
      if (any(I == 1L)) return(-fq * bq(I[-match(1L, I)]))
      m <- sum(I == 2L)
      S <- I[I != 2L]
      acc <- 0
      for (j in 0:(m - 1L)) {
        cj <- -(a$q - th[[1L]]) * (-1)^j * factorial(j) / th[[2L]]^(j + 1)
        acc <- acc + choose(m - 1L, j) * cj * fq * bq(c(rep(2L, m - 1L - j), S))
      }
      return(acc)
    }
    if (!any_pts) return(numeric(n))
    v <- fw * bell_f_ratio(params[I], ell)
    v[!live] <- 0
    sgn * ld_group_sum(as.numeric(v), ix, n)
  })
  stats::setNames(out, deriv_names(params, order))
}

#' Families With a Compiled CDF Quadrature
#'
#' @description
#' `compiled_quad_classes` names the classes whose compiled distribution
#' function integrates its shape derivatives, which [compiled_cdf_deriv_k()]
#' serves at orders three and four; `compiled_ls_classes` names those among
#' them that are location-scale in their first two parameters.
#'
#' @return Character vectors of S7 class names.
#'
#' @keywords internal
compiled_quad_classes <- c(
  "Gamma1Distrib", "Gamma2Distrib", "ChisqDistrib", "GenGamma1Distrib",
  "GenGamma2Distrib", "Beta1Distrib", "Beta2Distrib", "VonMises1Distrib",
  "VonMises2Distrib", "StudentT1Distrib", "StudentT2Distrib",
  "PseudoHuberDistrib", "PseudoHuber2Distrib", "SkewTDistrib")
#' @rdname compiled_quad_classes
compiled_ls_classes <- c("StudentT1Distrib", "StudentT2Distrib",
                         "PseudoHuberDistrib", "PseudoHuber2Distrib",
                         "SkewTDistrib")

#' Third and Fourth CDF Derivatives of a Continuous Family
#'
#' @description
#' Returns the derivatives of \eqn{F} of order three or four by
#' [compiled_cdf_deriv_k()] for the families it covers, and by the stencil of
#' [numerical_cdf_deriv_k()] otherwise.
#'
#' @inheritParams compiled_cdf_deriv_k
#'
#' @return A named list of numeric vectors, as [compiled_cdf_deriv_k()].
#'
#' @keywords internal
continuous_cdf_deriv_k <- function(distrib, q, theta, order) {
  if (attr(S7::S7_class(distrib), "name") %in% compiled_quad_classes) {
    compiled_cdf_deriv_k(distrib, q, theta, order)
  } else {
    numerical_cdf_deriv_k(distrib, q, theta, order)
  }
}
