
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
#' integrated as an exponential tail; in \eqn{\log y} the panels are at most
#' one unit wide, and they restart at the scale \eqn{\theta} where the
#' density's fast decay begins when that point lies more than one unit
#' above the center. For a small \eqn{k} part of the mass lies below the
#' smallest representable double, so on \eqn{(0, \infty)} the part below
#' \eqn{y_0 = \min(q, 10^{-290})} is integrated in closed form wherever it
#' holds at least \eqn{10^{-17}} of \eqn{F(q)}: there the log-density and
#' its parameter derivatives are linear in \eqn{L = \log y}, so they are
#' read at \eqn{10^{-290}} and 15 units of \eqn{L} below it and extended to
#' \eqn{L_0 = \log y_0}, and \eqn{\int_{-\infty}^{L_0} e^{c + k u} g(u)\, du =
#' f(y_0)\, y_0 \sum_j (-1)^j g^{(j)}(L_0)/k^{j+1}} ([compiled_cdf_head()]
#' at orders three and four). When the part of the support above
#' \eqn{q} holds less mass than the part below, the integral runs over the
#' upper part with the opposite sign, so that a derivative that is small
#' beside one is not a sum of terms of order one; the mass is the family's
#' closed distribution function where it has one, and the integral of the
#' density otherwise. The two beta families always integrate the lower part
#' when their second shape is below one. For the families that
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
#' \eqn{-\int_{q}^{hi} f B_I} on the side `cdf_rule_cpp()` chooses, by the
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
  head <- compiled_cdf_head(distrib, rule, th, order, params)
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
    res <- numeric(n)
    if (any_pts) {
      v <- fw * bell_f_ratio(params[I], ell)
      v[!live] <- 0
      res <- sgn * ld_group_sum(as.numeric(v), ix, n)
    }
    if (!is.null(head)) res[head$idx] <- res[head$idx] + head$integral(I)
    res
  })
  stats::setNames(out, deriv_names(params, order))
}

#' Power-Law Head of a Compiled CDF Derivative
#'
#' @description
#' The part of \eqn{\partial_I F(q) = \int_0^q f B_I} below the cut
#' \eqn{y_0} that `cdf_rule_cpp()` reports for a family on
#' \eqn{(0, \infty)} whose density behaves like \eqn{y^{\alpha - 1}} near
#' zero (a gamma or a generalized gamma with a small shape), where part of
#' the mass lies below the smallest representable double.
#'
#' @details
#' On \eqn{(0, y_0]}, with \eqn{L = \log y}, every derivative of \eqn{\ell}
#' in the parameters is linear in \eqn{L}, so it is read at
#' \eqn{10^{-290}} and 15 units of \eqn{L} below it and extended to any
#' \eqn{L}; \eqn{B_I} is then a polynomial in \eqn{L} of degree at most
#' four, and
#' \deqn{\int_{-\infty}^{L_0} e^{c + \alpha u} B_I(u)\, du =
#'   \frac{f(y_0)\, y_0}{\alpha} \int_0^\infty e^{-t}
#'   B_I(L_0 - t/\alpha)\, dt}
#' is taken by the three-point Gauss-Laguerre rule, exact for a polynomial
#' of degree five. The cut, \eqn{\alpha} and \eqn{f(y_0)\,y_0} come from
#' `cdf_rule_cpp()`, which uses the same head at orders one and two.
#'
#' @param distrib A family named in [compiled_cdf_deriv_k()].
#' @param rule The list `cdf_rule_cpp()` returns.
#' @param th The named list of parameter vectors, one value per point.
#' @param order The derivative order, 3 or 4.
#' @param params The parameter names.
#'
#' @return `NULL` when no point has a head; otherwise a list with `idx`, the
#'   points that have one, and `integral`, a function of a multi-index that
#'   returns the head's integral at those points.
#'
#' @keywords internal
compiled_cdf_head <- function(distrib, rule, th, order, params) {
  hi <- rule$head_idx
  if (!length(hi)) return(NULL)
  t <- c(0.41577455678347908, 2.2942803602790417, 6.2899450829374792)
  w <- c(0.71109300992917302, 0.27851773356924085, 0.010389256501586136)
  thh <- lapply(th, `[`, hi)
  e0 <- parent_ell(distrib, rule$head_y0, thh, order, params)
  e1 <- parent_ell(distrib, rule$head_y1, thh, order, params)
  nh <- length(hi)
  # node k of point i at position (k - 1) nh + i, as its distance in L
  # from the reading point: (L0 - La) - t_k / alpha
  off <- rep(rule$head_off, times = 3L) -
    rep(t, each = nh) / rep(rule$head_alpha, times = 3L)
  ell <- function(block) {
    a0 <- e0(block)
    rep(a0, times = 3L) + rep((a0 - e1(block)) / rule$head_dL, times = 3L) * off
  }
  scale <- exp(rule$head_logE) / rule$head_alpha
  list(idx = hi, integral = function(I) {
    b <- matrix(bell_f_ratio(params[I], ell), nrow = nh)
    scale * as.numeric(b %*% w)
  })
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
