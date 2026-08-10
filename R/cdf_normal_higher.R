#' @include cdf_survival_higher.R
NULL

# The families whose distribution function is a sum of normal tails.
#
# Both the inverse gaussian and the elastic net have
#
#   F = const + sum_k s_k exp(w_k) Phi(x_k),
#
# so every derivative is a Leibniz split between the weight and the tail, with
# a Faa di Bruno on each side. A family supplies, per term, its sign, the
# partial derivatives of the LOG weight and those of the argument, and gets
# four orders. The weight is never formed on its own: it is combined with the
# tail on the log scale, `exp(w + pnorm(x, log.p = TRUE))`, because the
# inverse gaussian's `exp(2/(phi mu))` overflows exactly where its `Phi(b)`
# underflows.

#' Derivatives of the Standard Normal Distribution Function
#'
#' @description
#' Returns the factor \eqn{h_m} in
#' \eqn{\Phi^{(m)}(x) = h_m(x)\varphi(x)} for \eqn{m = 1, \ldots, 4}, the
#' Hermite polynomials.
#'
#' @param x A numeric vector.
#' @param m The order, 1 to 4.
#'
#' @return A numeric vector.
#'
#' @keywords internal
phi_hermite <- function(x, m) {
  switch(m,
    rep_len(1, length(x)),
    -x,
    x^2 - 1,
    3 * x - x^3
  )
}

#' CDF Derivatives of a Sum of Weighted Normal Tails
#'
#' @description
#' Evaluates \eqn{\partial^{I}F} for every component of the requested order,
#' for \eqn{F = c_0 + \sum_k s_k e^{w_k}\Phi(x_k)}.
#'
#' @details
#' The Leibniz rule splits the positions of \eqn{I} between the weight and the
#' tail, \eqn{\partial^{S}e^{w} = e^{w}B_{S}(w)} being the complete Bell
#' polynomial in the partials of the log weight and \eqn{\partial^{T}\Phi(x)}
#' one Faa di Bruno pass over \eqn{x}. Nothing is transcribed: both sums run on
#' the package's own partition enumeration.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param order The derivative order, 1 to 4.
#' @param terms A list of terms, each a list with \code{sign}, \code{logw},
#'   \code{wderiv}, \code{x} and \code{xderiv}.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso \code{\link{bell_f_ratio}}
#' @keywords internal
phi_terms_cdf_deriv_k <- function(distrib, q, order, terms) {
  params <- distrib@params
  n <- length(q)

  one_term <- function(tm, I) {
    # exp(w) times the tail and times the density, both formed on the log
    # scale so that a large weight and a small tail cancel before either is
    # exponentiated
    ewPhi <- exp(tm$logw + stats::pnorm(tm$x, log.p = TRUE))
    ewphi <- exp(tm$logw + stats::dnorm(tm$x, log = TRUE))
    pos <- seq_along(I)
    acc <- rep(0, n)
    for (mask in seq_len(bitwShiftL(1L, length(pos))) - 1L) {
      take <- pos[bitwAnd(mask, bitwShiftL(1L, pos - 1L)) > 0L]
      rest <- setdiff(pos, take)
      wfac <- if (!length(take)) 1 else bell_f_ratio(params[I[take]], tm$wderiv)
      tail <- if (!length(rest)) {
        ewPhi
      } else {
        s <- rep(0, n)
        for (part in index_partitions(params[I[rest]])) {
          term <- ewphi * phi_hermite(tm$x, length(part))
          for (b in part) term <- term * tm$xderiv(b)
          s <- s + term
        }
        s
      }
      acc <- acc + wfac * tail
    }
    tm$sign * acc
  }

  idx <- deriv_indices(params, order)
  out <- lapply(idx, function(I) {
    acc <- rep(0, n)
    for (tm in terms) acc <- acc + one_term(tm, I)
    acc
  })
  stats::setNames(out, deriv_names(params, order))
}

#' Register the Four CDF Derivative Orders of a Normal-Tail Family
#'
#' @description
#' Turns a function returning the terms into the four methods.
#'
#' @param cls The S7 class.
#' @param term_fn A function of \code{(distrib, q, theta)} returning the term
#'   list \code{\link{phi_terms_cdf_deriv_k}} consumes.
#'
#' @return Invisibly \code{NULL}; called for the registration.
#'
#' @keywords internal
register_phi_terms_cdf <- function(cls, term_fn) {
  make <- function(o) {
    force(o)
    function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
      tms <- term_fn(distrib, q, theta)
      tabs <- lapply(seq_len(o),
                     function(k) phi_terms_cdf_deriv_k(distrib, q, k, tms))
      cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), tabs, o,
                  lower.tail, log)
    }
  }
  S7::method(distrib_grad_cdf, cls) <- make(1L)
  S7::method(distrib_hess_cdf, cls) <- make(2L)
  S7::method(distrib_deriv3_cdf, cls) <- make(3L)
  S7::method(distrib_deriv4_cdf, cls) <- make(4L)
  invisible(NULL)
}

#' Partial Derivatives of a Product of a Location Term and a Scale Term
#'
#' @description
#' Builds the evaluator for a quantity \eqn{U(\mu)V(\phi)}, whose mixed partial
#' is the product of the two one-variable derivatives.
#'
#' @param nm The two parameter names, in order.
#' @param uderiv A function of the order returning \eqn{U^{(j)}}.
#' @param vderiv A function of the order returning \eqn{V^{(k)}}.
#'
#' @return A function of a character vector of parameter names.
#'
#' @keywords internal
separable_deriv <- function(nm, uderiv, vderiv) {
  function(block) {
    uderiv(sum(block == nm[1L])) * vderiv(sum(block == nm[2L]))
  }
}

#' @title Inverse-Gaussian Log-CDF Derivatives
#' @name distrib_grad_cdf.InvGauss1Distrib
#' @description
#' Closed form at every order. The distribution function is
#' \eqn{\Phi(a) + e^{c}\Phi(b)} with
#' \eqn{a = (q/\mu - 1)/\sqrt{\phi q}}, \eqn{b = -(q/\mu + 1)/\sqrt{\phi q}}
#' and \eqn{c = 2/(\phi\mu)}. Each of the three is a product of a function of
#' the mean and a function of the dispersion, so its mixed partial derivatives
#' are products of one-variable ones, and the two terms are then a Leibniz
#' split between the weight and the tail. The weight and the tail are combined
#' on the log scale: \eqn{e^{c}} overflows exactly where \eqn{\Phi(b)}
#' underflows.
#' @param distrib An \code{InvGauss1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso \code{\link{invgauss1_distrib}}
#' @keywords internal
register_phi_terms_cdf(InvGauss1Distrib, function(distrib, q, theta) {
  mu <- theta[[1]]
  ph <- theta[[2]]
  nm <- distrib@params
  n <- length(q)
  rt <- sqrt(ph * q)
  a <- (q / mu - 1) / rt
  b <- -(q / mu + 1) / rt
  cval <- 2 / (ph * mu)

  # A = q/mu - 1 and C = q/mu + 1 differ only at order zero
  A <- function(j) if (j == 0L) q / mu - 1 else q * (-1)^j * factorial(j) / mu^(1 + j)
  C <- function(j) if (j == 0L) q / mu + 1 else q * (-1)^j * factorial(j) / mu^(1 + j)
  # B = (phi q)^(-1/2)
  B <- function(k) {
    coef <- if (k == 0L) 1 else prod(-0.5 - seq_len(k) + 1)
    coef * q^-0.5 * ph^(-0.5 - k)
  }
  # 1/mu and 1/phi, for c = 2/(phi mu)
  M <- function(j) (-1)^j * factorial(j) / mu^(1 + j)
  P <- function(k) (-1)^k * factorial(k) / ph^(1 + k)

  list(
    list(sign = 1, logw = rep(0, n), wderiv = function(block) rep(0, n),
         x = a, xderiv = separable_deriv(nm, A, B)),
    list(sign = 1, logw = rep_len(cval, n),
         wderiv = separable_deriv(nm, function(j) 2 * M(j), P),
         x = b,
         xderiv = separable_deriv(nm, function(j) -C(j), B))
  )
})


#' Derivatives of a Power of an Affine Argument
#'
#' @description
#' Returns \eqn{d^{k}u^{p}/dv^{k}} for \eqn{u = v} or \eqn{u = 1 - v}, the two
#' shapes the elastic net's scale and its argument are built from.
#'
#' @param u The base, already evaluated.
#' @param p The exponent.
#' @param k The derivative order.
#' @param inner The derivative of the base in the variable, 1 or -1.
#'
#' @return A numeric vector.
#'
#' @keywords internal
dpow_affine <- function(u, p, k, inner) {
  if (k == 0L) return(u^p)
  inner^k * prod(p - seq_len(k) + 1) * u^(p - k)
}

#' @title Elastic-Net Log-CDF Derivatives
#' @name distrib_grad_cdf.EnetDistrib
#' @description
#' Closed form at every order. Each half of the distribution function is a
#' truncated Gaussian, so with \eqn{z = q - \mu}, \eqn{s = \sqrt{c}} and
#' \eqn{x = a/\sqrt{c}} it is \eqn{e^{w}\Phi(X)} below the location and
#' \eqn{1 - e^{w}\Phi(X)} above it, for \eqn{X = \pm sz - x} and a weight
#' \eqn{w = -\log M(x) + x^{2}/2 + \mathrm{const}} written through the Mills
#' ratio the family already carries. Both \eqn{s} and \eqn{x} are products of
#' a function of \eqn{\lambda} and a function of \eqn{\alpha}, so their mixed
#' partial derivatives are products of one-variable ones.
#'
#' The location is the non-regular direction, as in the Laplace the family
#' contains: the second derivative in \eqn{\mu} carries a point mass at
#' \eqn{q = \mu}, and the formulas below hold on either side of it.
#' @param distrib An \code{EnetDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso \code{\link{enet_distrib}}
#' @keywords internal
register_phi_terms_cdf(EnetDistrib, function(distrib, q, theta) {
  p <- .enet_g_higher(.enet_parts(theta))
  nm <- distrib@params
  n <- length(q)
  lam <- p$lam
  al <- p$al
  x0 <- p$x
  s <- sqrt(p$c)
  z <- q - p$mu
  eps <- ifelse(z <= 0, 1, -1)

  # s = sqrt(lambda) sqrt(1 - alpha) and x = sqrt(lambda) alpha / sqrt(1 - alpha)
  sqrt_lam <- function(j) dpow_affine(lam, 0.5, j, 1)
  v_s <- function(k) dpow_affine(1 - al, 0.5, k, -1)
  v_x <- function(k) {
    # d^k of alpha (1 - alpha)^(-1/2), by the product rule. The second term is
    # guarded rather than multiplied by zero: R evaluates it either way, and
    # dpow_affine has no order below zero.
    out <- al * dpow_affine(1 - al, -0.5, k, -1)
    if (k > 0L) out <- out + k * dpow_affine(1 - al, -0.5, k - 1L, -1)
    out
  }
  # a block naming the location contributes nothing to either
  sep <- function(vfun) function(block) {
    if (any(block == nm[1L])) return(rep(0, n))
    rep_len(sqrt_lam(sum(block == nm[2L])) * vfun(sum(block == nm[3L])), n)
  }
  ds <- sep(v_s)
  dx <- sep(v_x)

  # w(x) = -log M(x) + x^2/2 + const, so its derivatives in x are written
  # through G = dlogM/dx, which the family already computes
  wx <- list(-p$g + x0, -p$dg + 1, -p$d2g, -p$d3g)
  wderiv <- function(block) {
    if (any(block == nm[1L])) return(rep(0, n))
    acc <- rep(0, n)
    for (part in index_partitions(block)) {
      term <- rep_len(wx[[length(part)]], n)
      for (b in part) term <- term * dx(b)
      acc <- acc + term
    }
    acc
  }

  # X = eps s (q - mu) - x, linear in the location
  xderiv <- function(block) {
    nmu <- sum(block == nm[1L])
    if (nmu >= 2L) return(rep(0, n))
    rest <- block[block != nm[1L]]
    if (nmu == 1L) {
      if (!length(rest)) return(-eps * rep_len(s, n))
      return(-eps * ds(rest))
    }
    eps * z * ds(block) - dx(block)
  }

  logw <- -.enet_logQ(x0) - log(2)
  list(list(sign = eps, logw = rep_len(logw, n),
            wderiv = wderiv, x = eps * s * z - x0, xderiv = xderiv))
})
