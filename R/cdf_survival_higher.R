#' @include cdf_mapped_higher.R
NULL

# The families whose survival function is an exponential of something
# elementary.
#
# With S = 1 - F written as exp(L), every derivative of F in the parameters is
# d^I F = -S B_I(L), the complete Bell polynomial in L's own partial
# derivatives -- the same identity the wrappers use, read on the survival
# function instead of on the density. So a family needs only to say what L is
# and what its partials are, and all four orders follow at once. Nothing here
# is transcribed from an expansion: `bell_f_ratio()` runs the partition sum.

#' CDF Derivatives From an Exponential Survival Function
#'
#' @description
#' Returns \eqn{\partial^{I}F} for every component of the requested order,
#' given \eqn{L = \log(1-F)} and a function evaluating its partial derivatives.
#'
#' @details
#' \eqn{S = e^{L}} gives \eqn{\partial^{I}S = S\,B_{I}}, with \eqn{B_{I}} the
#' complete Bell polynomial in the partials of \eqn{L}, and \eqn{F = 1 - S}
#' turns that into \eqn{\partial^{I}F = -S\,B_{I}}. The survival function is
#' evaluated as `exp(L)` rather than as `1 - F`, which keeps the far
#' tail from cancelling.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#' @param Lval The value of \eqn{L} at `q`.
#' @param Lderiv A function of a character vector of parameter names returning
#'   the corresponding partial derivative of \eqn{L}.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso [bell_f_ratio()]
#' @keywords internal
surv_cdf_deriv_k <- function(distrib, q, theta, order, Lval, Lderiv,
                             inside = NULL) {
  params <- distrib@params
  # below the support F is identically zero and so is every derivative; L is
  # still finite there and would otherwise produce a survival above one. A
  # family whose support depends on a parameter -- the generalized Pareto at a
  # negative shape -- says so itself rather than being read off the bounds.
  if (is.null(inside)) inside <- q > distrib@bounds[1L]
  S <- exp(Lval)
  idx <- deriv_indices(params, order)
  out <- lapply(idx, function(I) -S * bell_f_ratio(params[I], Lderiv) * inside)
  stats::setNames(out, deriv_names(params, order))
}

#' Register the Four CDF Derivative Orders of an Exponential Survival Family
#'
#' @description
#' Turns a function returning \eqn{L} and its partial-derivative evaluator into
#' the four methods.
#'
#' @param cls The S7 class.
#' @param pieces A function of `(distrib, q, theta)` returning a list with
#'   `Lval` and `Lderiv`.
#'
#' @return Invisibly `NULL`; called for the registration.
#'
#' @keywords internal
register_surv_cdf <- function(cls, pieces) {
  make <- function(o) {
    force(o)
    function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
      p <- pieces(distrib, q, theta)
      params <- distrib@params
      inside <- if (is.null(p$inside)) q > distrib@bounds[1L] else p$inside
      if (!lower.tail && log) {
        # log S is L, so the components ARE its partial derivatives. Taking
        # them the general way instead would divide by S = 1 - F read off the
        # natural scale, which is exactly zero past q/mu = 37 for an
        # exponential while L is finite to the end of the range.
        out <- lapply(deriv_indices(params, o),
                      function(I) rep_len(p$Lderiv(params[I]), length(q)) * inside)
        return(stats::setNames(out, deriv_names(params, o)))
      }
      tabs <- lapply(seq_len(o), function(k)
        surv_cdf_deriv_k(distrib, q, theta, k, p$Lval, p$Lderiv, inside))
      if (!lower.tail) {
        # d^I S = -d^I F, with S itself taken as exp(L)
        return(lapply(tabs[[o]], function(v) -v))
      }
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

#' @title Exponential Log-CDF Derivatives
#' @name distrib_grad_cdf.ExponentialDistrib
#' @description
#' Closed form at every order from the survival function
#' \eqn{S = \exp(-q/\mu)}, whose logarithm has the partial derivatives
#' \eqn{\partial^{j}L/\partial\mu^{j} = -q(-1)^{j}j!/\mu^{j+1}}.
#' @param distrib An `ExponentialDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing `mu`.
#' @param lower.tail Logical; if `TRUE` (default), the lower tail.
#' @param log Logical; if `TRUE` (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso [exponential_distrib()]
#' @keywords internal
register_surv_cdf(ExponentialDistrib, function(distrib, q, theta) {
  mu <- theta[[1]]
  list(
    Lval = -q / mu,
    Lderiv = function(block) {
      j <- length(block)
      -q * (-1)^j * factorial(j) / mu^(j + 1)
    }
  )
})

#' @title Weibull Log-CDF Derivatives
#' @name distrib_grad_cdf.Weibull1Distrib
#' @description
#' Closed form at every order from the survival function
#' \eqn{S = \exp\{-(q/\mu)^{\sigma}\}}. Writing \eqn{h = \sigma(\log q -
#' \log\mu)} the exponent is \eqn{L = -e^{h}}, so its partial derivatives are
#' \eqn{-e^{h}} times the complete Bell polynomial in the partials of \eqn{h},
#' and those are elementary: \eqn{\partial^{j}h/\partial\mu^{j} =
#' \sigma(-1)^{j}(j-1)!/\mu^{j}}, the same without the \eqn{\sigma} when one
#' index names the shape, and zero when two do.
#' @param distrib A `Weibull1Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing `mu` and `sigma`.
#' @param lower.tail Logical; if `TRUE` (default), the lower tail.
#' @param log Logical; if `TRUE` (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso [weibull1_distrib()]
#' @keywords internal
register_surv_cdf(Weibull1Distrib, function(distrib, q, theta) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  nm <- distrib@params
  h <- sigma * (log(q) - log(mu))
  # h is linear in the shape and a logarithm in the scale, so a block naming
  # the shape twice is zero and one naming it once drops the factor sigma
  hderiv <- function(block) {
    ns <- sum(block == nm[2L])
    j <- sum(block == nm[1L])
    if (ns >= 2L) return(0 * q)
    base <- if (j == 0L) rep_len(1, length(q)) else
      (-1)^j * factorial(j - 1L) / mu^j
    if (ns == 1L) {
      if (j == 0L) log(q) - log(mu) else base
    } else {
      sigma * base
    }
  }
  eh <- exp(h)
  list(
    Lval = -eh,
    Lderiv = function(block) -eh * bell_f_ratio(block, hderiv)
  )
})


# --- the generalized Pareto ------------------------------------------------
#
# Its survival function is exp(L) too, with L = -log1p(xi q/sigma)/xi, so the
# route above applies. What it needs is a form of L free of the 1/xi: writing
# u = xi q / sigma and
#
#   Lambda(u) = log1p(u)/u,        L = -(q/sigma) Lambda(u),
#
# every division by the shape disappears, the whole removable singularity
# sitting inside Lambda, which is analytic with Lambda(0) = 1. The exponential
# limit is then an ordinary point rather than a special case.

#' Derivatives of log1p(u)/u
#'
#' @description
#' Returns \eqn{\Lambda(u) = \log(1+u)/u} and its first four derivatives, one
#' vector per order.
#'
#' @details
#' Differentiating \eqn{u\Lambda = \log(1+u)} gives
#' \eqn{u\Lambda^{(r)} + r\Lambda^{(r-1)} = (-1)^{r-1}(r-1)!/(1+u)^{r}}, a
#' recursion that is exact away from the origin and useless at it: it divides
#' by \eqn{u} and subtracts two nearly equal quantities, and measured against
#' the Taylor series the fourth derivative is wrong by a factor of \eqn{10^{39}}
#' at \eqn{u = 10^{-14}}, by 1.7 at \eqn{10^{-4}} and by \eqn{3\times10^{-8}} at
#' \eqn{10^{-2}}. The series
#' \eqn{\Lambda^{(r)}(u) = \sum_{m\ge r}(-1)^{m}\frac{m!}{(m-r)!}
#' \frac{u^{m-r}}{m+1}} is used instead below \eqn{\lvert u\rvert = 1/2}, where
#' the two agree to \eqn{10^{-16}}, and its truncation is set so that the
#' \eqn{m^{4}} weight at the switch point stays under the rounding.
#'
#' @param u A numeric vector, greater than \eqn{-1}.
#'
#' @return A list of five numeric vectors, orders 0 to 4.
#'
#' @seealso [gpd_surv_pieces()]
#' @keywords internal
gpd_lambda_derivs <- function(u) {
  R <- 4L
  out <- lapply(seq_len(R + 1L), function(i) rep(NA_real_, length(u)))
  small <- abs(u) <= 0.5
  if (any(small)) {
    us <- u[small]
    for (r in 0:R) {
      acc <- rep(0, length(us))
      for (m in r:(r + 120L)) {
        acc <- acc + (-1)^m * exp(lfactorial(m) - lfactorial(m - r)) *
          us^(m - r) / (m + 1)
      }
      out[[r + 1L]][small] <- acc
    }
  }
  if (any(!small)) {
    ub <- u[!small]
    prev <- log1p(ub) / ub
    out[[1L]][!small] <- prev
    for (r in seq_len(R)) {
      prev <- ((-1)^(r - 1L) * factorial(r - 1L) / (1 + ub)^r - r * prev) / ub
      out[[r + 1L]][!small] <- prev
    }
  }
  out
}

#' The Exponential Survival Pieces of a Generalized Pareto
#'
#' @description
#' Returns \eqn{L = \log(1-F)} and an evaluator of its partial derivatives in
#' \eqn{(\sigma, \xi)}.
#'
#' @details
#' \eqn{L = -z\,\Lambda(u)} with \eqn{z = q/\sigma} and \eqn{u = \xi z}. The
#' scale enters \eqn{z} as a plain reciprocal and \eqn{u} is bilinear in the
#' shape and \eqn{z}, so a block naming the shape twice contributes nothing to
#' \eqn{u}; the partials of \eqn{\Lambda(u)} follow by Faa di Bruno over that,
#' and the product with \eqn{z} by Leibniz over the scale indices alone, the
#' shape not entering \eqn{z}.
#'
#' @param distrib A `GPDDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing `sigma` and `xi`.
#'
#' @return A list with `Lval`, `Lderiv` and `inside`.
#'
#' @seealso [gpd_lambda_derivs()], [register_surv_cdf()]
#' @keywords internal
gpd_surv_pieces <- function(distrib, q, theta) {
  sigma <- theta[[1]]
  xi <- theta[[2]]
  nm <- distrib@params
  n <- length(q)
  z <- q / sigma
  u <- xi * z
  # a negative shape bounds the support above, at u = -1; past it the survival
  # is zero and every derivative of F vanishes
  inside <- q > 0 & u > -1
  u_safe <- ifelse(inside, u, 0)
  z_safe <- ifelse(inside, z, 0)
  lam <- gpd_lambda_derivs(u_safe)

  # d^j z / d sigma^j
  dz <- lapply(0:4, function(j)
    rep_len(q * (-1)^j * factorial(j) / sigma^(1 + j), n))
  # d^{j sigma, k xi} u; u is bilinear, so two shape indices give zero
  du <- function(j, k) {
    if (k >= 2L) return(rep(0, n))
    if (k == 1L) return(dz[[j + 1L]])
    xi * dz[[j + 1L]]
  }
  # d^S Lambda(u) by Faa di Bruno over u
  dG <- function(j, k) {
    S <- c(rep("s", j), rep("x", k))
    if (!length(S)) return(lam[[1L]])
    acc <- 0
    for (part in index_partitions(S)) {
      term <- lam[[length(part) + 1L]]
      for (b in part) term <- term * du(sum(b == "s"), sum(b == "x"))
      acc <- acc + term
    }
    acc
  }

  list(
    Lval = -z_safe * lam[[1L]],
    inside = inside,
    Lderiv = function(block) {
      j <- sum(block == nm[1L])
      k <- sum(block == nm[2L])
      acc <- 0
      for (i in 0:j) acc <- acc + choose(j, i) * dz[[i + 1L]] * dG(j - i, k)
      -acc
    }
  )
}

#' @title Generalized Pareto Log-CDF Derivatives
#' @name distrib_grad_cdf.GPDDistrib
#' @description
#' Closed form at every order, from the survival function
#' \eqn{S = (1 + \xi q/\sigma)^{-1/\xi}}. Its logarithm is written
#' \eqn{-(q/\sigma)\Lambda(\xi q/\sigma)} with \eqn{\Lambda(u) = \log(1+u)/u},
#' which carries no division by the shape, so the exponential limit
#' \eqn{\xi \to 0} is an ordinary point of the formula rather than a branch.
#' @param distrib A `GPDDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing `sigma` and `xi`.
#' @param lower.tail Logical; if `TRUE` (default), the lower tail.
#' @param log Logical; if `TRUE` (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso [gpd_distrib()]
#' @keywords internal
register_surv_cdf(GPDDistrib, gpd_surv_pieces)
