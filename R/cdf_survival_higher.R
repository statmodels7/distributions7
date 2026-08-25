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
#' Returns \eqn{\partial^I F} for every component of the requested order, for a
#' family whose survival function is the exponential of something elementary.
#' Given \eqn{L = \log(1-F)} and a function evaluating its partial derivatives,
#' all four orders follow at once, so a family has only to say what \eqn{L} is.
#'
#' @details
#' # The identity
#'
#' \eqn{S = e^{L}} gives \eqn{\partial^I S = S\,B_I}, with \eqn{B_I} the
#' complete Bell polynomial in the partials of \eqn{L}, and \eqn{F = 1 - S}
#' turns that into \eqn{\partial^I F = -S\,B_I}. It is the same identity the
#' distribution wrappers use, read on the survival function instead of on the
#' density, and [bell_f_ratio()] runs the partition sum, so nothing here is
#' transcribed from an expansion.
#'
#' # Two things the arithmetic needs
#'
#' The survival function is evaluated as `exp(Lval)` and not as `1 - F`, which
#' keeps the far tail from canceling. And below the support \eqn{F} is
#' identically zero and so is every derivative, while \eqn{L} is still finite
#' there and would otherwise give a survival above one; the `inside` mask is
#' what suppresses that.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{S = 1 - F} the survival function,
#' \eqn{L = \log S} and \eqn{B_I} the complete Bell polynomial.
#'
#' @param distrib An object inheriting from `distrib`. Its `params` name and
#'   order the components, and its lower bound supplies the default mask.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The derivative order, 1 to 4.
#' @param Lval The value of \eqn{L} at `q`, a numeric vector.
#' @param Lderiv A function of a character vector of parameter names returning
#'   the corresponding partial derivative of \eqn{L}. An empty block is the
#'   zeroth order and is never asked for.
#' @param inside A logical vector saying which quantiles lie inside the
#'   support, or `NULL` (the default), which reads `q > distrib@bounds[1]`. A
#'   family whose support depends on a parameter, as the generalized Pareto's
#'   does at a negative shape, supplies its own; the fixed bounds cannot see
#'   it.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names], and exactly zero
#'   wherever `inside` is `FALSE`.
#'
#' @seealso [register_surv_cdf()], which turns a pieces function into the four
#'   methods; [bell_f_ratio()] for the partition sum;
#'   [gpd_surv_pieces()] for the most involved of the three families.
#'
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
#' Turns a function returning \eqn{L = \log(1-F)} and its partial-derivative
#' evaluator into the four S7 methods, so that a family states its survival
#' function once instead of four times. Three families are registered through
#' it: the exponential, the Weibull and the generalized Pareto.
#'
#' @details
#' All four orders are registered, [distrib_grad_cdf()] included, so these
#' families take the closed route from the first order up. Where the upper tail
#' is asked for on the natural scale the derivatives of \eqn{S} are returned
#' directly, \eqn{\partial^I S = -\partial^I F}, which is the one case where
#' the survival function is the quantity the construction produces first.
#'
#' `force(o)` inside the factory is what keeps the four registrations from
#' sharing one order.
#'
#' @param cls The S7 class to register on.
#' @param pieces A function of `(distrib, q, theta)` returning a list with
#'   `Lval` and `Lderiv`, and optionally `inside`, as
#'   [surv_cdf_deriv_k()] documents.
#'
#' @return Invisibly `NULL`. Called for the registration.
#'
#' @seealso [surv_cdf_deriv_k()], the body it registers;
#'   [distrib_grad_cdf.ExponentialDistrib()],
#'   [distrib_grad_cdf.Weibull1Distrib()] and
#'   [distrib_grad_cdf.GPDDistrib()], the three families.
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
#'
#' @description
#' Closed form at every order from one to four, from the survival function
#' \eqn{S = e^{-q/\mu}}. Its logarithm is \eqn{L = -q/\mu}, whose partial
#' derivatives are \eqn{\partial^j L/\partial\mu^j = -q(-1)^j j!/\mu^{j+1}},
#' and [surv_cdf_deriv_k()] turns those into the derivatives of \eqn{F}.
#'
#' @details
#' Against a product stencil on the same cdf: \eqn{2.7\times10^{-11}} at order
#' 1 and \eqn{2.7\times10^{-5}} at order 4. Below the support every derivative
#' is exactly zero, the mask in [surv_cdf_deriv_k()] suppressing the finite
#' value \eqn{L} would otherwise give there.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{F} the distribution function and
#' \eqn{S = 1 - F} the survival function.
#'
#' @param distrib An `ExponentialDistrib` object, from [exponential_distrib()].
#' @param q A numeric vector of quantiles. Values at or below zero give
#'   derivatives of exactly zero.
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per component of the order the
#'   generic asked for, which for a one-parameter family is one at every order.
#'
#' @seealso [surv_cdf_deriv_k()] for the identity;
#'   [distrib_grad_cdf.Weibull1Distrib()] and [distrib_grad_cdf.GPDDistrib()],
#'   the two families that contain this one; [exponential_distrib()].
#'
#' @examples
#' d <- exponential_distrib()
#' q <- c(0.5, 2, 5)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' fd <- numerical_cdf_deriv(d, q, list(mu = 3), order = 1)
#' max(abs(distrib_grad_cdf(d, q, list(mu = 3), log = FALSE)$mu / fd$mu - 1))
#'
#' # Exactly zero below the support.
#' distrib_grad_cdf(d, c(-1, 0.5), list(mu = 3), log = FALSE)$mu
#'
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
#'
#' @description
#' Closed form at every order from one to four, from the survival function
#' \eqn{S = \exp\{-(q/\mu)^{\sigma}\}}. Writing
#' \eqn{h = \sigma(\log q - \log\mu)} the exponent is \eqn{L = -e^{h}}, so its
#' partial derivatives are \eqn{-e^{h}} times the complete Bell polynomial in
#' the partials of \eqn{h}, and those are elementary.
#'
#' @details
#' # The partials of h
#'
#' \eqn{\partial^j h/\partial\mu^j = \sigma(-1)^j(j-1)!/\mu^j}, the same
#' without the factor \eqn{\sigma} when one index names the shape, and exactly
#' zero when two do: \eqn{h} is linear in the shape. That is what keeps the
#' expansion short at the higher orders.
#'
#' # An exact zero worth knowing about
#'
#' At \eqn{q = \mu} the exponent \eqn{h} vanishes and so does
#' \eqn{\partial h/\partial\sigma = \log q - \log\mu}, so the shape component
#' of the gradient is exactly zero there. A relative comparison against a
#' numerical derivative at that point measures nothing; an absolute one is what
#' to use, and it puts the closed route within \eqn{2.3\times10^{-11}} of a
#' central difference of the cdf.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the scale, \eqn{\sigma > 0} the shape,
#' \eqn{h = \sigma(\log q - \log\mu)}, \eqn{F} the distribution function and
#' \eqn{S = 1 - F} the survival function. The mean is
#' \eqn{\mu\,\Gamma(1+1/\sigma)}.
#'
#' @param distrib A `Weibull1Distrib` object, from [weibull1_distrib()].
#' @param q A numeric vector of quantiles. Values at or below zero give
#'   derivatives of exactly zero.
#' @param theta A named list with components `mu` (the scale, positive) and
#'   `sigma` (the shape, positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of the order the generic asked for,
#'   keyed as [`deriv_names(distrib@params, order)`][deriv_names]: two
#'   components for the gradient, three for the Hessian, four at order 3 and
#'   five at order 4.
#'
#' @seealso [surv_cdf_deriv_k()] for the identity;
#'   [distrib_grad_cdf.ExponentialDistrib()], the shape-1 case;
#'   [weibull1_distrib()].
#'
#' @examples
#' d <- weibull1_distrib()
#' q <- c(0.5, 2, 5)
#' th <- list(mu = 2, sigma = 3)
#'
#' # Against a central difference of the cdf, on an absolute scale.
#' fd <- numerical_cdf_deriv(d, q, th, order = 1)
#' max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) - unlist(fd)))
#'
#' # The shape component is exactly zero at q = mu.
#' distrib_grad_cdf(d, 2, th, log = FALSE)$sigma
#'
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
#' vector per order. The function is analytic at the origin, with
#' \eqn{\Lambda(0) = 1}, and it is the device that removes every division by
#' the generalized Pareto's shape from that family's survival function.
#'
#' @details
#' # Two routes, and where they change over
#'
#' Differentiating \eqn{u\Lambda = \log(1+u)} gives the recursion
#' \deqn{u\,\Lambda^{(r)} + r\,\Lambda^{(r-1)}
#'       = \frac{(-1)^{r-1}(r-1)!}{(1+u)^{r}},}
#' which is exact away from the origin and useless at it: it divides by
#' \eqn{u} and subtracts two nearly equal quantities. Measured against the
#' Taylor series, its fourth derivative is wrong by a factor of \eqn{10^{39}}
#' at \eqn{u = 10^{-14}}, by 1.7 at \eqn{10^{-4}} and by
#' \eqn{3\times10^{-8}} at \eqn{10^{-2}}.
#'
#' Below \eqn{|u| = 1/2} the series
#' \deqn{\Lambda^{(r)}(u) = \sum_{m \ge r} (-1)^{m}\frac{m!}{(m-r)!}
#'       \frac{u^{m-r}}{m+1}}
#' is used instead, where the two agree to \eqn{10^{-16}}. Its truncation is
#' set so that the \eqn{m^4} weight at the switch point stays under the
#' rounding.
#'
#' # Why the expression is arranged this way
#'
#' Differentiating \eqn{L = -\log(1+\xi q/\sigma)/\xi} directly gives terms in
#' \eqn{\xi^{-1-m}} that cancel only in the limit, which needs a guard and is
#' fragile whatever the guard. Writing \eqn{L = -(q/\sigma)\Lambda(u)} puts the
#' whole removable singularity inside one univariate function, and the
#' exponential limit \eqn{\xi \to 0} becomes an ordinary point of the formula.
#'
#' @section Notation:
#' \eqn{u = \xi q/\sigma} with \eqn{\xi} the shape and \eqn{\sigma} the scale.
#'
#' @param u A numeric vector, greater than \eqn{-1}. Values at or below
#'   \eqn{-1} are outside the support and are masked out by the caller before
#'   they reach here.
#'
#' @return A list of five numeric vectors the length of `u`, orders 0 to 4. At
#'   \eqn{u = 0} they are 1, \eqn{-1/2}, \eqn{2/3}, \eqn{-3/2} and
#'   \eqn{24/5}.
#'
#' @seealso [gpd_surv_pieces()], the one consumer;
#'   [distrib_grad_cdf.GPDDistrib()] for the family.
#'
#' @examples
#' # The limits at the origin, reached through the series branch.
#' vapply(distributions7:::gpd_lambda_derivs(1e-14), function(v) v[1],
#'        numeric(1))
#'
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
#' \eqn{(\sigma, \xi)}, in the form [register_surv_cdf()] wants. Writing
#' \eqn{L = -z\,\Lambda(u)} with \eqn{z = q/\sigma} and \eqn{u = \xi z} is what
#' keeps every division by the shape out of the expression.
#'
#' @details
#' # How the partials split
#'
#' The scale enters \eqn{z} as a plain reciprocal, and \eqn{u} is bilinear in
#' the shape and \eqn{z}, so a block naming the shape twice contributes nothing
#' to \eqn{u}. The partials of \eqn{\Lambda(u)} follow by Faa di Bruno over
#' that, and the product with \eqn{z} by Leibniz over the scale indices alone,
#' the shape not entering \eqn{z}.
#'
#' # The support
#'
#' A negative shape bounds the support above, at \eqn{u = -1}, so the mask is
#' `q > 0 & u > -1` and is supplied here, the family's fixed bounds being
#' unable to see it. Past the upper endpoint every derivative of \eqn{F} is
#' exactly zero.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale, \eqn{\xi} the shape of either sign,
#' \eqn{z = q/\sigma}, \eqn{u = \xi z} and \eqn{\Lambda(u) = \log(1+u)/u}.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `sigma` (positive) and `xi` (any
#'   real value), each a numeric vector of length 1 or `n`.
#'
#' @return A list with `Lval` (a numeric vector), `Lderiv` (a function of a
#'   block of parameter names) and `inside` (a logical vector).
#'
#' @seealso [gpd_lambda_derivs()] for the univariate function;
#'   [surv_cdf_deriv_k()] and [register_surv_cdf()];
#'   [distrib_grad_cdf.GPDDistrib()] for the family page.
#'
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
#'
#' @description
#' Closed form at every order from one to four, from the survival function
#' \eqn{S = (1 + \xi q/\sigma)^{-1/\xi}}. Its logarithm is written
#' \eqn{L = -(q/\sigma)\,\Lambda(\xi q/\sigma)} with
#' \eqn{\Lambda(u) = \log(1+u)/u}, which carries no division by the shape, so
#' the exponential limit \eqn{\xi \to 0} is an ordinary point of the formula
#' and not a branch.
#'
#' @details
#' # The support moves with the shape
#'
#' At \eqn{\xi \ge 0} the support is \eqn{(0, \infty)}; at \eqn{\xi < 0} it is
#' bounded above at \eqn{\sigma/|\xi|}, and past that endpoint every derivative
#' is exactly zero. The mask is computed in [gpd_surv_pieces()], the family's
#' fixed bounds being unable to record a support that moves with a parameter.
#'
#' # What it is worth, and the limit as a check
#'
#' Against a product stencil on the same cdf at \eqn{\sigma = 1},
#' \eqn{\xi = 0.3}: \eqn{5.2\times10^{-11}} at order 1,
#' \eqn{2.7\times10^{-7}} at order 2, \eqn{1.7\times10^{-5}} at order 3 and
#' \eqn{8.4\times10^{-4}} at order 4. At \eqn{\xi = 0} the scale component
#' equals the exponential family's to the last bit, and the fourth derivative
#' reads the same value at \eqn{\xi = 10^{-7}} and at \eqn{10^{-9}}, which is
#' what a removable singularity handled properly looks like.
#'
#' @section Notation:
#' \eqn{\sigma > 0} is the scale, \eqn{\xi} the shape of either sign,
#' \eqn{u = \xi q/\sigma}, \eqn{\Lambda(u) = \log(1+u)/u}, \eqn{F} the
#' distribution function and \eqn{S = 1 - F} the survival function.
#'
#' @param distrib A `GPDDistrib` object, from [gpd_distrib()].
#' @param q A numeric vector of quantiles. Values outside the support give
#'   derivatives of exactly zero.
#' @param theta A named list with components `sigma` (positive) and `xi` (any
#'   real value), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of the order the generic asked for,
#'   keyed as [`deriv_names(distrib@params, order)`][deriv_names]: two
#'   components for the gradient, three for the Hessian, four at order 3 and
#'   five at order 4.
#'
#' @seealso [gpd_surv_pieces()] and [gpd_lambda_derivs()] for the construction;
#'   [distrib_grad_cdf.ExponentialDistrib()], the \eqn{\xi = 0} case;
#'   [gpd_distrib()].
#'
#' @examples
#' d <- gpd_distrib()
#' q <- c(0.5, 2, 5)
#'
#' # At shape zero the scale component is the exponential family's.
#' rbind(gpd = distrib_grad_cdf(d, q, list(sigma = 3, xi = 0),
#'                              log = FALSE)$sigma,
#'       exponential = distrib_grad_cdf(exponential_distrib(), q,
#'                                      list(mu = 3), log = FALSE)$mu)
#'
#' # A negative shape bounds the support at sigma / |xi| = 2.
#' distrib_grad_cdf(d, c(1, 2, 3), list(sigma = 1, xi = -0.5),
#'                  log = FALSE)$sigma
#'
#' @keywords internal
register_surv_cdf(GPDDistrib, gpd_surv_pieces)
