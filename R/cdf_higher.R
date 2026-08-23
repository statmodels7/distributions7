#' @include cdf_derivatives.R partition_sums.R y_higher.R
NULL

# ===========================================================================
# Third and fourth derivatives of the distribution function.
#
# The governing identity is the one of cdf_derivatives.R at any order,
#
#   d^I F(q) / F(q) = E[ d^I f / f  |  Y <= q ],
#
# so the two routes are the same two: for a discrete family the conditional
# expectation is a finite sum and the identity is exact, and for a continuous
# one the cdf is analytic and differencing it is cheaper and more accurate
# than a semi-infinite quadrature.
#
# What changes at orders three and four is only that the quantity summed,
# d^I f / f, is a longer partition sum -- which bell_f_ratio() already
# computes -- and that the conversion to the log scale is the general
# moment-to-cumulant relation rather than the written-out second-order
# formula, which log_deriv() already computes. Both helpers are the wrappers'
# and are reused rather than copied.
# ===========================================================================

#' CDF Derivatives on the Requested Tail and Scale, at Any Order
#'
#' @description
#' The general form of [cdf_tail_scale()]: converts derivatives of
#' \eqn{F} of every order up to the one wanted into derivatives of whichever
#' tail was asked for, on the natural or the log scale.
#'
#' @details
#' Switching tail flips the sign, since \eqn{S = 1 - F}. Switching to the log
#' scale is the moment-to-cumulant relation
#' \eqn{d^I \log P = \sum_\pi (-1)^{|\pi|-1}(|\pi|-1)!
#' \prod_{B} (d^B P / P)}, which at second order is the familiar
#' \eqn{d^2 P/P - (dP/P)(dP/P)} and at third and fourth is what
#' [log_deriv()] sums. Only the ratios are needed, which is why
#' every order up to `order` has to be supplied.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param Fq The cdf evaluated at the quantile.
#' @param dF A list of length `order`; element `k` is the table of
#'   \eqn{k}-th derivatives of \eqn{F}, keyed as
#'   [`deriv_names(distrib@params, k)`][deriv_names].
#' @param order The derivative order wanted, 1 to 4.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of derivative component vectors of the requested
#'   order.
#'
#' @seealso [cdf_tail_scale()], [log_deriv()]
#' @keywords internal
cdf_scale_k <- function(distrib, Fq, dF, order, lower.tail, log) {
  params <- distrib@params
  P <- if (lower.tail) Fq else 1 - Fq
  sgn <- if (lower.tail) 1 else -1

  dP <- lapply(dF, function(tab) lapply(tab, function(v) sgn * v))
  if (!log) return(dP[[order]])

  ratio <- function(block) {
    dP[[length(block)]][[canon_key(block, params)]] / P
  }
  nm <- deriv_names(params, order)
  idx <- deriv_indices(params, order)
  stats::setNames(lapply(idx, function(I) log_deriv(params[I], ratio)), nm)
}

#' CDF Derivatives of a Discrete Distribution, at Any Order
#'
#' @description
#' The general form of [discrete_cdf_deriv()]: evaluates
#' \eqn{d^I F(q) = \sum_{y \le q} f(y)\,(d^I f/f)(y)} for any order up to
#' four.
#'
#' @details
#' The quantity summed is the complete Bell polynomial in the log-derivatives,
#' which [bell_f_ratio()] computes, so the order enters only through
#' how many of the family's derivative tables are fetched. At orders one and
#' two this reproduces the written-out \eqn{f g} and
#' \eqn{f(h + g_i g_j)} of [discrete_cdf_deriv()].
#'
#' As there, a test must not check this against the partial-expectation sum,
#' which is the same sum twice; finite differences of the cdf are the
#' independent reference.
#'
#' @param distrib An object inheriting from class `"discrete_distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso [discrete_cdf_deriv()], [bell_f_ratio()]
#' @keywords internal
discrete_cdf_deriv_k <- function(distrib, q, theta, order) {
  params <- distrib@params
  lo <- distrib@bounds[1]
  n <- length(q)
  nm <- deriv_names(params, order)
  idx <- deriv_indices(params, order)

  out <- stats::setNames(lapply(nm, function(i) numeric(n)), nm)
  for (k in seq_len(n)) {
    if (q[k] < lo) next
    th_k <- lapply(theta[seq_along(params)],
                   function(v) if (length(v) > 1) v[k] else v)
    grid <- seq(lo, floor(q[k] + 1e-9))
    fy <- distrib_pdf(distrib, grid, th_k)
    ell <- parent_ell(distrib, grid, th_k, order, params)
    for (m in seq_along(nm)) {
      out[[m]][k] <- sum(fy * bell_f_ratio(params[idx[[m]]], ell))
    }
  }
  out
}

#' Numerical CDF Derivatives of Any Order
#'
#' @description
#' One product stencil of the requested order applied to
#' [distrib_cdf()], which is analytic for every family in the
#' catalog.
#'
#' @details
#' A repeated parameter contributes the matching higher one-dimensional
#' factor and distinct parameters each contribute a central two-point factor,
#' so the whole thing is one stencil rather than a difference of a
#' difference. The step is \eqn{\varepsilon^{1/(k+2)}} scaled by the
#' parameter, which balances the \eqn{h^{2}} truncation against the
#' \eqn{\varepsilon/h^{k}} rounding, and is chosen per observation because
#' `theta` may vary by observation.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order The derivative order, 3 or 4.
#' @param h_rel Relative finite-difference step.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso [numerical_cdf_deriv()]
#' @keywords internal
numerical_cdf_deriv_k <- function(distrib, q, theta, order,
                                  h_rel = .Machine$double.eps^(1 / (order + 2))) {
  params <- distrib@params
  hs <- lapply(params, function(p) h_rel * pmax(1, abs(theta[[p]])))
  names(hs) <- params

  fac <- list(
    list(o = c(-1, 1), w = c(-0.5, 0.5)),
    list(o = c(-1, 0, 1), w = c(1, -2, 1)),
    list(o = c(-2, -1, 1, 2), w = c(-0.5, 1, -1, 0.5)),
    list(o = c(-2, -1, 0, 1, 2), w = c(1, -4, 6, -4, 1))
  )
  nm <- deriv_names(params, order)
  idx <- deriv_indices(params, order)

  stats::setNames(lapply(idx, function(I) {
    who <- params[I]
    tb <- table(who)
    ks <- names(tb)
    fs <- fac[as.integer(tb)]
    grid <- expand.grid(lapply(fs, function(f) seq_along(f$o)))
    acc <- 0
    for (r in seq_len(nrow(grid))) {
      t2 <- theta
      w <- 1
      for (j in seq_along(ks)) {
        pick <- grid[r, j]
        t2[[ks[j]]] <- theta[[ks[j]]] + fs[[j]]$o[pick] * hs[[ks[j]]]
        w <- w * fs[[j]]$w[pick]
      }
      acc <- acc + w * distrib_cdf(distrib, q, t2)
    }
    den <- 1
    for (j in seq_along(ks)) den <- den * hs[[ks[j]]]^as.integer(tb)[j]
    acc / den
  }), nm)
}

#' CDF Derivative Tables of Every Order Up To One
#'
#' @description
#' The derivatives of \eqn{F} of orders 1 to `order`, by whichever route
#' the class uses.
#'
#' @details
#' The conversion to the log scale needs every order below the one wanted, not
#' just the one wanted, because the moment-to-cumulant relation is a sum over
#' partitions of the multi-index and a partition into \eqn{k} blocks asks for
#' \eqn{k} lower-order ratios. Collecting them in one place keeps the choice of
#' route -- exact sum for a discrete family, one stencil on the analytic
#' distribution function for a continuous one -- in a single statement.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order The highest order wanted, 1 to 4.
#'
#' @return A list of length `order` of named derivative tables of
#'   \eqn{F}.
#'
#' @seealso [cdf_scale_k()], [discrete_cdf_deriv_k()],
#'   [numerical_cdf_deriv_k()]
#' @keywords internal
cdf_tables <- function(distrib, q, theta, order) {
  discrete <- S7::S7_inherits(distrib, discrete_distrib)
  lapply(seq_len(order), function(k) {
    if (discrete) {
      discrete_cdf_deriv_k(distrib, q, theta, k)
    } else if (k <= 2L) {
      numerical_cdf_deriv(distrib, q, theta, order = k)
    } else {
      numerical_cdf_deriv_k(distrib, q, theta, k)
    }
  })
}

#' Third and Fourth Derivatives of the Log Distribution Function
#'
#' @description
#' \eqn{\partial^{3}\log F(q)/\partial\theta_i\partial\theta_j\partial\theta_k}
#' and its fourth-order analogue, on either tail.
#'
#' @details
#' These complete the sequence begun by [distrib_grad_cdf()] and
#' [distrib_hess_cdf()]. What consumes them is truncation: with only
#' the first two orders available, [truncated()] pays one quadrature
#' per component at orders three and four, and with these it pays two calls on
#' the parent instead.
#'
#' A discrete family uses the exact finite sum and a continuous one one
#' product stencil on its analytic distribution function; see
#' [discrete_cdf_deriv_k()] and [numerical_cdf_deriv_k()].
#' A family with a closed form registers its own method, as at the orders
#' below.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#' @param ... Passed to methods.
#'
#' @return A named list of derivative component vectors, keyed as
#'   [`deriv_names(distrib@params, 3)`][deriv_names] or `4`.
#'
#' @seealso [distrib_hess_cdf()], [truncated()]
#'
#' @examples
#' distrib_deriv3_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1))
#'
#' @export
distrib_deriv3_cdf <- S7::new_generic(
  "distrib_deriv3_cdf", "distrib",
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @rdname distrib_deriv3_cdf
#'
#' @examples
#' distrib_deriv4_cdf(poisson_distrib(), 3, list(mu = 2))
#'
#' @export
distrib_deriv4_cdf <- S7::new_generic(
  "distrib_deriv4_cdf", "distrib",
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @title Default Third Log-CDF Derivatives
#' @name distrib_deriv3_cdf.distrib
#' @description The route of [cdf_tables()] for the class, put on
#'   the requested tail and scale.
#' @param distrib An object inheriting from class `"distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical.
#' @param log Logical.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @keywords internal
S7::method(distrib_deriv3_cdf, distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log = TRUE, ...) {
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta),
              cdf_tables(distrib, q, theta, 3L), 3L, lower.tail, log)
}

#' @rdname distrib_deriv3_cdf.distrib
#' @name distrib_deriv4_cdf.distrib
#' @return A named list of fourth-derivative components.
S7::method(distrib_deriv4_cdf, distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log = TRUE, ...) {
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta),
              cdf_tables(distrib, q, theta, 4L), 4L, lower.tail, log)
}

#' Location-Scale CDF Derivatives at Any Order
#'
#' @description
#' Closed-form derivatives of \eqn{F} of any order up to four, for a family
#' that is location-scale in its first two parameters.
#'
#' @details
#' With \eqn{z = (q-\mu)/\sigma} the distribution function is
#' \eqn{F(q) = F_0(z)}, so every derivative in \eqn{(\mu, \sigma)} is one
#' Faa di Bruno pass over that composition. The inner derivatives are
#' \eqn{F_0^{(m)}(z) = \sigma^{m} \partial^{m} F/\partial q^{m}}, and
#' \eqn{\partial^{m} F/\partial q^{m} = f(q) B_{m-1}}, the complete Bell
#' polynomial in the response derivatives of \eqn{\log f}. The map is
#' \eqn{\partial^{i+j} z/\partial\mu^{i}\partial\sigma^{j}}, which vanishes
#' for \eqn{i \ge 2} because \eqn{z} is linear in the location.
#'
#' This is what the response derivatives of order three and four are for:
#' with only [distrib_grad_y()] and [distrib_hess_y()]
#' the construction stops at second order, which is where
#' [loc_scale_cdf_deriv()] stops. At orders one and two the two
#' agree exactly.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, location first and scale second.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso [loc_scale_cdf_deriv()], [chain_assemble()]
#' @keywords internal
loc_scale_cdf_deriv_k <- function(distrib, q, theta, order) {
  s <- theta[[2]]
  z <- (q - theta[[1]]) / s
  f <- distrib_pdf(distrib, q, theta)

  # the response derivatives of log f, as many as the order needs
  ly <- list(distrib_grad_y(distrib, q, theta),
             distrib_hess_y(distrib, q, theta),
             if (order >= 3L) distrib_deriv3_y(distrib, q, theta) else NULL,
             if (order >= 4L) distrib_deriv4_y(distrib, q, theta) else NULL)
  ell <- function(block) ly[[length(block)]]

  # F0^(m)(z) = sigma^m f B_{m-1}
  D <- lapply(seq_len(order), function(m) {
    b <- if (m == 1L) 1 else bell_f_ratio(rep("y", m - 1L), ell)
    stats::setNames(list(s^m * f * b), paste(rep("z", m), collapse = "_"))
  })

  zmap <- list(z = list(
    "1" = -1 / s, "2" = -z / s,
    "1,2" = 1 / s^2, "2,2" = 2 * z / s^2,
    "1,2,2" = -2 / s^3, "2,2,2" = -6 * z / s^3,
    "1,2,2,2" = 6 / s^4, "2,2,2,2" = 24 * z / s^4
  ))
  chain_assemble(D, "z", zmap, distrib@params[1:2], order, length(q))
}

#' Location-Scale Third and Fourth Log-CDF Derivatives
#'
#' @description
#' The [distrib_deriv3_cdf()] and [distrib_deriv4_cdf()]
#' bodies shared by the location-scale families.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @seealso [loc_scale_cdf_deriv_k()]
#' @keywords internal
loc_scale_deriv_cdf_k <- function(order) {
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    dF <- lapply(seq_len(order),
                 function(k) loc_scale_cdf_deriv_k(distrib, q, theta, k))
    cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), dF, order,
                lower.tail, log)
  }
}

# the location-scale families, the same four the orders below cover
for (.cls in list(Gaussian1Distrib, LogisticDistrib, CauchyDistrib,
                  LaplaceDistrib)) {
  S7::method(distrib_deriv3_cdf, .cls) <- loc_scale_deriv_cdf_k(3L)
  S7::method(distrib_deriv4_cdf, .cls) <- loc_scale_deriv_cdf_k(4L)
}
rm(.cls)
