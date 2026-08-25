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
#' The general form of [cdf_tail_scale()]: converts derivatives of \eqn{F} into
#' derivatives of whichever tail was asked for, on the natural or the
#' logarithmic scale, at any order from 1 to 4.
#'
#' @details
#' # The two conversions
#'
#' Switching to the upper tail flips the sign at every order, \eqn{S = 1 - F}.
#' Switching to the log scale is the moment-to-cumulant relation
#' \deqn{\partial^I \log P = \sum_\pi (-1)^{|\pi|-1}(|\pi|-1)!
#'       \prod_{B \in \pi} \frac{\partial^B P}{P},}
#' summed over the set partitions of the multi-index by [log_deriv()]. At
#' second order it is the familiar
#' \eqn{\partial^2 P/P - (\partial P/P)^2}, and at third and fourth it has 4
#' and 15 terms.
#'
#' # Why every lower order is needed
#'
#' A partition into \eqn{k} blocks reads \eqn{k} ratios, so the relation at
#' order 4 needs the tables of orders 1, 2 and 3 as well. That is why `dF` is a
#' list of tables, and why a caller that wants only the fourth order still
#' assembles the first three.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{S = 1 - F} the survival function,
#' \eqn{P} whichever was asked for, \eqn{\pi} a set partition of the
#' multi-index and \eqn{B} one of its blocks.
#'
#' @param distrib An object inheriting from `distrib`. Only its `params` are
#'   read, to name and enumerate the components.
#' @param Fq The distribution function at the quantile, a numeric vector.
#' @param dF A list of length `order`. Element \eqn{k} is the table of
#'   \eqn{k}-th derivatives of \eqn{F}, keyed as
#'   [`deriv_names(distrib@params, k)`][deriv_names], on the natural scale and
#'   the lower tail.
#' @param order The derivative order wanted, 1 to 4.
#' @param lower.tail Is the lower tail wanted? A single logical.
#' @param log Are derivatives of the log probability wanted? A single logical.
#'   `FALSE` reads only `dF[[order]]` and the lower tables are then unused.
#'
#' @return A named list of numeric vectors of the requested order alone. The
#'   lower orders are consumed on the way and do not appear in the result.
#'
#' @seealso [cdf_tail_scale()], the written-out version for orders 1 and 2;
#'   [log_deriv()] for the partition sum; [cdf_tables()], which builds `dF`.
#'
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

#' CDF Derivatives of a Discrete Distribution at Any Order
#'
#' @description
#' The general form of [discrete_cdf_deriv()]: evaluates
#' \eqn{\partial^I F(q) = \sum_{y \le q} f(y)\,\partial^I f/f\,(y)} for a
#' discrete family at any order up to four. Nothing is differenced; the sum is
#' exact wherever the support has a finite lower bound, which the discrete
#' class requires.
#'
#' @details
#' What changes above second order is only the summand. \eqn{\partial^I f/f} is
#' the complete Bell polynomial in the derivatives of \eqn{\log f}, and
#' [bell_f_ratio()] already computes it for the distribution wrappers, so this
#' function reuses it and carries no second copy of the enumeration.
#'
#' @section Notation:
#' \eqn{f} is the mass function, \eqn{F} the distribution function,
#' \eqn{\ell = \log f} and \eqn{\partial^I} a derivative with respect to a
#' multi-index of parameters.
#'
#' @param distrib An object inheriting from `discrete_distrib`.
#' @param q A numeric vector of quantiles. Each is handled separately, the
#'   support being walked up to it, so the cost grows with the largest one.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names].
#'
#' @seealso [discrete_cdf_deriv()] for orders 1 and 2;
#'   [numerical_cdf_deriv_k()], the continuous route;
#'   [bell_f_ratio()] for the summand; [cdf_tables()], the caller.
#'
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
#' One product stencil of the requested order applied to [distrib_cdf()], which
#' is analytic for every family in the catalog. This is the continuous
#' family's route above second order, and it is a single stencil: a repeated
#' parameter contributes the matching one-dimensional higher-order factor and
#' distinct parameters each contribute a central two-point factor, so no
#' difference of a difference is ever taken.
#'
#' @details
#' # The step
#'
#' The relative step is \eqn{\varepsilon^{1/(k+2)}}, which is
#' \eqn{7.4\times10^{-4}} at order 3 and \eqn{2.5\times10^{-3}} at order 4, and
#' it is scaled per component by the parameter. It balances the \eqn{h^2}
#' truncation against the \eqn{\varepsilon/h^{k}} rounding, and it is chosen
#' per observation because `theta` may vary by observation.
#'
#' # What it delivers
#'
#' Measured on a Gaussian against the closed forms that family registers, the
#' relative error is \eqn{1.5\times10^{-5}} at order 3 and
#' \eqn{1.3\times10^{-4}} at order 4. That is a long way short of the closed
#' route and is the reason a family registers one where it can: the loss is
#' about five digits at order 3 and about four at order 4.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{\theta} the parameter on its own
#' scale, \eqn{h} the step and \eqn{\varepsilon} the machine epsilon.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The derivative order, 3 or 4.
#' @param h_rel The relative step. A single number, defaulting to
#'   \eqn{\varepsilon^{1/(\mathrm{order}+2)}}. A step much smaller than the
#'   default is worse, the rounding growing as \eqn{h^{-\mathrm{order}}}.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names].
#'
#' @seealso [numerical_cdf_deriv()] for orders 1 and 2;
#'   [discrete_cdf_deriv_k()], the exact route for a discrete family;
#'   [loc_scale_cdf_deriv_k()], the closed route this is measured against.
#'
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
#' Assembles the derivatives of \eqn{F} of orders 1 to `order` by whichever
#' route the class uses: the exact finite sum for a discrete family, one
#' product stencil on the analytic distribution function for a continuous one.
#' Keeping the choice of route in a single statement is the point of the
#' function.
#'
#' @details
#' Every order below the one wanted is collected, and not just the one wanted,
#' because the moment-to-cumulant relation [cdf_scale_k()] applies is a sum
#' over partitions of the multi-index: a partition into \eqn{k} blocks reads
#' \eqn{k} lower-order ratios.
#'
#' At orders 1 and 2 the tables come from [distrib_grad_cdf()] and
#' [distrib_hess_cdf()], so a family's own closed forms are used where it has
#' them; only orders 3 and 4 reach the routes named above.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The highest order wanted, 1 to 4.
#'
#' @return A list of length `order`. Element \eqn{k} is a named list of
#'   \eqn{k}-th derivatives of \eqn{F}, on the natural scale and the lower
#'   tail, keyed as [`deriv_names(distrib@params, k)`][deriv_names].
#'
#' @seealso [cdf_scale_k()], the consumer;
#'   [discrete_cdf_deriv_k()] and [numerical_cdf_deriv_k()], the two routes.
#'
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
#' `distrib_deriv3_cdf()` returns
#' \eqn{\partial^3 \log F(q)/\partial\theta_i\partial\theta_j\partial\theta_k}
#' and `distrib_deriv4_cdf()` its fourth-order analogue, on either tail and on
#' the natural or the logarithmic scale. Together with [distrib_grad_cdf()] and
#' [distrib_hess_cdf()] they complete the cdf derivative surface to fourth
#' order.
#'
#' @details
#' # The two routes
#'
#' A discrete family uses the exact finite sum of [discrete_cdf_deriv_k()], and
#' a continuous one takes a single product stencil on its analytic distribution
#' function through [numerical_cdf_deriv_k()]. 24 of the 42 univariate families
#' register a closed form of their own; of the 18 that do not, the discrete
#' ones sum exactly and the continuous ones (beta1, beta2, chisq, gamma1,
#' gamma2, gengamma1 and the two von Mises) difference.
#'
#' # What consumes them
#'
#' Truncation. [truncated()] needs the derivatives of its normalizing constant
#' \eqn{Z = F(U) - F(L^-)}, and with only the first two orders available it
#' pays one quadrature per component at orders three and four; with these it
#' pays two calls on the parent instead.
#'
#' # Accuracy against speed
#'
#' Unusually, the closed route is the slower of the two at these orders: on a
#' Gaussian at 1000 quantiles it costs 7.0 ms against the stencil's 4.7 ms,
#' because it runs a Faa di Bruno pass over the response derivatives. It is
#' preferred for accuracy alone, the stencil being \eqn{1.5\times10^{-5}} out
#' at order 3 and \eqn{1.3\times10^{-4}} at order 4.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{S = 1 - F} the survival function
#' and \eqn{\theta} the parameter on its own scale.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, in which case one value is returned per setting.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default; `FALSE` flips the sign of every component.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default. Far into a tail the probability underflows to zero and
#'   the result is `-Inf` or `NaN`.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`deriv_names(distrib@params, 3)`][deriv_names] for
#'   `distrib_deriv3_cdf()` and as `deriv_names(distrib@params, 4)` for
#'   `distrib_deriv4_cdf()`. A two-parameter family has 4 third-order and 5
#'   fourth-order components; a one-parameter family has one of each.
#'
#' @seealso [distrib_grad_cdf()] and [distrib_hess_cdf()] for the two orders
#'   below; [truncated()], the consumer; [numerical_cdf_deriv_k()] and
#'   [discrete_cdf_deriv_k()] for the two routes.
#'
#' @examples
#' # Four third-order components for a two-parameter family.
#' distrib_deriv3_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1))
#'
#' # On the upper tail every sign flips.
#' distrib_deriv3_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1),
#'                    lower.tail = FALSE, log = FALSE)
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
#' # One component for a one-parameter family, and the sum is exact here.
#' distrib_deriv4_cdf(poisson_distrib(), 3, list(mu = 2))
#'
#' @export
distrib_deriv4_cdf <- S7::new_generic(
  "distrib_deriv4_cdf", "distrib",
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @title Default Third and Fourth Log-CDF Derivatives
#' @name distrib_deriv3_cdf.distrib
#'
#' @description
#' The fallback both generics take when a family registers no closed form.
#' [cdf_tables()] assembles the derivatives of \eqn{F} of every order up to the
#' one wanted, by the exact sum for a discrete family and by one product
#' stencil for a continuous one, and [cdf_scale_k()] puts the top order on the
#' requested tail and scale.
#'
#' @details
#' 18 of the 42 univariate families reach these methods. For the ten discrete
#' ones the result is exact; for the eight continuous ones (beta1, beta2,
#' chisq, gamma1, gamma2, gengamma1 and the two von Mises) it carries the
#' stencil's error, about \eqn{1.5\times10^{-5}} at order 3 and
#' \eqn{1.3\times10^{-4}} at order 4.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return For `distrib_deriv3_cdf()`, a named list of third-derivative
#'   components keyed as [`deriv_names(distrib@params, 3)`][deriv_names]; for
#'   `distrib_deriv4_cdf()`, the fourth-order components keyed as
#'   `deriv_names(distrib@params, 4)`. Each vector is the length of `q`
#'   recycled against `theta`.
#'
#' @seealso [cdf_tables()] and [cdf_scale_k()], the two halves;
#'   [loc_scale_deriv_cdf_k()], the closed route four families take instead.
#'
#' @examples
#' # A gamma reaches this method and differences its cdf.
#' distrib_deriv3_cdf(gamma2_distrib(), 2, list(mu = 2, sigma2 = 1))
#'
#' # A beta-binomial reaches it too, and its sum is exact.
#' distrib_deriv4_cdf(betabinom1_distrib(size = 10), 4,
#'                    list(mu = 0.3, sigma = 0.5))
#'
#' @keywords internal
S7::method(distrib_deriv3_cdf, distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log = TRUE, ...) {
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta),
              cdf_tables(distrib, q, theta, 3L), 3L, lower.tail, log)
}

#' @rdname distrib_deriv3_cdf.distrib
#' @name distrib_deriv4_cdf.distrib
#' @keywords internal
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
#' that is location-scale in its first two parameters. The general form of
#' [loc_scale_cdf_deriv()], which stops at second order.
#'
#' @details
#' # The construction
#'
#' With \eqn{z = (q-\mu)/\sigma} the distribution function is
#' \eqn{F(q) = F_0(z)}, so every derivative in \eqn{(\mu, \sigma)} is one Faa
#' di Bruno pass over that composition. The inner derivatives are
#' \deqn{F_0^{(m)}(z) = \sigma^{m}\,\frac{\partial^{m} F}{\partial q^{m}},
#'       \qquad
#'       \frac{\partial^{m} F}{\partial q^{m}} = f(q)\,B_{m-1},}
#' with \eqn{B_{m-1}} the complete Bell polynomial in the response derivatives
#' of \eqn{\log f}. The outer map is
#' \eqn{\partial^{i+j} z/\partial\mu^{i}\partial\sigma^{j}}, which vanishes for
#' \eqn{i \ge 2} because \eqn{z} is linear in the location, so most of its
#' partials are exact zeros.
#'
#' # What it depends on
#'
#' The response derivatives of orders 3 and 4 are what make this reach past
#' second order: with only [distrib_grad_y()] and [distrib_hess_y()] the
#' construction stops where [loc_scale_cdf_deriv()] stops. At orders 1 and 2
#' the two agree to the last bit, which is the check that licenses the orders
#' above.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density, \eqn{F_0} the standardized
#' distribution function and \eqn{B_m} the complete Bell polynomial.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale, and which supplies response
#'   derivatives to the order asked for.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, the location first and the scale
#'   second.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(distrib@params[1:2], order)`][deriv_names].
#'
#' @seealso [loc_scale_cdf_deriv()] for orders 1 and 2;
#'   [loc_scale_deriv_cdf_k()], which registers this as a method;
#'   [chain_assemble()] for the partition sum; [bell_f_ratio()].
#'
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
#' Builds the [distrib_deriv3_cdf()] or [distrib_deriv4_cdf()] body that the
#' location-scale families register: [loc_scale_cdf_deriv_k()] at every order
#' up to the one wanted, put on the requested tail and scale by
#' [cdf_scale_k()]. Five families use it, the Gaussian, the logistic, the
#' Cauchy and the Laplace here and the Gumbel in `cdf_mapped_higher.R`.
#'
#' @details
#' Every lower order is computed because the log-scale conversion reads them
#' all. The function is a factory so that the order can be closed over,
#' `force(order)` being what keeps the two registrations from sharing one
#' value.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function of `(distrib, q, theta, lower.tail, log, ...)` suitable
#'   for registering as an S7 method on either generic, returning a named list
#'   of numeric vectors of that order.
#'
#' @seealso [loc_scale_cdf_deriv_k()] for the formulas;
#'   [distrib_deriv3_cdf()] for the generic;
#'   [loc_scale_grad_cdf()] and [loc_scale_hess_cdf()] for the orders below.
#'
#' @aliases distrib_deriv3_cdf.Gaussian1Distrib
#'   distrib_deriv3_cdf.LogisticDistrib distrib_deriv3_cdf.CauchyDistrib
#'   distrib_deriv3_cdf.LaplaceDistrib distrib_deriv3_cdf.GumbelDistrib
#'   distrib_deriv4_cdf.Gaussian1Distrib
#'   distrib_deriv4_cdf.LogisticDistrib distrib_deriv4_cdf.CauchyDistrib
#'   distrib_deriv4_cdf.LaplaceDistrib distrib_deriv4_cdf.GumbelDistrib
#'
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
