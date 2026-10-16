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
#' derivatives of whichever tail was requested, on the natural or the
#' logarithmic scale, at any order from 1 to 4.
#'
#' @details
#' # The two conversions
#'
#' Switching to the upper tail flips the sign at every order, \eqn{S = 1 - F};
#' the probability \eqn{S} by which the log scale divides is the family's own
#' upper tail when `Fq` comes from [tail_cdf()].
#' Switching to the log scale is the moment-to-cumulant relation
#' \deqn{\partial^I \log P = \sum_\pi (-1)^{|\pi|-1}(|\pi|-1)!
#'       \prod_{B \in \pi} \frac{\partial^B P}{P},}
#' summed over the set partitions of the multi-index by [log_deriv()]. At
#' second order it is the familiar
#' \eqn{\partial^2 P/P - (\partial P/P)^2}, and at third and fourth order it
#' has 5 and 15 terms.
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
#' \eqn{P} whichever was requested, \eqn{\pi} a set partition of the
#' multi-index and \eqn{B} one of its blocks.
#'
#' @param distrib An object inheriting from `distrib`. Only its `params` are
#'   read, to name and enumerate the components.
#' @param Fq The value of [tail_cdf()], the probability of the requested
#'   tail, or a plain numeric vector of distribution function values.
#' @param dF A list of length `order`. Element \eqn{k} is the table of
#'   \eqn{k}-th derivatives of \eqn{F}, keyed as
#'   [`deriv_names(distrib@params, k)`][deriv_names], on the natural scale and
#'   the lower tail.
#' @param order The derivative order wanted, 1 to 4.
#' @param lower.tail A single logical; `TRUE` selects the lower tail and `FALSE`
#'   the upper tail.
#' @param log A single logical; `TRUE` selects the derivatives of the log
#'   probability and `FALSE` those of the probability. `FALSE` reads only
#'   `dF[[order]]` and the lower tables are then unused.
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
  sgn <- if (lower.tail) 1 else -1
  # order two in the order of hess_names(), as every Hessian of the package
  hnames <- function(out) if (order == 2L) out[hess_names(params)] else out

  dP <- lapply(dF, function(tab) lapply(tab, function(v) sgn * v))
  if (!log) return(hnames(dP[[order]][deriv_names(params, order)]))

  P <- tail_probability(Fq, lower.tail)
  ratio <- function(block) {
    dP[[length(block)]][[canon_key(block, params)]] / P
  }
  nm <- deriv_names(params, order)
  idx <- deriv_indices(params, order)
  hnames(stats::setNames(lapply(idx, function(I) log_deriv(params[I], ratio)), nm))
}

#' @title The Probability of the Requested Tail
#' @name tail_probability
#'
#' @description
#' Internal helpers of [cdf_tail_scale()] and [cdf_scale_k()].
#' `tail_cdf()` returns the probability of the requested tail at `q`, from
#' the family's own distribution function for that tail (so that the upper
#' tail is not formed as \eqn{1 - F}, which loses its digits once \eqn{F} is
#' close to one), with the attribute `"tail"` set to `"lower"` or `"upper"`.
#' Called with the arguments that [log_tail_route()] has stored in
#' `log_tail_stash`, it returns the stored value.
#' `tail_probability()` returns the requested tail: the value itself when it
#' carries the attribute `"tail"` for that tail, and otherwise, the value being
#' read as the distribution function, \eqn{F} or \eqn{1 - F}.
#'
#' @param distrib A distribution object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameter values.
#' @param lower.tail A single logical; `TRUE` selects the lower tail and `FALSE`
#'   the upper tail.
#' @param Fq The value of `tail_cdf()`, or a plain numeric vector of
#'   distribution function values.
#'
#' @return A numeric vector of probabilities; for `tail_cdf()`, with the
#'   attribute `"tail"`.
#'
#' @examples
#' d <- logistic_distrib()
#' Fq <- distributions7:::tail_cdf(d, 40, list(mu = 0, sigma = 1), FALSE)
#' # The upper tail survives where 1 - F is zero.
#' c(from_family = distributions7:::tail_probability(Fq, FALSE),
#'   one_minus_F = 1 - plogis(40))
#'
#' @keywords internal
tail_cdf <- function(distrib, q, theta, lower.tail) {
  # the value log_tail_route() computed for these arguments
  st <- log_tail_stash
  if (!is.null(st$Fq) && identical(st$lower, lower.tail) &&
      identical(attr(S7::S7_class(distrib), "name"), st$cls) &&
      identical(q, st$q) && identical(theta[names(st$theta)], st$theta))
    return(st$Fq)
  structure(as.numeric(distrib_cdf(distrib, q, theta, lower.tail = lower.tail)),
            tail = if (lower.tail) "lower" else "upper")
}

#' @rdname tail_probability
tail_probability <- function(Fq, lower.tail) {
  v <- as.numeric(Fq)
  tl <- attr(Fq, "tail")
  if (is.null(tl)) tl <- "lower"
  if (identical(tl, if (lower.tail) "lower" else "upper")) v else 1 - v
}

#' CDF Derivatives of a Discrete Distribution at Any Order
#'
#' @description
#' The general form of [discrete_cdf_deriv()]: evaluates
#' \eqn{\partial^I F(q) = \sum_{y \le q} f(y)\,\partial^I f/f\,(y)} for a
#' discrete family at any order up to four. Nothing is differenced; the sum is
#' exact wherever the support has a finite lower bound, which the discrete class
#' requires. Where \eqn{P(Y > q)} is below 0.01 the sum is computed as minus the
#' sum above \eqn{q}, the derivatives of the mass summing to zero over the
#' support.
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
  # a constant that varies by observation (a binomial's size) is taken at
  # the row, as the parameters are
  rows <- distrib_const_rows(distrib) > 1L
  # the upper tail of every quantile, read once to choose the side of each sum
  surv <- distrib_cdf(distrib, q, theta[seq_along(params)], lower.tail = FALSE)
  surv <- rep_len(surv, n)
  for (k in seq_len(n)) {
    if (is.na(q[k])) {
      for (m in seq_along(nm)) out[[m]][k] <- NA_real_
      next
    }
    if (q[k] < lo) next
    th_k <- lapply(theta[seq_along(params)],
                   function(v) if (length(v) > 1) v[k] else v)
    dk <- if (rows) distrib_at_rows(distrib, k, n) else distrib
    term <- function(grid, m) {
      ell <- parent_ell(dk, grid, th_k, order, params)
      distrib_pdf(dk, grid, th_k) * bell_f_ratio(params[idx[[m]]], ell)
    }
    top <- floor(q[k] + 1e-9)
    # the lower sum loses about eps / S(q) of relative accuracy, so it is
    # replaced only where the upper tail is small
    if (isTRUE(surv[k] < 0.01)) {
      # the derivatives of the mass sum to zero over the support, so in the
      # upper tail the sum below q is minus the small sum above it
      for (m in seq_along(nm)) {
        out[[m]][k] <- -discrete_support_sum(function(K, I) term(K, m),
                                             top + 1, distrib@bounds[2], 1L)
      }
    } else {
      grid <- seq(lo, top)
      fy <- distrib_pdf(dk, grid, th_k)
      ell <- parent_ell(dk, grid, th_k, order, params)
      for (m in seq_along(nm)) {
        out[[m]][k] <- sum(fy * bell_f_ratio(params[idx[[m]]], ell))
      }
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
#' # Accuracy
#'
#' Against the closed forms of a location-scale family the stencil has a
#' relative error, normwise over each component, of the order of
#' \eqn{10^{-5}} at order 3 and \eqn{10^{-4}} at order 4. The closed forms are
#' accurate to rounding error, which is why a family registers one where it
#' can.
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

#' CDF Derivative Tables up to a Given Order
#'
#' @description
#' Assembles the derivatives of \eqn{F} of orders 1 to `order`. A discrete
#' family uses the exact finite sum at every order. A continuous family uses
#' its own [distrib_grad_cdf()] and [distrib_hess_cdf()] at orders 1 and 2
#' where it registers them, and [numerical_cdf_deriv()] otherwise; at orders
#' 3 and 4 it uses [continuous_cdf_deriv_k()], which is the compiled
#' quadrature for some families and one product stencil on the distribution
#' function for the others.
#'
#' @details
#' Every order below the one wanted is collected, and not just the one wanted,
#' because the moment-to-cumulant relation [cdf_scale_k()] applies is a sum
#' over partitions of the multi-index: a partition into \eqn{k} blocks reads
#' \eqn{k} lower-order ratios.
#'
#' Taking orders 1 and 2 from the family's own methods matters on the log
#' scale, where [cdf_scale_k()] combines them with the higher orders: a
#' differenced table at order 1 or 2 would limit the accuracy of an exact one
#' at order 3 or 4.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The highest order wanted, 1 to 4.
#' @param closed Logical of length 1. `TRUE`, the default, takes orders 1 and
#'   2 from the family's own [distrib_grad_cdf()] and [distrib_hess_cdf()]
#'   where [has_exact_cdf_deriv()] finds them; `FALSE` takes them from
#'   [numerical_cdf_deriv()]. A caller that is itself registered as one of
#'   those two methods passes `FALSE`, which keeps it from calling itself.
#'
#' @return A list of length `order`. Element \eqn{k} is a named list of
#'   \eqn{k}-th derivatives of \eqn{F}, on the natural scale and the lower
#'   tail, keyed as [`deriv_names(distrib@params, k)`][deriv_names].
#'
#' @seealso [cdf_scale_k()], the consumer;
#'   [discrete_cdf_deriv_k()] and [numerical_cdf_deriv_k()], the two routes.
#'
#' @keywords internal
cdf_tables <- function(distrib, q, theta, order, closed = TRUE) {
  discrete <- S7::S7_inherits(distrib, discrete_distrib)
  lapply(seq_len(order), function(k) {
    if (discrete) {
      discrete_cdf_deriv_k(distrib, q, theta, k)
    } else if (k <= 2L) {
      # a family's own first and second orders are exact where it has them;
      # a difference here would limit the third and fourth orders built on it
      if (closed && has_exact_cdf_deriv(distrib, k)) {
        own_cdf_deriv(distrib, q, theta, k)
      } else {
        numerical_cdf_deriv(distrib, q, theta, order = k)
      }
    } else {
      continuous_cdf_deriv_k(distrib, q, theta, k)
    }
  })
}

#' A Family's Own First or Second CDF Derivatives as a Table
#'
#' @description
#' Calls the family's [distrib_grad_cdf()] (`k = 1`) or [distrib_hess_cdf()]
#' (`k = 2`) on the lower tail and the natural scale, and orders the result as
#' [`deriv_names(distrib@params, k)`][deriv_names], the keying [cdf_tables()]
#' and [cdf_scale_k()] use.
#'
#' @param distrib An object inheriting from `continuous_distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param k The order, 1 or 2.
#'
#' @return A named list of numeric vectors, the derivatives of \eqn{F} of
#'   order `k`.
#'
#' @seealso [cdf_tables()] and [partial_loc_scale_deriv_cdf_k()], the two
#'   callers.
#'
#' @examples
#' d <- gamma2_distrib()
#' distributions7:::own_cdf_deriv(d, c(1, 2, 4), list(mu = 2, sigma2 = 1), 2L)
#'
#' @keywords internal
own_cdf_deriv <- function(distrib, q, theta, k) {
  gen <- if (k == 1L) distrib_grad_cdf else distrib_hess_cdf
  out <- gen(distrib, q, theta, lower.tail = TRUE, log = FALSE)
  out[deriv_names(distrib@params, k)]
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
#' A discrete family uses the exact finite sum of [discrete_cdf_deriv_k()]. A
#' continuous family whose compiled distribution function integrates its shape
#' derivatives (gamma1, gamma2, chisq, both generalized gammas, both betas,
#' both von Mises, both Student t, both pseudo-Huber and the skew t) takes the
#' exact integrals of [compiled_cdf_deriv_k()]; any other continuous family
#' without a closed form of its own takes a single product stencil on its
#' distribution function through [numerical_cdf_deriv_k()].
#'
#' # What consumes them
#'
#' Truncation. [truncated()] needs the derivatives of its normalizing constant
#' \eqn{Z = F(U) - F(L^-)}, and with only the first two orders available it
#' pays one quadrature per component at orders three and four; with these it
#' pays two calls on the parent instead.
#'
#' # Accuracy
#'
#' The closed route runs a Faa di Bruno pass over the response derivatives and
#' costs about as much as the stencil. It is used for its accuracy: the
#' stencil has a relative error of order \eqn{10^{-5}} at order 3 and
#' \eqn{10^{-4}} at order 4.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{S = 1 - F} the survival function
#' and \eqn{\theta} the parameter on its own scale.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, in which case one value is returned per setting.
#' @param lower.tail A single logical, `TRUE` by default; `TRUE` selects the
#'   lower tail and `FALSE` the upper tail. With `log = FALSE`, `FALSE` flips
#'   the sign of every component; with `log = TRUE` it returns the derivatives
#'   of \eqn{\log S}, which are not the negatives of the lower-tail ones.
#' @param log A single logical, `TRUE` by default; `TRUE` selects the
#'   derivatives of the log probability and `FALSE` those of the probability.
#'   Far into a tail the probability underflows to zero and the result is `-Inf`
#'   or `NaN`.
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
    args <- check_derivative_args(distrib, q, theta)
    q <- args$y
    theta <- args$theta
    # the log-scale derivatives of the families log_tail_route() lists
    r <- log_tail_route(distrib, q, theta, 3L, lower.tail, log, distrib_deriv3_cdf)
    if (!is.null(r)) return(r)
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
    args <- check_derivative_args(distrib, q, theta)
    q <- args$y
    theta <- args$theta
    # the log-scale derivatives of the families log_tail_route() lists
    r <- log_tail_route(distrib, q, theta, 4L, lower.tail, log, distrib_deriv4_cdf)
    if (!is.null(r)) return(r)
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
#' Among the shipped univariate families the discrete ones reach these
#' methods, and their result is exact. Every shipped continuous family
#' registers its own methods: closed forms, or the compiled quadratures of
#' [compiled_cdf_deriv_k()]. The continuous default, one product stencil on
#' the distribution function, serves families defined by the user and the
#' wrappers, and carries a relative error of order \eqn{10^{-5}} at order 3
#' and \eqn{10^{-4}} at order 4.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param lower.tail A single logical, `TRUE` by default; `TRUE` selects the
#'   lower tail and `FALSE` the upper tail.
#' @param log A single logical, `TRUE` by default; `TRUE` selects the
#'   derivatives of the log probability and `FALSE` those of the probability.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return For `distrib_deriv3_cdf()`, a named list of third-derivative
#'   components keyed as [`deriv_names(distrib@params, 3)`][deriv_names]; for
#'   `distrib_deriv4_cdf()`, the fourth-order components keyed as
#'   `deriv_names(distrib@params, 4)`. Each vector is the length of `q`
#'   recycled against `theta`.
#'
#' @seealso [cdf_tables()] and [cdf_scale_k()], the two halves;
#'   [loc_scale_deriv_cdf_k()], the closed route that five families take
#'   instead.
#'
#' @examples
#' # A truncated gamma reaches this method and differences its cdf.
#' distrib_deriv3_cdf(truncated(gamma2_distrib(), upper = 6), 2,
#'                    list(mu = 2, sigma2 = 1))
#'
#' # A beta-binomial reaches it too, and its sum is exact.
#' distrib_deriv4_cdf(betabinom1_distrib(size = 10), 4,
#'                    list(mu = 0.3, sigma = 0.5))
#'
#' @keywords internal
S7::method(distrib_deriv3_cdf, distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log = TRUE, ...) {
  cdf_scale_k(distrib, tail_cdf(distrib, q, theta, lower.tail),
              cdf_tables(distrib, q, theta, 3L), 3L, lower.tail, log)
}

#' @rdname distrib_deriv3_cdf.distrib
#' @name distrib_deriv4_cdf.distrib
#' @keywords internal
S7::method(distrib_deriv4_cdf, distrib) <- function(distrib, q, theta,
                                                    lower.tail = TRUE,
                                                    log = TRUE, ...) {
  cdf_scale_k(distrib, tail_cdf(distrib, q, theta, lower.tail),
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
#' the two agree to rounding error, which is the check that supports the
#' orders above.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density, \eqn{F_0} the standardized
#' distribution function and \eqn{B_m} the complete Bell polynomial.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale, and which supplies response
#'   derivatives to the requested order.
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
loc_scale_cdf_deriv_k <- function(distrib, q, theta, order, over_f = FALSE) {
  s <- theta[[2]]
  z <- (q - theta[[1]]) / s
  f <- if (over_f) 1 else distrib_pdf(distrib, q, theta)

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
#' [loc_scale_cdf_route()]. Four families use it, the Gaussian, the logistic
#' and the Laplace here and the Gumbel in `cdf_mapped_higher.R`; the Cauchy
#' has its own closed form, [cauchy_cdf_method()].
#'
#' @details
#' Every lower order is computed because the log-scale conversion reads them
#' all. The function is a factory so that the order is fixed in the closure
#' of each registered method.
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
#'   distrib_deriv3_cdf.LogisticDistrib
#'   distrib_deriv3_cdf.LaplaceDistrib distrib_deriv3_cdf.GumbelDistrib
#'   distrib_deriv4_cdf.Gaussian1Distrib
#'   distrib_deriv4_cdf.LogisticDistrib
#'   distrib_deriv4_cdf.LaplaceDistrib distrib_deriv4_cdf.GumbelDistrib
#'
#' @keywords internal
loc_scale_deriv_cdf_k <- function(order) {
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    loc_scale_cdf_route(distrib, q, theta, order, lower.tail, log)
  }
}

#' @title Location-Scale Distribution Function Derivatives on Either Tail
#' @name loc_scale_cdf_route
#' @description
#' The body shared by the location-scale families' [distrib_grad_cdf()],
#' [distrib_hess_cdf()], [distrib_deriv3_cdf()] and [distrib_deriv4_cdf()].
#' Every table of [loc_scale_cdf_deriv_k()] is the density times a factor, so
#' on the log scale the factors are multiplied by
#' \eqn{f/P = \exp(\log f - \log P)}, from the family's log-density and the
#' log probability of the requested tail: neither \eqn{f} nor \eqn{P} is
#' formed where it underflows, and the upper tail is the family's own.
#' @param distrib A location-scale distribution object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameter values.
#' @param order The derivative order, 1 to 4.
#' @param lower.tail,log Logical of length 1.
#' @return A named list of numeric vectors, keyed as [deriv_names()] at the
#'   order, or as [hess_names()] at order 2.
#' @examples
#' d <- logistic_distrib()
#' distributions7:::loc_scale_cdf_route(d, 40, list(mu = 0, sigma = 1), 1L,
#'                                      FALSE, TRUE)
#' @keywords internal
loc_scale_cdf_route <- function(distrib, q, theta, order, lower.tail, log) {
  tabs <- lapply(seq_len(order), function(k) {
    loc_scale_cdf_deriv_k(distrib, q, theta, k, over_f = TRUE)
  })
  lf <- distrib_pdf(distrib, q, theta, log = TRUE)
  if (!log) {
    f <- exp(lf)
    return(cdf_scale_k(distrib, NULL,
                       lapply(tabs, function(t) lapply(t, function(v) v * f)),
                       order, lower.tail, FALSE))
  }
  sgn <- if (lower.tail) 1 else -1
  rho <- sgn * exp(lf - distrib_cdf(distrib, q, theta, lower.tail = lower.tail,
                                    log.p = TRUE))
  cdf_scale_k(distrib, rep(1, length(rho)),
              lapply(tabs, function(t) lapply(t, function(v) v * rho)),
              order, TRUE, TRUE)
}

# the location-scale families, the same four the orders below cover
# the Cauchy has its own closed form, cauchy_cdf_method()
S7::method(distrib_deriv3_cdf, CauchyDistrib) <- cauchy_cdf_method(3L)
S7::method(distrib_deriv4_cdf, CauchyDistrib) <- cauchy_cdf_method(4L)

for (.cls in list(Gaussian1Distrib, LogisticDistrib, LaplaceDistrib)) {
  S7::method(distrib_deriv3_cdf, .cls) <- loc_scale_deriv_cdf_k(3L)
  S7::method(distrib_deriv4_cdf, .cls) <- loc_scale_deriv_cdf_k(4L)
}
rm(.cls)
