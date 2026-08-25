#' @include distrib.R generics.R utility_functions.R expected_derivatives.R link_scale.R partition_sums.R
NULL

# ===========================================================================
# Third and fourth derivatives for the wrappers.
#
# Every wrapper's log-likelihood has the same shape: the parent's log-density
# plus (or instead of) the logarithm of some quantity L that depends on the
# parameters.
#
#   zero-inflated, at y = 0     l_Y = log L0,      L0 = zi + (1-zi) f0
#   hurdle, at y > 0            l_Y = l(y) - log(1 - f0) + log(1 - za)
#   zero-adjusted continuous    l_Y = l_W(y) + log(1 - za)        (y != 0)
#   truncated                   l_T = l(y) - log Z
#
# So two general facts do all the work, and both are sums over set partitions
# of the multi-index -- the same enumeration the Bartlett identities use.
#
# (1) Derivatives of a density through its log. For a multi-index I,
#
#         d^I f / f  =  sum over partitions pi of I of  prod_{B in pi} l^(B),
#
#     the complete Bell polynomial in the parent's derivatives. This is the
#     lemma behind the Bartlett identities, applied in the other direction.
#
# (2) Derivatives of a logarithm. Inverting (1),
#
#         d^I log L  =  sum over pi  (-1)^{|pi|-1} (|pi|-1)!  prod_B (d^B L / L),
#
#     the moment-to-cumulant relation. Only the ratios d^B L / L are needed,
#     never L's derivatives on their own.
#
# Order 1 and 2 of every wrapper fall out of these as special cases, and agree
# with the closed forms written by hand in zero_inflated.R, zero_adjusted.R and
# truncated.R -- which is what makes those a check on this file and vice versa.
# ===========================================================================


#' Memoize a Ratio Function on Its Block
#'
#' @description
#' Caches a ratio function by the canonical key of its block.
#'
#' @details
#' A partition of a fourth-order index asks for the same small blocks many times
#' over. For the truncated wrapper each distinct block costs a quadrature, so
#' memoizing across the partition sum is the difference between one integration
#' per block and one per occurrence.
#'
#' @param f The ratio function to wrap.
#' @param params The parameter names, in declaration order.
#'
#' @return A function with the same signature as `f`, backed by a cache.
#'
#' @seealso [trunc_deriv_k()]
#' @keywords internal
memo_ratio <- function(f, params) {
  cache <- list()
  function(block) {
    k <- canon_key(block, params)
    if (is.null(cache[[k]])) cache[[k]] <<- f(block)
    cache[[k]]
  }
}

#' Multi-Indices of a Given Order, as Parameter Names
#'
#' @description
#' The multi-indices of a given order, in the order [deriv_names()]
#' lists them, expressed as parameter names.
#'
#' @details
#' A thin wrapper on [deriv_indices()]. It is deliberately not
#' `deriv_index_list()` from `link_scale.R`, whose order-2 case is
#' ordered for [hess_names()] -- diagonal first -- while
#' `deriv_names()` is lexicographic; pairing those would silently attach the
#' name `"mu_sigma"` to the index `(sigma, sigma)`. The orders actually
#' registered here are 3 and 4, where the two orderings agree; generating the
#' indices locally removes the mismatch instead of depending on it never being
#' reached.
#'
#' @param params A character vector of parameter names.
#' @param order The derivative order.
#'
#' @return A list of character vectors, each of length `order`.
#'
#' @seealso [deriv_indices()], [deriv_names()]
#' @keywords internal
order_indices <- function(params, order) {
  lapply(deriv_indices(params, order), function(i) params[i])
}

#' Derivatives of log(p) and log(1 - p)
#'
#' @description
#' The two elementary logarithmic derivatives the zero wrappers need:
#' \eqn{d^k \log p = (-1)^{k-1}(k-1)!/p^k} and
#' \eqn{d^k \log(1-p) = -(k-1)!/(1-p)^k}.
#'
#' @param p A numeric vector of probabilities.
#' @param k The derivative order.
#' @param complement Logical; `TRUE` for \eqn{\log(1-p)}.
#'
#' @return A numeric vector.
#'
#' @keywords internal
log_pow_deriv <- function(p, k, complement = FALSE) {
  if (complement) -gamma(k) / (1 - p)^k        # d^k log(1-p) = -(k-1)!/(1-p)^k
  else (-1)^(k - 1) * gamma(k) / p^k           # d^k log(p)   = (-1)^{k-1}(k-1)!/p^k
}

#' Split a Multi-Index Into Parent and Wrapper Parts
#'
#' @description
#' Separates a multi-index into the parent's parameters and the number of times
#' the wrapper's own parameter appears.
#'
#' @details
#' The wrapper's parameter is always the last one it declares, so the split needs
#' only its name.
#'
#' @param idx A character vector of parameter names, with repetition.
#' @param mix_name The wrapper parameter's name.
#'
#' @return A list with `theta`, the parent part, and `n_mix`, a count.
#'
#' @keywords internal
split_index <- function(idx, mix_name) {
  list(theta = idx[idx != mix_name], n_mix = sum(idx == mix_name))
}

#' Assemble One Order of a Wrapper's Derivatives
#'
#' @description
#' The shared skeleton: builds the named list of components of a given order by
#' calling `component` once per multi-index.
#'
#' @param distrib The wrapper distribution.
#' @param order The derivative order.
#' @param component A function of one multi-index returning that component.
#'
#' @return A named list of derivative component vectors.
#'
#' @keywords internal
assemble_deriv <- function(distrib, order, component) {
  params <- distrib@params
  nms <- deriv_names(params, order)
  idxs <- order_indices(params, order)
  stats::setNames(lapply(seq_along(nms), function(t) component(idxs[[t]])), nms)
}


# --------------------------------------------------------------------------
# Transformed distributions
#
# Nothing to do: the Jacobian does not depend on theta, so every derivative of
# the transformed log-likelihood is the parent's evaluated at x = g^-1(y). The
# methods exist only so that the numerical fallback is never reached.
# --------------------------------------------------------------------------

#' Derivatives of a Transformed Distribution
#'
#' @description
#' Builds the order-`k` derivative method for
#' [transformation()].
#'
#' @details
#' There is nothing to do. The Jacobian does not depend on \eqn{\theta}, so every
#' derivative of the transformed log-likelihood is the parent's evaluated at
#' \eqn{x = g^{-1}(y)}. The methods exist only so that the numerical fallback is
#' never reached.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @keywords internal
trans_deriv_k <- function(order) {
  function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
    x <- distrib@transformer@trans_inv(y)
    if (order == 3L) distrib_deriv3(distrib@parent_distrib, x, theta, expected = expected,
                                    approx = match.arg(approx), nsim = nsim)
    else distrib_deriv4(distrib@parent_distrib, x, theta, expected = expected,
                        approx = match.arg(approx), nsim = nsim)
  }
}


# --------------------------------------------------------------------------
# Zero-inflation
#
# At y > 0 the likelihood separates, log(1-zi) + l(y), so a mixed index gives
# zero and the pure ones are immediate. At y = 0 it is log L0 with
# L0 = zi + (1-zi) f0, which is *affine* in zi -- so a block containing two or
# more zi's contributes nothing, and the ratios d^B L0 / L0 are:
#
#   B all theta          w0 * (d^A f0 / f0)          w0 = (1-zi) f0 / L0
#   B with one zi        -f0 (d^A f0/f0) / L0,  or (1-f0)/L0 when A is empty
#   B with two or more   0
# --------------------------------------------------------------------------

#' Derivatives of a Zero-Inflated Distribution
#'
#' @description
#' Builds the order-`k` derivative method for
#' [zero_inflated()].
#'
#' @details
#' At \eqn{y > 0} the likelihood separates, \eqn{\log(1 - \zeta) + \ell(y)}, so a
#' mixed index gives zero and the pure ones are immediate. At \eqn{y = 0} it is
#' \eqn{\log L_0} with \eqn{L_0 = \zeta + (1-\zeta) f_0}, which is **affine
#' in \eqn{\zeta}** -- so any block containing two or more \eqn{\zeta}'s
#' contributes nothing, which is what keeps the partition sum small. Writing
#' \eqn{w_0 = (1-\zeta) f_0 / L_0}, the ratios are \eqn{w_0 (d^A f_0/f_0)} for a
#' block of parameters alone, and \eqn{-f_0 (d^A f_0/f_0)/L_0} for a block
#' carrying one \eqn{\zeta}.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @keywords internal
zi_deriv_k <- function(order) {
  function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
    if (expected) {
      return(expected_derivative(distrib, y, theta, order = order,
                                 approx = match.arg(approx), nsim = nsim))
    }
    parent <- distrib@parent_distrib
    params <- distrib@params
    p_names <- names(theta)[seq_len(distrib@n_params - 1L)]
    zi_name <- params[distrib@n_params]
    pars <- theta[seq_len(distrib@n_params - 1L)]
    zi <- theta[[distrib@n_params]]
    n <- length(y)

    f0 <- distrib_pdf(parent, 0, pars)
    L0 <- zi + (1 - zi) * f0
    w0 <- (1 - zi) * f0 / L0

    ell_0 <- parent_ell(parent, 0, pars, order, p_names)

    ratio <- memo_ratio(function(block) {
      sp <- split_index(block, zi_name)
      if (sp$n_mix >= 2L) return(0)
      if (sp$n_mix == 0L) return(w0 * bell_f_ratio(sp$theta, ell_0))
      if (length(sp$theta) == 0L) return((1 - f0) / L0)
      -f0 * bell_f_ratio(sp$theta, ell_0) / L0
    }, params)

    assemble_deriv(distrib, order, function(idx) {
      sp <- split_index(idx, zi_name)
      at_zero <- log_deriv(idx, ratio)
      at_pos <- if (sp$n_mix == 0L) {
        # y > 0 contributes log(1-zi) + l(y): the theta part is the parent's
        distrib_deriv_component(parent, y, pars, idx, p_names, order)
      } else if (length(sp$theta) == 0L) {
        log_pow_deriv(zi, order, complement = TRUE)
      } else {
        0
      }
      ifelse(y == 0, rep_len(at_zero, n), rep_len(at_pos, n))
    })
  }
}

#' One Component of the Parent's Derivative
#'
#' @description
#' Fetches a single component of the parent's derivative of a given order, by
#' canonical key.
#'
#' @param parent The parent distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param idx A character vector of parameter names, the multi-index.
#' @param params The parent's parameter names, in declaration order.
#' @param order The derivative order.
#'
#' @return A numeric vector.
#'
#' @keywords internal
distrib_deriv_component <- function(parent, y, theta, idx, params, order) {
  key <- canon_key(idx, params)
  switch(as.character(order),
    "3" = distrib_deriv3(parent, y, theta)[[key]],
    "4" = distrib_deriv4(parent, y, theta)[[key]]
  )
}


# --------------------------------------------------------------------------
# Zero-adjustment, discrete (hurdle)
#
# The likelihood separates completely, so every mixed index vanishes at every
# order. At y > 0 the theta part is the parent's derivative minus that of
# log(1 - f0), whose ratios are d^B(1-f0)/(1-f0) = -(d^B f0/f0) f0/(1-f0).
# --------------------------------------------------------------------------

#' Derivatives of a Zero-Adjusted Discrete Distribution
#'
#' @description
#' Builds the order-`k` derivative method for the hurdle form of
#' [zero_adjusted()].
#'
#' @details
#' The likelihood separates completely, so every mixed index vanishes at every
#' order. At \eqn{y > 0} the parameter part is the parent's derivative minus that
#' of \eqn{\log(1 - f_0)}, whose ratios are
#' \eqn{d^B(1-f_0)/(1-f_0) = -(d^B f_0/f_0) f_0/(1-f_0)}.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @keywords internal
za_disc_deriv_k <- function(order) {
  function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
    if (expected) {
      return(expected_derivative(distrib, y, theta, order = order,
                                 approx = match.arg(approx), nsim = nsim))
    }
    parent <- distrib@parent_distrib
    params <- distrib@params
    p_names <- names(theta)[seq_len(distrib@n_params - 1L)]
    za_name <- params[distrib@n_params]
    pars <- theta[seq_len(distrib@n_params - 1L)]
    za <- theta[[distrib@n_params]]
    n <- length(y)

    f0 <- distrib_pdf(parent, 0, pars)
    ell_0 <- parent_ell(parent, 0, pars, order, p_names)

    ratio <- memo_ratio(function(block) {
      -f0 * bell_f_ratio(block, ell_0) / (1 - f0)
    }, params)

    assemble_deriv(distrib, order, function(idx) {
      sp <- split_index(idx, za_name)
      if (sp$n_mix > 0L && length(sp$theta) > 0L) return(rep(0, n))

      if (length(sp$theta) == 0L) {
        # pure za: log(za) at zero, log(1-za) above it
        return(ifelse(y == 0,
                      log_pow_deriv(za, order, complement = FALSE),
                      log_pow_deriv(za, order, complement = TRUE)))
      }
      # pure theta: zero at y = 0, parent minus the truncation constant above
      at_pos <- distrib_deriv_component(parent, y, pars, idx, p_names, order) -
        rep_len(log_deriv(idx, ratio), n)
      ifelse(y == 0, 0, at_pos)
    })
  }
}


# --------------------------------------------------------------------------
# Zero-adjustment, continuous
#
# No truncation constant, so the theta part is simply the parent's, switched
# off at the atom.
# --------------------------------------------------------------------------

#' Derivatives of a Zero-Adjusted Continuous Distribution
#'
#' @description
#' Builds the order-`k` derivative method for the mixed form of
#' [zero_adjusted()].
#'
#' @details
#' A continuous parent puts no mass at zero, so there is no truncation constant
#' to correct for: the parameter part is simply the parent's derivative, switched
#' off at the atom.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @keywords internal
za_cont_deriv_k <- function(order) {
  function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
    if (expected) {
      return(expected_derivative(distrib, y, theta, order = order,
                                 approx = match.arg(approx), nsim = nsim))
    }
    parent <- distrib@parent_distrib
    params <- distrib@params
    p_names <- names(theta)[seq_len(distrib@n_params - 1L)]
    za_name <- params[distrib@n_params]
    pars <- theta[seq_len(distrib@n_params - 1L)]
    za <- theta[[distrib@n_params]]
    n <- length(y)

    assemble_deriv(distrib, order, function(idx) {
      sp <- split_index(idx, za_name)
      if (sp$n_mix > 0L && length(sp$theta) > 0L) return(rep(0, n))
      if (length(sp$theta) == 0L) {
        return(ifelse(y == 0,
                      log_pow_deriv(za, order, complement = FALSE),
                      log_pow_deriv(za, order, complement = TRUE)))
      }
      at_pos <- distrib_deriv_component(parent, y, pars, idx, p_names, order)
      ifelse(y == 0, 0, at_pos)
    })
  }
}


# --------------------------------------------------------------------------
# Truncation
#
# l_T = l - log Z, and the ratios are truncated expectations of the complete
# Bell polynomial:  d^B Z / Z = E_T[ d^B f / f ].  Each distinct block costs
# one quadrature or summation, which is why they are memoized.
# --------------------------------------------------------------------------

#' Derivatives of a Truncated Distribution
#'
#' @description
#' Builds the order-`k` derivative method for [truncated()].
#'
#' @details
#' Here \eqn{\ell_T = \ell - \log Z}, and the ratios are truncated expectations
#' of the same complete Bell quantity the other wrappers use,
#' \eqn{d^B Z / Z = \mathbb{E}_T[d^B f / f]}. Each distinct block costs one
#' quadrature or summation, which is why they are memoized across the partition
#' sum.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @seealso [memo_ratio()]
#' @keywords internal
trunc_deriv_k <- function(order) {
  function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
    if (expected) {
      return(expected_derivative(distrib, y, theta, order = order,
                                 approx = match.arg(approx), nsim = nsim))
    }
    parent <- distrib@parent_distrib
    params <- distrib@params
    n <- length(y)

    # d^B Z / Z. Where the parent has genuine cdf derivatives of the block's
    # order, Z's derivative is F(u) - F(l^-) differentiated, which is two
    # calls on the parent instead of one quadrature per component. The tables
    # are fetched once per order rather than once per block, and a NULL means
    # that order has no such route and the block falls back to the
    # expectation.
    dZ <- lapply(seq_len(order), function(k) trunc_mass_derivs(distrib, theta, k))
    Zval <- trunc_constants(distrib, theta)$Z
    ratio <- memo_ratio(function(block) {
      k <- length(block)
      if (!is.null(dZ[[k]])) return(dZ[[k]][[canon_key(block, params)]] / Zval)
      # expectation() calls its integrand as f(y = ., theta = .), by name
      expectation(distrib, function(y, theta) {
        bell_f_ratio(block, parent_ell(parent, y, theta, length(block), params))
      }, theta)
    }, params)

    assemble_deriv(distrib, order, function(idx) {
      own <- distrib_deriv_component(parent, y, theta, idx, params, order)
      rep_len(own - rep_len(log_deriv(idx, ratio), n), n)
    })
  }
}


# --------------------------------------------------------------------------
# Registration
# --------------------------------------------------------------------------

#' @title Transformed Third Derivatives
#' @name distrib_deriv3.TransformedDistrib
#'
#' @description
#' Exactly the parent's third derivatives, read at the preimage
#' \eqn{x = g^{-1}(y)}. A transformation of the response contributes
#' \eqn{\log|dg^{-1}/dy|} to the log-density, and that Jacobian does not
#' depend on \eqn{\theta}, so it vanishes from every derivative in
#' \eqn{\theta} and leaves the parent's untouched.
#'
#' Measured on `transformation(gaussian1_distrib(), exp_transform())`, the
#' whole component list is `all.equal` to the Gaussian's at `log(y)`, at this
#' order and at the fourth.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param y A numeric vector of observations on the transformed scale.
#' @param theta A named list of the **parent's** parameters; a transformation
#'   adds none.
#' @param expected Logical of length 1. Passed to the parent, so the parent
#'   decides whether an expectation is closed or approximated.
#' @param ... Passed to the parent's method, `approx` and `nsim` among them.
#'
#' @return A named list of third-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters, each a
#'   numeric vector of length `length(y)`.
#'
#' @examples
#' # A lognormal built by transforming a Gaussian.
#' logn <- transformation(gaussian1_distrib(), exp_transform())
#' y <- c(1.2, 2.5)
#' th <- list(mu = 0.3, sigma = 1.1)
#'
#' distrib_deriv3(logn, y, th)
#'
#' # It is the parent's answer at the preimage, component for component.
#' all.equal(distrib_deriv3(logn, y, th),
#'           distrib_deriv3(gaussian1_distrib(), log(y), th))
#'
#' @seealso [transformation()] for the wrapper;
#'   [distrib_deriv4.TransformedDistrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
S7::method(distrib_deriv3, TransformedDistrib) <- trans_deriv_k(3L)

#' @title Transformed Fourth Derivatives
#' @name distrib_deriv4.TransformedDistrib
#'
#' @description
#' Exactly the parent's fourth derivatives, read at the preimage
#' \eqn{x = g^{-1}(y)}, for the same reason as at third order: the Jacobian
#' \eqn{\log|dg^{-1}/dy|} carries no \eqn{\theta}, so it contributes nothing to
#' any derivative in \eqn{\theta} whatever the order.
#'
#' Nothing here is approximated and no partition sum is taken. The wrapper's
#' work is entirely in the argument, not in the derivative.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param y A numeric vector of observations on the transformed scale.
#' @param theta A named list of the **parent's** parameters; a transformation
#'   adds none.
#' @param expected Logical of length 1. Passed to the parent.
#' @param ... Passed to the parent's method.
#'
#' @return A named list of fourth-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters, each a
#'   numeric vector of length `length(y)`. A two-parameter parent gives five.
#'
#' @examples
#' logn <- transformation(gaussian1_distrib(), exp_transform())
#' y <- c(1.2, 2.5)
#' th <- list(mu = 0.3, sigma = 1.1)
#'
#' names(distrib_deriv4(logn, y, th))
#' all.equal(distrib_deriv4(logn, y, th),
#'           distrib_deriv4(gaussian1_distrib(), log(y), th))
#'
#' @seealso [transformation()] for the wrapper;
#'   [distrib_deriv3.TransformedDistrib()] for the order below;
#'   [distrib_deriv4()] for the generic.
S7::method(distrib_deriv4, TransformedDistrib) <- trans_deriv_k(4L)

#' @title Zero-Inflated Third Derivatives
#' @name distrib_deriv3.ZeroInflatedDistrib
#'
#' @description
#' Two regimes, and the observation decides which. At \eqn{y > 0} the
#' log-likelihood is \eqn{\log(1-\zeta) + \ell(y;\theta)}, which **separates**:
#' every mixed component in \eqn{\zeta} and a parent parameter is exactly zero.
#' At \eqn{y = 0} it is \eqn{\log L_0} with
#' \eqn{L_0 = \zeta + (1-\zeta) f(0;\theta)}, and the derivatives follow from
#' the moment-to-cumulant expansion over set partitions,
#' \deqn{d^I \log L_0 = \sum_{\pi} (-1)^{|\pi|-1}(|\pi|-1)!
#'       \prod_{B \in \pi} \frac{d^B L_0}{L_0}.}
#'
#' \eqn{L_0} is **affine** in \eqn{\zeta}, so any block of a partition naming
#' \eqn{\zeta} twice contributes nothing and the sum is shorter than it looks.
#' A whole component may still be non-zero: its partition into singletons
#' survives.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param y A numeric vector of counts.
#' @param theta A named list with the parent's parameters followed by `zi`, the
#'   inflation probability in \eqn{(0, 1)}.
#' @param expected Logical of length 1. `TRUE` takes the expectation, which for
#'   a lattice parent is an exact sum over the support.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of third-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters and `zi`,
#'   each a numeric vector of length `length(y)`.
#'
#' @examples
#' zi <- zero_inflated(poisson_distrib())
#' th <- list(mu = 3, zi = 0.2)
#'
#' # At a positive count the likelihood separates, so a mixed component is 0.
#' distrib_deriv3(zi, 2, th)[["mu_zi_zi"]]
#'
#' # At zero it is not: log L0 mixes the two.
#' round(unlist(distrib_deriv3(zi, 0, th)), 4)
#'
#' @seealso [zero_inflated()] for the wrapper and for why it cannot be
#'   stacked; [distrib_deriv4.ZeroInflatedDistrib()] for the order above;
#'   [log_deriv()], the partition sum used here.
S7::method(distrib_deriv3, ZeroInflatedDistrib) <- zi_deriv_k(3L)

#' @title Zero-Inflated Fourth Derivatives
#' @name distrib_deriv4.ZeroInflatedDistrib
#'
#' @description
#' The same two regimes as at third order, one order along. At \eqn{y > 0} the
#' log-likelihood separates into \eqn{\log(1-\zeta)} and the parent's, so a
#' component mixing \eqn{\zeta} with a parent parameter is exactly zero; at
#' \eqn{y = 0} the moment-to-cumulant sum over the partitions of a four-index
#' set gives \eqn{d^I \log L_0}, with every block naming \eqn{\zeta} twice
#' contributing nothing because \eqn{L_0} is affine in it.
#'
#' The partitions of four indices number fifteen against five for three, so
#' this order is the more expensive of the two by that ratio and by nothing
#' else: the same identity and the same ratios \eqn{d^B L_0 / L_0} are read.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param y A numeric vector of counts.
#' @param theta A named list with the parent's parameters followed by `zi`, the
#'   inflation probability in \eqn{(0, 1)}.
#' @param expected Logical of length 1. `TRUE` takes the expectation.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of fourth-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters and `zi`,
#'   each a numeric vector of length `length(y)`. A one-parameter parent gives
#'   five.
#'
#' @examples
#' zi <- zero_inflated(poisson_distrib())
#' th <- list(mu = 3, zi = 0.2)
#'
#' names(distrib_deriv4(zi, 0, th))
#' round(unlist(distrib_deriv4(zi, 0, th)), 4)
#'
#' # Separation at a positive count holds at this order too.
#' distrib_deriv4(zi, 2, th)[["mu_mu_zi_zi"]]
#'
#' @seealso [zero_inflated()] for the wrapper;
#'   [distrib_deriv3.ZeroInflatedDistrib()] for the order below;
#'   [log_deriv()], the partition sum used here.
S7::method(distrib_deriv4, ZeroInflatedDistrib) <- zi_deriv_k(4L)

#' @title Hurdle Third Derivatives
#' @name distrib_deriv3.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' The hurdle likelihood **factorizes** into a binary part and a positive part,
#' so at every order the mixed components in `za` and a parent parameter are
#' exactly zero, and the two halves are differentiated separately. The
#' \eqn{\theta} part is the parent's derivative less that of the truncation
#' constant \eqn{\log(1 - f_0)}, which is the mass the parent puts at zero and
#' the hurdle removes.
#'
#' This is the structural difference from [zero_inflated()], and it is one a
#' reader can check in a line: there the mixed block is not zero, because
#' inflation *adds* to the parent's mass at zero and no single zero can be
#' attributed to a mechanism.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of counts.
#' @param theta A named list with the parent's parameters followed by `za`, the
#'   probability of a zero, in \eqn{(0, 1)}.
#' @param expected Logical of length 1. `TRUE` takes the expectation, an exact
#'   sum over the support for a lattice parent.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of third-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters and `za`,
#'   each a numeric vector of length `length(y)`.
#'
#' @examples
#' h <- zero_adjusted(poisson_distrib())
#' th <- list(mu = 3, za = 0.25)
#'
#' round(unlist(distrib_deriv3(h, c(0, 2, 5), th)[["mu_mu_mu"]]), 6)
#'
#' # The mixed block is exactly zero at every observation, which the
#' # zero-inflated wrapper's is not.
#' distrib_deriv3(h, c(0, 2, 5), th)[["mu_mu_za"]]
#' distrib_deriv3(zero_inflated(poisson_distrib()), 0,
#'                list(mu = 3, zi = 0.25))[["mu_mu_zi"]]
#'
#' @seealso [zero_adjusted()] for the wrapper and the counting rule it
#'   enforces; [distrib_deriv4.ZeroAdjustedDiscreteDistrib()] for the order
#'   above; [distrib_deriv3.ZeroInflatedDistrib()] for the model this is not.
S7::method(distrib_deriv3, ZeroAdjustedDiscreteDistrib) <- za_disc_deriv_k(3L)

#' @title Hurdle Fourth Derivatives
#' @name distrib_deriv4.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' The factorization of the hurdle likelihood holds at every order, so this is
#' the third-order picture one step along: the components mixing `za` with a
#' parent parameter are exactly zero, the `za` part is the fourth derivative of
#' a binomial log-likelihood in one probability, and the \eqn{\theta} part is
#' the parent's fourth derivative less that of \eqn{\log(1 - f_0)}.
#'
#' The truncation constant is where the work is. Its fourth derivative comes
#' from the moment-to-cumulant expansion over the fifteen partitions of four
#' indices, read on the ratios \eqn{d^B f_0 / f_0} alone.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of counts.
#' @param theta A named list with the parent's parameters followed by `za`, the
#'   probability of a zero, in \eqn{(0, 1)}.
#' @param expected Logical of length 1. `TRUE` takes the expectation.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of fourth-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters and `za`,
#'   each a numeric vector of length `length(y)`.
#'
#' @examples
#' h <- zero_adjusted(poisson_distrib())
#' th <- list(mu = 3, za = 0.25)
#'
#' names(distrib_deriv4(h, c(0, 2, 5), th))
#'
#' # Every mixed component vanishes, at this order as at the one below.
#' vapply(c("mu_mu_mu_za", "mu_mu_za_za", "mu_za_za_za"),
#'        function(k) distrib_deriv4(h, 2, th)[[k]], numeric(1))
#'
#' @seealso [zero_adjusted()] for the wrapper;
#'   [distrib_deriv3.ZeroAdjustedDiscreteDistrib()] for the order below;
#'   [log_deriv()], the partition sum the truncation constant uses.
S7::method(distrib_deriv4, ZeroAdjustedDiscreteDistrib) <- za_disc_deriv_k(4L)

#' @title Zero-Adjusted Third Derivatives, Continuous Parent
#' @name distrib_deriv3.ZeroAdjustedContinuousDistrib
#'
#' @description
#' A continuous parent puts **no** mass at zero, so there is nothing to
#' truncate away and no normalizing constant to differentiate. Away from the
#' atom the log-likelihood is \eqn{\log(1 - \pi) + \ell(y;\theta)}: the
#' \eqn{\theta} derivatives are the parent's unchanged, and every component
#' mixing `za` with a parent parameter is exactly zero.
#'
#' The object is a **mixed** distribution, a density on the positive line plus
#' an atom at zero, and it declares that atom through [distrib_atoms()], so
#' [check_distrib()] tests it correctly: the density integrates to
#' \eqn{1 - \pi}, not to 1.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations, zero included.
#' @param theta A named list with the parent's parameters followed by `za`, the
#'   probability of the atom, in \eqn{(0, 1)}.
#' @param expected Logical of length 1. `TRUE` takes the expectation over the
#'   mixed law, atom included.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of third-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters and `za`,
#'   each a numeric vector of length `length(y)`. A two-parameter parent gives
#'   ten.
#'
#' @examples
#' za <- zero_adjusted(gamma2_distrib())
#' th <- list(mu = 2, sigma2 = 0.7, za = 0.25)
#'
#' names(distrib_deriv3(za, c(0, 1.5, 3), th))
#'
#' # Away from the atom the theta derivatives are the parent's, exactly.
#' all.equal(distrib_deriv3(za, c(1.5, 3), th)[["mu_mu_mu"]],
#'           distrib_deriv3(gamma2_distrib(), c(1.5, 3),
#'                          list(mu = 2, sigma2 = 0.7))[["mu_mu_mu"]])
#'
#' # And the mixed components vanish.
#' distrib_deriv3(za, c(0, 1.5, 3), th)[["mu_mu_za"]]
#'
#' @seealso [zero_adjusted()] for the wrapper;
#'   [distrib_atoms()], which declares the atom;
#'   [distrib_deriv4.ZeroAdjustedContinuousDistrib()] for the order above;
#'   [distrib_deriv3.ZeroAdjustedDiscreteDistrib()], where a truncation
#'   constant does appear.
S7::method(distrib_deriv3, ZeroAdjustedContinuousDistrib) <- za_cont_deriv_k(3L)

#' @title Zero-Adjusted Fourth Derivatives, Continuous Parent
#' @name distrib_deriv4.ZeroAdjustedContinuousDistrib
#'
#' @description
#' The same separation as at third order, one step along, and for the same
#' reason: a continuous parent has no mass at zero, so no normalizing constant
#' enters and the log-likelihood is \eqn{\log(1 - \pi) + \ell(y;\theta)} away
#' from the atom. The \eqn{\theta} components are the parent's fourth
#' derivatives unchanged, the pure `za` component is that of
#' \eqn{\log(1-\pi)}, and everything mixing the two is exactly zero.
#'
#' No partition sum is taken here at all, so this is the cheapest of the six
#' wrappers at this order.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations, zero included.
#' @param theta A named list with the parent's parameters followed by `za`, the
#'   probability of the atom, in \eqn{(0, 1)}.
#' @param expected Logical of length 1. `TRUE` takes the expectation over the
#'   mixed law.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of fourth-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters and `za`,
#'   each a numeric vector of length `length(y)`. A two-parameter parent gives
#'   fifteen.
#'
#' @examples
#' za <- zero_adjusted(gamma2_distrib())
#' th <- list(mu = 2, sigma2 = 0.7, za = 0.25)
#'
#' length(distrib_deriv4(za, c(0, 1.5, 3), th))
#' all.equal(distrib_deriv4(za, c(1.5, 3), th)[["mu_mu_mu_mu"]],
#'           distrib_deriv4(gamma2_distrib(), c(1.5, 3),
#'                          list(mu = 2, sigma2 = 0.7))[["mu_mu_mu_mu"]])
#'
#' @seealso [zero_adjusted()] for the wrapper;
#'   [distrib_deriv3.ZeroAdjustedContinuousDistrib()] for the order below;
#'   [distrib_deriv4.ZeroAdjustedDiscreteDistrib()], where a truncation
#'   constant does appear.
S7::method(distrib_deriv4, ZeroAdjustedContinuousDistrib) <- za_cont_deriv_k(4L)

#' @title Truncated Third Derivatives
#' @name distrib_deriv3.TruncatedContinuousDistrib
#'
#' @description
#' Truncation adds no parameter and one \eqn{\theta}-dependent normalizing
#' constant: \eqn{\ell_T = \ell - \log Z} with
#' \eqn{Z(\theta) = F(U) - F(L^-)}. The derivatives of \eqn{\log Z} follow from
#' the moment-to-cumulant expansion over set partitions, read on the truncated
#' expectations \eqn{\mathbb{E}_T[\partial^B f / f]}:
#' \deqn{d^I \log Z = \sum_{\pi} (-1)^{|\pi|-1}(|\pi|-1)!
#'       \prod_{B \in \pi} \mathbb{E}_T\!\left[\frac{\partial^B f}{f}\right].}
#'
#' Each **distinct block** costs one quadrature, memoized across the partition
#' sum, which is why a truncated derivative is far dearer than its parent's.
#' Where the parent has genuinely closed cdf derivatives, or is a lattice
#' family whose cdf derivatives are an exact sum, that route replaces the
#' quadrature and is about three and a half times faster.
#'
#' @param distrib A `TruncatedContinuousDistrib` object, from [truncated()].
#' @param y A numeric vector of observations inside the truncation interval.
#' @param theta A named list of the **parent's** parameters; truncation adds
#'   none, the endpoints being constants.
#' @param expected Logical of length 1. `TRUE` takes the expectation under the
#'   truncated law.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of third-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters, each a
#'   numeric vector of length `length(y)`.
#'
#' @examples
#' tr <- truncated(gaussian1_distrib(), lower = -1, upper = 3)
#' th <- list(mu = 0.3, sigma = 1.1)
#'
#' round(unlist(distrib_deriv3(tr, c(0, 1, 2), th)[["mu_mu_mu"]]), 6)
#'
#' # log Z is what separates it from the parent: the parent's third derivative
#' # in mu is identically zero and the truncated one is not.
#' c(truncated = distrib_deriv3(tr, 1, th)[["mu_mu_mu"]],
#'   parent = distrib_deriv3(gaussian1_distrib(), 1, th)[["mu_mu_mu"]])
#'
#' @seealso [truncated()] for the wrapper and the endpoint convention;
#'   [distrib_deriv4.TruncatedContinuousDistrib()] for the order above;
#'   [distrib_grad_cdf()], the cheaper route where the parent has one.
S7::method(distrib_deriv3, TruncatedContinuousDistrib) <- trunc_deriv_k(3L)

#' @title Truncated Fourth Derivatives
#' @name distrib_deriv4.TruncatedContinuousDistrib
#'
#' @description
#' The third-order construction one step along: \eqn{\ell_T = \ell - \log Z},
#' and \eqn{d^I \log Z} comes from the moment-to-cumulant expansion over the
#' fifteen partitions of four indices, each block an expectation
#' \eqn{\mathbb{E}_T[\partial^B f / f]} under the truncated law.
#'
#' The cost is set by the number of **distinct** blocks, not by the number of
#' partitions: a block met twice is computed once and looked up. At four
#' indices over two parameters that is five distinct blocks against fifteen
#' partitions.
#'
#' @param distrib A `TruncatedContinuousDistrib` object, from [truncated()].
#' @param y A numeric vector of observations inside the truncation interval.
#' @param theta A named list of the **parent's** parameters; truncation adds
#'   none.
#' @param expected Logical of length 1. `TRUE` takes the expectation under the
#'   truncated law.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of fourth-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters, each a
#'   numeric vector of length `length(y)`. A two-parameter parent gives five.
#'
#' @examples
#' tr <- truncated(gaussian1_distrib(), lower = -1, upper = 3)
#' th <- list(mu = 0.3, sigma = 1.1)
#'
#' names(distrib_deriv4(tr, c(0, 1, 2), th))
#' round(unlist(distrib_deriv4(tr, 1, th)), 5)
#'
#' @seealso [truncated()] for the wrapper;
#'   [distrib_deriv3.TruncatedContinuousDistrib()] for the order below;
#'   [log_deriv()], the partition sum both orders use.
S7::method(distrib_deriv4, TruncatedContinuousDistrib) <- trunc_deriv_k(4L)

#' @title Truncated Third Derivatives, Discrete Parent
#' @name distrib_deriv3.TruncatedDiscreteDistrib
#'
#' @description
#' The same construction as for a continuous parent, \eqn{\ell_T = \ell -
#' \log Z} with \eqn{d^I \log Z} from the moment-to-cumulant expansion over
#' set partitions, with every expectation taken by **summation over the
#' retained support** in place of a quadrature.
#'
#' That makes it both exact and cheap. The support of a truncated lattice
#' family is finite whenever both endpoints are, so each block is a finite sum
#' of terms the parent already computes, and the memoization across the
#' partition sum applies as before.
#'
#' @param distrib A `TruncatedDiscreteDistrib` object, from [truncated()].
#' @param y A numeric vector of counts inside the truncation interval. Both
#'   endpoints are **included**.
#' @param theta A named list of the **parent's** parameters; truncation adds
#'   none.
#' @param expected Logical of length 1. `TRUE` takes the expectation under the
#'   truncated law, itself an exact sum.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of third-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters, each a
#'   numeric vector of length `length(y)`. A one-parameter parent gives one.
#'
#' @examples
#' tr <- truncated(poisson_distrib(), lower = 1, upper = 8)
#' th <- list(mu = 3)
#'
#' distrib_deriv3(tr, c(2, 4, 6), th)
#'
#' # The parent's third derivative in mu is 2 y / mu^3 - ... and the truncated
#' # one differs by the derivative of log Z, which is the whole wrapper.
#' c(truncated = distrib_deriv3(tr, 4, th)[["mu_mu_mu"]],
#'   parent = distrib_deriv3(poisson_distrib(), 4, th)[["mu_mu_mu"]])
#'
#' @seealso [truncated()] for the wrapper and the endpoint convention;
#'   [distrib_deriv4.TruncatedDiscreteDistrib()] for the order above;
#'   [distrib_deriv3.TruncatedContinuousDistrib()], where the same blocks cost
#'   a quadrature apiece.
S7::method(distrib_deriv3, TruncatedDiscreteDistrib) <- trunc_deriv_k(3L)

#' @title Truncated Fourth Derivatives, Discrete Parent
#' @name distrib_deriv4.TruncatedDiscreteDistrib
#'
#' @description
#' The third-order construction one step along, with the expectations again
#' taken by summation over the retained support: \eqn{\ell_T = \ell - \log Z},
#' and \eqn{d^I \log Z} from the moment-to-cumulant expansion over the fifteen
#' partitions of four indices.
#'
#' Because every block is a finite sum, this order costs what the order below
#' costs times the ratio of distinct blocks, and nothing is approximated
#' anywhere in it.
#'
#' @param distrib A `TruncatedDiscreteDistrib` object, from [truncated()].
#' @param y A numeric vector of counts inside the truncation interval. Both
#'   endpoints are **included**.
#' @param theta A named list of the **parent's** parameters; truncation adds
#'   none.
#' @param expected Logical of length 1. `TRUE` takes the expectation under the
#'   truncated law.
#' @param ... Passed on, `approx` and `nsim` among them.
#'
#' @return A named list of fourth-derivative components, keyed
#'   lexicographically by [deriv_names()] on the parent's parameters, each a
#'   numeric vector of length `length(y)`.
#'
#' @examples
#' tr <- truncated(poisson_distrib(), lower = 1, upper = 8)
#' th <- list(mu = 3)
#'
#' distrib_deriv4(tr, c(2, 4, 6), th)
#'
#' # A truncated lattice family with both endpoints finite has a finite
#' # support, so every expectation behind this is an exact sum.
#' sum(distrib_pdf(tr, 1:8, th))
#'
#' @seealso [truncated()] for the wrapper;
#'   [distrib_deriv3.TruncatedDiscreteDistrib()] for the order below;
#'   [log_deriv()], the partition sum both orders use.
S7::method(distrib_deriv4, TruncatedDiscreteDistrib) <- trunc_deriv_k(4L)

