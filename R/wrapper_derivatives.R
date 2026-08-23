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
#' \eqn{\log L_0} with \eqn{L_0 = \zeta + (1-\zeta) f_0}, which is \strong{affine
#' in \eqn{\zeta}} -- so any block containing two or more \eqn{\zeta}'s
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
#' @description
#' Exactly the parent's, evaluated at \eqn{x = g^{-1}(y)}: the Jacobian does not
#' depend on \eqn{\theta}, so it leaves every derivative in \eqn{\theta} untouched.
#' @param distrib A `TransformedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [transformation()]
S7::method(distrib_deriv3, TransformedDistrib) <- trans_deriv_k(3L)

#' @title Transformed Fourth Derivatives
#' @name distrib_deriv4.TransformedDistrib
#' @description As [`the third()`][distrib_deriv3.TransformedDistrib], at fourth order.
#' @param distrib A `TransformedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [transformation()]
S7::method(distrib_deriv4, TransformedDistrib) <- trans_deriv_k(4L)

#' @title Zero-Inflated Third Derivatives
#' @name distrib_deriv3.ZeroInflatedDistrib
#' @description
#' At \eqn{y > 0} the likelihood separates; at \eqn{y = 0} it is \eqn{\log L_0}
#' with \eqn{L_0} affine in \eqn{\zeta}, so the derivatives follow from the
#' moment-to-cumulant expansion over set partitions.
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [zero_inflated()]
S7::method(distrib_deriv3, ZeroInflatedDistrib) <- zi_deriv_k(3L)

#' @title Zero-Inflated Fourth Derivatives
#' @name distrib_deriv4.ZeroInflatedDistrib
#' @description As [`the third()`][distrib_deriv3.ZeroInflatedDistrib], at fourth order.
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [zero_inflated()]
S7::method(distrib_deriv4, ZeroInflatedDistrib) <- zi_deriv_k(4L)

#' @title Hurdle Third Derivatives
#' @name distrib_deriv3.ZeroAdjustedDiscreteDistrib
#' @description
#' The likelihood separates, so mixed components vanish at every order; the
#' \eqn{\theta} part is the parent's derivative less that of the truncation
#' constant \eqn{\log(1-f_0)}.
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `za`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [zero_adjusted()]
S7::method(distrib_deriv3, ZeroAdjustedDiscreteDistrib) <- za_disc_deriv_k(3L)

#' @title Hurdle Fourth Derivatives
#' @name distrib_deriv4.ZeroAdjustedDiscreteDistrib
#' @description As [`the third()`][distrib_deriv3.ZeroAdjustedDiscreteDistrib], at fourth order.
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `za`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [zero_adjusted()]
S7::method(distrib_deriv4, ZeroAdjustedDiscreteDistrib) <- za_disc_deriv_k(4L)

#' @title Zero-Adjusted Continuous Third Derivatives
#' @name distrib_deriv3.ZeroAdjustedContinuousDistrib
#' @description
#' There is no truncation constant, so away from the atom the \eqn{\theta}
#' derivatives are the parent's and the mixed ones vanish.
#' @param distrib A `ZeroAdjustedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `za`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [zero_adjusted()]
S7::method(distrib_deriv3, ZeroAdjustedContinuousDistrib) <- za_cont_deriv_k(3L)

#' @title Zero-Adjusted Continuous Fourth Derivatives
#' @name distrib_deriv4.ZeroAdjustedContinuousDistrib
#' @description As [`the third()`][distrib_deriv3.ZeroAdjustedContinuousDistrib], at fourth order.
#' @param distrib A `ZeroAdjustedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `za`.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [zero_adjusted()]
S7::method(distrib_deriv4, ZeroAdjustedContinuousDistrib) <- za_cont_deriv_k(4L)

#' @title Truncated Third Derivatives (Continuous)
#' @name distrib_deriv3.TruncatedContinuousDistrib
#' @description
#' \eqn{\ell_T = \ell - \log Z}, and the derivatives of \eqn{\log Z} follow from
#' the truncated expectations \eqn{\mathbb{E}_T[\partial^B f / f]} through the
#' moment-to-cumulant expansion. Each distinct block costs one quadrature.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [truncated()]
S7::method(distrib_deriv3, TruncatedContinuousDistrib) <- trunc_deriv_k(3L)

#' @title Truncated Fourth Derivatives (Continuous)
#' @name distrib_deriv4.TruncatedContinuousDistrib
#' @description As [`the third()`][distrib_deriv3.TruncatedContinuousDistrib], at fourth order.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [truncated()]
S7::method(distrib_deriv4, TruncatedContinuousDistrib) <- trunc_deriv_k(4L)

#' @title Truncated Third Derivatives (Discrete)
#' @name distrib_deriv3.TruncatedDiscreteDistrib
#' @description As the continuous case, with the expectations taken by summation.
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [truncated()]
S7::method(distrib_deriv3, TruncatedDiscreteDistrib) <- trunc_deriv_k(3L)

#' @title Truncated Fourth Derivatives (Discrete)
#' @name distrib_deriv4.TruncatedDiscreteDistrib
#' @description As [`the third()`][distrib_deriv3.TruncatedDiscreteDistrib], at fourth order.
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if `TRUE`, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso [truncated()]
S7::method(distrib_deriv4, TruncatedDiscreteDistrib) <- trunc_deriv_k(4L)

