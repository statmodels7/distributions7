#' @include distrib.R generics.R utility_functions.R expected_derivatives.R link_scale.R

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

# Set partitions of a multi-index, given as a character vector of parameter
# names with repetition. Returns a list of partitions, each a list of blocks.
index_partitions <- function(idx) {
  lapply(set_partitions(length(idx)), function(p) lapply(p, function(b) idx[b]))
}

# The canonical component name of a block: parameters in the order the
# distribution declares them, joined by "_", exactly as deriv_names() builds it.
canon_key <- function(block, params) {
  paste(block[order(match(block, params))], collapse = "_")
}

# (1) d^I f / f from the parent's log-derivatives.
bell_f_ratio <- function(idx, ell) {
  total <- 0
  for (p in index_partitions(idx)) {
    term <- 1
    for (b in p) term <- term * ell(b)
    total <- total + term
  }
  total
}

# (2) d^I log L from the ratios d^B L / L.
log_deriv <- function(idx, ratio) {
  total <- 0
  for (p in index_partitions(idx)) {
    k <- length(p)
    term <- (-1)^(k - 1) * gamma(k)          # (-1)^{k-1} (k-1)!
    for (b in p) term <- term * ratio(b)
    total <- total + term
  }
  total
}

# The parent's derivative components to a given order, and a lookup keyed by
# block. Orders are fetched once per call, not once per block.
parent_ell <- function(parent, y, theta, max_order, params) {
  d <- vector("list", max_order)
  d[[1]] <- distrib_gradient(parent, y, theta)
  if (max_order >= 2) d[[2]] <- distrib_hessian(parent, y, theta)
  if (max_order >= 3) d[[3]] <- distrib_deriv3(parent, y, theta)
  if (max_order >= 4) d[[4]] <- distrib_deriv4(parent, y, theta)
  function(block) d[[length(block)]][[canon_key(block, params)]]
}

# Memoise a ratio function on the canonical key of its block: a partition of a
# fourth-order index asks for the same small blocks many times over, and for the
# truncated wrapper each distinct block costs a quadrature.
memo_ratio <- function(f, params) {
  cache <- list()
  function(block) {
    k <- canon_key(block, params)
    if (is.null(cache[[k]])) cache[[k]] <<- f(block)
    cache[[k]]
  }
}

# The multi-indices of a given order, in the order deriv_names() lists them.
#
# Built here rather than taken from deriv_index_list(), whose order-2 case is
# deliberately ordered for hess_names() -- diagonal first -- while deriv_names()
# is lexicographic. Pairing the two would silently attach the name "mu_sigma" to
# the index (sigma, sigma). The orders actually registered are 3 and 4, where the
# two agree, but a mismatch that only bites if someone reuses this helper is
# exactly the kind worth removing rather than commenting on.
order_indices <- function(params, order) {
  p <- length(params)
  idx <- as.matrix(do.call(expand.grid, rep(list(seq_len(p)), order)))
  idx <- idx[, rev(seq_len(order)), drop = FALSE]
  idx <- idx[apply(idx, 1L, function(r) all(diff(r) >= 0)), , drop = FALSE]
  lapply(seq_len(nrow(idx)), function(k) params[idx[k, ]])
}

# d^k log(x) / dx^k evaluated through a chain of the form log(a + b*p): the two
# cases the wrappers need are log(pi) and log(1 - pi).
log_pow_deriv <- function(p, k, complement = FALSE) {
  if (complement) -gamma(k) / (1 - p)^k        # d^k log(1-p) = -(k-1)!/(1-p)^k
  else (-1)^(k - 1) * gamma(k) / p^k           # d^k log(p)   = (-1)^{k-1}(k-1)!/p^k
}

# Split a multi-index into its parent part and the count of the wrapper's own
# parameter, whose name is always the last one the wrapper declares.
split_index <- function(idx, mix_name) {
  list(theta = idx[idx != mix_name], n_mix = sum(idx == mix_name))
}

# Shared skeleton: assemble one named list of components of the given order.
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

# One component of the parent's derivative of the given order, by canonical key.
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
# one quadrature or summation, which is why they are memoised.
# --------------------------------------------------------------------------

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

    # expectation() calls its integrand as f(y = ., theta = .), by name
    ratio <- memo_ratio(function(block) {
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
#' @param distrib A \code{TransformedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{transformation}}
S7::method(distrib_deriv3, TransformedDistrib) <- trans_deriv_k(3L)

#' @title Transformed Fourth Derivatives
#' @name distrib_deriv4.TransformedDistrib
#' @description As \code{\link[=distrib_deriv3.TransformedDistrib]{the third}}, at fourth order.
#' @param distrib A \code{TransformedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{transformation}}
S7::method(distrib_deriv4, TransformedDistrib) <- trans_deriv_k(4L)

#' @title Zero-Inflated Third Derivatives
#' @name distrib_deriv3.ZeroInflatedDistrib
#' @description
#' At \eqn{y > 0} the likelihood separates; at \eqn{y = 0} it is \eqn{\log L_0}
#' with \eqn{L_0} affine in \eqn{\zeta}, so the derivatives follow from the
#' moment-to-cumulant expansion over set partitions.
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{zero_inflated}}
S7::method(distrib_deriv3, ZeroInflatedDistrib) <- zi_deriv_k(3L)

#' @title Zero-Inflated Fourth Derivatives
#' @name distrib_deriv4.ZeroInflatedDistrib
#' @description As \code{\link[=distrib_deriv3.ZeroInflatedDistrib]{the third}}, at fourth order.
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{zero_inflated}}
S7::method(distrib_deriv4, ZeroInflatedDistrib) <- zi_deriv_k(4L)

#' @title Hurdle Third Derivatives
#' @name distrib_deriv3.ZeroAdjustedDiscreteDistrib
#' @description
#' The likelihood separates, so mixed components vanish at every order; the
#' \eqn{\theta} part is the parent's derivative less that of the truncation
#' constant \eqn{\log(1-f_0)}.
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_deriv3, ZeroAdjustedDiscreteDistrib) <- za_disc_deriv_k(3L)

#' @title Hurdle Fourth Derivatives
#' @name distrib_deriv4.ZeroAdjustedDiscreteDistrib
#' @description As \code{\link[=distrib_deriv3.ZeroAdjustedDiscreteDistrib]{the third}}, at fourth order.
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_deriv4, ZeroAdjustedDiscreteDistrib) <- za_disc_deriv_k(4L)

#' @title Zero-Adjusted Continuous Third Derivatives
#' @name distrib_deriv3.ZeroAdjustedContinuousDistrib
#' @description
#' There is no truncation constant, so away from the atom the \eqn{\theta}
#' derivatives are the parent's and the mixed ones vanish.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_deriv3, ZeroAdjustedContinuousDistrib) <- za_cont_deriv_k(3L)

#' @title Zero-Adjusted Continuous Fourth Derivatives
#' @name distrib_deriv4.ZeroAdjustedContinuousDistrib
#' @description As \code{\link[=distrib_deriv3.ZeroAdjustedContinuousDistrib]{the third}}, at fourth order.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_deriv4, ZeroAdjustedContinuousDistrib) <- za_cont_deriv_k(4L)

#' @title Truncated Third Derivatives (Continuous)
#' @name distrib_deriv3.TruncatedContinuousDistrib
#' @description
#' \eqn{\ell_T = \ell - \log Z}, and the derivatives of \eqn{\log Z} follow from
#' the truncated expectations \eqn{\mathbb{E}_T[\partial^B f / f]} through the
#' moment-to-cumulant expansion. Each distinct block costs one quadrature.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{truncated}}
S7::method(distrib_deriv3, TruncatedContinuousDistrib) <- trunc_deriv_k(3L)

#' @title Truncated Fourth Derivatives (Continuous)
#' @name distrib_deriv4.TruncatedContinuousDistrib
#' @description As \code{\link[=distrib_deriv3.TruncatedContinuousDistrib]{the third}}, at fourth order.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{truncated}}
S7::method(distrib_deriv4, TruncatedContinuousDistrib) <- trunc_deriv_k(4L)

#' @title Truncated Third Derivatives (Discrete)
#' @name distrib_deriv3.TruncatedDiscreteDistrib
#' @description As the continuous case, with the expectations taken by summation.
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{truncated}}
S7::method(distrib_deriv3, TruncatedDiscreteDistrib) <- trunc_deriv_k(3L)

#' @title Truncated Fourth Derivatives (Discrete)
#' @name distrib_deriv4.TruncatedDiscreteDistrib
#' @description As \code{\link[=distrib_deriv3.TruncatedDiscreteDistrib]{the third}}, at fourth order.
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @return A named list of derivative components.
#' @seealso \code{\link{truncated}}
S7::method(distrib_deriv4, TruncatedDiscreteDistrib) <- trunc_deriv_k(4L)
