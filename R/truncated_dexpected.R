#' @include truncated.R dexpected_families.R wrapper_derivatives.R
NULL

#' The Derivatives of a Truncated Family's Expected Information
#'
#' @description
#' Computes the first (`k = 1`) or second (`k = 2`) derivatives of the
#' expected information of a [truncated()] family on the parameter scale, as
#' expectations under the truncated law taken by [expectation()], the sums
#' over the support or the quadratures that give the expected information
#' itself.
#'
#' @details
#' With \eqn{\ell} the truncated log-likelihood, \eqn{\ell - \log Z}, and
#' \eqn{G = \ell_i \ell_j}, the second Bartlett identity gives
#' \eqn{E[\ell_{ij}] = -E[G]}, and differentiating under the expectation,
#' \deqn{\partial_c E[G] = E[G_c + G\ell_c],}
#' \deqn{\partial_{cd} E[G] = E[G_{cd} + G_c\ell_d + G_d\ell_c +
#'   G(\ell_c\ell_d + \ell_{cd})],}
#' with \eqn{G_c = \ell_{ic}\ell_j + \ell_i\ell_{jc}} and
#' \eqn{G_{cd} = \ell_{icd}\ell_j + \ell_{ic}\ell_{jd} + \ell_{id}\ell_{jc} +
#' \ell_i\ell_{jcd}}. The truncation points do not depend on the parameters,
#' so no boundary term arises. The truncated derivatives are the parent's less
#' those of \eqn{\log Z}, which do not depend on \eqn{y}: they are computed
#' once per parameter combination by `trunc_logz_derivs()` and passed to the
#' integrand, which reads only the parent's derivatives at the nodes.
#' Registered through [register_dexpected()].
#'
#' @param d A truncated family, from [truncated()].
#' @param y The response, read for its length.
#' @param th The parameters.
#' @param k The order, `1L` or `2L`.
#' @param t Unused; the thread count of the registration.
#'
#' @return A named list keyed as [dexpected_names()] (`k = 1`) or
#'   [d2expected_names()] (`k = 2`), each component of length `length(y)`.
#'
#' @seealso [distrib_expected_hessian.TruncatedContinuousDistrib()]
#'
#' @keywords internal
trunc_dexpected <- function(d, y, th, k, t = 1L) {
  n <- length(y)
  parent <- d@parent_distrib
  P <- d@params
  p <- length(P)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  if (k == 1L) {
    parts <- trunc_route_parts(d, y, th, "dinfo")
    if (!is.null(parts)) {
      # d_c of -(S_ab / Z - m_a m_b), with d_c m_a = M_ac - m_a m_c
      Z <- parts$Z
      m <- lapply(parts$Zi, function(v) v / Z)
      M <- function(a, c) parts$Zij[[hess_pair_name(P, a, c)]] / Z
      dm <- function(a, c) M(a, c) - m[[P[a]]] * m[[P[c]]]
      out <- list()
      for (r in seq_len(nrow(pairs))) {
        a <- pairs[r, 1L]; b <- pairs[r, 2L]
        S <- parts$S[[hess_pair_name(P, a, b)]]
        for (c in seq_len(p)) {
          key <- dexpected_key(P, a, b, c)
          out[[key]] <- rep_len(-(parts$D[[key]] / Z - S * parts$Zi[[P[c]]] / (Z * Z) -
                                    (dm(a, c) * m[[P[b]]] + m[[P[a]]] * dm(b, c))), n)
        }
      }
      return(out[dexpected_names(P)])
    }
  }
  hk <- function(a, b) hess_pair_name(P, a, b)
  tk <- function(a, b, c) paste(P[sort(c(a, b, c))], collapse = "_")
  # the derivatives of log Z, passed to the integrand through expectation()'s
  # dots under names that cannot collide with a parameter's
  lz <- trunc_logz_derivs(d, th, k + 1L)
  names(lz) <- paste0(".lz_", names(lz))
  E <- function(f) rep_len(-do.call(expectation, c(list(d, f, th), lz)), n)
  # the truncated derivatives at the nodes, from the parent's
  parts <- function(y, theta, z, order) {
    gp <- distrib_gradient(parent, y, theta)
    g <- lapply(seq_len(p), function(i) gp[[P[i]]] - z[[paste0(".lz_", P[i])]])
    hp <- distrib_hessian(parent, y, theta)
    h <- function(a, b) hp[[hk(a, b)]] - z[[paste0(".lz_", hk(a, b))]]
    out <- list(g = g, h = h)
    if (order >= 2L) {
      tp <- distrib_deriv3(parent, y, theta)
      out$t3 <- function(a, b, c) tp[[tk(a, b, c)]] - z[[paste0(".lz_", tk(a, b, c))]]
    }
    out
  }
  out <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    if (k == 1L) {
      for (c in seq_len(p)) {
        out[[dexpected_key(P, a, b, c)]] <- E(function(y, theta, ...) {
          q <- parts(y, theta, list(...), 1L)
          g <- q$g; h <- q$h
          h(a, c) * g[[b]] + g[[a]] * h(b, c) + g[[a]] * g[[b]] * g[[c]]
        })
      }
      next
    }
    for (r2 in seq_len(nrow(pairs))) {
      c <- pairs[r2, 1L]; dd <- pairs[r2, 2L]
      out[[d2expected_key(P, a, b, c, dd)]] <- E(function(y, theta, ...) {
        q <- parts(y, theta, list(...), 2L)
        g <- q$g; h <- q$h; t3 <- q$t3
        G <- g[[a]] * g[[b]]
        Gc <- h(a, c) * g[[b]] + g[[a]] * h(b, c)
        Gd <- h(a, dd) * g[[b]] + g[[a]] * h(b, dd)
        Gcd <- t3(a, c, dd) * g[[b]] + h(a, c) * h(b, dd) + h(a, dd) * h(b, c) +
          g[[a]] * t3(b, c, dd)
        Gcd + Gc * g[[dd]] + Gd * g[[c]] + G * (g[[c]] * g[[dd]] + h(c, dd))
      })
    }
  }
  out[if (k == 1L) dexpected_names(P) else d2expected_names(P)]
}

#' @rdname trunc_dexpected
#' @details
#' `trunc_logz_derivs()` returns the derivatives of \eqn{\log Z} to order
#' `order`, keyed by [deriv_names()] at each order (the pairs as
#' [hess_names()] spells them), from the ratios \eqn{d^B Z / Z} by
#' [log_deriv()]: a ratio is read from the parent's distribution function
#' where its derivatives of that order exist, and is otherwise the
#' expectation of the complete Bell polynomial under the truncated law.
#' @param theta The parameters.
#' @param order The highest order, 1 to 4.
trunc_logz_derivs <- function(d, theta, order) {
  parent <- d@parent_distrib
  P <- d@params
  dZ <- lapply(seq_len(order), function(k) trunc_mass_derivs(d, theta, k))
  Zval <- trunc_constants(d, theta)$Z
  ratio <- memo_ratio(function(block) {
    k <- length(block)
    if (!is.null(dZ[[k]])) return(dZ[[k]][[canon_key(block, P)]] / Zval)
    expectation(d, function(y, theta) {
      bell_f_ratio(block, parent_ell(parent, y, theta, length(block), P))
    }, theta)
  }, P)
  out <- list()
  for (o in seq_len(order)) {
    for (idx in deriv_indices(P, o)) {
      out[[canon_key(P[idx], P)]] <- log_deriv(P[idx], ratio)
    }
  }
  out
}

register_dexpected(TruncatedContinuousDistrib, trunc_dexpected)
register_dexpected(TruncatedDiscreteDistrib, trunc_dexpected)
