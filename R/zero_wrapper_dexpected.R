#' @include zero_inflated.R zero_adjusted.R dexpected_hessian.R
NULL

#' The Derivatives of the Zero Wrappers' Expected Information
#'
#' @description
#' Compute the first derivatives of the expected information of
#' [zero_inflated()] and [zero_adjusted()] families in closed form, from the
#' parent's expected information and its derivatives and from the parent's
#' mass, score and Hessian at zero, on the parameter scale.
#'
#' @details
#' With \eqn{f = f(0)}, \eqn{s} and \eqn{H} the parent's score and Hessian at
#' zero, \eqn{E} and \eqn{D} the parent's expected information and its
#' derivatives, \eqn{\zeta} the zero-inflation probability and
#' \eqn{L = \zeta + (1-\zeta) f}, the zero-inflated expected information is
#' \deqn{E^{ZI}_{ij} = (1-\zeta) E_{ij} + \zeta(1-\zeta) f s_i s_j / L,}
#' \deqn{E^{ZI}_{i\zeta} = -f s_i / L, \qquad
#'   E^{ZI}_{\zeta\zeta} = -(1-f)^2/L - (1-f)/(1-\zeta),}
#' the observed Hessian at zero cancelling from the first, and the
#' derivatives follow with \eqn{\partial_c f = f s_c} and
#' \eqn{\partial_c s_i = H_{ic}}. For the zero-adjusted discrete family, with
#' \eqn{q = 1 - f} and \eqn{\pi} the probability of zero,
#' \deqn{E^{ZA}_{ij} = (1-\pi)\,[E_{ij}/q + f s_i s_j / q^2], \qquad
#'   E^{ZA}_{\pi\pi} = -1/(\pi(1-\pi)),}
#' the mixed block being zero; for the continuous one
#' \eqn{E^{ZA}_{ij} = (1-\pi) E_{ij}}.
#'
#' @param distrib A family from [zero_inflated()] or [zero_adjusted()].
#' @param y The response, read for its length.
#' @param theta The parameters, the parent's followed by the mixing
#'   probability.
#' @param threads The thread count passed to the parent's methods.
#'
#' @return A named list keyed as [dexpected_names()], each component of
#'   length `length(y)`, on the parameter scale.
#'
#' @seealso [distrib_dexpected_hessian()]
#'
#' @keywords internal
zero_wrapper_dexpected <- function(distrib, y, theta, threads = 1L) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  z <- pars$mix
  n <- length(y)
  P <- parent@params
  p <- length(P)
  allp <- distrib@params
  continuous <- S7::S7_inherits(distrib, ZeroAdjustedContinuousDistrib)
  inflated <- S7::S7_inherits(distrib, ZeroInflatedDistrib)
  E <- distrib_expected_hessian(parent, y, pars$orig, threads = threads)
  D <- distrib_dexpected_hessian(parent, y, pars$orig, threads = threads)
  Ep <- function(a, b) rep_len(E[[hess_pair_name(P, a, b)]], n)
  Dp <- function(a, b, c) rep_len(D[[dexpected_key(P, a, b, c)]], n)
  if (!continuous) {
    f0 <- exp(distrib_pdf(parent, 0, pars$orig, log = TRUE))
    sg <- distrib_gradient(parent, 0, pars$orig)
    hs <- distrib_hessian(parent, 0, pars$orig)
    s <- function(a) rep_len(sg[[P[a]]], n)
    H <- function(a, b) rep_len(hs[[hess_pair_name(P, a, b)]], n)
  }
  z <- rep_len(z, n)
  out <- list()
  pairs <- which(upper.tri(diag(p + 1L), diag = TRUE), arr.ind = TRUE)
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1L]; j <- pairs[r, 2L]
    for (c in seq_len(p + 1L)) {
      key <- dexpected_key(allp, i, j, c)
      out[[key]] <- if (inflated) {
        zi_dexpected_entry(i, j, c, p, z, f0, s, H, Ep, Dp)
      } else if (continuous) {
        za_cont_dexpected_entry(i, j, c, p, z, Ep, Dp, n)
      } else {
        za_disc_dexpected_entry(i, j, c, p, z, f0, s, H, Ep, Dp, n)
      }
    }
  }
  out[dexpected_names(allp)]
}

#' @rdname zero_wrapper_dexpected
#' @details
#' Each entry is computed by `zi_dexpected_entry()`,
#' `za_disc_dexpected_entry()` or `za_cont_dexpected_entry()` with the
#' expressions of the scalar registry (`pt_wrappers.h`), so that the two
#' routes agree to the last bit.
#' @param i,j,c The indices of one entry: the pair and the parameter
#'   differentiated, the parent's parameters first and the probability last.
#' @param p The number of the parent's parameters.
#' @param zi,za The wrapper's probability, one value per observation.
#' @param f0 The parent's mass at zero.
#' @param s,H Functions of parameter indices returning the parent's score and
#'   Hessian at zero.
#' @param Ep,Dp Functions of parameter indices returning the parent's expected
#'   information and its derivatives.
#' @param n The number of observations.
zi_dexpected_entry <- function(i, j, c, p, zi, f0, s, H, Ep, Dp) {
  l0 <- zi + (1 - zi) * f0
  m <- p + 1L
  if (i <= p && j <= p && c <= p) {
    return((1 - zi) * Dp(i, j, c) +
             zi * (1 - zi) * f0 * (H(i, c) * s(j) + s(i) * H(j, c)) / l0 +
             zi^2 * (1 - zi) * f0 * s(i) * s(j) * s(c) / l0^2)
  }
  if (i <= p && j <= p) {
    return(-Ep(i, j) + f0 * s(i) * s(j) *
             ((1 - 2 * zi) / l0 - zi * (1 - zi) * (1 - f0) / l0^2))
  }
  if (i <= p && j == m) {
    if (c <= p) {
      return(-(f0 / l0) * (s(c) * s(i) + H(i, c)) +
               (1 - zi) * f0^2 * s(i) * s(c) / l0^2)
    }
    return(f0 * s(i) * (1 - f0) / l0^2)
  }
  if (c <= p) {
    return(f0 * s(c) * (2 * (1 - f0) / l0 + (1 - zi) * (1 - f0)^2 / l0^2 +
                          1 / (1 - zi)))
  }
  (1 - f0)^3 / l0^2 - (1 - f0) / (1 - zi)^2
}

#' @rdname zero_wrapper_dexpected
za_disc_dexpected_entry <- function(i, j, c, p, za, f0, s, H, Ep, Dp, n) {
  q <- 1 - f0
  m <- p + 1L
  if (i <= p && j <= p && c <= p) {
    return((1 - za) * (Dp(i, j, c) / q + Ep(i, j) * f0 * s(c) / q^2 +
                         f0 * (s(c) * s(i) * s(j) + H(i, c) * s(j) +
                                 s(i) * H(j, c)) / q^2 +
                         2 * f0^2 * s(i) * s(j) * s(c) / q^3))
  }
  if (i <= p && j <= p) return(-(Ep(i, j) / q + f0 * s(i) * s(j) / q^2))
  if (i <= p || c <= p) return(rep(0, n))
  (1 - 2 * za) / (za * (1 - za))^2
}

#' @rdname zero_wrapper_dexpected
za_cont_dexpected_entry <- function(i, j, c, p, za, Ep, Dp, n) {
  m <- p + 1L
  if (i <= p && j <= p && c <= p) return((1 - za) * Dp(i, j, c))
  if (i <= p && j <= p) return(-Ep(i, j))
  if (i <= p || c <= p) return(rep(0, n))
  (1 - 2 * za) / (za * (1 - za))^2
}

for (.zw_cls in list(ZeroInflatedDistrib, ZeroAdjustedDiscreteDistrib,
                     ZeroAdjustedContinuousDistrib)) {
  S7::method(distrib_dexpected_hessian, .zw_cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
      threads = 1L) {
    dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                       function(k) zero_wrapper_dexpected(distrib, y, theta,
                                                          threads))
  }
}
rm(.zw_cls)
