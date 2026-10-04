#' @include zero_inflated.R zero_adjusted.R dexpected_hessian.R
#' @include dexpected_families.R
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

#' The Second Derivatives of the Zero Wrappers' Expected Information
#'
#' @description
#' Compute the second derivatives of the expected information of
#' [zero_inflated()] and [zero_adjusted()] families in closed form, on the
#' parameter scale, by differentiating the expressions of
#' [zero_wrapper_dexpected()] once more.
#'
#' @details
#' The parent enters through its expected information and its first and
#' second derivatives, and, for a discrete parent, through its mass
#' \eqn{f}, score \eqn{s}, Hessian \eqn{H} and third derivatives \eqn{T} at
#' zero, with \eqn{\partial_c f = f s_c}, \eqn{\partial_c s_i = H_{ic}} and
#' \eqn{\partial_c H_{ij} = T_{ijc}}. Each entry of the zero-inflated
#' information is a product or a quotient of these quantities and of the
#' probability, and its second derivative is formed by Leibniz's rule; a
#' quotient \eqn{N/L} uses
#' \deqn{\partial_{cd}(N/L) = \frac{N_{cd}}{L} - \frac{N_c L_d + N_d L_c}{L^2}
#'   - \frac{N L_{cd}}{L^2} + \frac{2 N L_c L_d}{L^3}.}
#' The zero-adjusted discrete entries \eqn{(1-\pi)[E_{ij}/q + f s_i s_j/q^2]}
#' are differentiated through \eqn{1/q} and \eqn{1/q^2}, and the continuous
#' ones are \eqn{(1-\pi)} times the parent's.
#'
#' @inheritParams zero_wrapper_dexpected
#'
#' @return A named list keyed as [d2expected_names()], each component of
#'   length `length(y)`, on the parameter scale.
#'
#' @seealso [distrib_d2expected_hessian()]
#'
#' @keywords internal
zero_wrapper_d2expected <- function(distrib, y, theta, threads = 1L) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  n <- length(y)
  P <- parent@params
  p <- length(P)
  m <- p + 1L
  allp <- distrib@params
  continuous <- S7::S7_inherits(distrib, ZeroAdjustedContinuousDistrib)
  inflated <- S7::S7_inherits(distrib, ZeroInflatedDistrib)
  E <- distrib_expected_hessian(parent, y, pars$orig, threads = threads)
  D <- distrib_dexpected_hessian(parent, y, pars$orig, threads = threads)
  D2 <- distrib_d2expected_hessian(parent, y, pars$orig, threads = threads)
  # the parent's expected information and its derivatives; a derivative in
  # the probability (index m) is zero
  e0 <- function(a, b) rep_len(E[[hess_pair_name(P, a, b)]], n)
  e1 <- function(a, b, c) {
    if (c == m) return(0)
    rep_len(D[[dexpected_key(P, a, b, c)]], n)
  }
  e2 <- function(a, b, c, d) {
    if (c == m || d == m) return(0)
    rep_len(D2[[d2expected_key(P, a, b, c, d)]], n)
  }
  z <- rep_len(pars$mix, n)
  tm <- list(e0 = e0, e1 = e1, e2 = e2, z = z, m = m, n = n)
  if (!continuous) {
    f0 <- exp(distrib_pdf(parent, 0, pars$orig, log = TRUE))
    sg <- distrib_gradient(parent, 0, pars$orig)
    hs <- distrib_hessian(parent, 0, pars$orig)
    d3 <- distrib_deriv3(parent, 0, pars$orig)
    tm$f <- rep_len(f0, n)
    tm$s <- function(a) rep_len(sg[[P[a]]], n)
    tm$H <- function(a, b) rep_len(hs[[hess_pair_name(P, a, b)]], n)
    tm$T <- function(a, b, c) rep_len(d3[[paste(P[sort(c(a, b, c))], collapse = "_")]], n)
  }
  entry <- if (inflated) {
    zi_d2expected_entry
  } else if (continuous) {
    za_cont_d2expected_entry
  } else {
    za_disc_d2expected_entry
  }
  out <- list()
  pairs <- which(upper.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1L]; j <- pairs[r, 2L]
    for (r2 in seq_len(nrow(pairs))) {
      c <- pairs[r2, 1L]; d <- pairs[r2, 2L]
      out[[d2expected_key(allp, i, j, c, d)]] <- rep_len(entry(i, j, c, d, tm), n)
    }
  }
  out[d2expected_names(allp)]
}

#' @rdname zero_wrapper_d2expected
#' @details
#' Each entry is computed by `zi_d2expected_entry()`,
#' `za_disc_d2expected_entry()` or `za_cont_d2expected_entry()`, and
#' `d2_quotient()` gives the second derivative of a quotient.
#' @param i,j,c,d The indices of one entry: the pair \eqn{(i,j)} and the
#'   pair \eqn{(c,d)} differentiated, the parent's parameters first and the
#'   probability last.
#' @param tm A list of the parent's quantities: functions `e0`, `e1`, `e2`
#'   of parameter indices returning the expected information and its first and
#'   second derivatives, the probability `z`, the index `m` of the
#'   probability, the length `n`, and, for a discrete parent, the mass `f` and
#'   functions `s`, `H`, `T` returning the score, Hessian and third
#'   derivatives at zero.
#' @param N,Nc,Nd,Ncd,L,Lc,Ld,Lcd A numerator and a denominator with their
#'   first derivatives in \eqn{c} and \eqn{d} and their second derivative.
d2_quotient <- function(N, Nc, Nd, Ncd, L, Lc, Ld, Lcd) {
  Ncd / L - (Nc * Ld + Nd * Lc) / L^2 - N * Lcd / L^2 + 2 * N * Lc * Ld / L^3
}

#' @rdname zero_wrapper_d2expected
#' @details
#' `zw_disc_parts()` returns the quantities shared by the discrete entries,
#' as functions of parameter indices: \eqn{f s_i s_j} and \eqn{f s_i} with
#' their first and second derivatives, which vanish in the probability, and
#' \eqn{q = 1 - f} with its derivatives.
zw_disc_parts <- function(tm) {
  f <- tm$f; s <- tm$s; H <- tm$H; T <- tm$T; m <- tm$m
  list(
    K0 = function(i, j) f * s(i) * s(j),
    K1 = function(i, j, c) {
      if (c == m) return(0)
      f * (s(c) * s(i) * s(j) + H(i, c) * s(j) + s(i) * H(j, c))
    },
    K2 = function(i, j, c, d) {
      if (c == m || d == m) return(0)
      f * (s(d) * (s(c) * s(i) * s(j) + H(i, c) * s(j) + s(i) * H(j, c)) +
             H(c, d) * s(i) * s(j) + s(c) * H(i, d) * s(j) +
             s(c) * s(i) * H(j, d) + T(i, c, d) * s(j) + H(i, c) * H(j, d) +
             H(i, d) * H(j, c) + s(i) * T(j, c, d))
    },
    M0 = function(i) f * s(i),
    M1 = function(i, c) if (c == m) 0 else f * (s(c) * s(i) + H(i, c)),
    M2 = function(i, c, d) {
      if (c == m || d == m) return(0)
      f * (s(d) * (s(c) * s(i) + H(i, c)) + H(c, d) * s(i) + s(c) * H(i, d) +
             T(i, c, d))
    },
    q0 = 1 - f,
    q1 = function(c) if (c == m) 0 else -f * s(c),
    q2 = function(c, d) {
      if (c == m || d == m) return(0)
      -f * (s(c) * s(d) + H(c, d))
    }
  )
}

#' @rdname zero_wrapper_d2expected
zi_d2expected_entry <- function(i, j, c, d, tm) {
  z <- tm$z; m <- tm$m
  k <- zw_disc_parts(tm)
  f <- tm$f; s <- tm$s; H <- tm$H
  # L = zeta + (1 - zeta) f
  L <- z + (1 - z) * f
  L1 <- function(a) if (a == m) 1 - f else (1 - z) * f * s(a)
  L2 <- if (c == m && d == m) {
    0
  } else if (c == m) {
    -f * s(d)
  } else if (d == m) {
    -f * s(c)
  } else {
    (1 - z) * f * (s(c) * s(d) + H(c, d))
  }
  if (i < m && j < m) {
    # (1 - zeta) E_ij + w G, G = K / L, w = zeta (1 - zeta)
    u1 <- function(a) if (a == m) -1 else 0
    w <- z * (1 - z)
    w1 <- function(a) if (a == m) 1 - 2 * z else 0
    w2 <- if (c == m && d == m) -2 else 0
    K0 <- k$K0(i, j)
    G1 <- function(a) k$K1(i, j, a) / L - K0 * L1(a) / L^2
    G2 <- d2_quotient(K0, k$K1(i, j, c), k$K1(i, j, d), k$K2(i, j, c, d),
                      L, L1(c), L1(d), L2)
    return(u1(c) * tm$e1(i, j, d) + u1(d) * tm$e1(i, j, c) +
             (1 - z) * tm$e2(i, j, c, d) + w2 * K0 / L + w1(c) * G1(d) +
             w1(d) * G1(c) + w * G2)
  }
  if (i < m) {
    # -f s_i / L
    return(-d2_quotient(k$M0(i), k$M1(i, c), k$M1(i, d), k$M2(i, c, d),
                        L, L1(c), L1(d), L2))
  }
  # -(1 - f)^2 / L - (1 - f) / (1 - zeta)
  q <- k$q0
  N2 <- 2 * (k$q1(c) * k$q1(d) + q * k$q2(c, d))
  v <- 1 / (1 - z)
  v1 <- function(a) if (a == m) 1 / (1 - z)^2 else 0
  v2 <- if (c == m && d == m) 2 / (1 - z)^3 else 0
  -d2_quotient(q^2, 2 * q * k$q1(c), 2 * q * k$q1(d), N2, L, L1(c), L1(d), L2) -
    (k$q2(c, d) * v + k$q1(c) * v1(d) + k$q1(d) * v1(c) + q * v2)
}

#' @rdname zero_wrapper_d2expected
za_disc_d2expected_entry <- function(i, j, c, d, tm) {
  z <- tm$z; m <- tm$m
  if (i == m && j == m) {
    return(if (c == m && d == m) -2 / z^3 - 2 / (1 - z)^3 else 0)
  }
  if (j == m || (c == m && d == m)) return(0)
  k <- zw_disc_parts(tm)
  # A = E_ij / q + K / q^2, which does not depend on the probability
  q <- k$q0
  r1 <- function(a) -k$q1(a) / q^2
  t1 <- function(a) -2 * k$q1(a) / q^3
  A1 <- function(a) {
    tm$e1(i, j, a) / q + tm$e0(i, j) * r1(a) + k$K1(i, j, a) / q^2 +
      k$K0(i, j) * t1(a)
  }
  if (c == m) return(-A1(d))
  if (d == m) return(-A1(c))
  r2 <- -k$q2(c, d) / q^2 + 2 * k$q1(c) * k$q1(d) / q^3
  t2 <- -2 * k$q2(c, d) / q^3 + 6 * k$q1(c) * k$q1(d) / q^4
  A2 <- tm$e2(i, j, c, d) / q + tm$e1(i, j, c) * r1(d) +
    tm$e1(i, j, d) * r1(c) + tm$e0(i, j) * r2 + k$K2(i, j, c, d) / q^2 +
    k$K1(i, j, c) * t1(d) + k$K1(i, j, d) * t1(c) + k$K0(i, j) * t2
  (1 - z) * A2
}

#' @rdname zero_wrapper_d2expected
za_cont_d2expected_entry <- function(i, j, c, d, tm) {
  z <- tm$z; m <- tm$m
  if (i == m && j == m) {
    return(if (c == m && d == m) -2 / z^3 - 2 / (1 - z)^3 else 0)
  }
  if (j == m || (c == m && d == m)) return(0)
  if (c == m) return(-tm$e1(i, j, d))
  if (d == m) return(-tm$e1(i, j, c))
  (1 - z) * tm$e2(i, j, c, d)
}

for (.zw_cls in list(ZeroInflatedDistrib, ZeroAdjustedDiscreteDistrib,
                     ZeroAdjustedContinuousDistrib)) {
  register_dexpected(.zw_cls, function(d, y, th, k, t) {
    if (k == 1L) zero_wrapper_dexpected(d, y, th, t)
    else zero_wrapper_d2expected(d, y, th, t)
  })
}
rm(.zw_cls)
