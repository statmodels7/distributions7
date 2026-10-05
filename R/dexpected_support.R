#' @include dexpected_hessian.R dexpected_families.R negbin1_distrib.R
#' @include betabinom1_distrib.R betabinom2_distrib.R
NULL

#' @title Derivatives of the Expected Information, Sums Over the Support
#' @name distrib_dexpected_hessian.NegBin1Distrib
#' @aliases distrib_d2expected_hessian.NegBin1Distrib
#'   distrib_dexpected_hessian.BetaBinom1Distrib
#'   distrib_d2expected_hessian.BetaBinom1Distrib
#'   distrib_dexpected_hessian.BetaBinom2Distrib
#'   distrib_d2expected_hessian.BetaBinom2Distrib
#' @description
#' The first and second derivatives of the expected information in the
#' parameters for three families whose expected information is itself a sum
#' over the support.
#'
#' @details
#' **NB1.** The expected information is \eqn{-A}, \eqn{W/\theta^2} and
#' \eqn{-W/(\mu\theta)}, with \eqn{A} and \eqn{W} the two sums over the support
#' of [distrib_expected_hessian.NegBin1Distrib()]. Their derivatives are sums
#' over the same mass with the score of the mass beside the summand,
#' \eqn{\partial_c F = \mathbb{E}[f_c + f s_c]} and
#' \eqn{\partial_{ce} F = \mathbb{E}[f_{ce} + f_c s_e + f_e s_c + f(s_c s_e +
#' h_{ce})]}, taken in \eqn{(\mu, \theta)} by a compiled kernel, with the
#' scores \eqn{s} and second derivatives \eqn{h} of the log-mass written so
#' that no term of order \eqn{\mu/\theta^2} is formed. Measured against exact
#' sums at 50 digits, the error relative to the largest component of each
#' order is about 2e-10 at \eqn{\mu = 100}, \eqn{\theta = 0.05}, and grows
#' toward the Poisson limit, to 4e-08 at \eqn{\theta = 0.005} and 5e-06 at
#' \eqn{\theta = 0.001} (\eqn{\mu = 10}), the fourth derivative in \eqn{\theta}
#' being the worst component.
#'
#' **Beta-binomial.** The support \eqn{\{0, \dots, n\}} is finite, so the
#' identities
#' \deqn{\partial_c \mathbb{E}[\ell_{ab}] = \mathbb{E}[\ell_{abc} + \ell_{ab}\ell_c]}
#' and its second-order counterpart,
#' \deqn{\partial_{cd} \mathbb{E}[\ell_{ab}] = \mathbb{E}[\ell_{abcd} + \ell_{abd}\ell_c
#'   + \ell_{abc}\ell_d + \ell_{ab}\ell_{cd} + \ell_{ab}\ell_c\ell_d],}
#' are exact finite sums against the mass. They are taken in the shapes by a
#' compiled kernel, where every derivative of the log-mass is a polygamma
#' difference at an integer shift and is therefore a sum of negative powers,
#' accumulated beside the mass with no polygamma called;
#' [betabinom1_distrib()] reads the same sums through its map to the shapes
#' ([betabinom1_dexpected()]). [support_dexpected()] computes the same
#' identities from the family's own observed derivatives, sharing no
#' arithmetic with the kernel, and the two agree to about \eqn{10^{-12}}
#' relative; the kernel is 5 to 15 times faster at \eqn{n = 4000}.
#'
#' @param distrib A distribution object of one of the classes above.
#' @param y A numeric vector of observations, read for its length.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads A single positive integer, passed to the kernels.
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#' @seealso [distrib_dexpected_hessian()], [distrib_d2expected_hessian()]
#' @keywords internal
NULL

register_dexpected(NegBin1Distrib, function(d, y, th, k, t)
  (if (k == 1L) negbin1_dexpected1_cpp(y, th[[1]], th[[2]], t)
    else negbin1_dexpected2_cpp(y, th[[1]], th[[2]], t)))


#' Derivatives of the Expected Information as Exact Sums Over a Finite Support
#'
#' @description
#' Sums the moment identities of [distrib_dexpected_hessian()] and
#' [distrib_d2expected_hessian()] over the support \eqn{\{0, \dots, N\}}
#' against the mass, reading the family's own observed derivatives.
#'
#' @details
#' The support being finite, each sum is the expectation exactly: no
#' quadrature, no truncation. Every observation's parameters are crossed with
#' the \eqn{N + 1} support points and the derivative generics are called once
#' on the whole grid.
#'
#' @param distrib A discrete distribution with a finite support.
#' @param y The response, read for its length.
#' @param theta A named list of parameters, each of length one or `length(y)`.
#' @param order `1L` or `2L`.
#' @param N The largest point of the support.
#' @param threads Passed to the derivative generics.
#'
#' @return A named list on the parameter scale, keyed as [dexpected_names()]
#'   at order 1 and as [dexpected_names()] followed by [d2expected_names()]
#'   at order 2, so that [dexpected_analytic()] reads both orders from one
#'   call.
#'
#' @seealso [distrib_dexpected_hessian.BetaBinom1Distrib]
#'
#' @keywords internal
support_dexpected <- function(distrib, y, theta, order, N, threads = 1L) {
  P <- distrib@params
  np <- length(P)
  n <- length(y)
  m <- N + 1L
  yy <- rep(0:N, each = n)
  tl <- lapply(theta, function(v) rep(rep_len(v, n), times = m))
  pm <- distrib_pdf(distrib, yy, tl)
  g <- distrib_gradient(distrib, yy, tl, threads = threads)
  h <- distrib_hessian(distrib, yy, tl, threads = threads)
  d3 <- distrib_deriv3(distrib, yy, tl, threads = threads)
  d4 <- if (order == 2L) distrib_deriv4(distrib, yy, tl, threads = threads)
  srt <- function(ix) paste(P[sort(ix)], collapse = "_")
  E <- function(v) rowSums(matrix(pm * v, n, m))
  pairs <- which(upper.tri(diag(np), diag = TRUE), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L] != pairs[, 2L]), , drop = FALSE]
  out <- list()
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1L]; j <- pairs[r, 2L]
    hab <- h[[hess_pair_name(P, i, j)]]
    # the first order is returned at order 2 as well, being cheap once the
    # family's derivatives over the support are in hand
    for (k in seq_len(np)) {
      out[[dexpected_key(P, i, j, k)]] <- E(d3[[srt(c(i, j, k))]] + hab * g[[P[k]]])
    }
    if (order == 1L) next
    for (s in seq_len(nrow(pairs))) {
      k <- pairs[s, 1L]; l <- pairs[s, 2L]
      out[[d2expected_key(P, i, j, k, l)]] <- E(
        d4[[srt(c(i, j, k, l))]] +
          d3[[srt(c(i, j, l))]] * g[[P[k]]] + d3[[srt(c(i, j, k))]] * g[[P[l]]] +
          hab * h[[hess_pair_name(P, k, l)]] + hab * g[[P[k]]] * g[[P[l]]])
    }
  }
  out
}

register_dexpected(BetaBinom2Distrib, function(d, y, th, k, t) {
  if (k == 1L) betabinom_shapes_dexpected1_cpp(y, th[[1]], th[[2]], d@size, t)
  else betabinom_shapes_dexpected2_cpp(y, th[[1]], th[[2]], d@size, t)
})
register_dexpected(BetaBinom1Distrib, function(d, y, th, k, t)
  betabinom1_dexpected(y, th[[1]], th[[2]], d@size, k, t))

#' The Beta-Binomial's Expected-Information Derivatives in Mean and Dispersion
#'
#' @description
#' [betabinom1_distrib()] is the shapes \eqn{a = \mu/\sigma} and
#' \eqn{b = (1-\mu)/\sigma} under another name, so its derivatives are the
#' shapes' ones, summed exactly over the support by the compiled kernel, and
#' carried across with the map's partials written out (at order 1 by the
#' compiled kernel `betabinom1_dexpected1_cpp`, at order 2 by
#' [dexpected_chain()]):
#' \eqn{a_\mu = 1/\sigma}, \eqn{a_\sigma = -\mu/\sigma^2},
#' \eqn{a_{\mu\sigma} = -1/\sigma^2}, \eqn{a_{\sigma\sigma} = 2\mu/\sigma^3},
#' \eqn{a_{\mu\sigma\sigma} = 2/\sigma^3},
#' \eqn{a_{\sigma\sigma\sigma} = -6\mu/\sigma^4}, and for \eqn{b} the same
#' with \eqn{\mu} replaced by \eqn{1-\mu} and every derivative in \eqn{\mu}
#' changing sign; every derivative taking \eqn{\mu} twice is zero.
#'
#' @param y The response, read for its length.
#' @param mu,sigma The mean proportion and the dispersion.
#' @param size The number of trials.
#' @param order `1L` or `2L`.
#' @param threads Passed to the kernel.
#'
#' @return A named list on the parameter scale, keyed as [dexpected_names()]
#'   at order 1 and as [d2expected_names()] at order 2.
#'
#' @seealso [support_dexpected()], [dexpected_chain()]
#'
#' @keywords internal
betabinom1_dexpected <- function(y, mu, sigma, size, order, threads = 1L) {
  if (order == 1L) {
    return(betabinom1_dexpected1_cpp(y, mu, sigma, size, threads))
  }
  n <- length(y)
  s <- sigma; m1 <- 1 - mu
  a <- mu / s; b <- m1 / s
  kp <- betabinom_shapes_dexpected1_cpp(y, a, b, size, threads)
  Pp <- c("alpha", "beta"); Pn <- c("mu", "sigma")
  # the derivative of order k of the information reads the map to order k + 1
  maps <- md_betabinom1(list(mu, sigma), order + 1L)
  E <- kp[hess_names(Pp)]
  d1 <- kp[dexpected_names(Pp)]
  if (order == 1L) return(dexpected_chain(Pp, Pn, E, d1, NULL, maps, 1L, n))
  d2 <- betabinom_shapes_dexpected2_cpp(y, a, b, size, threads)
  dexpected_chain(Pp, Pn, E, d1, d2, maps, 2L, n)
}
