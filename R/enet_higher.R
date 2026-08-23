#' @include enet_distrib.R reparametrize.R
NULL

#' Higher Derivatives of the Elastic Net's Log Mills Ratio
#'
#' @description
#' \eqn{G''} and \eqn{G'''}, obtained by differentiating the identity
#' \eqn{G' = 1 + xG - G^2} rather than the ratio itself.
#'
#' @details
#' \deqn{G'' = G + xG' - 2GG', \qquad
#'   G''' = 2G' + xG'' - 2(G')^{2} - 2GG''.}
#' Written this way each order is a polynomial in \eqn{x}, \eqn{G} and the
#' orders below, so the cancellation that afflicts \eqn{G = x - 1/M} for
#' large \eqn{x} is confined to \eqn{G} itself, where `.enet_G`
#' already switches to an asymptotic series.
#'
#' @param p The value of `.enet_parts`.
#'
#' @return `p` with `d2g` and `d3g` added.
#'
#' @keywords internal
.enet_g_higher <- function(p) {
  g <- p$g
  dg <- p$dg
  d2g <- g + p$x * dg - 2 * g * dg
  p$d2g <- d2g
  p$d3g <- 2 * dg + p$x * d2g - 2 * dg^2 - 2 * g * d2g
  p
}

#' The Elastic Net in Its Two Rates
#'
#' @description
#' The log-density's derivatives with respect to \eqn{(\mu, a, c)} -- the
#' location and the two rates -- up to the requested order.
#'
#' @details
#' In these coordinates the log-density is
#' \eqn{-a|z| - cz^{2}/2 - \log 2 - L(x) + \tfrac{1}{2}\log c} with
#' \eqn{z = y - \mu}, \eqn{x = ac^{-1/2}} and \eqn{L = \log M}, so the data
#' term is quadratic in \eqn{z} and linear in each rate. Every derivative
#' of order three or more is therefore a derivative of the normalizing
#' constant, apart from
#' \eqn{\partial^{3}\ell/\partial\mu^{2}\partial c = -1}, and the
#' normalizer is one pass of [chain_assemble()] over
#' \eqn{L(x(a,c))}: the inner derivatives are \eqn{-G, -G', -G'', -G'''}
#' and the map has \eqn{\partial^{2}x/\partial a^{2} = 0}, which leaves a
#' short table.
#'
#' The location is not differentiable at \eqn{z = 0}, as in the Laplace,
#' and the sign function is what the first derivative carries there.
#'
#' @param y A numeric vector of observations.
#' @param p The value of `.enet_parts`, extended by
#'   `.enet_g_higher`.
#' @param order The derivative order, 1 to 4.
#'
#' @return A list of length `order`; element `k` is the table of
#'   order-`k` derivatives, keyed as
#'   [`deriv_names(c("mu", "a", "c"), k)`][deriv_names].
#'
#' @keywords internal
.enet_ac_derivs <- function(y, p, order) {
  cc <- p$c
  z <- y - p$mu
  s <- sign(z)
  one <- rep(1, length(z))
  zero <- rep(0, length(z))

  # d^j/dc^j of c^{-1/2}, j = 0..4
  u <- c(cc^-0.5, -0.5 * cc^-1.5, 0.75 * cc^-2.5,
         -1.875 * cc^-3.5, 6.5625 * cc^-4.5)
  # x = a c^{-1/2}: a key with two or more a's is an exact zero, which the
  # table says by leaving it out
  xmap <- list("1" = u[1], "2" = p$a * u[2], "1,2" = u[2],
               "2,2" = p$a * u[3], "1,2,2" = u[3],
               "2,2,2" = p$a * u[4], "1,2,2,2" = u[4],
               "2,2,2,2" = p$a * u[5])
  LD <- list(list(x = -p$g), list(x_x = -p$dg),
             list(x_x_x = -p$d2g), list(x_x_x_x = -p$d3g))
  N <- lapply(seq_len(order), function(k) {
    chain_assemble(LD[seq_len(k)], "x", list(x = xmap), c("a", "c"), k, 1L)
  })
  # and the + (1/2) log c of the normalizer, which is pure in c
  lg <- c(0.5 / cc, -0.5 / cc^2, 1 / cc^3, -3 / cc^4)
  for (k in seq_len(order)) {
    key <- paste(rep("c", k), collapse = "_")
    N[[k]][[key]] <- N[[k]][[key]] + lg[k]
  }

  out <- lapply(seq_len(order), function(k) {
    nm <- deriv_names(c("mu", "a", "c"), k)
    stats::setNames(lapply(nm, function(i) zero), nm)
  })
  for (k in seq_len(order)) {
    for (nm in names(N[[k]])) {
      out[[k]][[nm]] <- out[[k]][[nm]] + N[[k]][[nm]] * one
    }
  }

  out[[1]][["mu"]] <- p$a * s + cc * z
  out[[1]][["a"]] <- out[[1]][["a"]] - abs(z)
  out[[1]][["c"]] <- out[[1]][["c"]] - z^2 / 2
  if (order >= 2) {
    out[[2]][["mu_mu"]] <- -cc * one
    out[[2]][["mu_a"]] <- s
    out[[2]][["mu_c"]] <- z
  }
  if (order >= 3) out[[3]][["mu_mu_c"]] <- -one
  # nothing at order four touches mu: the data term is quadratic in z and
  # linear in each rate
  out
}

#' The Map From the Rates to the Elastic Net's Parameters
#'
#' @description
#' The partial tables of \eqn{(\mu, a, c)} as functions of
#' \eqn{(\mu, \lambda, \alpha)}, with \eqn{a = \lambda\alpha} and
#' \eqn{c = \lambda(1-\alpha)}.
#'
#' @details
#' The map is bilinear, so its second partials are the constants
#' \eqn{\pm 1} and its third and higher vanish. Those second partials are
#' what carries a low-order derivative in the rates up to a high-order one
#' in the parameters, which is why the third derivatives in
#' \eqn{(\lambda, \alpha)} still depend on the data although the third
#' derivatives in the rates do not.
#'
#' @param p The value of `.enet_parts`.
#'
#' @return A list of keyed tables, one per rate coordinate.
#'
#' @keywords internal
.enet_rate_maps <- function(p) {
  list(
    mu = list("1" = 1),
    a  = list("2" = p$al, "3" = p$lam, "2,3" = 1),
    c  = list("2" = 1 - p$al, "3" = -p$lam, "2,3" = -1)
  )
}

#' The Elastic Net's Derivatives of a Given Order
#'
#' @description
#' Assembles the log-density's derivatives in \eqn{(\mu, a, c)} and carries
#' them onto \eqn{(\mu, \lambda, \alpha)}.
#'
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `lambda` and `alpha`.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of components.
#'
#' @keywords internal
.enet_chain <- function(y, theta, order) {
  p <- .enet_g_higher(.enet_parts(theta))
  chain_assemble(.enet_ac_derivs(y, p, order), c("mu", "a", "c"),
                 .enet_rate_maps(p), c("mu", "lambda", "alpha"),
                 order, length(y))
}

#' @title Elastic-Net Third and Fourth Derivatives
#' @name distrib_deriv3.EnetDistrib
#'
#' @description
#' The third and fourth derivatives of the log-density, assembled in the
#' two rates and carried onto \eqn{(\lambda, \alpha)} by the bilinear map.
#'
#' @details
#' Written in \eqn{(\mu, a, c)} the log-density is quadratic in \eqn{z} and
#' linear in each rate, so at these orders only the normalizing constant
#' contributes, apart from
#' \eqn{\partial^{3}\ell/\partial\mu^{2}\partial c = -1}; see
#' `.enet_ac_derivs`. Running the same assembly at orders one and two
#' reproduces the hand-written score and Hessian, which is what licenses it
#' at the orders where there is nothing to compare against.
#'
#' The expected derivatives have no closed form here and go through
#' [expected_derivative()].
#'
#' @param distrib An `EnetDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `lambda` and `alpha`.
#' @param expected Logical; whether to return the expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the
#'   generic.
#' @param approx How the expectation is approximated.
#' @param nsim Monte Carlo sample size, used when `approx = "mc"`.
#' @param ... Unused.
#'
#' @return A named list of third-derivative components.
#'
#' @seealso [enet_distrib()]
S7::method(distrib_deriv3, EnetDistrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  .enet_chain(y, theta, 3L)
}

#' @rdname distrib_deriv3.EnetDistrib
#' @name distrib_deriv4.EnetDistrib
#' @return A named list of fourth-derivative components.
S7::method(distrib_deriv4, EnetDistrib) <- function(distrib, y, theta,
                                                    expected = FALSE,
                                                    scale = c("parameter", "link"),
                                                    approx = c("integrate", "bartlett", "mc", "opg"),
                                                    nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  .enet_chain(y, theta, 4L)
}
