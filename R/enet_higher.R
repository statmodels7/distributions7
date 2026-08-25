#' @include enet_distrib.R reparametrize.R
NULL

#' @title Higher Derivatives of the Elastic Net's Log Mills Ratio
#'
#' @description
#' Adds \eqn{G''} and \eqn{G'''} to the list of pieces, obtained by
#' differentiating the identity \eqn{G' = 1 + xG - G^2} instead of the ratio
#' itself:
#' \deqn{G'' = G + xG' - 2GG', \qquad
#'   G''' = 2G' + xG'' - 2(G')^{2} - 2GG''.}
#'
#' @details
#' Written this way each order is a polynomial in \eqn{x}, \eqn{G} and the
#' orders below, so the cancellation that afflicts \eqn{G = x - 1/M} for large
#' \eqn{x} is confined to \eqn{G} itself, where `.enet_G()` already switches to
#' an asymptotic series past \eqn{|x| = 10^3}. Nothing above first order
#' introduces a new cancellation.
#'
#' @param p The value of `.enet_parts()`: a list carrying `mu`, `lam`, `al`,
#'   `a`, `c`, `x`, `g` and `dg`.
#'
#' @return `p` with two components added: `d2g` for \eqn{G''} and `d3g` for
#'   \eqn{G'''}, each of the length `g` has.
#'
#' @section Notation:
#' \eqn{x = a/\sqrt c}, \eqn{M} the Mills ratio and
#' \eqn{G = \mathrm{d}\log M/\mathrm{d}x}.
#'
#' @seealso [enet_distrib()] for the family and
#'   [distrib_deriv3.EnetDistrib()] for the order these serve.
#'
#' @examples
#' p <- distributions7:::.enet_g_higher(
#'   distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5)))
#' c(x = p$x, G = p$g, dG = p$dg, d2G = p$d2g, d3G = p$d3g)
#'
#' # G'' against a central difference of the identity for G'.
#' gp <- function(x) { g <- distributions7:::.enet_G(x); 1 + x * g - g^2 }
#' eps <- 1e-5
#' c(analytic = p$d2g, numeric = (gp(p$x + eps) - gp(p$x - eps)) / (2 * eps))
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

#' @title The Elastic Net in Its Two Rates
#'
#' @description
#' Returns the log-density's derivatives with respect to the location and the
#' two rates, \eqn{(\mu, a, c)}, up to the requested order. It is the inner
#' half of [.enet_chain()]; the outer half carries the result onto
#' \eqn{(\mu, \lambda, \alpha)}.
#'
#' @details
#' # Why these coordinates
#'
#' In them the log-density is
#' \eqn{-a|z| - cz^{2}/2 - \log 2 - L(x) + \tfrac{1}{2}\log c} with
#' \eqn{z = y - \mu}, \eqn{x = ac^{-1/2}} and \eqn{L = \log M}, so the data
#' term is quadratic in \eqn{z} and linear in each rate. Every derivative of
#' order three or more is therefore a derivative of the **normalizing
#' constant**, apart from \eqn{\partial^{3}\ell/\partial\mu^{2}\partial c = -1}.
#'
#' The normalizer is one pass of [chain_assemble()] over \eqn{L(x(a,c))}: the
#' inner derivatives are \eqn{-G, -G', -G'', -G'''} and the map has
#' \eqn{\partial^{2}x/\partial a^{2} = 0}, so a key with two or more \eqn{a}'s
#' is an exact zero and the table leaves it out. The remaining
#' \eqn{\tfrac12\log c} is pure in \eqn{c} and is added component by component.
#'
#' # The kink
#'
#' The location is not differentiable at \eqn{z = 0}, as in the Laplace, and
#' the sign function is what the first derivative carries there. Nothing at
#' order four touches \eqn{\mu} at all.
#'
#' @param y A numeric vector of observations.
#' @param p The value of `.enet_parts()`, extended by [.enet_g_higher()] when
#'   `order` is 3 or 4.
#' @param order A single integer, 1 to 4: the highest derivative order wanted.
#'
#' @return A list of length `order`; element `k` is a named list of the
#'   order-`k` derivatives, keyed as
#'   [`deriv_names(c("mu", "a", "c"), k)`][deriv_names] names them, each entry
#'   a numeric vector of the length of `y`.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha} is the Laplace rate, \eqn{c = \lambda(1-\alpha)} the
#' Gaussian one, \eqn{x = a/\sqrt c}, \eqn{L = \log M} with \eqn{M} the Mills
#' ratio, and \eqn{z = y - \mu}.
#'
#' @seealso [.enet_chain()], which calls this and carries the result onward;
#'   [chain_assemble()] for the partition sum; and [enet_distrib()] for the
#'   family.
#'
#' @examples
#' p <- distributions7:::.enet_g_higher(
#'   distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5)))
#' y <- c(-1.5, 0.4, 2.1)
#' out <- distributions7:::.enet_ac_derivs(y, p, 3L)
#' vapply(out, length, 0L)          # 3, 6 and 10 components
#'
#' # The one third-order component that carries data.
#' out[[3]][["mu_mu_c"]]
#'
#' # Everything else at order three is free of y, the data term being
#' # quadratic in z and linear in each rate.
#' vapply(out[[3]][setdiff(names(out[[3]]), "mu_mu_c")],
#'        function(v) diff(range(v)), 0)
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

#' @title The Map From the Rates to the Elastic Net's Parameters
#'
#' @description
#' Returns the partial derivatives of \eqn{(\mu, a, c)} as functions of
#' \eqn{(\mu, \lambda, \alpha)}, with \eqn{a = \lambda\alpha} and
#' \eqn{c = \lambda(1-\alpha)}, in the keyed-table form [chain_assemble()]
#' consumes.
#'
#' @details
#' The map is bilinear, so its second partials are the constants \eqn{\pm 1}
#' and its third and higher vanish; the table says so by leaving them out.
#'
#' Those second partials are what carries a low-order derivative in the rates
#' up to a high-order one in the parameters. That is why the third derivatives
#' in \eqn{(\lambda, \alpha)} still depend on the data although the third
#' derivatives in the rates do not.
#'
#' @param p The value of `.enet_parts()`, read for `al` and `lam`.
#'
#' @return A named list of three keyed tables, one per rate coordinate, with
#'   keys `"1"`, `"2"`, `"3"` for the first partials in \eqn{\mu},
#'   \eqn{\lambda}, \eqn{\alpha} and `"2,3"` for the one non-zero second
#'   partial.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, and the numeric keys
#' index \eqn{(\mu, \lambda, \alpha)} in that order.
#'
#' @seealso [.enet_chain()], its only caller, and [chain_assemble()] for the
#'   partition sum that consumes it.
#'
#' @examples
#' p <- distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5))
#' distributions7:::.enet_rate_maps(p)
#'
#' # The two rates add up to lambda whatever alpha is, so their first
#' # partials in lambda sum to one and those in alpha cancel.
#' m <- distributions7:::.enet_rate_maps(p)
#' c(in_lambda = m$a[["2"]] + m$c[["2"]],
#'   in_alpha = m$a[["3"]] + m$c[["3"]])
#'
#' @keywords internal
.enet_rate_maps <- function(p) {
  list(
    mu = list("1" = 1),
    a  = list("2" = p$al, "3" = p$lam, "2,3" = 1),
    c  = list("2" = 1 - p$al, "3" = -p$lam, "2,3" = -1)
  )
}

#' @title The Elastic Net's Derivatives of a Given Order
#'
#' @description
#' Assembles the log-density's derivatives in \eqn{(\mu, a, c)} through
#' [.enet_ac_derivs()] and carries them onto \eqn{(\mu, \lambda, \alpha)}
#' through [chain_assemble()] and the bilinear table of [.enet_rate_maps()].
#' It is what [distrib_deriv3.EnetDistrib()] and
#' [distrib_deriv4.EnetDistrib()] return.
#'
#' @details
#' The licence for orders three and four is that the **same assembly** run at
#' orders one and two reproduces the hand-written score and observed Hessian:
#' measured, exactly 0 at order one and \eqn{1.1\times10^{-16}} at order two.
#' Against an independent route, one product stencil on the analytic
#' log-density, the pure-\eqn{\lambda} components agree to
#' \eqn{2\times10^{-4}} at order three and \eqn{2\times10^{-3}} at order four,
#' which is the stencil's own accuracy at those orders.
#'
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, read
#'   positionally by `.enet_parts()`.
#' @param order A single integer, 1 to 4: the derivative order wanted.
#'
#' @return A named list of the distinct order-`order` components in
#'   \eqn{(\mu, \lambda, \alpha)}, named as [deriv_names()] names them, each a
#'   numeric vector of the length of `y`.
#'
#' @seealso [.enet_ac_derivs()] for the inner half, [.enet_rate_maps()] for
#'   the map, [chain_assemble()] for the partition sum, and
#'   [distrib_deriv3.EnetDistrib()] for the method that calls this.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # At order one it reproduces the hand-written score exactly, which is what
#' # licenses the orders where there is nothing to compare against.
#' g <- distrib_gradient(d, y, th)
#' a1 <- distributions7:::.enet_chain(y, th, 1L)
#' max(abs(unlist(a1[names(g)]) - unlist(g)))
#'
#' # And at order two, the hand-written Hessian.
#' h <- distrib_hessian(d, y, th)
#' a2 <- distributions7:::.enet_chain(y, th, 2L)
#' max(abs(unlist(a2[names(h)]) - unlist(h)))
#'
#' @keywords internal
.enet_chain <- function(y, theta, order) {
  p <- .enet_g_higher(.enet_parts(theta))
  chain_assemble(.enet_ac_derivs(y, p, order), c("mu", "a", "c"),
                 .enet_rate_maps(p), c("mu", "lambda", "alpha"),
                 order, length(y))
}

#' @title Elastic-Net Third Derivatives
#' @name distrib_deriv3.EnetDistrib
#'
#' @description
#' Computes the ten third derivatives of the log-density in closed form,
#' through [.enet_chain()]. Written in \eqn{(\mu, a, c)} the log-density is
#' quadratic in \eqn{z = y-\mu} and linear in each rate, so at this order only
#' the normalizing constant contributes, apart from
#' \eqn{\partial^{3}\ell/\partial\mu^{2}\partial c = -1}. The bilinear map to
#' \eqn{(\lambda, \alpha)} then puts the data back in: `mu_mu_lambda` and
#' `mu_mu_alpha` are not constant across observations even though their
#' preimages are.
#'
#' The licence for this order is that the same assembly at orders one and two
#' reproduces the hand-written score and Hessian, measured at exactly 0 and
#' \eqn{1.1\times10^{-16}}.
#'
#' With `expected = TRUE` the value is an expectation, and there it is **not**
#' closed form: the call routes to [expected_derivative()], so `approx` and
#' `nsim` are read.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length matters.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of ten numeric vectors, from `mu_mu_mu` to
#'   `alpha_alpha_alpha` as [deriv_names()] names them.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, \eqn{z = y - \mu}, and
#' \eqn{G = \mathrm{d}\log M/\mathrm{d}x} with \eqn{M} the Mills ratio.
#'
#' @seealso [distrib_hessian.EnetDistrib()] for the order below,
#'   [distrib_deriv4.EnetDistrib()] for the order above, [.enet_chain()] for
#'   the assembly, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # Against a central difference of the analytic Hessian.
#' eps <- 1e-5
#' rbind(analytic = d3$mu_mu_lambda,
#'       numeric = (distrib_hessian(d, y, list(mu = 0, lambda = 2 + eps,
#'                                             alpha = 0.5))$mu_mu -
#'                  distrib_hessian(d, y, list(mu = 0, lambda = 2 - eps,
#'                                             alpha = 0.5))$mu_mu) / (2 * eps))
#'
#' # The log-density is quadratic in the response, so the third derivative in
#' # the location alone is exactly zero.
#' d3$mu_mu_mu
#'
#' # The pure-lambda component against one stencil on the log-density, a
#' # route that shares no algebra with the assembly.
#' ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, lambda = v,
#'                                              alpha = 0.5), log = TRUE))
#' c(ours = sum(d3$lambda_lambda_lambda),
#'   stencil = numericals7::fd_derivative(ld, 2, 3L, h = 0.02))
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

#' @title Elastic-Net Fourth Derivatives
#' @name distrib_deriv4.EnetDistrib
#'
#' @description
#' Computes the fifteen fourth derivatives of the log-density in closed form,
#' through [.enet_chain()], in the notation of
#' [distrib_deriv3.EnetDistrib()]. At this order **nothing touches the
#' location**: the data term is quadratic in \eqn{z} and linear in each rate,
#' so every component naming \eqn{\mu} more than twice is exactly zero, and
#' what remains is a derivative of the normalizing constant carried through the
#' bilinear map.
#'
#' With `expected = TRUE` the value is an expectation and is not closed form,
#' as at third order.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length matters.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of fifteen numeric vectors, from `mu_mu_mu_mu` to
#'   `alpha_alpha_alpha_alpha`.
#'
#' @seealso [distrib_deriv3.EnetDistrib()] for the order below,
#'   [.enet_chain()] for the assembly, and [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' d4 <- distrib_deriv4(d, y, th)
#' length(d4)
#'
#' # Nothing at this order touches the location.
#' d4$mu_mu_mu_mu
#'
#' # Against a central difference of the third order.
#' eps <- 1e-5
#' rbind(analytic = d4$mu_mu_lambda_lambda,
#'       numeric = (distrib_deriv3(d, y, list(mu = 0, lambda = 2 + eps,
#'                                            alpha = 0.5))$mu_mu_lambda -
#'                  distrib_deriv3(d, y, list(mu = 0, lambda = 2 - eps,
#'                                            alpha = 0.5))$mu_mu_lambda) /
#'                 (2 * eps))
#'
#' # The pure-lambda component against one stencil on the log-density.
#' ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, lambda = v,
#'                                              alpha = 0.5), log = TRUE))
#' c(ours = sum(d4$lambda_lambda_lambda_lambda),
#'   stencil = numericals7::fd_derivative(ld, 2, 4L, h = 0.05))
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
