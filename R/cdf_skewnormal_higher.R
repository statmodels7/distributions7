#' @include cdf_normal_higher.R
NULL

# The skew normal's distribution function, closed in the shape as well.
#
# F = Phi(z) - 2T(z, alpha) with Owen's T, whose two partial derivatives are
# elementary:
#
#   dT/dh     = -phi(h) (Phi(a h) - 1/2)
#   dT/dalpha =  phi(h) phi(a h) / (1 + a^2)
#
# so that dF/dz = 2 phi(z) Phi(alpha z), which is the density and is the check
# that the first identity is the right one, and
# dF/dalpha = -2 phi(z) phi(alpha z)/(1 + alpha^2). Everything above those is a
# product of normal densities, Hermite polynomials and a rational function of
# the shape, so all four orders close. The integral in Owen's T never has to be
# differentiated: it is differentiated away at the first order.

#' Derivatives of the Reciprocal of One Plus a Square
#'
#' @description
#' Returns \eqn{d^{j}/da^{j}\,(1+a^{2})^{-1}} for \eqn{j = 0, \ldots, 3}, the
#' factor the skew normal's shape derivative carries.
#'
#' @param a A numeric vector.
#' @param j The order, 0 to 3.
#'
#' @return A numeric vector.
#'
#' @keywords internal
recip_1p_sq <- function(a, j) {
  v <- 1 + a^2
  switch(j + 1L,
    1 / v,
    -2 * a / v^2,
    (6 * a^2 - 2) / v^3,
    24 * a * (1 - a^2) / v^4
  )
}

#' The Skew Normal's Distribution Function in Standard Coordinates
#'
#' @description
#' Returns \eqn{\partial^{i}_{z}\partial^{j}_{\alpha}G} for
#' \eqn{G(z, \alpha) = \Phi(z) - 2T(z, \alpha)}, over the pairs with
#' \eqn{1 \le i + j \le} \code{order}, as a list indexed by \code{i + 1} then
#' \code{j + 1}.
#'
#' @details
#' The first derivatives are \eqn{G_{z} = 2\varphi(z)\Phi(u)} and
#' \eqn{G_{\alpha} = -2\varphi(z)\varphi(u)R} with \eqn{u = \alpha z} and
#' \eqn{R = (1+\alpha^{2})^{-1}}; everything above them is one Leibniz rule
#' over those factors, with \eqn{\Phi(u)} and \eqn{\varphi(u)} differentiated
#' by Faa di Bruno over \eqn{u}, which is bilinear and so has only three
#' non-zero partial derivatives.
#'
#' @param z The standardized quantile.
#' @param al The shape.
#' @param order The highest total order, 1 to 4.
#'
#' @return A nested list, \code{[[i + 1]][[j + 1]]}.
#'
#' @seealso \code{\link{skewnormal1_distrib}}
#' @keywords internal
sn_cdf_std_derivs <- function(z, al, order) {
  n <- length(z)
  u <- al * z
  p <- stats::dnorm(z)
  du_ <- stats::dnorm(u)

  # p^(k) = phi^(k)(z) = Phi^(k+1)(z), so the Hermite index is one above the
  # derivative order; at k = 0 that polynomial is 1 and the case needs no
  # branch of its own
  pk <- function(k) phi_hermite(z, k + 1L) * p
  # Phi^(m)(u) and phi^(m)(u); phi^(m) = Phi^(m+1)
  Phim <- function(m) if (m == 0L) stats::pnorm(u) else phi_hermite(u, m) * du_
  phim <- function(m) phi_hermite(u, m + 1L) * du_

  # u = alpha z is bilinear, so only three of its partials are non-zero
  du <- function(nz, na) {
    if (nz == 1L && na == 0L) return(rep_len(al, n))
    if (nz == 0L && na == 1L) return(rep_len(z, n))
    if (nz == 1L && na == 1L) return(rep(1, n))
    rep(0, n)
  }
  # d^{i,j} Phi(u)
  dP <- function(i, j) {
    S <- c(rep("z", i), rep("a", j))
    if (!length(S)) return(Phim(0L))
    acc <- rep(0, n)
    for (part in index_partitions(S)) {
      term <- Phim(length(part))
      for (b in part) term <- term * du(sum(b == "z"), sum(b == "a"))
      acc <- acc + term
    }
    acc
  }

  out <- lapply(seq_len(order + 1L), function(i) vector("list", order + 1L))
  for (i in 0:order) {
    for (j in 0:(order - i)) {
      if (i + j == 0L) next
      out[[i + 1L]][[j + 1L]] <- if (i >= 1L) {
        # d^{i-1,j} of 2 p Phi(u), p depending on z alone
        acc <- rep(0, n)
        for (k in 0:(i - 1L)) {
          acc <- acc + choose(i - 1L, k) * pk(k) * dP(i - 1L - k, j)
        }
        2 * acc
      } else {
        # d^{j-1}_alpha of -2 p phi(u) R, and d^m_alpha phi(u) = phi^(m)(u) z^m
        acc <- rep(0, n)
        for (m in 0:(j - 1L)) {
          acc <- acc + choose(j - 1L, m) * phim(m) * z^m *
            recip_1p_sq(al, j - 1L - m)
        }
        -2 * p * acc
      }
    }
  }
  out
}

#' @title Skew Normal Log-CDF Derivatives
#' @name distrib_grad_cdf.SkewNormal1Distrib
#' @description
#' Closed form at every order, in the shape as well as in the location and the
#' scale. With \eqn{z = (q-\mu)/\sigma} the distribution function is
#' \eqn{\Phi(z) - 2T(z, \alpha)}, and Owen's \eqn{T} has elementary partial
#' derivatives, so the integral in its definition is differentiated away at the
#' first order and never has to be differentiated again. The location and the
#' scale then enter only through \eqn{z}, by the same chain rule the other
#' location-scale families use.
#' @param distrib A \code{SkewNormal1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{sigma} and \code{alpha}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso \code{\link{skewnormal1_distrib}}
#' @keywords internal
sn_cdf_deriv_k <- function(distrib, q, theta, order) {
  params <- distrib@params
  s <- theta[[2]]
  al <- theta[[3]]
  z <- (q - theta[[1]]) / s
  n <- length(q)
  G <- sn_cdf_std_derivs(z, al, order)

  # z is linear in the location and a reciprocal in the scale; a block naming
  # the location twice is therefore zero, and the table is the one the other
  # location-scale families are chained through
  zmap <- list(
    "1" = -1 / s, "2" = -z / s,
    "1,2" = 1 / s^2, "2,2" = 2 * z / s^2,
    "1,2,2" = -2 / s^3, "2,2,2" = -6 * z / s^3,
    "1,2,2,2" = 6 / s^4, "2,2,2,2" = 24 * z / s^4
  )
  dz <- function(block) {
    v <- zmap[[paste(sort(match(block, params[1:2])), collapse = ",")]]
    if (is.null(v)) rep(0, n) else rep_len(v, n)
  }

  idx <- deriv_indices(params, order)
  out <- lapply(idx, function(I) {
    nms <- params[I]
    j <- sum(nms == params[3L])
    rest <- nms[nms != params[3L]]
    if (!length(rest)) return(rep_len(G[[1L]][[j + 1L]], n))
    acc <- rep(0, n)
    for (part in index_partitions(rest)) {
      term <- rep_len(G[[length(part) + 1L]][[j + 1L]], n)
      for (b in part) term <- term * dz(b)
      acc <- acc + term
    }
    acc
  })
  stats::setNames(out, deriv_names(params, order))
}

local({
  make <- function(o) {
    force(o)
    function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
      tabs <- lapply(seq_len(o),
                     function(k) sn_cdf_deriv_k(distrib, q, theta, k))
      cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), tabs, o,
                  lower.tail, log)
    }
  }
  S7::method(distrib_grad_cdf, SkewNormal1Distrib) <- make(1L)
  S7::method(distrib_hess_cdf, SkewNormal1Distrib) <- make(2L)
  S7::method(distrib_deriv3_cdf, SkewNormal1Distrib) <- make(3L)
  S7::method(distrib_deriv4_cdf, SkewNormal1Distrib) <- make(4L)
})

# the second parametrization is the first at its direct parameters, and its
# gate now passes at every order
register_mapped_cdf_k(SkewNormal2Distrib, skewnormal1_distrib,
                      sn2_theta, md_skewnormal2, orders = 1:4)
