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
#' Returns \eqn{\partial^j (1+a^2)^{-1}/\partial a^j} for \eqn{j = 0, \ldots, 3}:
#' \eqn{(1+a^2)^{-1}}, \eqn{-2a(1+a^2)^{-2}}, \eqn{(6a^2-2)(1+a^2)^{-3}} and
#' \eqn{24a(1-a^2)(1+a^2)^{-4}}. This is the factor the skew normal's shape
#' derivative carries, and writing it out saves a Faa di Bruno pass over a
#' quantity whose derivatives are four short expressions.
#'
#' @param a A numeric vector, the shape.
#' @param j The order, 0 to 3. Any other value returns `NULL`, `switch()`
#'   falling through; three is the highest the fourth-order cdf derivative
#'   needs.
#'
#' @return A numeric vector the length of `a`.
#'
#' @seealso [sn_cdf_std_derivs()], the one consumer.
#'
#' @examples
#' # The first derivative at a = 2 is -2a / (1 + a^2)^2 = -0.16.
#' distributions7:::recip_1p_sq(2, 1)
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
#' Returns \eqn{\partial^i_z \partial^j_\alpha G} for
#' \eqn{G(z, \alpha) = \Phi(z) - 2T(z, \alpha)} over the pairs with
#' \eqn{1 \le i + j \le} `order`, with \eqn{T} Owen's T. This is the whole of
#' the skew normal's cdf derivative surface in standard coordinates; the
#' location and the scale are chained on afterwards by [sn_cdf_deriv_k()].
#'
#' @details
#' # Why the integral never has to be differentiated
#'
#' Owen's T is defined by an integral, and both of its partial derivatives are
#' elementary:
#' \deqn{\frac{\partial T}{\partial h} = -\varphi(h)\left\{\Phi(ah) -
#'       \tfrac12\right\}, \qquad
#'       \frac{\partial T}{\partial a} = \frac{\varphi(h)\,\varphi(ah)}{1+a^2}.}
#' The integral is therefore differentiated away at the first order and never
#' has to be differentiated again. A quantity defined by an integral is not
#' thereby a quantity whose derivatives need one.
#'
#' # The first derivatives, and the check on them
#'
#' \eqn{G_z = 2\varphi(z)\Phi(u)} and
#' \eqn{G_\alpha = -2\varphi(z)\varphi(u)R} with \eqn{u = \alpha z} and
#' \eqn{R = (1+\alpha^2)^{-1}}. The first is the skew normal's own density,
#' which confirms the identity above; the second is checked against the
#' gradient the family reports.
#'
#' Everything above the first order is one Leibniz rule over those factors,
#' with \eqn{\Phi(u)} and \eqn{\varphi(u)} differentiated by Faa di Bruno over
#' \eqn{u}. That map is bilinear, so it has only three non-zero partials and
#' the expansion stays short.
#'
#' @section Notation:
#' \eqn{z} is the standardized quantile, \eqn{\alpha} the shape,
#' \eqn{u = \alpha z}, \eqn{\Phi} and \eqn{\varphi} the standard normal
#' distribution and density, and \eqn{T} Owen's T.
#'
#' @param z The standardized quantile, a numeric vector.
#' @param al The shape, a numeric vector of any sign recyclable against `z`.
#' @param order The highest total order wanted, 1 to 4.
#'
#' @return A nested list indexed `[[i + 1]][[j + 1]]`, holding
#'   \eqn{\partial^i_z\partial^j_\alpha G} as a numeric vector. Entries with
#'   \eqn{i + j} above `order` are absent.
#'
#' @seealso [sn_cdf_deriv_k()], which chains this onto the location and the
#'   scale; [recip_1p_sq()] for the shape factor;
#'   [numericals7::owen_t()] for the function itself.
#'
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

#' Assemble the Skew Normal's CDF Derivatives of a Given Order
#'
#' @description
#' Chains the standard-coordinate table of [sn_cdf_std_derivs()] through
#' \eqn{z = (q-\mu)/\sigma}, the shape passing straight through as the second
#' index of that table. Together the two functions give the family closed cdf
#' derivatives at all four orders in all three parameters.
#'
#' @details
#' The map is the one the other location-scale families are chained through:
#' \eqn{z} is linear in the location and a reciprocal in the scale, so a block
#' naming the location twice contributes an exact zero and only eight partials
#' of the map are non-zero up to order four. The shape does not enter \eqn{z}
#' at all, which is why it can be carried as an index of the inner table rather
#' than through the chain.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\alpha} the
#' shape, \eqn{z = (q-\mu)/\sigma} and \eqn{F} the distribution function.
#'
#' @param distrib A `SkewNormal1Distrib` object, whose `params` name and order
#'   the components.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive) and
#'   `alpha` (any sign).
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names].
#'
#' @seealso [sn_cdf_std_derivs()] for the inner table;
#'   [distrib_grad_cdf.SkewNormal1Distrib()], the family page;
#'   [loc_scale_cdf_deriv_k()] for the same chain without a shape.
#'
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

#' @title Skew Normal Log-CDF Derivatives
#' @name distrib_grad_cdf.SkewNormal1Distrib
#' @aliases distrib_hess_cdf.SkewNormal1Distrib
#'   distrib_deriv3_cdf.SkewNormal1Distrib
#'   distrib_deriv4_cdf.SkewNormal1Distrib
#'
#' @description
#' Closed form at every order from one to four, in the shape as well as in the
#' location and the scale. With \eqn{z = (q-\mu)/\sigma} the distribution
#' function is \eqn{\Phi(z) - 2T(z,\alpha)}, and Owen's \eqn{T} has elementary
#' partial derivatives in both arguments, so the integral in its definition is
#' differentiated away at the first order and never has to be differentiated
#' again. The location and the scale then enter only through \eqn{z}, by the
#' same chain rule the other location-scale families use.
#'
#' @details
#' # What it is worth
#'
#' Against a product stencil on the same cdf, at \eqn{\mu = 0.3},
#' \eqn{\sigma = 1.2}, \eqn{\alpha = 2}: \eqn{3.6\times10^{-9}} at order 1,
#' \eqn{2.1\times10^{-6}} at order 2, \eqn{9.5\times10^{-5}} at order 3 and
#' \eqn{3.8\times10^{-4}} at order 4.
#'
#' # The family this leaves behind
#'
#' The skew normal used to be grouped with the Student t and the pseudo-Huber,
#' whose shape components are differenced. It is not any more, and the reason
#' is the elementary partials above: only the skew t still differences a shape
#' here, its degrees of freedom entering through a Student t distribution
#' function.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\alpha} the
#' shape, \eqn{z = (q-\mu)/\sigma}, \eqn{\Phi} the standard normal distribution
#' function and \eqn{T} Owen's T.
#'
#' @param distrib A `SkewNormal1Distrib` object, from [skewnormal1_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive) and
#'   `alpha` (any sign), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of the order the generic asked for,
#'   keyed as [`deriv_names(distrib@params, order)`][deriv_names]: three
#'   components for the gradient, six for the Hessian, ten at order 3 and
#'   fifteen at order 4.
#'
#' @seealso [sn_cdf_std_derivs()] and [sn_cdf_deriv_k()] for the construction;
#'   [distrib_grad_cdf.SkewTDistrib()], which does difference its shape;
#'   [skewnormal1_distrib()].
#'
#' @examples
#' d <- skewnormal1_distrib()
#' q <- c(-1, 0.5, 2)
#' th <- list(mu = 0.3, sigma = 1.2, alpha = 2)
#'
#' # The location component is exact, the density itself.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' # At shape zero the family is Gaussian, and so are its cdf derivatives.
#' g0 <- distrib_grad_cdf(d, q, list(mu = 0.3, sigma = 1.2, alpha = 0),
#'                        log = FALSE)
#' gg <- distrib_grad_cdf(gaussian1_distrib(), q, list(mu = 0.3, sigma = 1.2),
#'                        log = FALSE)
#' max(abs(g0$mu - gg$mu))
#'
#' # Three, six, ten and fifteen components as the order rises.
#' lengths(list(distrib_grad_cdf(d, q, th),
#'              distrib_hess_cdf(d, q, th),
#'              distrib_deriv3_cdf(d, q, th),
#'              distrib_deriv4_cdf(d, q, th)))
#'
#' @keywords internal
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
