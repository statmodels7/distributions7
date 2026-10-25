#' @include expected_loc_scale.R d34expected_analytic.R skewnormal1_distrib.R skewt_distrib.R pseudohuber_distrib.R pseudohuber2_distrib.R skewnormal2_distrib.R
NULL

# The third and fourth derivatives of the expected information of the
# location-scale families skewnormal1, skewt, pseudohuber, pseudohuber2 and
# skewnormal2, by the quadrature of loc_scale_expected().
#
# The expectation is taken as E_ab = -E[l_a l_b], with l the log-density and
# l_a its derivative in the parameter a. This form is equal to E[l_ab] for
# every theta, and it needs the derivatives of l only up to order four for the
# third derivative of E_ab and up to order five for the fourth. The form
# E[l_ab] would need order five and six.
#
# For G = -l_a l_b and a multiset S of parameter indices, differentiation under
# the integral gives d_S E[G] = E[(1/f) d_S(G f)], and the Leibniz rule over
# the positions of S turns this into
#   d_S E[G] = E[ sum_{T subset of S} d_T G * B_(S \ T) ],
# where B_U = (1/f) d_U f is the sum over the set partitions of U of the
# products of the derivatives of l over the blocks, and
#   d_T G = - sum_{T1 subset of T} l_(a u T1) * l_(b u (T \ T1)).
# The subsets are over positions, so repeated indices count as they should.
# Only the components of the requested order are formed.
#
# E_ab does not depend on the location and scales as sigma^-k with k the
# number of indices among a and b on the location or the scale, so only the
# derivatives in the shape parameters are integrated, at sigma = 1.
#
# The nodes and weights are those of loc_scale_expected(), with the weights of
# the skew t increased by loc_scale_weights() for the tail past the last node.

#' The Third or Fourth Derivative of a Location-Scale Expected Information
#'
#' @description
#' Computes the components of [distrib_d3expected_hessian()] or
#' [distrib_d4expected_hessian()] on the parameter scale for a family whose
#' first two parameters are a location and a scale, by the quadrature of
#' [loc_scale_expected()] applied to \eqn{E_{ab} = -E[\ell_a \ell_b]}.
#'
#' @details
#' Differentiating under the integral, \eqn{\partial_S E[G] =
#' E[\sum_{T \subseteq S} \partial_T G \, B_{S \setminus T}]} with
#' \eqn{B_U = f^{-1} \partial_U f}, a sum over the set partitions of \eqn{U}.
#' The form \eqn{-E[\ell_a \ell_b]} reads the derivatives of the log-density
#' up to order five at the fourth derivative, where \eqn{E[\ell_{ab}]} would
#' read order six. A component with an index on the location is zero, and
#' \eqn{j} indices on the scale contribute the factor
#' \eqn{(-1)^j (k)_j \sigma^{-(k+j)}}, with \eqn{(k)_j} the rising factorial
#' and \eqn{k} the number of indices among \eqn{a, b} on the location or the
#' scale. The quadrature runs once per distinct value of the shapes.
#'
#' @param distrib A distribution object of a location-scale family.
#' @param theta A named list of parameters.
#' @param order `3L` or `4L`.
#' @param n The number of observations.
#' @param threads The thread count passed to the family's derivatives.
#'
#' @return A named list keyed as [d3expected_names()] or [d4expected_names()].
#'
#' @keywords internal
loc_scale_d34expected <- function(distrib, theta, order, n, threads = 1L) {
  if (!order %in% c(3L, 4L)) {
    stop("loc_scale_d34expected() takes order 3 or 4.", call. = FALSE)
  }
  params <- distrib@params
  p <- length(params)
  rule <- loc_scale_nodes(distrib)
  N <- length(rule$x)

  # the distinct shapes, and which one each observation carries
  shp <- params[-(1:2)]
  S <- vapply(shp, function(nm) rep_len(theta[[nm]], n), numeric(n))
  if (!is.matrix(S)) S <- matrix(S, nrow = n)
  key <- if (length(shp)) do.call(paste, c(as.data.frame(S), sep = "\r")) else rep("", n)
  first <- !duplicated(key)
  U <- S[first, , drop = FALSE]
  idx <- match(key, key[first])
  m <- nrow(U)

  # every node against every distinct shape, at location 0 and scale 1
  yy <- rep(rule$x, m) + rep(loc_scale_center(distrib, U), each = N)
  th <- c(list(rep(0, N * m), rep(1, N * m)),
          lapply(seq_along(shp), function(j) rep(U[, j], each = N)))
  names(th) <- params
  fw <- loc_scale_pdf(distrib, yy, th) * loc_scale_weights(distrib, U, rule)
  live <- fw != 0

  # the derivatives of the log-density at the nodes, parameter scale; the
  # fifth order is needed only for the fourth derivative of E_ab
  g <- distrib_gradient(distrib, yy, th, threads = threads)
  H <- distrib_hessian(distrib, yy, th, threads = threads)
  T3 <- distrib_deriv3(distrib, yy, th, threads = threads)
  T4 <- distrib_deriv4(distrib, yy, th, threads = threads)
  T5 <- if (order == 4L) distrib_deriv5(distrib, yy, th, threads = threads)

  nm_k <- function(v) paste(params[sort(v)], collapse = "_")
  # l_v for a multiset v of parameter indices, of length one to five
  Lk <- function(v) {
    v <- sort(v)
    switch(length(v),
           g[[v]],
           H[[hess_pair_name(params, v[1L], v[2L])]],
           T3[[nm_k(v)]],
           T4[[nm_k(v)]],
           T5[[nm_k(v)]])
  }

  memo <- new.env(hash = TRUE, parent = emptyenv())
  memo_get <- function(key, make) {
    val <- memo[[key]]
    if (is.null(val)) {
      val <- make()
      assign(key, val, envir = memo)
    }
    val
  }
  # the subsets of k positions, as logical masks
  masks <- function(k) {
    if (k == 0L) return(list(logical(0)))
    lapply(seq_len(2^k) - 1L, function(i) {
      as.logical(bitwAnd(as.integer(i), as.integer(2^(seq_len(k) - 1L))))
    })
  }
  # B_U, the sum over the set partitions of U of the products of l over the
  # blocks; the block of the first position is taken with each subset of the rest
  B <- function(u) {
    if (!length(u)) return(1)
    u <- sort(u)
    memo_get(paste0("B:", paste(u, collapse = ",")), function() {
      first <- u[1L]
      rest <- u[-1L]
      acc <- 0
      for (mk in masks(length(rest))) {
        acc <- acc + Lk(c(first, rest[mk])) * B(rest[!mk])
      }
      acc
    })
  }
  # d_T G_ab for G = -l_a l_b, by the Leibniz rule over the positions of T
  G <- function(a, b, T) {
    T <- sort(T)
    memo_get(paste0("G:", a, ",", b, ":", paste(T, collapse = ",")), function() {
      acc <- 0
      for (mk in masks(length(T))) {
        acc <- acc + Lk(c(a, T[mk])) * Lk(c(b, T[!mk]))
      }
      -acc
    })
  }
  # the integrand of d_S E_ab, before the quadrature
  Q <- function(a, b, Sv) {
    acc <- 0
    for (mk in masks(length(Sv))) {
      acc <- acc + G(a, b, Sv[mk]) * B(Sv[!mk])
    }
    acc
  }
  integrate_rows <- function(v) {
    v[!live] <- 0
    out <- colSums(matrix(v * fw, N, m))
    out[!is.finite(out)] <- NA_real_
    out
  }

  sig <- rep_len(theta[[params[2L]]], n)
  # E_ab does not depend on the location, and E_ab(mu, sigma, s) = sigma^-k
  # c_ab(s) with k the number of indices among a and b on the location or the
  # scale. So a derivative with an index on the location is zero, and one with
  # j indices on the scale is (-1)^j (k)_j sigma^-(k + j) times the derivative
  # of c_ab in the shape indices alone, which the quadrature takes at sigma = 1.
  pr <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  tuples <- deriv_indices(params, order)
  nmo <- if (order == 3L) d3expected_names(params) else d4expected_names(params)
  out <- stats::setNames(vector("list", length(nmo)), nmo)
  for (r in seq_len(nrow(pr))) {
    a <- pr[r, 1L]
    b <- pr[r, 2L]
    k_ab <- sum(c(a, b) <= 2L)
    for (Sv in tuples) {
      key <- if (order == 3L) {
        d3expected_key(params, a, b, Sv[1L], Sv[2L], Sv[3L])
      } else {
        d4expected_key(params, a, b, Sv[1L], Sv[2L], Sv[3L], Sv[4L])
      }
      if (any(Sv == 1L)) {
        out[[key]] <- numeric(n)
        next
      }
      j <- sum(Sv == 2L)
      shape_idx <- Sv[Sv > 2L]
      rising <- prod(k_ab + seq_len(j) - 1)
      val <- (-1)^j * rising * integrate_rows(Q(a, b, shape_idx))[idx] /
        sig^(k_ab + j)
      out[[key]] <- val
    }
  }
  out[nmo]
}

register_d34expected(SkewNormal1Distrib, function(d, y, th, k, t)
  loc_scale_d34expected(d, th, k, length(y), t))

register_d34expected(SkewTDistrib, function(d, y, th, k, t)
  loc_scale_d34expected(d, th, k, length(y), t))

register_d34expected(PseudoHuberDistrib, function(d, y, th, k, t)
  loc_scale_d34expected(d, th, k, length(y), t))

register_d34expected(PseudoHuber2Distrib, function(d, y, th, k, t)
  loc_scale_d34expected(d, th, k, length(y), t))

#' The Third or Fourth Derivative of the Centered Skew Normal's Expected
#' Information
#'
#' @description
#' Computes the components by [loc_scale_d34expected()] where
#' \eqn{|\gamma_1|} is at least the threshold of the series that orders zero
#' to two use near zero, and by one stencil on the second derivative of the
#' expected information, as the default method does, at the observations
#' below it.
#'
#' @details
#' Near \eqn{\gamma_1 = 0} the derivatives of the log-density in
#' \eqn{\gamma_1} diverge in the representation that the quadrature uses,
#' while the second derivative of the expected information comes from a series
#' that is exact there.
#'
#' @param d A `skewnormal2` distribution object.
#' @param y A numeric vector of observations, read for its length.
#' @param th A named list of parameters.
#' @param k `3L` or `4L`.
#' @param t The thread count.
#'
#' @return A named list keyed as [d3expected_names()] or [d4expected_names()],
#'   on the parameter scale.
#'
#' @keywords internal
sn2_d34expected <- function(d, y, th, k, t) {
  th <- align_theta(d, th)
  n <- length(y)
  small <- rep_len(abs(th[[3L]]) < sn2_ge(), n)
  if (!any(small)) return(loc_scale_d34expected(d, th, k, n, t))
  out <- numerical_dexpected_higher(d, y, th, "parameter", k)
  if (all(small)) return(out)
  big <- which(!small)
  sub <- lapply(th, function(v) if (length(v) == 1L) v else v[big])
  qd <- loc_scale_d34expected(d, sub, k, length(big), t)
  for (nm in names(out)) {
    v <- rep_len(out[[nm]], n)
    v[big] <- qd[[nm]]
    out[[nm]] <- v
  }
  out
}

register_d34expected(SkewNormal2Distrib, sn2_d34expected)
