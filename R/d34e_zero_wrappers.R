#' @include zero_inflated.R zero_adjusted.R zero_wrapper_dexpected.R
#' @include d34expected_analytic.R
NULL

#' Derivative Algebra on Index Multisets
#'
#' @description
#' Small operations on quantities represented as functions of a sorted integer
#' vector `S`, a multiset of parameter indices, returning the derivative of the
#' quantity along `S` (one value per observation). `zw34_mul()` forms the
#' product by Leibniz's rule over the sub-multisets of `S`, `zw34_recip()` the
#' reciprocal from the identity \eqn{A\,(1/A) = 1} differentiated along `S`,
#' and the remaining functions are the sum, the difference and the constants.
#' They serve [zero_wrapper_d34expected()].
#'
#' @param f A function of `S` to be memoized.
#' @param A,B Quantities, as functions of a sorted integer vector.
#' @param x A numeric vector, the value of a constant quantity.
#' @param k A non-negative integer, the length of `S`.
#'
#' @return `zw34_memo()`, `zw34_mul()`, `zw34_recip()`, `zw34_add()`,
#'   `zw34_neg()` and `zw34_const()` return a function of `S`;
#'   `zw34_subsets()` returns a list of logical vectors of length `k`, one per
#'   sub-multiset.
#'
#' @keywords internal
zw34_memo <- function(f) {
  cache <- new.env(parent = emptyenv())
  function(S) {
    key <- paste0("s", paste(S, collapse = ","))
    v <- cache[[key]]
    if (is.null(v)) {
      v <- f(S)
      cache[[key]] <- v
    }
    v
  }
}

#' @rdname zw34_memo
zw34_subsets <- function(k) {
  if (k == 0L) return(list(logical(0)))
  lapply(0:(2^k - 1), function(b) bitwAnd(b, 2^(seq_len(k) - 1)) > 0)
}

#' @rdname zw34_memo
zw34_mul <- function(A, B) {
  zw34_memo(function(S) {
    acc <- 0
    for (U in zw34_subsets(length(S))) acc <- acc + A(S[U]) * B(S[!U])
    acc
  })
}

#' @rdname zw34_memo
zw34_recip <- function(A) {
  G <- zw34_memo(function(S) {
    if (length(S) == 0L) return(1 / A(S))
    acc <- 0
    for (U in zw34_subsets(length(S))[-1L]) acc <- acc + A(S[U]) * G(S[!U])
    -acc * G(integer(0))
  })
  G
}

#' @rdname zw34_memo
zw34_add <- function(A, B) function(S) A(S) + B(S)

#' @rdname zw34_memo
zw34_neg <- function(A) function(S) -A(S)

#' @rdname zw34_memo
zw34_const <- function(x) function(S) if (length(S) == 0L) x else 0

#' The Third and Fourth Derivatives of the Zero Wrappers' Expected Information
#'
#' @description
#' Compute the third or fourth derivative of the expected Hessian of
#' [zero_inflated()] and [zero_adjusted()] families, on the parameter scale,
#' by differentiating the expressions of [zero_wrapper_dexpected()] along each
#' multiset of parameters. The parent enters through its expected Hessian and
#' its derivatives up to the same order, and, for a discrete parent, through
#' the log-mass at zero and its derivatives to order five.
#'
#' @details
#' The entries are those of [zero_wrapper_dexpected()]. With \eqn{f} the
#' parent's mass at zero, \eqn{s_i} its score there, \eqn{E_{ij}} the parent's
#' expected Hessian, \eqn{\zeta} the probability and
#' \eqn{L = \zeta + (1-\zeta) f}, the zero-inflated entries are
#' \eqn{(1-\zeta)E_{ij} + \zeta(1-\zeta) f s_i s_j / L},
#' \eqn{-f s_i / L} and \eqn{-(1-f)^2/L - (1-f)/(1-\zeta)}, and the
#' zero-adjusted discrete ones are \eqn{(1-\pi)[E_{ij}/q + f s_i s_j/q^2]}
#' with \eqn{q = 1-f}. The derivative of \eqn{f} along a multiset of
#' parameters is formed from the log-mass derivatives by
#' \eqn{\partial_{aR} f = \sum_{U \subseteq R} \partial_U f\,
#' \partial_{a R\setminus U}\log f}, products by Leibniz's rule and
#' reciprocals from \eqn{A\,(1/A) = 1}, through the operations of
#' [zw34_memo()]. The orders of the parent's expected Hessian derivatives are
#' read from its own methods, closed form where the parent has them and the
#' default stencil otherwise. The order-3 call reads the parent to order 3 and
#' the order-4 call to order 4.
#'
#' @param distrib A family from [zero_inflated()] or [zero_adjusted()].
#' @param y The response, read for its length.
#' @param theta The parameters, the parent's followed by the mixing
#'   probability.
#' @param order `3L` or `4L`.
#' @param threads The thread count passed to the parent's methods.
#'
#' @return A named list keyed as [d3expected_names()] or [d4expected_names()],
#'   each component of length `length(y)`, on the parameter scale.
#'
#' @seealso [zero_wrapper_dexpected()], [distrib_d3expected_hessian()]
#'
#' @keywords internal
zero_wrapper_d34expected <- function(distrib, y, theta, order, threads = 1L) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  n <- length(y)
  P <- parent@params
  p <- length(P)
  m <- p + 1L
  allp <- distrib@params
  continuous <- S7::S7_inherits(distrib, ZeroAdjustedContinuousDistrib)
  inflated <- S7::S7_inherits(distrib, ZeroInflatedDistrib)
  has_m <- function(S) any(S == m)

  # the parent's expected information and its derivatives, by order
  PD <- list(
    distrib_expected_hessian(parent, y, pars$orig, threads = threads),
    distrib_dexpected_hessian(parent, y, pars$orig, threads = threads),
    distrib_d2expected_hessian(parent, y, pars$orig, threads = threads),
    distrib_d3expected_hessian(parent, y, pars$orig, threads = threads))
  if (order == 4L) {
    PD[[5L]] <- distrib_d4expected_hessian(parent, y, pars$orig,
                                           threads = threads)
  }
  pkey <- function(a, b, S) {
    switch(length(S) + 1L,
      hess_pair_name(P, a, b),
      dexpected_key(P, a, b, S[1L]),
      d2expected_key(P, a, b, S[1L], S[2L]),
      d3expected_key(P, a, b, S[1L], S[2L], S[3L]),
      d4expected_key(P, a, b, S[1L], S[2L], S[3L], S[4L]))
  }
  Epar <- function(a, b) {
    function(S) {
      if (length(S) > 0L && has_m(S)) return(0)
      rep_len(PD[[length(S) + 1L]][[pkey(a, b, S)]], n)
    }
  }

  zv <- function(S) {
    if (length(S) == 0L) rep_len(pars$mix, n)
    else if (length(S) == 1L && S == m) 1
    else 0
  }
  one_z <- function(S) {
    if (length(S) == 0L) 1 - rep_len(pars$mix, n)
    else if (length(S) == 1L && S == m) -1
    else 0
  }

  if (!continuous) {
    # the parent's log-mass at zero and its derivatives to order five
    lf0 <- distrib_pdf(parent, 0, pars$orig, log = TRUE)
    LD <- list(distrib_gradient(parent, 0, pars$orig),
               distrib_hessian(parent, 0, pars$orig),
               distrib_deriv3(parent, 0, pars$orig),
               distrib_deriv4(parent, 0, pars$orig))
    if (order == 4L) LD[[5L]] <- distrib_deriv5(parent, 0, pars$orig)
    lder <- function(S) {
      k <- length(S)
      if (k == 0L) return(rep_len(lf0, n))
      if (has_m(S)) return(0)
      S <- sort(S)
      key <- if (k == 1L) P[S] else if (k == 2L) hess_pair_name(P, S[1L], S[2L])
             else paste(P[S], collapse = "_")
      rep_len(LD[[k]][[key]], n)
    }
    f <- zw34_memo(function(S) {
      if (length(S) == 0L) return(rep_len(exp(lf0), n))
      if (has_m(S)) return(0)
      a <- S[1L]
      R <- S[-1L]
      acc <- 0
      for (U in zw34_subsets(length(R))) {
        acc <- acc + f(R[U]) * lder(sort(c(a, R[!U])))
      }
      acc
    })
    q <- function(S) {
      if (length(S) == 0L) rep_len(-expm1(lf0), n) else -f(S)
    }
    sc <- function(i) function(S) lder(sort(c(i, S)))
    Lz <- zw34_add(zv, zw34_mul(one_z, f))
  }

  if (!continuous) {
    rL <- zw34_recip(Lz)
    rq <- zw34_recip(q)
    rq2 <- zw34_mul(rq, rq)
    rz1 <- zw34_recip(one_z)
    q2 <- zw34_mul(q, q)
  }
  cache <- new.env(parent = emptyenv())
  # the quantity of the pair (i, j) as a function of the multiset
  pair_fun <- function(i, j) {
    key <- paste(i, j)
    if (!is.null(cache[[key]])) return(cache[[key]])
    h <- if (continuous) {
      if (i == m && j == m) {
        function(S) {
          if (any(S != m)) return(0)
          k <- length(S)
          z <- rep_len(pars$mix, n)
          -(((-1)^k * factorial(k) * z^(-k - 1)) +
              factorial(k) * (1 - z)^(-k - 1))
        }
      } else if (i == m || j == m) {
        function(S) 0
      } else {
        zw34_mul(one_z, Epar(i, j))
      }
    } else if (inflated) {
      if (i < m && j < m) {
        K <- zw34_mul(zw34_mul(f, sc(i)), sc(j))
        w <- zw34_mul(zv, one_z)
        zw34_add(zw34_mul(one_z, Epar(i, j)),
                 zw34_mul(w, zw34_mul(K, rL)))
      } else if (i < m) {
        zw34_neg(zw34_mul(zw34_mul(f, sc(i)), rL))
      } else {
        zw34_add(zw34_neg(zw34_mul(q2, rL)), zw34_neg(zw34_mul(q, rz1)))
      }
    } else {
      if (i == m && j == m) {
        function(S) {
          if (any(S != m)) return(0)
          k <- length(S)
          z <- rep_len(pars$mix, n)
          -(((-1)^k * factorial(k) * z^(-k - 1)) +
              factorial(k) * (1 - z)^(-k - 1))
        }
      } else if (i == m || j == m) {
        function(S) 0
      } else {
        K <- zw34_mul(zw34_mul(f, sc(i)), sc(j))
        zw34_mul(one_z, zw34_add(zw34_mul(Epar(i, j), rq),
                                 zw34_mul(K, rq2)))
      }
    }
    cache[[key]] <- h
    h
  }

  tuples <- deriv_indices(allp, order)
  nm <- if (order == 3L) d3expected_names(allp) else d4expected_names(allp)
  pairs <- which(upper.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  out <- stats::setNames(vector("list", length(nm)), nm)
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1L]
    j <- pairs[r, 2L]
    h <- pair_fun(i, j)
    for (S in tuples) {
      key <- if (order == 3L) {
        d3expected_key(allp, i, j, S[1L], S[2L], S[3L])
      } else {
        d4expected_key(allp, i, j, S[1L], S[2L], S[3L], S[4L])
      }
      out[[key]] <- rep_len(h(S), n)
    }
  }
  out
}

for (.zw_cls in list(ZeroInflatedDistrib, ZeroAdjustedDiscreteDistrib,
                     ZeroAdjustedContinuousDistrib)) {
  # The third and fourth derivatives of a zero wrapper's expected information,
  # from the parent's derivatives to the same order.
  register_d34expected(.zw_cls, function(d, y, th, k, t)
    zero_wrapper_d34expected(d, y, th, k, t))
}
rm(.zw_cls)
