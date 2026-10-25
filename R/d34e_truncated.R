#' @include truncated.R truncated_dexpected.R wrapper_derivatives.R d34expected_analytic.R
NULL

# The third and fourth derivatives of the expected information of a truncated
# family, on the parameter scale, as expectations under the truncated law. The
# construction is the one of trunc_dexpected() carried to orders three and
# four: E_ab = -E[l_a l_b], differentiated under the expectation.

#' Set Partitions of a Vector
#'
#' @description
#' Lists the partitions of the elements of a vector, each a list of blocks.
#'
#' @param v A vector, treated by position.
#'
#' @return A list of partitions; each partition is a list of vectors.
#'
#' @keywords internal
d34e_set_partitions <- function(v) {
  if (length(v) == 0L) return(list(list()))
  first <- v[1L]
  out <- list()
  for (p in d34e_set_partitions(v[-1L])) {
    out[[length(out) + 1L]] <- c(list(first), p)
    for (j in seq_along(p)) {
      q <- p
      q[[j]] <- c(first, p[[j]])
      out[[length(out) + 1L]] <- q
    }
  }
  out
}

#' The Terms of a Derivative of the Expectation of a Product of Two Scores
#'
#' @description
#' Expands \eqn{\partial_S E[\ell_a \ell_b]} as the expectation of a polynomial
#' in the derivatives of the log-density \eqn{\ell}.
#'
#' @details
#' With \eqn{f = e^{\ell}} the derivative of the density in the parameters of
#' a multiset \eqn{R} is \eqn{f} times the complete Bell polynomial
#' \eqn{Y_R} of the derivatives of \eqn{\ell}, the sum over the set partitions
#' of \eqn{R} of the products of \eqn{\ell^{(B)}} over the blocks. Then
#' \deqn{\partial_S E[\ell_a\ell_b] = E\Big[\sum_{T \subseteq S}
#'   \partial_T(\ell_a\ell_b)\, Y_{S \setminus T}\Big],\qquad
#'   \partial_T(\ell_a\ell_b) = \sum_{T = T_1 \sqcup T_2}
#'   \ell_{aT_1}\,\ell_{bT_2},}
#' where the sums run over the subsets of the positions of \eqn{S}. The terms
#' are collected by their sorted factors. The highest derivative of
#' \eqn{\ell} that appears is of order \eqn{|S| + 1}.
#'
#' @param P The parameter names.
#' @param a,b The indices of the pair.
#' @param S The integer indices of the derivative, with repetition.
#'
#' @return A list with `coef`, the coefficient of each term, `factors`, a list
#'   of the keys of each term's factors, and `blocks`, a named list giving the
#'   parameter names of each key.
#'
#' @keywords internal
d34e_plan <- function(P, a, b, S) {
  k <- length(S)
  key_of <- function(ix) canon_key(P[ix], P)
  coef <- numeric(0)
  fac <- list()
  blocks <- list()
  add <- function(parts) {
    keys <- vapply(parts, key_of, character(1))
    for (i in seq_along(parts)) {
      if (is.null(blocks[[keys[i]]])) {
        blocks[[keys[i]]] <<- P[parts[[i]]][order(match(P[parts[[i]]], P))]
      }
    }
    keys <- sort(keys)
    nm <- paste(keys, collapse = "|")
    if (is.na(coef[nm])) {
      coef[nm] <<- 0
      fac[[nm]] <<- keys
    }
    coef[nm] <<- coef[nm] + 1
  }
  pos <- seq_len(k)
  for (mask in 0:(2^k - 1)) {
    inT <- bitwAnd(mask, 2^(pos - 1)) > 0
    Tt <- S[inT]
    R <- S[!inT]
    nT <- length(Tt)
    for (s in 0:(2^nT - 1)) {
      side <- bitwAnd(s, 2^seq_len(nT)/2) > 0
      T1 <- Tt[side]
      T2 <- Tt[!side]
      for (pt in d34e_set_partitions(R)) {
        add(c(list(c(a, T1)), list(c(b, T2)), pt))
      }
    }
  }
  list(coef = unname(coef), factors = unname(fac), blocks = blocks)
}

#' Evaluate the Terms of a Derivative Plan
#'
#' @description
#' Sums the products of the factors of a plan from [d34e_plan()].
#'
#' @param plan A plan.
#' @param L A named list of the derivatives of the log-density, keyed as the
#'   factors of the plan.
#'
#' @return A numeric vector.
#'
#' @keywords internal
d34e_eval_plan <- function(plan, L) {
  total <- 0
  for (i in seq_along(plan$coef)) {
    term <- plan$coef[i]
    for (key in plan$factors[[i]]) term <- term * L[[key]]
    total <- total + term
  }
  total
}

#' The Third and Fourth Derivatives of a Truncated Family's Expected Information
#'
#' @description
#' Computes the derivatives of order `k` (3 or 4) of the expected Hessian of a
#' [truncated()] family on the parameter scale, as expectations under the
#' truncated law taken by [expectation()].
#'
#' @details
#' With \eqn{\ell = \ell_P - \log Z} the truncated log-likelihood and
#' \eqn{E_{ab} = -E[\ell_a\ell_b]}, the component for a multiset \eqn{S} is
#' minus the expectation of the polynomial of [d34e_plan()], which reads the
#' derivatives of \eqn{\ell} to order \eqn{k + 1}. The truncation points do
#' not depend on the parameters, so there is no boundary term. Those of
#' \eqn{\ell} are the parent's less those of \eqn{\log Z}, computed once per
#' parameter combination by `trunc_logz_derivs()` and passed to the integrand.
#' The derivatives of \eqn{\log Z} to order four are read from the parent's
#' distribution function where [has_exact_cdf_deriv()] holds, and the
#' parent's distribution function has no derivatives of order five, so each
#' ratio \eqn{d^B Z/Z} of that order is the expectation, under the truncated
#' law, of the complete Bell polynomial of the parent's log-density
#' derivatives (a quadrature or a sum), which is accurate but not a closed
#' form. Order five enters the fourth derivative only. The form
#' \eqn{E[\ell_{ab}]} would need the sixth derivative of the parent.
#'
#' @param d A truncated family, from [truncated()].
#' @param y The response, read for its length.
#' @param th The parameters.
#' @param k The order, `3L` or `4L`.
#' @param t Unused; it is the thread count of the registration.
#'
#' @return A list of components keyed as [d3expected_names()] or
#'   [d4expected_names()], each of length `length(y)`.
#'
#' @seealso [trunc_dexpected()]
#'
#' @keywords internal
trunc_d34expected <- function(d, y, th, k, t = 1L) {
  n <- length(y)
  parent <- d@parent_distrib
  P <- d@params
  p <- length(P)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  lz <- trunc_logz_derivs(d, th, k + 1L)
  tuples <- deriv_indices(P, k)
  keyfun <- if (k == 3L) d3expected_key else d4expected_key
  out <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (S in tuples) {
      plan <- d34e_plan(P, a, b, S)
      bl <- plan$blocks
      zs <- lz[names(bl)]
      names(zs) <- paste0(".lz_", names(bl))
      f <- function(y, theta, ...) {
        z <- list(...)
        ellp <- parent_ell(parent, y, theta, k + 1L, P)
        L <- lapply(names(bl), function(key) ellp(bl[[key]]) - z[[paste0(".lz_", key)]])
        names(L) <- names(bl)
        d34e_eval_plan(plan, L)
      }
      val <- -do.call(expectation, c(list(d, f, th), zs))
      out[[do.call(keyfun, as.list(c(list(P, a, b), S)))]] <- rep_len(val, n)
    }
  }
  out[if (k == 3L) d3expected_names(P) else d4expected_names(P)]
}

register_d34expected(TruncatedContinuousDistrib, trunc_d34expected)
register_d34expected(TruncatedDiscreteDistrib, trunc_d34expected)
