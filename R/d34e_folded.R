#' @include folded.R folded_expected.R d34e_truncated.R d34expected_analytic.R
NULL

# The third and fourth derivatives of the expected information of a folded
# family, on the parameter scale, as expectations under the folded law by the
# quadrature of folded_expected().

#' The Third and Fourth Derivatives of a Folded Family's Expected Information
#'
#' @description
#' Computes the derivatives of order `k` (3 or 4) of the expected Hessian of a
#' [folded()] family on the parameter scale, as integrals over \eqn{y > 0}
#' against the folded density, on the nodes of `fold_rule_cpp()` that
#' [folded_expected()] uses.
#'
#' @details
#' With \eqn{\ell} the folded log-density and \eqn{E_{ab} = -E[\ell_a\ell_b]},
#' the component for a multiset \eqn{S} is minus the integral of the
#' polynomial of [d34e_plan()] against the folded density. The polynomial
#' reads the folded derivatives of \eqn{\ell} to order \eqn{k + 1}: those of
#' orders one to four are the family's own methods, and those of order five
#' are formed by `fold_deriv_k(5L)` from the parent's fifth derivatives, which
#' the fourth derivative of the information needs. The form
#' \eqn{E[\ell_{ab}]} would need the sixth derivative of the parent.
#'
#' Where the parent's mass beyond the fold is below the square of the machine
#' epsilon, the components are the parent's own derivatives of the expected
#' information, as in [folded_expected()].
#'
#' A parent with a kink at its center (a parameter that is not smooth) moves
#' the point where the folded score jumps, and the integral then gains
#' boundary terms that involve the mixed derivatives of \eqn{\ell} in the
#' response and the parameters beyond the first (the package has
#' `distrib_cross_y()` only for one response derivative and one parameter).
#' These parents, and a parent with no scalar route, are served by
#' [numerical_dexpected_higher()], one stencil on the second derivative.
#'
#' @param d A folded family, from [folded()].
#' @param y The response, read for its length.
#' @param th The parameters.
#' @param k The order, `3L` or `4L`.
#' @param t Unused; it is the thread count of the registration.
#'
#' @return A list of components keyed as [d3expected_names()] or
#'   [d4expected_names()], each of length `length(y)`.
#'
#' @seealso [folded_expected()], [trunc_d34expected()]
#'
#' @keywords internal
folded_d34expected <- function(d, y, th, k, t = 1L) {
  n <- length(y)
  parent <- d@parent_distrib
  P <- d@params
  p <- length(P)
  keys <- if (k == 3L) d3expected_names(P) else d4expected_names(P)
  keyfun <- if (k == 3L) d3expected_key else d4expected_key
  inner <- distrib_scalar_route(parent)
  cs <- NULL
  if (!is.null(inner)) {
    thl <- lapply(align_theta(d, th)[P], rep_len, length.out = n)
    M <- do.call(cbind, c(thl, lapply(inner$constants, rep_len, length.out = n)))
    if (!is.matrix(M)) M <- matrix(M, nrow = n)
    key <- do.call(paste, c(as.data.frame(M), sep = "\r"))
    first <- !duplicated(key)
    U <- M[first, , drop = FALSE]
    idx <- match(key, key[first])
    cs <- d7_center_scale_probe(inner$name, U)
    if (anyNA(cs)) cs <- NULL
  }
  if (is.null(cs) || !all(param_smoothness(d)[P])) {
    return(unclass(numerical_dexpected_higher(d, y, th, "parameter", k)))
  }
  m <- nrow(U)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  tuples <- deriv_indices(P, k)
  plans <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (S in tuples) {
      plans[[do.call(keyfun, as.list(c(list(P, a, b), S)))]] <- d34e_plan(P, a, b, S)
    }
  }
  deriv5 <- fold_deriv_k(5L)
  vals <- matrix(NA_real_, m, length(keys), dimnames = list(NULL, keys))
  for (u in seq_len(m)) {
    thp <- stats::setNames(as.list(U[u, seq_len(p)]), P)
    reflected <- distrib_cdf(parent, 0, thp, lower.tail = cs[u, 1] > 0)
    if (is.finite(reflected) && reflected < .Machine$double.eps^2) {
      ref <- if (k == 3L) distrib_d3expected_hessian(parent, 1, thp)
             else distrib_d4expected_hessian(parent, 1, thp)
      vals[u, ] <- vapply(keys, function(key) ref[[key]][1L], numeric(1))
      next
    }
    rule <- fold_rule_cpp(cs[u, 1], cs[u, 2])
    yv <- rule$y
    thu <- lapply(seq_len(p), function(j) rep(U[u, j], length(yv)))
    names(thu) <- P
    fw <- distrib_pdf(d, yv, thu) * rule$w
    live <- fw != 0
    tab <- list(distrib_gradient(d, yv, thu), distrib_hessian(d, yv, thu),
                distrib_deriv3(d, yv, thu))
    if (k >= 3L) tab[[4L]] <- distrib_deriv4(d, yv, thu)
    if (k >= 4L) tab[[5L]] <- deriv5(d, yv, thu)
    for (key in keys) {
      plan <- plans[[key]]
      L <- lapply(plan$blocks, function(bl) tab[[length(bl)]][[canon_key(bl, P)]])
      v <- d34e_eval_plan(plan, L)
      v[!live] <- 0
      out <- sum(v * fw)
      vals[u, key] <- if (is.finite(out)) -out else NA_real_
    }
  }
  out <- lapply(keys, function(key) unname(vals[idx, key]))
  names(out) <- keys
  out
}

register_d34expected(FoldedDistrib, folded_d34expected)
