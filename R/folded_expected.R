#' @include folded.R expected_loc_scale.R scalar_route.R dexpected_hessian.R
NULL

#' The Expected Information of a Folded Family
#'
#' @description
#' Computes the expected information of a [folded()] family (`order = 0`)
#' or its first derivatives (`order = 1`), on the parameter scale, as
#' expectations under the folded law.
#'
#' @details
#' An expectation under the folded law is an integral over \eqn{y > 0}
#' against the folded density \eqn{f(y) + f(-y)}. The integrand turns
#' sharply near zero, where the weight
#' \eqn{f(y)/(f(y)+f(-y))} goes from one half to one over a width of order
#' \eqn{s^2/|c|}, \eqn{c} and \eqn{s} a center and a scale of the parent
#' read from the scalar registry. The integral is taken over \eqn{y > 0} by a
#' tanh-sinh rule on \eqn{[0, |c|]}, whose nodes cluster at both ends, and
#' the rule of [loc_scale_rule()] scaled by \eqn{s} on \eqn{[|c|, \infty)};
#' the split at \eqn{|c|} also takes in a parent with a kink at its center.
#' The nodes are the scalar registry's. The expectations are formed through
#' the second Bartlett identity,
#' \deqn{E[\ell_{ij}] = -E[\ell_i \ell_j], \qquad
#'   \partial_c E[\ell_{ij}] = -E[\ell_{ic}\ell_j + \ell_i\ell_{jc} +
#'   \ell_i\ell_j\ell_c],}
#' so that the folded score and Hessian are all it reads. The sums are
#' accumulated as the scalar registry accumulates them. A parent that the
#' registry does not cover, or that is not a location-scale family on the
#' real line, is integrated by [expected_derivative()] with
#' `approx = "integrate"` instead.
#'
#' @param distrib A [folded()] family.
#' @param y The response, read for its length.
#' @param theta The parameters.
#' @param order `0L` for the expected information, `1L` for its first
#'   derivatives.
#'
#' @return A named list keyed as [hess_names()] (`order = 0`) or
#'   [dexpected_names()] (`order = 1`), each component of length
#'   `length(y)`.
#'
#' @seealso [distrib_expected_hessian()], [distrib_dexpected_hessian()]
#'
#' @keywords internal
folded_expected <- function(distrib, y, theta, order) {
  n <- length(y)
  parent <- distrib@parent_distrib
  P <- distrib@params
  p <- length(P)
  inner <- distrib_scalar_route(parent)
  cs <- NULL
  if (!is.null(inner) && !grepl("|", inner$name, fixed = TRUE)) {
    th <- lapply(align_theta(distrib, theta)[P], rep_len, length.out = n)
    M <- do.call(cbind, c(th, lapply(inner$constants, rep_len, length.out = n)))
    if (!is.matrix(M)) M <- matrix(M, nrow = n)
    key <- do.call(paste, c(as.data.frame(M), sep = "\r"))
    first <- !duplicated(key)
    U <- M[first, , drop = FALSE]
    idx <- match(key, key[first])
    cs <- d7_center_scale_probe(inner$name, U)
    if (anyNA(cs)) cs <- NULL
  }
  if (is.null(cs)) {
    if (order == 0L) {
      return(expected_derivative(distrib, y, theta, order = 2L,
                                 approx = "integrate"))
    }
    return(numerical_dexpected_hessian(distrib, y, theta, "parameter",
                                       "integrate", 10000))
  }
  m <- nrow(U)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  keys <- if (order == 0L) hess_names(P) else dexpected_names(P)
  vals <- matrix(NA_real_, m, length(keys), dimnames = list(NULL, keys))
  for (u in seq_len(m)) {
    rule <- fold_rule_cpp(cs[u, 1], cs[u, 2])
    yv <- rule$y
    thu <- lapply(seq_len(p), function(j) rep(U[u, j], length(yv)))
    names(thu) <- P
    fw <- distrib_pdf(distrib, yv, thu) * rule$w
    live <- fw != 0
    g <- distrib_gradient(distrib, yv, thu)
    h <- if (order >= 1L) distrib_hessian(distrib, yv, thu)
    total <- function(v) {
      v[!live] <- 0
      out <- sum(v * fw)
      if (is.finite(out)) -out else NA_real_
    }
    for (r in seq_len(nrow(pairs))) {
      a <- pairs[r, 1L]; b <- pairs[r, 2L]
      ga <- g[[P[a]]]; gb <- g[[P[b]]]
      if (order == 0L) {
        vals[u, hess_pair_name(P, a, b)] <- total(ga * gb)
        next
      }
      for (c in seq_len(p)) {
        gc <- g[[P[c]]]
        hac <- h[[hess_pair_name(P, a, c)]]
        hbc <- h[[hess_pair_name(P, b, c)]]
        vals[u, dexpected_key(P, a, b, c)] <- total(hac * gb + ga * hbc + ga * gb * gc)
      }
    }
  }
  out <- lapply(keys, function(k) unname(vals[idx, k]))
  names(out) <- keys
  out
}

S7::method(distrib_expected_hessian, FoldedDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  folded_expected(distrib, y, theta, 0L)
}

S7::method(distrib_dexpected_hessian, FoldedDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                     function(k) folded_expected(distrib, y, theta, 1L))
}

S7::method(expected_hessian_costly, FoldedDistrib) <- function(x, ...) TRUE

#' The Scalar Route of a Folded Family
#'
#' @description
#' Returns the route of a [folded()] family, `"FoldedDistrib|<inner name>"`
#' with the inner family's constants, or `NULL` when the parent has no route
#' or is not a location-scale family on the real line, the expected
#' information's quadrature reading the parent's center and scale.
#'
#' @param distrib A folded family.
#' @param ... Unused.
#'
#' @return A list with components `name` and `constants`, or `NULL`.
#'
#' @keywords internal
folded_scalar_route <- function(distrib, ...) {
  inner <- distrib_scalar_route(distrib@parent_distrib)
  if (is.null(inner) || grepl("|", inner$name, fixed = TRUE)) return(NULL)
  np <- length(distrib@params) + length(inner$constants)
  if (anyNA(d7_center_scale_probe(inner$name, matrix(1, 1, np)))) return(NULL)
  list(name = paste0("FoldedDistrib|", inner$name), constants = inner$constants)
}
S7::method(distrib_scalar_route, FoldedDistrib) <- folded_scalar_route
