#' @include truncated.R scalar_route.R
NULL

#' The Scalar Route of a Truncated Family
#'
#' @description
#' Returns the route of a [truncated()] family with a discrete parent,
#' `"TruncatedDiscreteDistrib|<inner name>"`, with the inner family's
#' constants followed by the truncation points `lower` and `upper`, as
#' [distrib_scalar_route()] describes, or `NULL` when the parent has no
#' route.
#'
#' @param distrib A truncated family.
#' @param ... Unused.
#'
#' @return A list with components `name` and `constants`, or `NULL`.
#'
#' @keywords internal
truncated_scalar_route <- function(distrib, ...) {
  inner <- distrib_scalar_route(distrib@parent_distrib)
  if (is.null(inner)) return(NULL)
  list(name = paste0("TruncatedDiscreteDistrib|", inner$name),
       constants = c(inner$constants,
                     list(lower = distrib@lower, upper = distrib@upper)))
}
S7::method(distrib_scalar_route, TruncatedDiscreteDistrib) <- truncated_scalar_route

#' The Sums of a Truncated Family over Support Points
#'
#' @description
#' Computes, for a [truncated()] family with a discrete parent that has a
#' scalar route, the retained mass \eqn{Z} and the sums that the truncated
#' density, score, Hessian, expected information and its first derivatives
#' read, over the support points that the compiled registry uses.
#'
#' @details
#' With \eqn{f}, \eqn{s_i} and \eqn{\ell_{ij}} the parent's mass, score and
#' Hessian at a support point, the sums are
#' \deqn{Z = \sum f, \quad Z_i = \sum f s_i, \quad
#'   Z_{ij} = \sum f (\ell_{ij} + s_i s_j), \quad S_{ij} = \sum f s_i s_j,}
#' \deqn{D_{ijc} = \sum f (s_i s_j s_c + \ell_{ic} s_j + s_i \ell_{jc}),}
#' so that \eqn{E_T[s_i s_j] = S_{ij}/Z} and \eqn{D_{ijc}} is the derivative
#' of \eqn{\sum f s_i s_j} in parameter \eqn{c}. `trunc_rule_cpp()` returns,
#' for each observation, the points the sums run over: the retained points
#' when they are finitely many, and otherwise a series over the retained
#' points stopped by a rule on the masses alone. For the density, the score
#' and the Hessian on an unbounded support, the points removed below `lower`
#' are used instead when their mass is at most one half, and the retained
#' sums are \eqn{1 - \sum f}, \eqn{-\sum f s_i} and
#' \eqn{-\sum f(\ell_{ij} + s_i s_j)}. The sums accumulate in long double in
#' increasing support point, as the registry's do, so that each diagonal
#' quantity is the registry's to the bit.
#'
#' @param distrib A truncated family.
#' @param y The observations, read for their number.
#' @param theta A named list of parameter values.
#' @param what One of `"z"`, `"grad"`, `"hess"`, `"info"` and `"dinfo"`,
#'   each adding sums to those of the previous one.
#'
#' @return `NULL` when the family has no route with a discrete parent;
#'   otherwise a list with `Z`, `Zi` (named by parameter), and, as `what`
#'   requires, `Zij` and `S` (named as [hess_names()]) and `D` (named as
#'   [dexpected_names()]), each a vector over the observations.
#'
#' @keywords internal
trunc_route_parts <- function(distrib, y, theta, what) {
  if (!S7::S7_inherits(distrib, TruncatedDiscreteDistrib)) return(NULL)
  route <- distrib_scalar_route(distrib)
  if (is.null(route)) return(NULL)
  lev <- match(what, c("z", "grad", "hess", "info", "dinfo"))
  parent <- distrib@parent_distrib
  P <- distrib@params
  p <- length(P)
  n <- max(length(y), lengths(theta))
  th <- lapply(align_theta(distrib, theta)[P], rep_len, length.out = n)
  tm <- do.call(cbind, c(unname(th), lapply(unname(route$constants), rep_len,
                                            length.out = n)))
  if (!is.matrix(tm)) tm <- matrix(tm, nrow = n)
  rule <- trunc_rule_cpp(route$name, tm, lev >= 4L)
  if (anyNA(rule$y1)) {
    stop(sprintf(paste0(
      "The sum over the retained support of '%s' did not reach its stopping ",
      "rule within %g terms."), distrib@distrib_name, 1e7), call. = FALSE)
  }
  len <- pmax(rule$y1 - rule$y0 + 1, 0)
  idx <- rep.int(seq_len(n), len)
  ys <- unlist(lapply(seq_len(n), function(i) {
    if (len[i] > 0) seq(rule$y0[i], rule$y1[i]) else numeric(0)
  }))
  thp <- lapply(th, `[`, idx)
  sumg <- function(v) ld_group_sum(as.numeric(v), idx, n)
  comp <- rule$branch == 1L
  any_pts <- length(ys) > 0L
  f <- if (any_pts) exp(distrib_pdf(parent, ys, thp, log = TRUE)) else numeric(0)
  g <- if (any_pts && lev >= 2L) distrib_gradient(parent, ys, thp)
  h <- if (any_pts && lev >= 3L) distrib_hessian(parent, ys, thp)
  hp <- function(a, b) h[[hess_pair_name(P, a, b)]]

  Fs <- sumg(f)
  Z <- ifelse(comp, 1 - Fs, Fs)
  if (any(!is.finite(Z)) || any(Z <= 0)) {
    stop(sprintf(paste0(
      "The truncation interval [%s, %s] carries no probability under these ",
      "parameter values (computed mass %s). A truncated distribution is not ",
      "defined there."
    ), format(distrib@lower), format(distrib@upper), format(min(Z))),
    call. = FALSE)
  }
  out <- list(Z = Z)
  if (lev < 2L) return(out)
  out$Zi <- stats::setNames(lapply(seq_len(p), function(i) {
    G <- sumg(f * g[[P[i]]])
    ifelse(comp, -G, G)
  }), P)
  if (lev < 3L) return(out)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  hn <- vapply(seq_len(nrow(pairs)), function(r)
    hess_pair_name(P, pairs[r, 1L], pairs[r, 2L]), "")
  out$Zij <- stats::setNames(lapply(seq_len(nrow(pairs)), function(r) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    H <- sumg(f * (hp(a, b) + g[[P[a]]] * g[[P[b]]]))
    ifelse(comp, -H, H)
  }), hn)
  if (lev < 4L) return(out)
  out$S <- stats::setNames(lapply(seq_len(nrow(pairs)), function(r) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    sumg(f * (g[[P[a]]] * g[[P[b]]]))
  }), hn)
  if (lev < 5L) return(out)
  D <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (c in seq_len(p)) {
      key <- dexpected_key(P, a, b, c)
      D[[key]] <- sumg(f * (g[[P[a]]] * g[[P[b]]] * g[[P[c]]] + hp(a, c) * g[[P[b]]] +
                              g[[P[a]]] * hp(b, c)))
    }
  }
  out$D <- D
  out
}
