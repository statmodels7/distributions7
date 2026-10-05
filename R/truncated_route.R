#' @include truncated.R scalar_route.R
NULL

#' The Scalar Route of a Truncated Family
#'
#' @description
#' Returns the route of a [truncated()] family,
#' `"TruncatedDiscreteDistrib|<inner name>"` or
#' `"TruncatedContinuousDistrib|<inner name>"`, with the inner family's
#' constants followed by the truncation points `lower` and `upper`, as
#' [distrib_scalar_route()] describes, or `NULL` when the parent has no
#' route. A continuous parent needs a compiled distribution function, and
#' the only wrapper it may carry is [fixed()].
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
  cls <- attr(S7::S7_class(distrib), "name")
  name <- paste0(cls, "|", inner$name)
  if (d7_scalar_thread_safe_probe(name) < 0L) return(NULL)
  list(name = name,
       constants = c(inner$constants,
                     list(lower = distrib@lower, upper = distrib@upper)))
}
S7::method(distrib_scalar_route, TruncatedDiscreteDistrib) <- truncated_scalar_route
S7::method(distrib_scalar_route, TruncatedContinuousDistrib) <- truncated_scalar_route

#' The Retained Mass of a Truncated Family and the Sums of Its Information
#'
#' @description
#' Computes, for a [truncated()] family whose parent has a scalar route, the
#' retained mass \eqn{Z}, its derivatives, and the sums that the truncated
#' expected information and its first derivatives read, by the same
#' arithmetic as the compiled registry.
#'
#' @details
#' With \eqn{f}, \eqn{s_i} and \eqn{\ell_{ij}} the parent's density, score
#' and Hessian, the quantities are
#' \deqn{Z = \int f, \quad Z_i = \int f s_i, \quad
#'   Z_{ij} = \int f (\ell_{ij} + s_i s_j), \quad S_{ij} = \int f s_i s_j,}
#' \deqn{D_{ijc} = \int f (s_i s_j s_c + \ell_{ic} s_j + s_i \ell_{jc}),}
#' over the retained set, so that \eqn{E_T[s_i s_j] = S_{ij}/Z} and
#' \eqn{D_{ijc}} is the derivative of \eqn{\int f s_i s_j} in parameter
#' \eqn{c}.
#'
#' For a discrete parent the integrals are sums over support points.
#' `trunc_rule_cpp()` returns, for each observation, the points they run
#' over: the retained points when they are finitely many, and otherwise a
#' series over the retained points stopped by a rule on the masses alone.
#' For the density, the score and the Hessian on an unbounded support, the
#' points removed below `lower` are used instead when their mass is at most
#' one half, and the retained sums are \eqn{1 - \sum f}, \eqn{-\sum f s_i}
#' and \eqn{-\sum f(\ell_{ij} + s_i s_j)}.
#'
#' For a continuous parent, \eqn{Z}, \eqn{Z_i} and \eqn{Z_{ij}} are the
#' differences of the compiled distribution function and its derivatives at
#' the two truncation points (`trunc_cont_ends_cpp()`), \eqn{Z} in the
#' survival function when the lower point lies above the median; \eqn{S}
#' and \eqn{D} are taken by the rule of `trunc_cont_rule_cpp()`, whose nodes
#' depend on the interval and on the parent's center, scale and kinks only.
#'
#' The sums accumulate in long double in the registry's order, so that each
#' diagonal quantity is the registry's to the bit.
#'
#' @param distrib A truncated family.
#' @param y The observations, read for their number.
#' @param theta A named list of parameter values.
#' @param what One of `"z"`, `"grad"`, `"hess"`, `"info"` and `"dinfo"`,
#'   each adding quantities to those of the previous one.
#' @param rows `NULL`, or the observations to compute; the function sets it
#'   when it computes each distinct parameter vector once.
#'
#' @return `NULL` when the family has no route; otherwise a list with `Z`,
#'   `Zi` (named by parameter), and, as `what` requires, `Zij` and `S` (named
#'   as [hess_names()]) and `D` (named as [dexpected_names()]), each a vector
#'   over the observations.
#'
#' @keywords internal
trunc_route_parts <- function(distrib, y, theta, what, rows = NULL) {
  if (!is_truncated(distrib)) return(NULL)
  route <- distrib_scalar_route(distrib)
  if (is.null(route)) return(NULL)
  disc <- S7::S7_inherits(distrib, TruncatedDiscreteDistrib)
  lev <- match(what, c("z", "grad", "hess", "info", "dinfo"))
  parent <- distrib@parent_distrib
  P <- distrib@params
  p <- length(P)
  n <- max(length(y), lengths(theta))
  th <- lapply(align_theta(distrib, theta)[P], rep_len, length.out = n)
  tm <- do.call(cbind, c(unname(th), lapply(unname(route$constants), rep_len,
                                            length.out = n)))
  if (!is.matrix(tm)) tm <- matrix(tm, nrow = n)
  # the quantities depend on the parameters alone, so rows that repeat a
  # parameter vector are computed once: an intercept-only fit rebuilt the
  # quadrature rule at every observation for one vector. The rows are keyed
  # by their exact bits, and a sum over the support or a rule (discrete, or
  # orders 4 and 5) is worth the keying; the closed-form ends take the
  # shortcut only where every column is constant
  if (is.null(rows) && n > 1L) {
    first <- if (all(tm == rep(tm[1L, ], each = n), na.rm = FALSE) %in% TRUE) {
      rep.int(1L, n)
    } else if (disc || lev >= 4L) {
      key <- do.call(paste, c(lapply(seq_len(ncol(tm)), function(j)
        sprintf("%a", tm[, j])), sep = "|"))
      match(key, key)
    }
    if (!is.null(first) && anyDuplicated(first)) {
      u <- which(first == seq_len(n))
      out <- trunc_route_parts(distrib, y, theta, what, rows = u)
      at <- match(first, u)
      return(rapply(out, function(v) v[at], how = "list"))
    }
  }
  n0 <- n
  if (!is.null(rows)) {
    tm <- tm[rows, , drop = FALSE]
    th <- lapply(th, `[`, rows)
    n <- length(rows)
  }
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L] != pairs[, 2L], pairs[, 1L], pairs[, 2L]), ,
                 drop = FALSE]
  hn <- vapply(seq_len(nrow(pairs)), function(r)
    hess_pair_name(P, pairs[r, 1L], pairs[r, 2L]), "")

  no_mass <- function(Z) {
    if (any(!is.finite(Z)) || any(Z <= 0)) {
      stop(sprintf(paste0(
        "The truncation interval [%s, %s] carries no probability under these ",
        "parameter values (computed mass %s). A truncated distribution is not ",
        "defined there."
      ), format(distrib@lower), format(distrib@upper), format(min(Z))),
      call. = FALSE)
    }
  }

  if (disc) {
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
    comp <- rule$branch == 1L
    fw <- function(lp) exp(lp)
  } else {
    ends <- trunc_cont_ends_cpp(route$name, tm, min(lev - 1L, 2L))
    no_mass(ends$z)
    out <- list(Z = ends$z)
    if (lev >= 2L) out$Zi <- stats::setNames(lapply(seq_len(p), function(i) ends$g[, i]), P)
    if (lev >= 3L) out$Zij <- stats::setNames(lapply(seq_len(nrow(pairs)),
                                                     function(r) ends$h[, r]), hn)
    if (lev < 4L) return(out)
    rule <- trunc_cont_rule_cpp(route$name, tm)
    idx <- rule$idx
    ys <- rule$y
    fw <- function(lp) exp(lp) * rule$w
    # a row where the rule could not be placed has no nodes, and its sums
    # are not available, as the registry reports them
    unplaced <- tabulate(idx, n) == 0L
  }

  thp <- lapply(th, `[`, idx)
  # the parent is evaluated at points that belong to the rows idx, and a
  # constant that varies by observation (a binomial's size) is taken there
  # too: recycled against the points, it gave a mass of 0.912 where it is
  # 0.922, and read past its end once the rows were subset
  parent <- distrib_at_rows(parent, if (is.null(rows)) idx else rows[idx], n0)
  sumg <- function(v) {
    r <- ld_group_sum(as.numeric(v), idx, n)
    if (!disc) r[unplaced] <- NaN
    r
  }
  any_pts <- length(ys) > 0L
  # a beta parent is evaluated from log y and log(1 - y), which the rule
  # hands over where y itself rounds to 1
  logs <- if (!disc && any_pts && !is.null(rule$ly)) {
    beta_logs_parts_cpp(sub("^[^|]*[|]", "", route$name), rule$ly, rule$l1y,
                        tm[idx, , drop = FALSE])
  }
  f <- if (!any_pts) numeric(0) else if (!is.null(logs)) fw(logs$logpdf) else
    fw(distrib_pdf(parent, ys, thp, log = TRUE))
  live <- f != 0
  # a node whose weight is zero contributes nothing, whatever its integrand
  sumf <- if (disc) sumg else function(v) { v[!live] <- 0; sumg(v) }
  g <- if (any_pts && lev >= 2L) {
    if (!is.null(logs)) stats::setNames(lapply(seq_len(p), function(i) logs$g[, i]), P)
    else distrib_gradient(parent, ys, thp)
  }
  h <- if (any_pts && lev >= 3L) {
    if (!is.null(logs)) stats::setNames(lapply(seq_len(nrow(pairs)),
                                               function(r) logs$h[, r]), hn)
    else distrib_hessian(parent, ys, thp)
  }
  hp <- function(a, b) h[[hess_pair_name(P, a, b)]]

  if (disc) {
    Fs <- sumg(f)
    Z <- ifelse(comp, 1 - Fs, Fs)
    no_mass(Z)
    out <- list(Z = Z)
    if (lev < 2L) return(out)
    out$Zi <- stats::setNames(lapply(seq_len(p), function(i) {
      G <- sumg(f * g[[P[i]]])
      ifelse(comp, -G, G)
    }), P)
    if (lev < 3L) return(out)
    out$Zij <- stats::setNames(lapply(seq_len(nrow(pairs)), function(r) {
      a <- pairs[r, 1L]; b <- pairs[r, 2L]
      H <- sumg(f * (hp(a, b) + g[[P[a]]] * g[[P[b]]]))
      ifelse(comp, -H, H)
    }), hn)
    if (lev < 4L) return(out)
  }
  out$S <- stats::setNames(lapply(seq_len(nrow(pairs)), function(r) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    sumf(f * (g[[P[a]]] * g[[P[b]]]))
  }), hn)
  if (lev < 5L) return(out)
  D <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (c in seq_len(p)) {
      D[[dexpected_key(P, a, b, c)]] <- sumf(f * (g[[P[a]]] * g[[P[b]]] * g[[P[c]]] +
                                                    hp(a, c) * g[[P[b]]] +
                                                    g[[P[a]]] * hp(b, c)))
    }
  }
  out$D <- D
  out
}

#' A Family at Selected Observations
#'
#' @description
#' Returns `distrib` with every constant that varies by observation (a
#' binomial's or a beta-binomial's `size`), its parents' included, taken at
#' the observations `idx`, so that the family can be evaluated at points that
#' belong to those observations.
#'
#' @param distrib A univariate family.
#' @param idx Integer indices of observations, one per point.
#' @param n The number of observations the constants are recycled to.
#'
#' @return The family, with its varying constants of length `length(idx)`.
#'
#' @keywords internal
distrib_at_rows <- function(distrib, idx, n) {
  if (S7::prop_exists(distrib, "parent_distrib")) {
    distrib@parent_distrib <- distrib_at_rows(distrib@parent_distrib, idx, n)
  }
  if (S7::prop_exists(distrib, "size") && length(distrib@size) > 1L) {
    distrib@size <- rep_len(distrib@size, n)[idx]
  }
  distrib
}
