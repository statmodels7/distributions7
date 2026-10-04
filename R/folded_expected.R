#' @include folded.R expected_loc_scale.R scalar_route.R dexpected_hessian.R
NULL

#' The Expected Information of a Folded Family
#'
#' @description
#' Computes the expected information of a [folded()] family (`order = 0`),
#' its first derivatives (`order = 1`) or its first and second derivatives
#' (`order = 2`), on the parameter scale, as expectations under the folded
#' law.
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
#' so that the folded score and Hessian are all it reads at orders 0 and 1.
#' The second derivatives differentiate \eqn{-E[G]}, \eqn{G = \ell_i\ell_j},
#' twice under the integral,
#' \deqn{\partial_{cd} E[G] = E[G_{cd} + G_c\ell_d + G_d\ell_c +
#'   G(\ell_c\ell_d + \ell_{cd})],}
#' and read the folded third derivatives as well. At order 2 a parent
#' without a route signals an error. The sums are
#' accumulated as the scalar registry accumulates them. A parent that the
#' registry does not cover, or that is not a location-scale family on the
#' real line, is integrated by [expected_derivative()] with
#' `approx = "integrate"` instead.
#'
#' @param distrib A [folded()] family.
#' @param y The response, read for its length.
#' @param theta The parameters.
#' @param order `0L` for the expected information, `1L` for its first
#'   derivatives, `2L` for its first and second derivatives.
#'
#' @return A named list keyed as [hess_names()] (`order = 0`),
#'   [dexpected_names()] (`order = 1`), or [dexpected_names()] followed by
#'   [d2expected_names()] (`order = 2`), each component of length
#'   `length(y)`.
#'
#' @seealso [distrib_expected_hessian()], [distrib_dexpected_hessian()],
#'   [distrib_d2expected_hessian()]
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
    if (order == 2L) {
      stop(sprintf(paste0(
        "'%s' has no analytic second derivative of its expected information:\n",
        "  its parent has no scalar route with a center and a scale."),
        distrib@distrib_name), call. = FALSE)
    }
    if (order == 0L) {
      return(expected_derivative(distrib, y, theta, order = 2L,
                                 approx = "integrate"))
    }
    return(numerical_dexpected_hessian(distrib, y, theta, "parameter",
                                       "integrate", 10000))
  }
  m <- nrow(U)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  keys <- switch(order + 1L, hess_names(P), dexpected_names(P),
                 c(dexpected_names(P), d2expected_names(P)))
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
    t3 <- if (order >= 2L) distrib_deriv3(distrib, yv, thu)
    l3 <- function(a, b, c) t3[[paste(P[sort(c(a, b, c))], collapse = "_")]]
    total <- function(v, extra = 0) {
      v[!live] <- 0
      out <- sum(v * fw) + extra
      if (is.finite(out)) -out else NA_real_
    }
    B <- if (order >= 1L) fold_kink_parts(distrib, U[u, seq_len(p)], cs[u, 1], order)
    J <- function(v) if (is.null(B)) 0 else B$L * (v[1L] - v[2L])
    ap <- if (is.null(B)) numeric(p) else B$ap
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
        jb <- if (ap[c] == 0) 0 else J(B$g[[a]] * B$g[[b]]) * ap[c]
        vals[u, dexpected_key(P, a, b, c)] <- total(hac * gb + ga * hbc + ga * gb * gc, jb)
      }
      if (order < 2L) next
      # d_cd E[l_a l_b] = E[G_cd + G_c s_d + G_d s_c + G (s_c s_d + l_cd)],
      # G = l_a l_b
      for (r2 in seq_len(nrow(pairs))) {
        c <- pairs[r2, 1L]; d <- pairs[r2, 2L]
        gc <- g[[P[c]]]; gd <- g[[P[d]]]
        hac <- h[[hess_pair_name(P, a, c)]]; hbc <- h[[hess_pair_name(P, b, c)]]
        had <- h[[hess_pair_name(P, a, d)]]; hbd <- h[[hess_pair_name(P, b, d)]]
        G <- ga * gb
        Gc <- hac * gb + ga * hbc
        Gd <- had * gb + ga * hbd
        Gcd <- l3(a, c, d) * gb + hac * hbd + had * hbc + ga * l3(b, c, d)
        jb <- 0
        if (ap[c] != 0 || ap[d] != 0) {
          # the boundary terms at the kink, read one-sidedly
          bg <- B$g; bh <- function(i, j) B$h[[hess_pair_name(P, i, j)]]
          Gb <- bg[[a]] * bg[[b]]
          Kc <- bh(a, c) * bg[[b]] + bg[[a]] * bh(b, c) + Gb * bg[[c]]
          Kd <- bh(a, d) * bg[[b]] + bg[[a]] * bh(b, d) + Gb * bg[[d]]
          Y <- B$liy[[a]] * bg[[b]] + bg[[a]] * B$liy[[b]] + Gb * B$ly
          jb <- J(Kc) * ap[d] + J(Kd) * ap[c] + J(Y) * ap[c] * ap[d]
        }
        vals[u, d2expected_key(P, a, b, c, d)] <-
          total(Gcd + Gc * gd + Gd * gc + G * (gc * gd + h[[hess_pair_name(P, c, d)]]), jb)
      }
    }
  }
  out <- lapply(keys, function(k) unname(vals[idx, k]))
  names(out) <- keys
  out
}

#' @rdname folded_expected
#' @details
#' A parent with a kink at its center \eqn{c} (a parameter that is not
#' smooth, as in [laplace_distrib()], [laplace2_distrib()] and
#' [enet_distrib()]) gives a folded score that jumps at \eqn{y^* = |c|},
#' a point that moves with \eqn{c}. With \eqn{F = L\,G} the integrand,
#' \eqn{J(X) = X(y^{*-}) - X(y^{*+})} and \eqn{a_k = \partial y^*/\partial
#' \theta_k = \mathrm{sign}(c)} for the center and zero otherwise,
#' \deqn{\partial_k \int F = \int \partial_k F + J(F)\,a_k,}
#' \deqn{\partial_{kl} \int F = \int \partial_{kl} F + J(\partial_k F)\,a_l
#'   + J(\partial_l F)\,a_k + J(\partial_y F)\,a_k a_l,}
#' the one-sided values read four units in the last place either side of
#' \eqn{y^*}. `fold_kink_parts()` returns those values: the folded density at
#' \eqn{y^*}, the folded score and Hessian, and, at order 2, the response
#' derivative \eqn{\ell_y} and the mixed derivatives \eqn{\ell_{iy}} from
#' the parent's at the two preimages, or `NULL` where the parent has no kink
#' or \eqn{c = 0}.
#' @param th One row of parameters, in the order of the family's.
#' @param center The parent's center at that row.
fold_kink_parts <- function(distrib, th, center, order) {
  P <- distrib@params
  kp <- which(!distrib@params_smooth[P])
  if (length(kp) != 1L || center == 0 || th[[kp]] != center) return(NULL)
  parent <- distrib@parent_distrib
  ys <- abs(center)
  yb <- ys * (1 + c(-4, 4) * .Machine$double.eps)
  thb <- lapply(seq_along(P), function(j) rep(th[[j]], 2L))
  names(thb) <- P
  ap <- numeric(length(P))
  ap[kp] <- sign(center)
  g <- distrib_gradient(distrib, yb, thb)
  out <- list(L = fold_parts(parent, ys, lapply(thb, `[`, 1L))$L,
              g = unname(g[P]), h = distrib_hessian(distrib, yb, thb), ap = ap)
  if (order < 2L) return(out)
  w <- fold_parts(parent, yb, thb)$w
  gp <- distrib_gradient(parent, yb, thb)
  gm <- distrib_gradient(parent, -yb, thb)
  yp <- distrib_grad_y(parent, yb, thb)
  ym <- distrib_grad_y(parent, -yb, thb)
  cp <- distrib_cross_y(parent, yb, thb)
  cm <- distrib_cross_y(parent, -yb, thb)
  ly <- w * yp - (1 - w) * ym
  out$ly <- ly
  out$liy <- lapply(seq_along(P), function(i) {
    w * (cp[[P[i]]] + gp[[P[i]]] * yp) -
      (1 - w) * (cm[[P[i]]] + gm[[P[i]]] * ym) - out$g[[i]] * ly
  })
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

S7::method(distrib_d2expected_hessian, FoldedDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  dexpected_analytic(distrib, y, theta, match.arg(scale), 2L, threads,
                     function(k) folded_expected(distrib, y, theta, k))
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
