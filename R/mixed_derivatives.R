#' Mixed Derivatives by One Stencil on the Highest Analytical Quantity
#'
#' @description
#' Computes the derivatives of the log-density of order `ay` in the response
#' and `ctheta` in the parameters, for a continuous family that has no method
#' of its own for them. Among the quantities the family implements itself
#' (the log-density, its response derivatives, its score and Hessian in the
#' parameters, and the mixed derivatives of lower order), the one that leaves
#' the smallest order to difference is taken, and each component is a single
#' tensor-product central difference of it, never a difference of a
#' difference.
#'
#' @details
#' A candidate quantity of order \eqn{(a, c)} in the response and the
#' parameters can serve where \eqn{a \le} `ay` and \eqn{c \le} `ctheta`, and it
#' leaves \eqn{r = (\mathrm{ay} - a) + (\mathrm{ctheta} - c)} orders to the
#' stencil. The log-density always serves; another quantity serves only where
#' the family's class registers its method, as [owns_method()] reads it, so
#' that a numerical fallback is never differenced again. The candidate with
#' the smallest \eqn{r} is used, and among equal ones the one with more
#' derivatives in the parameters.
#'
#' The stencil is the tensor product of one-dimensional central stencils of
#' [numericals7::fd_weights()] at second-order accuracy, one in the response
#' and one in each parameter that still has to be differenced, with the
#' parameters that occur most often taken by the analytical quantity first,
#' as [tensor_derivatives()] does. The step is
#' \eqn{h = \varepsilon^{1/(r+2)}\max(1, |x|)} in every direction, the
#' response step clamped inside the support by [fd_steps_y()] and a
#' parameter's step inside its bounds. Values shared by several components
#' are computed once.
#'
#' For [transformation()], which carries the parent's score and Hessian in
#' the parameters and no response derivative, the fourth-order
#' \eqn{\partial^4\ell/\partial y^2\,\partial\theta^2} is a second difference
#' in the response of the analytical Hessian: \eqn{r = 2}, where the previous
#' fallback differenced a difference in the same parameter.
#'
#' @param distrib A `continuous_distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param ay The order in the response, 1 to 3.
#' @param ctheta The order in the parameters, 1 or 2.
#'
#' @return A named list: one vector per parameter, keyed by `distrib@params`,
#'   where `ctheta` is 1; one vector per unordered pair, keyed as
#'   [`hess_names(distrib@params)`][hess_names], where it is 2.
#'
#' @seealso [tensor_derivatives()] for the derivatives in the parameters alone,
#'   [distrib_cross_y()], [distrib_cross2_y()], [distrib_cross3_y()],
#'   [distrib_grad_y_hess()] and [distrib_hess_y_hess()], whose continuous
#'   fallbacks call this.
#' @keywords internal
#'
#' @examples
#' # transformation() has no response derivative of its own: the mixed
#' # fourth order of a log-gamma is a stencil on its analytical Hessian.
#' d <- fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)
#' b <- c(-0.4, 0.3)
#' k <- 1 / 0.26
#' got <- distributions7:::mixed_tensor_derivatives(d, b, list(sigma2 = 0.26), 2, 2)
#' rbind(got$sigma2_sigma2, -2 * k^3 * exp(b))
mixed_tensor_derivatives <- function(distrib, y, theta, ay, ctheta) {
  ay <- as.integer(ay)
  ctheta <- as.integer(ctheta)
  params <- distrib@params
  p <- length(params)
  prs <- hess_pairs(params)
  pair_name <- function(i, j) {
    hit <- vapply(prs, function(pr) setequal(pr, params[c(i, j)]) &&
                    (i != j || pr[1] == pr[2]), logical(1))
    names(prs)[which(hit)[1]]
  }

  # the candidates, in (order in y, order in theta)
  cands <- list(
    list(gen = NULL, a = 0L, c = 0L),
    list(gen = distrib_grad_y, a = 1L, c = 0L),
    list(gen = distrib_hess_y, a = 2L, c = 0L),
    list(gen = distrib_deriv3_y, a = 3L, c = 0L),
    list(gen = distrib_deriv4_y, a = 4L, c = 0L),
    list(gen = distrib_gradient, a = 0L, c = 1L),
    list(gen = distrib_hessian, a = 0L, c = 2L),
    list(gen = distrib_cross_y, a = 1L, c = 1L),
    list(gen = distrib_cross2_y, a = 2L, c = 1L),
    list(gen = distrib_cross3_y, a = 3L, c = 1L),
    list(gen = distrib_grad_y_hess, a = 1L, c = 2L))
  usable <- Filter(function(b) {
    b$a <= ay && b$c <= ctheta && !(b$a == ay && b$c == ctheta) &&
      (is.null(b$gen) || owns_method(distrib, b$gen))
  }, cands)
  r_of <- vapply(usable, function(b) (ay - b$a) + (ctheta - b$c), integer(1))
  best <- which(r_of == min(r_of))
  base <- usable[[best[which.max(vapply(usable[best], function(b) b$c,
                                        integer(1)))]]]
  r_tot <- min(r_of)
  h_rel <- .Machine$double.eps^(1 / (r_tot + 2))

  ry <- ay - base$a
  hy <- if (ry > 0L) fd_steps_y(y, distrib@bounds, h_rel, "clamp", order = ry)
  hth <- lapply(seq_len(p), function(j) {
    x <- theta[[params[j]]]
    hj <- h_rel * pmax(1, abs(x))
    reach <- numericals7::fd_offsets(max(1L, ctheta - base$c), 2L)$reach
    b <- distrib@params_bounds[[params[j]]]
    if (!is.null(b)) {
      if (is.finite(b[1])) hj <- pmin(hj, 0.49 * (x - b[1]) / reach)
      if (is.finite(b[2])) hj <- pmin(hj, 0.49 * (b[2] - x) / reach)
    }
    hj
  })

  cache <- list()
  value_at <- function(oy, oth) {
    key <- paste(c(oy, oth), collapse = ",")
    v <- cache[[key]]
    if (is.null(v)) {
      yy <- if (oy != 0L) y + oy * hy else y
      th <- theta
      for (j in which(oth != 0L)) {
        th[[params[j]]] <- theta[[params[j]]] + oth[j] * hth[[j]]
      }
      v <- if (is.null(base$gen)) {
        list(distrib_pdf(distrib, yy, th, log = TRUE))
      } else {
        out <- base$gen(distrib, yy, th)
        if (is.list(out)) out else list(out)
      }
      cache[[key]] <<- v
    }
    v
  }
  stencil <- function(o) {
    off <- numericals7::fd_offsets(o, 2L)$central
    w <- numericals7::fd_weights(off, o)
    keep <- w != 0
    list(off = off[keep], w = w[keep])
  }

  one <- function(idx) {
    # idx: the parameters of the component, a multiset of positions
    cnt <- tabulate(idx, p)
    b <- integer(p)
    for (s in seq_len(base$c)) {
      j <- which.max(cnt - b)
      b[j] <- b[j] + 1L
    }
    rem <- cnt - b
    bj <- rep(seq_len(p), b)
    comp <- if (base$c == 0L) 1L else if (base$c == 1L) params[bj] else
      pair_name(bj[1], bj[2])
    dirs <- which(rem > 0L)
    st <- c(if (ry > 0L) list(stencil(ry)), lapply(dirs, function(j) stencil(rem[j])))
    grid <- as.matrix(expand.grid(lapply(st, function(s) seq_along(s$off))))
    acc <- 0
    for (g in seq_len(nrow(grid))) {
      oy <- 0L
      oth <- integer(p)
      wt <- 1
      q0 <- 0L
      if (ry > 0L) {
        oy <- st[[1]]$off[grid[g, 1]]
        wt <- st[[1]]$w[grid[g, 1]]
        q0 <- 1L
      }
      for (q in seq_along(dirs)) {
        oth[dirs[q]] <- st[[q0 + q]]$off[grid[g, q0 + q]]
        wt <- wt * st[[q0 + q]]$w[grid[g, q0 + q]]
      }
      acc <- acc + wt * value_at(oy, oth)[[comp]]
    }
    den <- if (ry > 0L) hy^ry else 1
    for (j in dirs) den <- den * hth[[j]]^rem[j]
    acc / den
  }

  if (ctheta == 1L) {
    out <- lapply(seq_len(p), function(j) one(j))
    names(out) <- params
  } else {
    out <- lapply(prs, function(pr) one(match(pr, params)))
    names(out) <- names(prs)
  }
  out
}
