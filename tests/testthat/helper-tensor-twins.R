# The R implementations of the tensor stencils as they stood in
# distributions7 0.78.1, before they became adapters on
# numericals7::fd_tensor(). They are kept as twins: at accuracy 2 the
# adapters evaluate the same nodes with the same steps away from a bound, so
# the two must agree to rounding (test-tensor-twins.R). twin_numerical_deriv5()
# takes the fourth order it differences as an argument, the generic having
# moved to accuracy 4.

twin_tensor_derivatives <- function(distrib, y, theta, order, base, rel_step = NULL,
                               skip = NULL) {
  params <- distrib@params
  p <- length(params)
  nms <- deriv_names(params, order)
  idx_of <- deriv_indices(params, order)
  r_tot <- order - base
  if (r_tot < 1L) stop("'order' must exceed 'base'.", call. = FALSE)
  if (is.null(rel_step)) rel_step <- .Machine$double.eps^(1 / (r_tot + 2))

  # one step per parameter, clamped for the widest stencil of this order
  reach <- numericals7::fd_offsets(r_tot, 2L)$half_width
  h <- lapply(seq_len(p), function(j) {
    x <- theta[[params[j]]]
    hj <- rel_step * pmax(1, abs(x))
    b <- distrib@params_bounds[[params[j]]]
    if (!is.null(b)) {
      if (is.finite(b[1])) hj <- pmin(hj, 0.49 * (x - b[1]) / reach)
      if (is.finite(b[2])) hj <- pmin(hj, 0.49 * (b[2] - x) / reach)
    }
    hj
  })

  gen <- switch(as.character(base), "1" = distrib_gradient,
                "2" = distrib_hessian, "3" = distrib_deriv3, NULL)
  cache <- list()
  value_at <- function(off) {
    key <- paste(off, collapse = ",")
    v <- cache[[key]]
    if (is.null(v)) {
      th <- theta
      for (j in which(off != 0L)) {
        th[[params[j]]] <- theta[[params[j]]] + off[j] * h[[j]]
      }
      v <- if (base == 0L) list(distrib_pdf(distrib, y, th, log = TRUE)) else
        gen(distrib, y, th)
      cache[[key]] <<- v
    }
    v
  }

  out <- vector("list", length(nms))
  names(out) <- nms
  for (t in seq_along(nms)) {
    if (nms[t] %in% skip) next
    cnt <- tabulate(idx_of[[t]], p)
    b <- integer(p)
    for (s in seq_len(base)) {
      j <- which.max(cnt - b)
      b[j] <- b[j] + 1L
    }
    rem <- cnt - b
    comp <- if (base == 0L) 1L else
      paste(params[rep(seq_len(p), b)], collapse = "_")
    dirs <- which(rem > 0L)
    st <- lapply(dirs, function(j) {
      off <- numericals7::fd_offsets(rem[j], 2L)$central
      w <- numericals7::fd_weights(off, rem[j])
      keep <- w != 0
      list(off = off[keep], w = w[keep])
    })
    grid <- as.matrix(expand.grid(lapply(st, function(s) seq_along(s$off))))
    acc <- 0
    for (g in seq_len(nrow(grid))) {
      off <- integer(p)
      wt <- 1
      for (q in seq_along(dirs)) {
        off[dirs[q]] <- st[[q]]$off[grid[g, q]]
        wt <- wt * st[[q]]$w[grid[g, q]]
      }
      acc <- acc + wt * value_at(off)[[comp]]
    }
    den <- 1
    for (j in dirs) den <- den * h[[j]]^rem[j]
    out[[nms[t]]] <- acc / den
  }
  out
}


twin_mixed_tensor_derivatives <- function(distrib, y, theta, ay, ctheta) {
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
  rel_step <- .Machine$double.eps^(1 / (r_tot + 2))

  ry <- ay - base$a
  hy <- if (ry > 0L) fd_steps_y(y, distrib@bounds, rel_step, "clamp", order = ry)
  hth <- lapply(seq_len(p), function(j) {
    x <- theta[[params[j]]]
    hj <- rel_step * pmax(1, abs(x))
    reach <- numericals7::fd_offsets(max(1L, ctheta - base$c), 2L)$half_width
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


twin_numerical_deriv5 <- function(distrib, y, theta, scale = c("parameter", "link"),
                             accuracy = 2L, deriv4 = NULL) {
  scale <- match.arg(scale)
  params <- distrib@params
  p <- length(params)
  nms <- deriv_names(params, 5)
  idx_of <- deriv_indices(params, 5)

  # The stencil is numericals7's, nodes and weights together. A zero weight is
  # dropped rather than evaluated: at the default accuracy the centre node
  # carries one, so the rule costs two calls per parameter.
  s <- numericals7::fd_offsets(1L, accuracy = accuracy)$central
  w <- numericals7::fd_weights(s, 1L)
  nodes <- which(w != 0)

  # The point the difference is taken at, and the map back to a theta the
  # fourth-order method can be called with. On the link scale the two differ;
  # on the parameter scale the map is the identity.
  if (scale == "link") {
    links <- distrib@link_params
    x0 <- lapply(params, function(nm)
      linkfunctions7::linkfun(links[[nm]], theta[[nm]]))
    bnds <- lapply(params, function(nm)
      linkfunctions7::eta_bounds(links[[nm]]))
    to_theta <- function(x) {
      th <- theta
      for (i in seq_len(p)) {
        th[[params[i]]] <- linkfunctions7::linkinv(links[[params[i]]], x[[i]])
      }
      th
    }
  } else {
    x0 <- lapply(params, function(nm) theta[[nm]])
    bnds <- lapply(params, function(nm) distrib@params_bounds[[nm]])
    to_theta <- function(x) {
      th <- theta
      for (i in seq_len(p)) th[[params[i]]] <- x[[i]]
      th
    }
  }

  # deriv4 replaces the fourth order differenced, so that the route of
  # 0.78.1 can be rebuilt on its own fallback
  D4 <- if (is.null(deriv4)) {
    function(th) distrib_deriv4(distrib, y, th, scale = scale)
  } else deriv4

  out <- vector("list", length(nms))
  names(out) <- nms

  # The component of a non-decreasing five-tuple is the four-tuple left when
  # its last index is dropped, differentiated in the parameter that index
  # names. Grouping by that index is what makes the whole order cost one
  # stencil per parameter rather than one per component.
  last <- vapply(idx_of, function(r) r[5L], integer(1))
  key4 <- vapply(idx_of, function(r) paste(params[r[1:4]], collapse = "_"),
                 character(1))

  for (m in seq_len(p)) {
    who <- which(last == m)
    if (!length(who)) next

    h <- numericals7::fd_step(x0[[m]], 1L, accuracy = accuracy,
                              bounds = bnds[[m]])

    acc <- vector("list", length(who))
    for (j in nodes) {
      x <- x0
      x[[m]] <- x0[[m]] + s[j] * h
      D <- D4(to_theta(x))
      for (t in seq_along(who)) {
        term <- w[j] * D[[key4[who[t]]]]
        acc[[t]] <- if (is.null(acc[[t]])) term else acc[[t]] + term
      }
    }

    for (t in seq_along(who)) out[[nms[who[t]]]] <- acc[[t]] / h
  }

  out
}
