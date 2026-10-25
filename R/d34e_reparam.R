#' @include reparametrize.R reparam_maps.R dexpected_reparam.R dexpected_families.R d34expected_analytic.R d34expected_hessian.R
NULL

#' Does a Reparametrization Carry the Map's Closed Derivatives?
#'
#' @description
#' Returns `TRUE` when the map's partial derivatives of a reparametrized family
#' were supplied through `map_derivs`, and `FALSE` when they come from the
#' finite-difference tables of [reparam_stencil_derivs()].
#'
#' @details
#' The numerical tables lose accuracy quickly with the order, so a chain of the
#' third or fourth order that reads them is no better than the stencil it would
#' replace. The two routes are told apart by the body of the `reparam_derivs`
#' function, which is the one of [reparam_stencil_derivs()] in the numerical
#' case.
#'
#' @param distrib A reparametrized family.
#'
#' @return A single logical.
#'
#' @keywords internal
reparam_has_map_derivs <- function(distrib) {
  ref <- reparam_stencil_derivs(function(x) x, "a", "b")
  !identical(body(distrib@reparam_derivs), body(ref))
}

#' Index Tables for the Chain Through a Reparametrization
#'
#' @description
#' Rewrites the keyed partial tables of a map as lists indexed by an integer
#' code of the sorted tuple of new coordinates, so that the chain formulas read
#' a partial without composing a key.
#'
#' @param maps One keyed table per parent parameter, as [reparam_tables()]
#'   returns them.
#' @param p,q The number of parent parameters and of new ones.
#' @param rmax The highest order of partial in `maps`.
#'
#' @return A list with `enc`, the code of a tuple (the empty tuple is 1),
#'   `tup_mat`, the sorted tuples of one length as the rows of a matrix, and
#'   `Jc`, the partial of parent parameter `i` along a tuple of new
#'   coordinates, `0` where the table has no such key.
#'
#' @keywords internal
reparam_map_index <- function(maps, p, q, rmax) {
  base <- q + 1
  enc <- function(idx) {
    if (!length(idx)) return(1)
    1 + sum(sort(idx) * base^(seq_along(idx) - 1L))
  }
  qs <- seq_len(q)
  tup_mat <- function(r) {
    g <- as.matrix(expand.grid(rep(list(qs), r)))
    if (r > 1L) g <- t(apply(g, 1L, sort))
    unique(g)
  }
  Jt <- lapply(seq_len(p), function(i) {
    tab <- vector("list", base^rmax)
    for (r in seq_len(rmax)) {
      g <- tup_mat(r)
      for (s in seq_len(nrow(g))) {
        v <- maps[[i]][[paste(g[s, ], collapse = ",")]]
        tab[[enc(g[s, ])]] <- if (is.null(v)) 0 else v
      }
    }
    tab
  })
  list(enc = enc, tup_mat = tup_mat, base = base,
       Jc = function(i, idx) Jt[[i]][[enc(idx)]])
}

#' The Parent's Derivatives of the Expected Information as Nested Lists
#'
#' @description
#' Reads the keyed derivatives of the parent's expected information into lists
#' indexed by the parent's parameter positions, so that the chain formulas
#' address \eqn{E_{ij,k\ldots}} by integers.
#'
#' @param Pp The parent's parameter names.
#' @param src A named list keyed as [hess_names()], [dexpected_names()],
#'   [d2expected_names()], [d3expected_names()] or [d4expected_names()].
#' @param r The number of derivative indices, 0 to 4.
#'
#' @return A list nested `2 + r` deep over `seq_along(Pp)`.
#'
#' @keywords internal
parent_expected_nested <- function(Pp, src, r) {
  seqp <- seq_along(Pp)
  key <- switch(r + 1L,
    function(i, j) hess_pair_name(Pp, i, j),
    function(i, j, k) dexpected_key(Pp, i, j, k),
    function(i, j, k, l) d2expected_key(Pp, i, j, k, l),
    function(i, j, k, l, m) d3expected_key(Pp, i, j, k, l, m),
    function(i, j, k, l, m, n) d4expected_key(Pp, i, j, k, l, m, n))
  rec <- function(idx) {
    if (length(idx) == r + 2L) return(src[[do.call(key, as.list(idx))]])
    lapply(seqp, function(v) rec(c(idx, v)))
  }
  rec(integer(0))
}

#' Assemble a Derivative of the Expected Information From the Factors
#'
#' @description
#' Applies Leibniz's rule to \eqn{\tilde E_{ab} = \sum_{ij} G_{ij} J_{ia}
#' J_{jb}} along every sorted tuple of new coordinates of the given length: the
#' positions of the tuple are assigned in all \eqn{3^r} ways to the factors
#' \eqn{G_{ij}}, \eqn{J_{ia}} and \eqn{J_{jb}}.
#'
#' @param Gtab Nested list over the parent pairs, each holding the derivatives
#'   of \eqn{G_{ij}} indexed by the code of the tuple.
#' @param ix The index tables from [reparam_map_index()].
#' @param Pn The new parameter names.
#' @param p The number of parent parameters.
#' @param r The order, `3L` or `4L`.
#' @param n The number of observations, to which every component is recycled.
#'
#' @return A named list keyed as [d3expected_names()] or [d4expected_names()]
#'   over `Pn`.
#'
#' @keywords internal
chain_leibniz <- function(Gtab, ix, Pn, p, r, n) {
  q <- length(Pn)
  seqp <- seq_len(p)
  Jc <- ix$Jc; enc <- ix$enc
  is0 <- function(x) length(x) == 1L && x == 0
  asg <- as.matrix(expand.grid(rep(list(0:2), r)))
  pairs <- which(upper.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  tup <- ix$tup_mat(r)
  out <- list()
  for (s in seq_len(nrow(pairs))) {
    a <- pairs[s, 1L]; b <- pairs[s, 2L]
    for (u in seq_len(nrow(tup))) {
      S <- tup[u, ]
      acc <- 0
      for (v in seq_len(nrow(asg))) {
        w <- asg[v, ]
        kG <- enc(S[w == 0L])
        iA <- c(a, S[w == 1L]); iB <- c(b, S[w == 2L])
        for (i in seqp) {
          ja <- Jc(i, iA)
          if (is0(ja)) next
          for (j in seqp) {
            jb <- Jc(j, iB)
            if (is0(jb)) next
            acc <- acc + Gtab[[i]][[j]][[kG]] * ja * jb
          }
        }
      }
      key <- if (r == 3L) d3expected_key(Pn, a, b, S[1L], S[2L], S[3L]) else
        d4expected_key(Pn, a, b, S[1L], S[2L], S[3L], S[4L])
      out[[key]] <- acc
    }
  }
  nm <- if (r == 3L) d3expected_names(Pn) else d4expected_names(Pn)
  lapply(out[nm], function(v) rep_len(v, max(n, length(v))))
}

#' The Third Derivative of the Expected Information Through a Reparametrization
#'
#' @description
#' Computes the third derivative of the expected information in the new
#' coordinates \eqn{\phi} from the parent's expected information and its first
#' three derivatives, and the map's partial derivatives up to the fourth order.
#'
#' @details
#' With \eqn{J_{ia} = \partial\theta_i/\partial\phi_a}, the expected
#' information transforms as \eqn{\tilde E_{ab} = \sum_{ij} G_{ij} J_{ia}
#' J_{jb}}, where \eqn{G_{ij}(\phi) = E_{ij}(\theta(\phi))}, and there is no
#' term in the map's second derivative because the score has mean zero. For a
#' multiset \eqn{S = \{c, d, e\}} of new coordinates, Leibniz's rule over the
#' three factors gives
#' \deqn{\partial_S \tilde E_{ab} = \sum_{ij} \sum_{S = S_G \cup S_A \cup S_B}
#'   \partial_{S_G} G_{ij}\; \theta_{i,\,a S_A}\; \theta_{j,\,b S_B},}
#' the sum running over the \eqn{3^3} assignments of the three positions to the
#' factors, with \eqn{\theta_{i,\,a S_A}} the partial derivative of
#' \eqn{\theta_i} in \eqn{\phi_a} and the coordinates of \eqn{S_A}. The
#' derivatives of \eqn{G_{ij}} are those of a composition of the parent's
#' \eqn{E_{ij}} with the map, written out to third order:
#' \deqn{\partial_c G = \sum_k E_{,k}\,\theta_{k,c},}
#' \deqn{\partial_{cd} G = \sum_{kl} E_{,kl}\,\theta_{k,c}\theta_{l,d}
#'   + \sum_k E_{,k}\,\theta_{k,cd},}
#' \deqn{\partial_{cde} G = \sum_{klm} E_{,klm}\,\theta_{k,c}\theta_{l,d}
#'   \theta_{m,e} + \sum_{kl} E_{,kl}\big(\theta_{k,cd}\theta_{l,e}
#'   + \theta_{k,ce}\theta_{l,d} + \theta_{k,de}\theta_{l,c}\big)
#'   + \sum_k E_{,k}\,\theta_{k,cde},}
#' where \eqn{E_{,k}}, \eqn{E_{,kl}} and \eqn{E_{,klm}} are the parent's first,
#' second and third derivatives of \eqn{E_{ij}}. The map is read up to the
#' fourth order, \eqn{\theta_{i,\,a S_A}} with \eqn{|S_A| = 3} being the
#' highest partial.
#'
#' @param Pp,Pn The parent's parameter names and the new ones.
#' @param E,d1,d2,d3 Named lists keyed as [hess_names()], [dexpected_names()],
#'   [d2expected_names()] and [d3expected_names()] on the parent's parameters.
#' @param maps One keyed table per parent parameter, as [reparam_tables()]
#'   returns them, to the order 4.
#' @param n The number of observations, to which every component is recycled.
#'
#' @return A named list keyed as [d3expected_names()] on the new parameters.
#'
#' @seealso [dexpected_chain()] for the orders one and two.
#'
#' @keywords internal
dexpected_chain3 <- function(Pp, Pn, E, d1, d2, d3, maps, n) {
  p <- length(Pp); q <- length(Pn); seqp <- seq_len(p)
  ix <- reparam_map_index(maps, p, q, 4L)
  Jc <- ix$Jc; enc <- ix$enc
  Epl <- parent_expected_nested(Pp, E, 0L)
  dEl <- parent_expected_nested(Pp, d1, 1L)
  d2El <- parent_expected_nested(Pp, d2, 2L)
  d3El <- parent_expected_nested(Pp, d3, 3L)
  Gtab <- lapply(seqp, function(i) vector("list", p))
  t1 <- ix$tup_mat(1L); t2 <- ix$tup_mat(2L); t3 <- ix$tup_mat(3L)
  for (i in seqp) for (j in i:p) {
    tab <- vector("list", ix$base^3)
    tab[[1L]] <- Epl[[i]][[j]]
    for (u in seq_len(nrow(t1))) {
      S <- t1[u, ]
      acc <- 0
      for (k in seqp) acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S[1L])
      tab[[enc(S)]] <- acc
    }
    for (u in seq_len(nrow(t2))) {
      S <- t2[u, ]
      acc <- 0
      for (k in seqp) {
        acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S)
        for (l in seqp) {
          acc <- acc + d2El[[i]][[j]][[k]][[l]] * Jc(k, S[1L]) * Jc(l, S[2L])
        }
      }
      tab[[enc(S)]] <- acc
    }
    for (u in seq_len(nrow(t3))) {
      S <- t3[u, ]
      acc <- 0
      for (k in seqp) {
        acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S)
        for (l in seqp) {
          acc <- acc + d2El[[i]][[j]][[k]][[l]] *
            (Jc(k, S[c(1L, 2L)]) * Jc(l, S[3L]) +
               Jc(k, S[c(1L, 3L)]) * Jc(l, S[2L]) +
               Jc(k, S[c(2L, 3L)]) * Jc(l, S[1L]))
          for (m in seqp) {
            acc <- acc + d3El[[i]][[j]][[k]][[l]][[m]] *
              Jc(k, S[1L]) * Jc(l, S[2L]) * Jc(m, S[3L])
          }
        }
      }
      tab[[enc(S)]] <- acc
    }
    Gtab[[i]][[j]] <- tab
    Gtab[[j]][[i]] <- tab
  }
  chain_leibniz(Gtab, ix, Pn, p, 3L, n)
}

#' The Fourth Derivative of the Expected Information Through a Reparametrization
#'
#' @description
#' Computes the fourth derivative of the expected information in the new
#' coordinates \eqn{\phi} from the parent's expected information and its first
#' four derivatives, and the map's partial derivatives up to the fifth order.
#'
#' @details
#' The construction is that of [dexpected_chain3()] with four positions. The
#' assignments of the positions to \eqn{G_{ij}}, \eqn{J_{ia}} and \eqn{J_{jb}}
#' number \eqn{3^4}, and the derivative of \eqn{G_{ij}} along
#' \eqn{S = \{c, d, e, f\}} is the sum over the fifteen set partitions of
#' \eqn{S}, written by block sizes:
#' \deqn{\partial_{S} G = \sum_k E_{,k}\theta_{k,S}
#'   + \sum_{kl} E_{,kl}\Big(\sum_{x \in S} \theta_{k,S\setminus x}\theta_{l,x}
#'   + \theta_{k,cd}\theta_{l,ef} + \theta_{k,ce}\theta_{l,df}
#'   + \theta_{k,cf}\theta_{l,de}\Big)}
#' \deqn{\quad + \sum_{klm} E_{,klm} \sum_{\{x,y\} \subset S}
#'   \theta_{k,xy}\,\theta_{l,z}\,\theta_{m,w}
#'   + \sum_{klmn} E_{,klmn}\,\theta_{k,c}\theta_{l,d}\theta_{m,e}\theta_{n,f},}
#' with \eqn{z, w} the two positions outside the pair \eqn{\{x, y\}}. The
#' factor \eqn{J_{ia}} is differentiated up to four times, so the map is read
#' to the fifth order.
#'
#' @param Pp,Pn The parent's parameter names and the new ones.
#' @param E,d1,d2,d3,d4 Named lists keyed as [hess_names()],
#'   [dexpected_names()], [d2expected_names()], [d3expected_names()] and
#'   [d4expected_names()] on the parent's parameters.
#' @param maps One keyed table per parent parameter to the order 5.
#' @param n The number of observations, to which every component is recycled.
#'
#' @return A named list keyed as [d4expected_names()] on the new parameters.
#'
#' @keywords internal
dexpected_chain4 <- function(Pp, Pn, E, d1, d2, d3, d4, maps, n) {
  p <- length(Pp); q <- length(Pn); seqp <- seq_len(p)
  ix <- reparam_map_index(maps, p, q, 5L)
  Jc <- ix$Jc; enc <- ix$enc
  is0 <- function(x) length(x) == 1L && x == 0
  Epl <- parent_expected_nested(Pp, E, 0L)
  dEl <- parent_expected_nested(Pp, d1, 1L)
  d2El <- parent_expected_nested(Pp, d2, 2L)
  d3El <- parent_expected_nested(Pp, d3, 3L)
  d4El <- parent_expected_nested(Pp, d4, 4L)
  Gtab <- lapply(seqp, function(i) vector("list", p))
  tabs <- lapply(1:4, ix$tup_mat)
  pair_pos <- list(c(1L, 2L), c(1L, 3L), c(1L, 4L), c(2L, 3L), c(2L, 4L),
                   c(3L, 4L))
  for (i in seqp) for (j in i:p) {
    tab <- vector("list", ix$base^4)
    tab[[1L]] <- Epl[[i]][[j]]
    for (u in seq_len(nrow(tabs[[1L]]))) {
      S <- tabs[[1L]][u, ]
      acc <- 0
      for (k in seqp) acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S[1L])
      tab[[enc(S)]] <- acc
    }
    for (u in seq_len(nrow(tabs[[2L]]))) {
      S <- tabs[[2L]][u, ]
      acc <- 0
      for (k in seqp) {
        acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S)
        for (l in seqp) {
          acc <- acc + d2El[[i]][[j]][[k]][[l]] * Jc(k, S[1L]) * Jc(l, S[2L])
        }
      }
      tab[[enc(S)]] <- acc
    }
    for (u in seq_len(nrow(tabs[[3L]]))) {
      S <- tabs[[3L]][u, ]
      acc <- 0
      for (k in seqp) {
        acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S)
        for (l in seqp) {
          acc <- acc + d2El[[i]][[j]][[k]][[l]] *
            (Jc(k, S[c(1L, 2L)]) * Jc(l, S[3L]) +
               Jc(k, S[c(1L, 3L)]) * Jc(l, S[2L]) +
               Jc(k, S[c(2L, 3L)]) * Jc(l, S[1L]))
          for (m in seqp) {
            acc <- acc + d3El[[i]][[j]][[k]][[l]][[m]] *
              Jc(k, S[1L]) * Jc(l, S[2L]) * Jc(m, S[3L])
          }
        }
      }
      tab[[enc(S)]] <- acc
    }
    for (u in seq_len(nrow(tabs[[4L]]))) {
      S <- tabs[[4L]][u, ]
      acc <- 0
      for (k in seqp) {
        acc <- acc + dEl[[i]][[j]][[k]] * Jc(k, S)
        j1k <- Jc(k, S[1L])
        for (l in seqp) {
          t2 <- 0
          for (x in 1:4) t2 <- t2 + Jc(k, S[-x]) * Jc(l, S[x])
          t2 <- t2 + Jc(k, S[c(1L, 2L)]) * Jc(l, S[c(3L, 4L)]) +
            Jc(k, S[c(1L, 3L)]) * Jc(l, S[c(2L, 4L)]) +
            Jc(k, S[c(1L, 4L)]) * Jc(l, S[c(2L, 3L)])
          acc <- acc + d2El[[i]][[j]][[k]][[l]] * t2
          for (m in seqp) {
            t3 <- 0
            for (pp in pair_pos) {
              rest <- S[-pp]
              t3 <- t3 + Jc(k, S[pp]) * Jc(l, rest[1L]) * Jc(m, rest[2L])
            }
            acc <- acc + d3El[[i]][[j]][[k]][[l]][[m]] * t3
            if (is0(j1k)) next
            for (nn in seqp) {
              acc <- acc + d4El[[i]][[j]][[k]][[l]][[m]][[nn]] * j1k *
                Jc(l, S[2L]) * Jc(m, S[3L]) * Jc(nn, S[4L])
            }
          }
        }
      }
      tab[[enc(S)]] <- acc
    }
    Gtab[[i]][[j]] <- tab
    Gtab[[j]][[i]] <- tab
  }
  chain_leibniz(Gtab, ix, Pn, p, 4L, n)
}

#' The Map's Partial Derivatives of the Fifth Order
#'
#' @description
#' Adds to the keyed tables of a reparametrization the partials of the fifth
#' order, each the derivative along its largest coordinate of the partial of
#' the fourth order below it.
#'
#' @details
#' The contract of `map_derivs` stops at the fourth order, and the fourth
#' derivative of the expected information reads the fifth. A partial of the
#' fifth order is the derivative, along one new coordinate \eqn{\phi_c}, of an
#' analytic partial of the fourth order. It is formed by a central difference
#' of the tables at \eqn{\phi_c \pm h} and at \eqn{\pm h/2}, extrapolated as
#' \eqn{(4D(h/2) - D(h))/3}, with \eqn{h} the relative step \eqn{2\cdot 10^{-4}}
#' of the coordinate kept inside its bounds by [fd_steps()]. The tables are
#' closed forms evaluated at four shifted points per coordinate, with no
#' evaluation of the parent, and the difference is taken once on an analytic
#' quantity.
#'
#' @param distrib A reparametrized family whose `reparam_derivs` supplies the
#'   fourth order.
#' @param theta A named list of the new parameters.
#' @param maps The tables to the fourth order, from [reparam_tables()].
#'
#' @return `maps` with the partials of the fifth order added under their keys.
#'
#' @keywords internal
reparam_map_order5 <- function(distrib, theta, maps) {
  params <- distrib@params
  q <- length(params); p <- length(maps)
  g <- as.matrix(expand.grid(rep(list(seq_len(q)), 5L)))
  g <- unique(t(apply(g, 1L, sort)))
  for (cc in seq_len(q)) {
    nm <- params[cc]
    h <- fd_steps(theta[[nm]], distrib@params_bounds[[nm]], 2e-4, order = 1L)
    tab_at <- function(s) {
      th <- theta
      th[[nm]] <- theta[[nm]] + s
      reparam_tables(distrib, th, 4L)
    }
    up1 <- tab_at(h); dn1 <- tab_at(-h)
    up2 <- tab_at(h / 2); dn2 <- tab_at(-h / 2)
    rows <- which(g[, 5L] == cc)
    for (r in rows) {
      key4 <- paste(g[r, 1:4], collapse = ",")
      key5 <- paste(g[r, ], collapse = ",")
      for (i in seq_len(p)) {
        val <- function(tb) {
          v <- tb[[i]][[key4]]
          if (is.null(v)) 0 else v
        }
        d1 <- (val(up1) - val(dn1)) / (2 * h)
        d2 <- (val(up2) - val(dn2)) / h
        v5 <- (4 * d2 - d1) / 3
        if (any(v5 != 0)) maps[[i]][[key5]] <- v5
      }
    }
  }
  maps
}

#' The Third Derivative of the Expected Information of a Reparametrized Family
#'
#' @description
#' Reads the parent's expected information and its first three derivatives at
#' the mapped parameters and the map's tables to the fourth order, and passes
#' them to [dexpected_chain3()].
#'
#' @param distrib A reparametrized distribution.
#' @param y The response, read for its length.
#' @param theta A named list of the new parameters.
#' @param threads The thread count passed to the parent's methods.
#'
#' @return A named list on the parameter scale, keyed as [d3expected_names()].
#'
#' @keywords internal
reparam_d3expected <- function(distrib, y, theta, threads = 1L) {
  parent <- distrib@parent_distrib
  th_par <- reparam_theta(distrib, theta)
  maps <- reparam_tables(distrib, theta, 4L)
  E <- distrib_expected_hessian(parent, y, th_par, scale = "parameter",
                                threads = threads)
  d1 <- distrib_dexpected_hessian(parent, y, th_par, scale = "parameter",
                                  threads = threads)
  d2 <- distrib_d2expected_hessian(parent, y, th_par, scale = "parameter",
                                   threads = threads)
  d3 <- distrib_d3expected_hessian(parent, y, th_par, scale = "parameter",
                                   threads = threads)
  dexpected_chain3(parent@params, distrib@params, E, d1, d2, d3, maps,
                   length(y))
}

#' The Fourth Derivative of the Expected Information of a Reparametrized Family
#'
#' @description
#' Reads the parent's expected information and its first four derivatives at
#' the mapped parameters and the map's tables to the fifth order, and passes
#' them to [dexpected_chain4()].
#'
#' @details
#' The tables of the fifth order are formed by [reparam_map_order5()] from the
#' fourth-order tables, which `map_derivs` supplies in closed form.
#'
#' @param distrib A reparametrized distribution.
#' @param y The response, read for its length.
#' @param theta A named list of the new parameters.
#' @param threads The thread count passed to the parent's methods.
#'
#' @return A named list on the parameter scale, keyed as [d4expected_names()].
#'
#' @keywords internal
reparam_d4expected <- function(distrib, y, theta, threads = 1L) {
  parent <- distrib@parent_distrib
  theta <- align_theta(distrib, theta)[distrib@params]
  th_par <- reparam_theta(distrib, theta)
  maps <- reparam_map_order5(distrib, theta, reparam_tables(distrib, theta, 4L))
  E <- distrib_expected_hessian(parent, y, th_par, scale = "parameter",
                                threads = threads)
  d1 <- distrib_dexpected_hessian(parent, y, th_par, scale = "parameter",
                                  threads = threads)
  d2 <- distrib_d2expected_hessian(parent, y, th_par, scale = "parameter",
                                   threads = threads)
  d3 <- distrib_d3expected_hessian(parent, y, th_par, scale = "parameter",
                                   threads = threads)
  d4 <- distrib_d4expected_hessian(parent, y, th_par, scale = "parameter",
                                   threads = threads)
  dexpected_chain4(parent@params, distrib@params, E, d1, d2, d3, d4, maps,
                   length(y))
}

#' The Third and Fourth Derivatives of a Reparametrized Family's Expected Information
#'
#' @description
#' Returns the third (`k = 3`) or fourth (`k = 4`) derivatives of the expected
#' information of a reparametrized family on the parameter scale, from the
#' parent's closed derivatives and the map's tables. The function is registered
#' through [register_d34expected()], which carries the derivatives to the link
#' scale with the wrapper's links.
#'
#' @details
#' The third derivative is the chain of [dexpected_chain3()] and the fourth that
#' of [dexpected_chain4()]. They are exact wherever the parent's derivatives to
#' that order are, and the map's partials to the fourth order are, the
#' partials of the fifth order being formed from the fourth by
#' [reparam_map_order5()]. A family built without `map_derivs` has numerical
#' tables that lose their leading digits from the third order, and both orders
#' then come from [numerical_dexpected_higher()], the single stencil on the
#' second derivative.
#'
#' @param d A reparametrized family, from [reparametrize()].
#' @param y The response, read for its length.
#' @param th The new parameters.
#' @param k The order, `3L` or `4L`.
#' @param t The thread count passed to the parent's methods.
#'
#' @return A named list keyed as [d3expected_names()] (`k = 3`) or
#'   [d4expected_names()] (`k = 4`).
#'
#' @keywords internal
reparam_d34expected <- function(d, y, th, k, t) {
  if (!reparam_has_map_derivs(d)) {
    return(numerical_dexpected_higher(d, y, th, "parameter", k))
  }
  if (k == 3L) reparam_d3expected(d, y, th, t) else
    reparam_d4expected(d, y, th, t)
}
register_d34expected(ReparamContinuousDistrib, reparam_d34expected)
register_d34expected(ReparamDiscreteDistrib, reparam_d34expected)
