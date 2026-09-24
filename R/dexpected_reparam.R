#' @include dexpected_hessian.R reparametrize.R dexpected_families.R
#' @include student_t1_distrib.R gengamma1_distrib.R
NULL

#' @title Derivatives of the Expected Information, Student t and Generalized Gamma
#' @name distrib_dexpected_hessian.StudentT1Distrib
#' @aliases distrib_d2expected_hessian.StudentT1Distrib
#'   distrib_dexpected_hessian.GenGamma1Distrib
#'   distrib_d2expected_hessian.GenGamma1Distrib
#' @description
#' The first and second derivatives of the expected information in the
#' parameters, from compiled kernels.
#'
#' @details
#' For the Student t, with \eqn{w_1 = 1 + 1/\nu} and \eqn{w_3 = 1 + 3/\nu},
#' \eqn{\mathbb{E}[\ell_{\mu\mu}] = -(w_1/w_3)/\sigma^2},
#' \eqn{\mathbb{E}[\ell_{\sigma\sigma}] = -(2/w_3)/\sigma^2},
#' \eqn{\mathbb{E}[\ell_{\sigma\nu}] = 2\nu^{-2}/(w_1 w_3 \sigma)} and
#' \deqn{\mathbb{E}[\ell_{\nu\nu}] = \tfrac14\{\psi'((\nu+1)/2) - \psi'(\nu/2)\}
#'   + \frac{\nu+5}{2\nu(\nu+1)(\nu+3)} = -\frac{7}{2\nu^4} + \frac{13}{\nu^5} - \dots,}
#' whose terms cancel from order \eqn{\nu^{-2}}; its derivatives in \eqn{\nu}
#' are taken from the asymptotic series above \eqn{\nu = 30}, the series'
#' coefficients being exact integers from the duplication identity
#' \eqn{\psi'((\nu+1)/2) - \psi'(\nu/2) = 4\psi'(\nu) - 2\psi'(\nu/2)}.
#'
#' For the generalized gamma by scale \eqn{a} and shapes \eqn{d, p}, each
#' component is \eqn{a^{\alpha} p^{\beta} F(k)} with \eqn{k = d/p} and
#' \eqn{F} a combination of \eqn{\psi(k)}, \eqn{\psi(k+1)} and their
#' derivatives, differentiated through \eqn{\partial k/\partial d = 1/p} and
#' \eqn{\partial k/\partial p = -k/p}.
#'
#' @param distrib A distribution object of one of the classes above.
#' @param y A numeric vector of observations, read for its length.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads A single positive integer, the kernel's thread count.
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#' @seealso [distrib_dexpected_hessian()], [distrib_d2expected_hessian()]
#' @keywords internal
NULL

register_dexpected(StudentT1Distrib, function(d, y, th, k, t)
  student_t1_dexpected_cpp(y, th[[1]], th[[2]], th[[3]], k, t))
register_dexpected(GenGamma1Distrib, function(d, y, th, k, t)
  gengamma1_dexpected_cpp(y, th[[1]], th[[2]], th[[3]], k, t))


#' The Expected Information's Derivatives Through a Reparametrization
#'
#' @description
#' Carries the parent's \eqn{\partial E/\partial\theta} and
#' \eqn{\partial^2 E/\partial\theta^2} into the new coordinates \eqn{\phi}.
#'
#' @details
#' The expected information transforms as a tensor, with no term in the map's
#' second derivative because the score has mean zero:
#' \eqn{\tilde E_{ab} = \sum_{ij} E_{ij} J_{ia} J_{jb}}, with
#' \eqn{J_{ia} = \partial\theta_i/\partial\phi_a}. Writing \eqn{J_{ia,c}} and
#' \eqn{J_{ia,cd}} for the map's second and third derivatives,
#' \deqn{\partial_c\tilde E_{ab} = \sum_{ij}\Big[\sum_k \partial_k E_{ij}
#'   J_{kc}\Big] J_{ia}J_{jb} + \sum_{ij} E_{ij}(J_{ia,c}J_{jb} + J_{ia}J_{jb,c}),}
#' and \eqn{\partial_{cd}\tilde E_{ab}} follows by Leibniz's rule once more,
#' reading the parent's second derivative and the map's third. The map's
#' partials are the keyed tables of [reparam_tables()], hand-written for the
#' shipped reparametrizations, so every term is analytic wherever the
#' parent's derivatives are.
#'
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param order `1L` or `2L`.
#' @param threads The thread count passed to the parent.
#'
#' @return A named list on the parameter scale, keyed as [dexpected_names()]
#'   or [d2expected_names()].
#'
#' @seealso [reparam_chain()], [dexpected_link()]
#'
#' @keywords internal
reparam_dexpected <- function(distrib, y, theta, order, threads = 1L) {
  parent <- distrib@parent_distrib
  Pp <- parent@params
  Pn <- distrib@params
  th_par <- reparam_theta(distrib, theta)
  maps <- reparam_tables(distrib, theta)
  E <- distrib_expected_hessian(parent, y, th_par, scale = "parameter",
                                threads = threads)
  d1 <- distrib_dexpected_hessian(parent, y, th_par, scale = "parameter",
                                  threads = threads)
  d2 <- if (order == 2L) {
    distrib_d2expected_hessian(parent, y, th_par, scale = "parameter",
                               threads = threads)
  }
  dexpected_chain(Pp, Pn, E, d1, d2, maps, order, length(y))
}

#' The Chain Rule for the Expected Information's Derivatives
#'
#' @description
#' The arithmetic of [reparam_dexpected()], given the parent's quantities
#' rather than the parent: \eqn{E}, \eqn{\partial E} and, at order 2,
#' \eqn{\partial^2 E} on the parent's scale, and the map's partials as keyed
#' tables. It serves a reparametrized family and any family written as a map
#' of another's coordinates, such as [betabinom1_distrib()] over the shapes.
#'
#' @param Pp,Pn The parent's parameter names and the new ones.
#' @param E,d1,d2 Named lists keyed as [hess_names()], [dexpected_names()]
#'   and [d2expected_names()] on the parent's parameters; `d2` is read only at
#'   order 2.
#' @param maps One keyed table per parent parameter, as [reparam_tables()]
#'   returns them: the key is the sorted, comma-joined indices of the new
#'   parameters differentiated, and a missing key is an exact zero.
#' @param order `1L` or `2L`.
#' @param n The number of observations, to which every component is recycled.
#'
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()]
#'   on the new parameters.
#'
#' @seealso [reparam_dexpected()]
#'
#' @keywords internal
dexpected_chain <- function(Pp, Pn, E, d1, d2, maps, order, n) {
  p <- length(Pp); q <- length(Pn)
  # every lookup is resolved once, into lists indexed by integers: composing
  # the keys inside the loops below was nearly all of this function's cost
  Jv <- function(i, ...) {
    v <- maps[[i]][[paste(sort(c(...)), collapse = ",")]]
    if (is.null(v)) 0 else v
  }
  J1 <- lapply(seq_len(p), function(i) lapply(seq_len(q), function(a) Jv(i, a)))
  J2 <- lapply(seq_len(p), function(i) lapply(seq_len(q), function(a)
    lapply(seq_len(q), function(c) Jv(i, a, c))))
  J3 <- if (order == 2L) {
    lapply(seq_len(p), function(i) lapply(seq_len(q), function(a)
      lapply(seq_len(q), function(c) lapply(seq_len(q), function(d)
        Jv(i, a, c, d)))))
  }
  Epl <- lapply(seq_len(p), function(i) lapply(seq_len(p), function(j)
    E[[hess_pair_name(Pp, i, j)]]))
  dEl <- lapply(seq_len(p), function(i) lapply(seq_len(p), function(j)
    lapply(seq_len(p), function(k) d1[[dexpected_key(Pp, i, j, k)]])))
  d2El <- if (order == 2L) {
    lapply(seq_len(p), function(i) lapply(seq_len(p), function(j)
      lapply(seq_len(p), function(k) lapply(seq_len(p), function(l)
        d2[[d2expected_key(Pp, i, j, k, l)]]))))
  }
  # the parent's first derivative contracted along phi_c: sum_k dE_ijk J_kc
  dEcl <- lapply(seq_len(p), function(i) lapply(seq_len(p), function(j)
    lapply(seq_len(q), function(c) Reduce(`+`, lapply(seq_len(p), function(k)
      dEl[[i]][[j]][[k]] * J1[[k]][[c]])))))
  pairs <- which(upper.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L] != pairs[, 2L]), , drop = FALSE]
  out <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    if (order == 1L) {
      for (c in seq_len(q)) {
        acc <- 0
        for (i in seq_len(p)) for (j in seq_len(p)) {
          acc <- acc + dEcl[[i]][[j]][[c]] * J1[[i]][[a]] * J1[[j]][[b]] +
            Epl[[i]][[j]] * (J2[[i]][[a]][[c]] * J1[[j]][[b]] +
                               J1[[i]][[a]] * J2[[j]][[b]][[c]])
        }
        out[[dexpected_key(Pn, a, b, c)]] <- acc
      }
      next
    }
    for (s in seq_len(nrow(pairs))) {
      c <- pairs[s, 1L]; d <- pairs[s, 2L]
      acc <- 0
      for (i in seq_len(p)) for (j in seq_len(p)) {
        u <- J1[[i]][[a]] * J1[[j]][[b]]
        # the parent's second derivative along (c, d) and its first along the
        # map's second derivative
        t2 <- 0
        for (k in seq_len(p)) {
          t2 <- t2 + dEl[[i]][[j]][[k]] * J2[[k]][[c]][[d]]
          for (l in seq_len(p)) {
            t2 <- t2 + d2El[[i]][[j]][[k]][[l]] * J1[[k]][[c]] * J1[[l]][[d]]
          }
        }
        acc <- acc + t2 * u +
          dEcl[[i]][[j]][[c]] * (J2[[i]][[a]][[d]] * J1[[j]][[b]] +
                                   J1[[i]][[a]] * J2[[j]][[b]][[d]]) +
          dEcl[[i]][[j]][[d]] * (J2[[i]][[a]][[c]] * J1[[j]][[b]] +
                                   J1[[i]][[a]] * J2[[j]][[b]][[c]]) +
          Epl[[i]][[j]] * (J3[[i]][[a]][[c]][[d]] * J1[[j]][[b]] +
                             J2[[i]][[a]][[c]] * J2[[j]][[b]][[d]] +
                             J2[[i]][[a]][[d]] * J2[[j]][[b]][[c]] +
                             J1[[i]][[a]] * J3[[j]][[b]][[c]][[d]])
      }
      out[[d2expected_key(Pn, a, b, c, d)]] <- acc
    }
  }
  lapply(out, function(v) rep_len(v, max(n, length(v))))
}

#' A Method Body for the Reparametrized Families
#'
#' @description
#' Returns the method of [distrib_dexpected_hessian()] (order 1) or
#' [distrib_d2expected_hessian()] (order 2) registered on the reparametrized
#' classes, reading [reparam_dexpected()].
#'
#' @param order `1L` or `2L`.
#'
#' @return A function with the generics' signature.
#'
#' @keywords internal
reparam_dexpected_method <- function(order) {
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
           threads = 1L) {
    scale <- match.arg(scale)
    dexpected_analytic(distrib, y, theta, scale, order, threads,
                       function(k) reparam_dexpected(distrib, y, theta, k, threads))
  }
}

#' @title Derivatives of the Expected Information, Reparametrized Families
#' @name distrib_dexpected_hessian.ReparamContinuousDistrib
#' @aliases distrib_d2expected_hessian.ReparamContinuousDistrib
#'   distrib_dexpected_hessian.ReparamDiscreteDistrib
#'   distrib_d2expected_hessian.ReparamDiscreteDistrib
#' @description
#' The parent's derivatives carried through the map by
#' [reparam_dexpected()]. A parent without an analytic second derivative
#' leaves the reparametrized family without one, the error naming the parent.
#' @param distrib A reparametrized distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the new parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads Passed to the parent.
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#' @keywords internal
NULL

S7::method(distrib_dexpected_hessian, ReparamContinuousDistrib) <- reparam_dexpected_method(1L)
S7::method(distrib_d2expected_hessian, ReparamContinuousDistrib) <- reparam_dexpected_method(2L)
S7::method(distrib_dexpected_hessian, ReparamDiscreteDistrib) <- reparam_dexpected_method(1L)
S7::method(distrib_d2expected_hessian, ReparamDiscreteDistrib) <- reparam_dexpected_method(2L)
