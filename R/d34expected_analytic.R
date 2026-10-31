#' @include d34expected_hessian.R dexpected_hessian.R link_scale.R
NULL

# The third and fourth derivatives of the expected information from a family's
# own compiled kernels. A kernel returns the parameter-scale components of ONE
# order, in the order of d3expected_names() or d4expected_names(); the link
# scale is carried across once, here, from the orders below it.

#' Register a Family's Third and Fourth Derivatives of the Expected Information
#'
#' @description
#' Registers methods of [distrib_d3expected_hessian()] and
#' [distrib_d4expected_hessian()] for one class, both reading one kernel on
#' the parameter scale through [d34expected_analytic()].
#'
#' @param cls An S7 class.
#' @param kernel A function of `(distrib, y, theta, order, threads)` returning,
#'   for `order` 3 or 4, the parameter-scale components in the order of
#'   [d3expected_names()] or [d4expected_names()]; names are set by the caller.
#'
#' @return `NULL`, invisibly; called for the registration.
#'
#' @keywords internal
register_d34expected <- function(cls, kernel) {
  S7::method(distrib_d3expected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      method = NULL, ...,
      threads = 1L) {
    d34expected_analytic(distrib, y, theta, match.arg(scale), 3L, threads,
                         function(k) kernel(distrib, y, theta, k, threads))
  }
  S7::method(distrib_d4expected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      method = NULL, ...,
      threads = 1L) {
    d34expected_analytic(distrib, y, theta, match.arg(scale), 4L, threads,
                         function(k) kernel(distrib, y, theta, k, threads))
  }
  invisible(NULL)
}


#' The Analytic Third and Fourth Derivatives of the Expected Information
#'
#' @description
#' Reads a family's kernel at order 3 or 4 on the parameter scale and carries
#' the result onto the link scale where it is requested, from the expected
#' information and its derivatives of every lower order.
#'
#' @param distrib A distribution object.
#' @param y,theta As the generic takes them.
#' @param scale `"parameter"` or `"link"`.
#' @param order `3L` or `4L`.
#' @param threads The thread count passed to the kernels.
#' @param kern A function of the order returning the kernel's list.
#'
#' @return A named list keyed as [d3expected_names()] or [d4expected_names()].
#'
#' @seealso [dexpected_link_higher()]
#'
#' @keywords internal
d34expected_analytic <- function(distrib, y, theta, scale, order, threads,
                                 kern) {
  params <- distrib@params
  nm <- if (order == 3L) d3expected_names(params) else d4expected_names(params)
  top <- kern(order)
  if (length(top) != length(nm)) {
    stop(sprintf("The kernel of order %d returned %d components, not %d.",
                 order, length(top), length(nm)), call. = FALSE)
  }
  top <- stats::setNames(unclass(top), nm)
  if (identical(scale, "parameter")) return(top)
  nat <- list(
    distrib_expected_hessian(distrib, y, theta, scale = "parameter",
                             threads = threads),
    distrib_dexpected_hessian(distrib, y, theta, scale = "parameter",
                              threads = threads),
    distrib_d2expected_hessian(distrib, y, theta, scale = "parameter",
                               threads = threads))
  if (order == 4L) {
    nat[[4L]] <- distrib_d3expected_hessian(distrib, y, theta,
                                            scale = "parameter",
                                            threads = threads)
  }
  nat[[order + 1L]] <- top
  h <- inverse_link_derivs(distrib, theta, order + 1L)
  dexpected_link_higher(params, nat, h, order)
}


#' The Derivatives of the Link-Scale Expected Information, Any Order
#'
#' @description
#' Applies Leibniz's rule and the diagonal Faa di Bruno formula to
#' \eqn{F_{ab}(\eta) = E_{ab}(\theta(\eta))\,h_a'(\eta_a)\,h_b'(\eta_b)}, at
#' orders one to five.
#'
#' @details
#' For a multiset \eqn{S} of coordinates, with \eqn{G_{ab} = E_{ab}\circ\theta}
#' and \eqn{u_{ab} = h_a' h_b'},
#' \deqn{\partial_S F_{ab} = \sum_{T \subseteq S} \binom{S}{T}\,
#'   \partial_T G_{ab}\;\partial_{S \setminus T} u_{ab},}
#' where the sum runs over the sub-multisets of \eqn{S} and the coefficient is
#' the product over the distinct indices of the binomial coefficients of their
#' multiplicities. Each parameter moves with its own coordinate only, so
#' \deqn{\partial_T G_{ab} = \sum_{j} \partial^{J} E_{ab}
#'   \prod_i B_{m_i, j_i}(h_i', h_i'', \dots),}
#' with \eqn{m_i} the multiplicity of index \eqn{i} in \eqn{T}, \eqn{1 \le j_i
#' \le m_i}, \eqn{J} the multiset holding each \eqn{i} \eqn{j_i} times, and
#' \eqn{B} the partial Bell polynomials of [bell_partial()]. The factor
#' \eqn{\partial_R u_{ab}} vanishes unless \eqn{R} names only \eqn{a} and
#' \eqn{b}; then it is \eqn{h_a^{(1+r_a)} h_b^{(1+r_b)}} for \eqn{a \ne b} and
#' \eqn{\sum_k \binom{r}{k} h_a^{(1+k)} h_a^{(1+r-k)}} for \eqn{a = b}.
#'
#' @param params The parameter names.
#' @param nat A list: `nat[[1]]` the expected information on the parameter
#'   scale, keyed as [hess_names()], and `nat[[m + 1]]` its derivatives of
#'   order `m`, keyed as [dexpected_names()], [d2expected_names()],
#'   [d3expected_names()], [d4expected_names()] and [d5expected_names()].
#' @param h The inverse link's derivatives, from [inverse_link_derivs()], to
#'   order `order + 1`.
#' @param order `1L` to `5L`.
#'
#' @return A named list on the link scale, keyed as the order's names.
#'
#' @keywords internal
dexpected_link_higher <- function(params, nat, h, order) {
  p <- length(params)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  tuples <- if (order == 1L) as.list(seq_len(p)) else
    if (order == 2L) lapply(seq_len(nrow(pairs)), function(r) pairs[r, ]) else
      deriv_indices(params, order)
  # the key of the derivative of E_ab along the sorted multiset T
  ekey <- function(a, b, T) {
    switch(length(T) + 1L,
      hess_pair_name(params, a, b),
      dexpected_key(params, a, b, T[1L]),
      d2expected_key(params, a, b, T[1L], T[2L]),
      d3expected_key(params, a, b, T[1L], T[2L], T[3L]),
      d4expected_key(params, a, b, T[1L], T[2L], T[3L], T[4L]),
      d5expected_key(params, a, b, T[1L], T[2L], T[3L], T[4L], T[5L]))
  }
  hk <- function(i, k) h[[i]][[k]]
  # d_T G_ab, T given by its multiplicities over the indices
  dG <- function(a, b, m) {
    if (sum(m) == 0L) return(nat[[1L]][[ekey(a, b, integer(0))]])
    idx <- which(m > 0L)
    grid <- as.matrix(expand.grid(lapply(m[idx], seq_len)))
    acc <- 0
    for (g in seq_len(nrow(grid))) {
      j <- unname(grid[g, ])
      J <- rep(idx, j)
      coef <- 1
      for (s in seq_along(idx)) {
        coef <- coef * bell_partial(m[idx[s]], j[s], h[[idx[s]]])
      }
      acc <- acc + nat[[length(J) + 1L]][[ekey(a, b, J)]] * coef
    }
    acc
  }
  # d_R u_ab, R given by its multiplicities
  du <- function(a, b, r) {
    if (any(r[-unique(c(a, b))] > 0L)) return(0)
    if (a != b) return(hk(a, 1L + r[a]) * hk(b, 1L + r[b]))
    ra <- r[a]
    acc <- 0
    for (k in 0:ra) acc <- acc + choose(ra, k) * hk(a, 1L + k) * hk(a, 1L + ra - k)
    acc
  }
  nm <- switch(order, dexpected_names(params), d2expected_names(params),
               d3expected_names(params), d4expected_names(params),
               d5expected_names(params))
  out <- stats::setNames(vector("list", length(nm)), nm)
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (S in tuples) {
      mS <- tabulate(S, p)
      sub <- as.matrix(expand.grid(lapply(mS, function(v) 0:v)))
      acc <- 0
      for (s in seq_len(nrow(sub))) {
        mT <- unname(sub[s, ])
        coef <- prod(choose(mS, mT))
        uR <- du(a, b, mS - mT)
        if (identical(uR, 0)) next
        acc <- acc + coef * dG(a, b, mT) * uR
      }
      out[[ekey(a, b, sort(S))]] <- acc
    }
  }
  out
}
