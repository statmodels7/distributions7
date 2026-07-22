#' @include distrib.R utility_functions.R
NULL

#' Derivatives on the Link (Real) Scale
#'
#' @name link_scale_derivatives
#'
#' @description
#' The derivative generics of the package (\code{\link{distrib_gradient}},
#' \code{\link{distrib_hessian}}, \code{\link{distrib_expected_hessian}},
#' \code{\link{distrib_deriv3}}, \code{\link{distrib_deriv4}}) accept a
#' \code{scale} argument selecting the parameterization the derivatives are taken
#' with respect to:
#' \itemize{
#'   \item \code{scale = "parameter"} (default): derivatives with respect to the
#'     parameters \eqn{\theta} on their natural, possibly constrained scale.
#'   \item \code{scale = "link"}: derivatives with respect to the unconstrained
#'     linear predictors \eqn{\eta_i = g_i(\theta_i)}, where \eqn{g_i} is the link
#'     stored in \code{distrib@link_params}. This is the scale on which
#'     optimization is normally carried out, since \eqn{\eta \in \mathbb{R}^p}.
#' }
#'
#' @details
#' Write \eqn{h_i = g_i^{-1}}, so that \eqn{\theta_i = h_i(\eta_i)}. Because each
#' parameter carries its own link, the Jacobian \eqn{\partial\theta/\partial\eta}
#' is \strong{diagonal} and the multivariate Faa di Bruno formula factorises. For a
#' derivative whose multi-index involves the distinct parameters
#' \eqn{p_1,\dots,p_r} with multiplicities \eqn{m_1,\dots,m_r},
#' \deqn{\frac{\partial^k \ell}{\partial \eta_{p_1}^{m_1} \cdots \partial \eta_{p_r}^{m_r}}
#'       = \sum_{j_1=1}^{m_1} \cdots \sum_{j_r=1}^{m_r}
#'         \ell_{\,p_1^{j_1} \cdots p_r^{j_r}} \prod_{t=1}^{r} B_{m_t, j_t}\!\left(h_{p_t}', h_{p_t}'', \dots\right)}
#' where \eqn{\ell_{\cdot}} are the parameter-scale derivatives and \eqn{B_{m,j}}
#' are the partial (incomplete) Bell polynomials. The low-order cases are
#' \deqn{\frac{\partial \ell}{\partial \eta_a} = \ell_a h', \qquad
#'       \frac{\partial^2 \ell}{\partial \eta_a^2} = \ell_{aa} (h')^2 + \ell_a h''}
#' \deqn{\frac{\partial^3 \ell}{\partial \eta_a^3} = \ell_{aaa}(h')^3 + 3\ell_{aa} h' h'' + \ell_a h'''}
#' \deqn{\frac{\partial^4 \ell}{\partial \eta_a^4} = \ell_{aaaa}(h')^4 + 6\ell_{aaa}(h')^2 h'' + \ell_{aa}\left(4h'h''' + 3(h'')^2\right) + \ell_a h''''}
#' while mixed derivatives across \emph{different} parameters simply multiply the
#' corresponding factors (e.g.
#' \eqn{\partial^2\ell/\partial\eta_a\partial\eta_b = \ell_{ab} h_a' h_b'}).
#'
#' For the \strong{expected} Hessian the first-order term vanishes because the
#' score has zero expectation, so the information transforms as the congruence
#' \eqn{\mathrm{diag}(h')\,\mathbb{E}[H]\,\mathrm{diag}(h')}.
#'
#' The inverse-link derivatives \eqn{h', h'', h''', h''''} are obtained from
#' \code{\link[linkfunctions7]{linkinvderiv}}, so link-scale derivatives are
#' available up to order 4 for every link in \pkg{linkfunctions7}.
#'
#' @seealso \code{\link{distrib_gradient}}, \code{\link{distrib_hessian}},
#'   \code{\link{distrib_deriv3}}, \code{\link{distrib_deriv4}}
NULL

# Derivatives on the link (real) scale.
#
# Each parameter carries its own link, eta_i = g_i(theta_i) and theta_i =
# h_i(eta_i) with h_i = g_i^{-1}. Because the map is applied coordinate-wise the
# Jacobian is diagonal, so the multivariate Faa di Bruno formula factorises: for
# a multi-index with distinct parameters p_1..p_r of multiplicities m_1..m_r,
#
#   d^k L / d eta_{p_1}^{m_1} ... d eta_{p_r}^{m_r}
#     = sum_{j_1=1}^{m_1} ... sum_{j_r=1}^{m_r}
#         l_{p_1^{j_1} ... p_r^{j_r}} * prod_t B_{m_t, j_t}(h_{p_t}', h_{p_t}'', ...)
#
# where l_{...} are the parameter-scale derivatives and B_{m,j} are the partial
# (incomplete) Bell polynomials.

# Partial Bell polynomials B_{m,j} for m <= 4. `h` is a list with h[[k]] the
# k-th derivative of the inverse link evaluated at eta (numeric vectors).
bell_partial <- function(m, j, h) {
  if (m == 1L) {
    return(h[[1]])                                   # B_{1,1}
  }
  if (m == 2L) {
    return(switch(j,
      h[[2]],                                        # B_{2,1}
      h[[1]]^2                                       # B_{2,2}
    ))
  }
  if (m == 3L) {
    return(switch(j,
      h[[3]],                                        # B_{3,1}
      3 * h[[1]] * h[[2]],                           # B_{3,2}
      h[[1]]^3                                       # B_{3,3}
    ))
  }
  if (m == 4L) {
    return(switch(j,
      h[[4]],                                        # B_{4,1}
      4 * h[[1]] * h[[3]] + 3 * h[[2]]^2,            # B_{4,2}
      6 * h[[1]]^2 * h[[2]],                         # B_{4,3}
      h[[1]]^4                                       # B_{4,4}
    ))
  }
  stop("Link-scale derivatives are implemented up to order 4.", call. = FALSE)
}

# Canonical (non-decreasing) index tuples for a derivative order, in exactly the
# same output order as the corresponding name helper: parameter order for
# order 1, hess_names() for order 2, deriv_names() for order >= 3.
deriv_index_list <- function(p, order) {
  if (order == 1L) {
    return(lapply(seq_len(p), function(i) i))
  }
  if (order == 2L) {
    out <- lapply(seq_len(p), function(i) c(i, i))          # diagonal first
    if (p >= 2L) {
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) out[[length(out) + 1L]] <- c(i, j)
      }
    }
    return(out)
  }
  idx <- as.matrix(do.call(expand.grid, rep(list(seq_len(p)), order)))
  idx <- idx[, rev(seq_len(order)), drop = FALSE]
  idx <- idx[apply(idx, 1L, function(r) all(diff(r) >= 0)), , drop = FALSE]
  lapply(seq_len(nrow(idx)), function(k) as.integer(idx[k, ]))
}

# Inverse-link derivatives h^(1..order) for every parameter, evaluated at
# eta = g(theta). Returns a list indexed by parameter position.
inverse_link_derivs <- function(distrib, theta, order) {
  params <- distrib@params
  links <- distrib@link_params   # one property read, not one per parameter

  # The order-specific generics are called directly. linkinvderiv() is a
  # convenience router: it dispatches on the link, then calls dlinkinv() and
  # friends, which dispatch again. That second dispatch is a third of the cost of
  # the call, and this function makes one call per parameter per order -- twelve
  # of them for a three-parameter distribution at fourth order. Every link in
  # linkfunctions7 implements the order-specific generics; they are what
  # linkinvderiv() itself reaches for.
  lapply(seq_along(params), function(i) {
    lk <- links[[params[i]]]
    e <- linkfunctions7::linkfun(lk, theta[[i]])
    switch(as.integer(order),
      list(linkfunctions7::dlinkinv(lk, e)),
      list(linkfunctions7::dlinkinv(lk, e),
           linkfunctions7::d2linkinv(lk, e)),
      list(linkfunctions7::dlinkinv(lk, e),
           linkfunctions7::d2linkinv(lk, e),
           linkfunctions7::d3linkinv(lk, e)),
      list(linkfunctions7::dlinkinv(lk, e),
           linkfunctions7::d2linkinv(lk, e),
           linkfunctions7::d3linkinv(lk, e),
           linkfunctions7::d4linkinv(lk, e))
    )
  })
}

# Parameter-scale derivatives of every order strictly below `order`, in the
# layout expected by to_link_scale(). For expected derivatives the first-order
# slot is zero, since the score has zero expectation.
link_scale_lower_orders <- function(distrib, y, theta, expected, order) {
  params <- distrib@params
  nat <- vector("list", order - 1L)

  nat[[1]] <- if (isTRUE(expected)) {
    stats::setNames(lapply(params, function(p) rep(0, length(y))), params)
  } else {
    distrib_gradient(distrib, y, theta)
  }

  if (order >= 3L) {
    nat[[2]] <- if (isTRUE(expected)) {
      distrib_expected_hessian(distrib, y, theta)
    } else {
      distrib_hessian(distrib, y, theta)
    }
  }
  if (order >= 4L) {
    nat[[3]] <- distrib_deriv3(distrib, y, theta, expected = expected)
  }

  nat
}

# Convert parameter-scale derivatives to link-scale ones.
#
# `nat` is a list with nat[[m]] the named list of parameter-scale derivatives of
# order m, for m = 1..order (nat[[1]] may be a list of zeros when the relevant
# first-order term vanishes, as for expected derivatives where E[score] = 0).
to_link_scale <- function(distrib, theta, nat, order) {
  params <- distrib@params
  p <- length(params)
  h <- inverse_link_derivs(distrib, theta, order)

  # At first order the chain rule is a diagonal rescaling, dl/deta_i = (dl/dtheta_i)
  # h_i'. Sending that through the general assembly below costs a few hundred
  # microseconds to perform a multiplication, and it is the order evaluated most
  # often, once per scoring iteration.
  if (order == 1L) {
    out <- lapply(seq_len(p), function(i) nat[[1L]][[params[i]]] * h[[i]][[1L]])
    names(out) <- params
    return(out)
  }

  terms <- deriv_index_list(p, order)
  out <- vector("list", length(terms))
  nms <- character(length(terms))

  for (t in seq_along(terms)) {
    idx <- terms[[t]]
    nms[t] <- paste(params[idx], collapse = "_")

    uniq <- unique(idx)
    mult <- tabulate(match(idx, uniq), length(uniq))

    # Nested sum over j_t = 1..m_t. The combinations are enumerated by decoding a
    # counter in mixed radix rather than built with expand.grid, which costs
    # sixty microseconds a term here -- more than everything else in this loop
    # put together, and paid once per component on every call.
    n_comb <- prod(mult)
    radix <- c(1L, cumprod(mult)[-length(mult)])
    acc <- 0

    for (r in seq_len(n_comb)) {
      j <- as.integer(((r - 1L) %/% radix) %% mult + 1L)
      nat_idx <- sort(rep(uniq, times = j))
      key <- paste(params[nat_idx], collapse = "_")
      term <- nat[[length(nat_idx)]][[key]]
      if (is.null(term)) {
        stop(sprintf("Missing parameter-scale derivative component '%s'.", key), call. = FALSE)
      }
      coef <- 1
      for (s in seq_along(uniq)) {
        coef <- coef * bell_partial(mult[s], j[s], h[[uniq[s]]])
      }
      acc <- acc + term * coef
    }

    out[[t]] <- acc
  }

  names(out) <- nms
  out
}
