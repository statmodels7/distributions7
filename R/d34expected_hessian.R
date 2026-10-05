#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R dexpected_hessian.R
NULL

# The third and fourth derivatives of the expected information in the
# parameters.
#
# The consumer is a score-driven filter whose score is scaled by a power of the
# expected information: the k-th derivative of the scaled score asks for the
# k-th derivative of the information, and the exact derivatives of a fit reach
# order four (the outer criterion's Hessian). No family writes these orders
# out. Every family that carries an analytic SECOND derivative gets them from
# one stencil on it -- a first difference at order three, one second-order
# stencil at order four -- which is a single difference of an analytic
# quantity, the licence the first derivative's default already uses. A family
# without an analytic second derivative has no third or fourth: differencing a
# first derivative twice, or a numerical one at all, is the nested differencing
# the package forbids.
#
# Measured against symbolic derivatives (base R D()) of the link-scale
# information of poisson, gaussian1, negbin2 and beta1: order three to
# 1.3e-11 .. 2.5e-10 relative, order four to 3.4e-09 .. 3.3e-07.

#' The Third Derivative of the Expected Information
#'
#' @description
#' \eqn{\partial^3\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d\,
#' \partial\theta_e}, one component per pair \eqn{(a,b)} and per unordered
#' triple \eqn{(c,d,e)}.
#'
#' @details
#' The components are symmetric in \eqn{(a,b)} and in \eqn{(c,d,e)}
#' separately, and are keyed by [d3expected_names()].
#'
#' The default method takes one central difference of
#' [distrib_d2expected_hessian()] along each parameter, so it requires a family
#' that registers an analytic second derivative of its expected information
#' and signals an error otherwise. The step is the cube root of machine epsilon
#' relative to the coordinate, kept inside a finite bound on the parameter
#' scale by [fd_steps()].
#'
#' On `scale = "link"` the difference is taken along the free scale of the
#' parameter, of the link-scale second derivative, so the chain rule is the
#' one [dexpected_link()] already applies.
#'
#' @param distrib A distribution object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, each of length 1 or
#'   `length(y)`.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Accepted for symmetry with
#'   [distrib_dexpected_hessian()]; no method reads them.
#' @param ... Passed to methods.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`d3expected_names(distrib@params)`][d3expected_names].
#'
#' @examples
#' d <- gaussian1_distrib()
#' str(distrib_d3expected_hessian(d, 0, list(mu = 0, sigma = 1)))
#'
#' @seealso [distrib_d2expected_hessian()], [distrib_d4expected_hessian()],
#'   [d3expected_names()]
#'
#' @export
distrib_d3expected_hessian <- S7::new_generic(
  "distrib_d3expected_hessian", "distrib",
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000,
           ...) {
    args <- check_derivative_args(distrib, y, theta)
    y <- args$y
    theta <- args$theta
    S7::S7_dispatch()
  })


#' @title Default Third Derivative of the Expected Information
#' @name distrib_d3expected_hessian.distrib
#' @description
#' One central difference of [distrib_d2expected_hessian()] per parameter;
#' an error where the family has no analytic second derivative.
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @return A named list keyed as [d3expected_names()].
#' @keywords internal
S7::method(distrib_d3expected_hessian, distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  numerical_dexpected_higher(distrib, y, theta, match.arg(scale), 3L)
}


#' The Fourth Derivative of the Expected Information
#'
#' @description
#' \eqn{\partial^4\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d\,
#' \partial\theta_e\,\partial\theta_f}, one component per pair \eqn{(a,b)} and
#' per unordered quadruple \eqn{(c,d,e,f)}.
#'
#' @details
#' The components are symmetric in \eqn{(a,b)} and in \eqn{(c,d,e,f)}
#' separately, and are keyed by [d4expected_names()].
#'
#' The default method applies one second-order central stencil to
#' [distrib_d2expected_hessian()]: the three-point stencil along one parameter
#' for a repeated pair \eqn{(e,e)}, the four-point mixed stencil for two
#' distinct parameters \eqn{(e,f)}. It requires a family that registers an
#' analytic second derivative of its expected information and signals an error
#' otherwise. The step is the fourth root of machine epsilon relative to the
#' coordinate, kept inside a finite bound on the parameter scale by
#' [fd_steps()].
#'
#' @inheritParams distrib_d3expected_hessian
#'
#' @return A named list of numeric vectors, keyed as
#'   [`d4expected_names(distrib@params)`][d4expected_names].
#'
#' @examples
#' d <- poisson_distrib()
#' str(distrib_d4expected_hessian(d, 2, list(mu = 3), scale = "link"))
#'
#' @seealso [distrib_d2expected_hessian()], [distrib_d3expected_hessian()],
#'   [d4expected_names()]
#'
#' @export
distrib_d4expected_hessian <- S7::new_generic(
  "distrib_d4expected_hessian", "distrib",
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000,
           ...) {
    args <- check_derivative_args(distrib, y, theta)
    y <- args$y
    theta <- args$theta
    S7::S7_dispatch()
  })


#' @title Default Fourth Derivative of the Expected Information
#' @name distrib_d4expected_hessian.distrib
#' @description
#' One second-order stencil on [distrib_d2expected_hessian()]; an error where
#' the family has no analytic second derivative.
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @return A named list keyed as [d4expected_names()].
#' @keywords internal
S7::method(distrib_d4expected_hessian, distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  numerical_dexpected_higher(distrib, y, theta, match.arg(scale), 4L)
}


#' The Names of the Expected Information's Third Derivative
#'
#' @description
#' One key per pair \eqn{(a,b)}, spelled as [hess_names()] spells it, and per
#' unordered triple, spelled as [deriv_names()] spells it at order three,
#' joined `ab` first.
#'
#' @param params A character vector of parameter names, in the family's order.
#'
#' @return A character vector, `length(hess_names(params)) *
#'   length(deriv_names(params, 3))` long.
#'
#' @examples
#' d3expected_names(c("mu", "sigma"))
#'
#' @seealso [d3expected_key()], [d2expected_names()]
#'
#' @export
d3expected_names <- function(params) {
  as.vector(t(outer(hess_names(params), deriv_names(params, 3L), paste,
                    sep = "_")))
}


#' The Names of the Expected Information's Fourth Derivative
#'
#' @description
#' One key per pair \eqn{(a,b)}, spelled as [hess_names()] spells it, and per
#' unordered quadruple, spelled as [deriv_names()] spells it at order four,
#' joined `ab` first.
#'
#' @param params A character vector of parameter names, in the family's order.
#'
#' @return A character vector, `length(hess_names(params)) *
#'   length(deriv_names(params, 4))` long.
#'
#' @examples
#' d4expected_names(c("mu", "sigma"))
#'
#' @seealso [d4expected_key()], [d3expected_names()]
#'
#' @export
d4expected_names <- function(params) {
  as.vector(t(outer(hess_names(params), deriv_names(params, 4L), paste,
                    sep = "_")))
}


#' The Key of One Component of the Expected Information's Third Derivative
#'
#' @description
#' The name under which [distrib_d3expected_hessian()] returns
#' \eqn{\partial^3\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d\,
#' \partial\theta_e}.
#'
#' @param params A character vector of parameter names, in the family's order.
#' @param a,b Indices of the information's pair; their order does not matter.
#' @param c,d,e Indices of the parameters differentiated in; their order does
#'   not matter either.
#'
#' @return A single string.
#'
#' @examples
#' d3expected_key(c("mu", "sigma"), 1, 1, 2, 1, 2)
#'
#' @seealso [d3expected_names()]
#'
#' @export
d3expected_key <- function(params, a, b, c, d, e) {
  paste0(hess_pair_name(params, a, b), "_",
         paste(params[sort(c(c, d, e))], collapse = "_"))
}


#' The Key of One Component of the Expected Information's Fourth Derivative
#'
#' @description
#' The name under which [distrib_d4expected_hessian()] returns
#' \eqn{\partial^4\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d\,
#' \partial\theta_e\,\partial\theta_f}.
#'
#' @param params A character vector of parameter names, in the family's order.
#' @param a,b Indices of the information's pair; their order does not matter.
#' @param c,d,e,f Indices of the parameters differentiated in; their order
#'   does not matter either.
#'
#' @return A single string.
#'
#' @examples
#' d4expected_key(c("mu", "sigma"), 1, 1, 2, 2, 2, 1)
#'
#' @seealso [d4expected_names()]
#'
#' @export
d4expected_key <- function(params, a, b, c, d, e, f) {
  paste0(hess_pair_name(params, a, b), "_",
         paste(params[sort(c(c, d, e, f))], collapse = "_"))
}


#' Does a Family Register an Analytic Second Derivative of Its Expected
#' Information?
#'
#' @description
#' `TRUE` where the [distrib_d2expected_hessian()] method that dispatches is
#' registered on a class other than the base class, whose method signals an
#' error.
#'
#' @param x A distribution object.
#'
#' @return A single logical.
#'
#' @keywords internal
has_d2expected_hessian <- function(x) {
  m <- tryCatch(S7::method(distrib_d2expected_hessian, S7::S7_class(x)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  !is_class(attr(m, "signature")[[1]], distrib)
}


#' The Third or Fourth Derivative of the Expected Information by One Stencil
#'
#' @description
#' The default route behind [distrib_d3expected_hessian()] and
#' [distrib_d4expected_hessian()]: one central stencil on the family's
#' analytic second derivative of its expected information.
#'
#' @details
#' A component of order three, keyed by the pair \eqn{(a,b)} and the sorted
#' triple \eqn{(c,d,e)}, is the first difference along \eqn{e} of the second
#' derivative in \eqn{(c,d)}. A component of order four, keyed by the sorted
#' quadruple \eqn{(c,d,e,f)}, is the second-order stencil in \eqn{(e,f)} of the
#' second derivative in \eqn{(c,d)}. Each shifted evaluation is computed once
#' and shared by every component that reads it: \eqn{2p} evaluations at order
#' three and \eqn{1 + 2p^2} at order four.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param order `3L` or `4L`.
#'
#' @return A named list keyed as [d3expected_names()] or [d4expected_names()].
#'
#' @keywords internal
numerical_dexpected_higher <- function(distrib, y, theta, scale, order) {
  if (!has_d2expected_hessian(distrib)) {
    stop(sprintf(paste0(
      "'%s' has no analytic second derivative of its expected information,\n",
      "  and its derivative of order %d is one stencil on that second\n",
      "  derivative."), distrib@distrib_name, order), call. = FALSE)
  }
  params <- distrib@params
  p <- length(params)
  theta <- align_theta(distrib, theta)
  link <- identical(scale, "link")
  # the stencil's own order: a first difference at order three, a second-order
  # stencil at order four, each at the step that balances its truncation
  # against rounding
  st <- order - 2L
  h_rel <- .Machine$double.eps^(1 / (st + 2))
  eta <- if (link) {
    stats::setNames(lapply(params, function(q)
      linkfunctions7::linkfun(distrib@link_params[[q]], theta[[q]])), params)
  } else {
    NULL
  }
  step <- stats::setNames(lapply(params, function(q) {
    if (link) {
      h_rel * pmax(1, abs(eta[[q]]))
    } else {
      fd_steps(theta[[q]], distrib@params_bounds[[q]], h_rel, order = st)
    }
  }), params)
  # the second derivative at theta moved by `k` steps along each parameter
  D2 <- function(k) {
    th <- theta
    for (j in which(k != 0)) {
      q <- params[j]
      if (link) {
        th[[q]] <- linkfunctions7::linkinv(distrib@link_params[[q]],
                                           eta[[q]] + k[j] * step[[q]])
      } else {
        th[[q]] <- theta[[q]] + k[j] * step[[q]]
      }
    }
    distrib_d2expected_hessian(distrib, y, th, scale = scale)
  }
  unit <- function(j, s) { k <- integer(p); k[j] <- s; k }
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  tuples <- deriv_indices(params, order)
  out_names <- if (order == 3L) d3expected_names(params) else
    d4expected_names(params)
  out <- stats::setNames(vector("list", length(out_names)), out_names)

  if (order == 3L) {
    # one pair of evaluations per direction, shared by every triple ending in it
    diffs <- lapply(seq_len(p), function(j) {
      up <- D2(unit(j, 1L)); dn <- D2(unit(j, -1L))
      list(up = up, dn = dn, h = step[[params[j]]])
    })
    for (r in seq_len(nrow(pairs))) {
      a <- pairs[r, 1L]; b <- pairs[r, 2L]
      for (tu in tuples) {
        e <- tu[3L]
        k2 <- d2expected_key(params, a, b, tu[1L], tu[2L])
        dj <- diffs[[e]]
        out[[d3expected_key(params, a, b, tu[1L], tu[2L], e)]] <-
          (dj$up[[k2]] - dj$dn[[k2]]) / (2 * dj$h)
      }
    }
    return(out)
  }

  # order four: the centre once, the three-point stencil along each parameter,
  # the four-point mixed stencil for each distinct pair
  centre <- D2(integer(p))
  axis <- lapply(seq_len(p), function(j) {
    list(up = D2(unit(j, 1L)), dn = D2(unit(j, -1L)))
  })
  mixed <- list()
  stencil <- function(e, f, key) {
    he <- step[[params[e]]]
    if (e == f) {
      return((axis[[e]]$up[[key]] - 2 * centre[[key]] + axis[[e]]$dn[[key]]) /
               he^2)
    }
    id <- paste(e, f)
    if (is.null(mixed[[id]])) {
      k <- function(se, sf) { v <- integer(p); v[e] <- se; v[f] <- sf; v }
      mixed[[id]] <<- list(pp = D2(k(1L, 1L)), pm = D2(k(1L, -1L)),
                           mp = D2(k(-1L, 1L)), mm = D2(k(-1L, -1L)))
    }
    m <- mixed[[id]]
    hf <- step[[params[f]]]
    (m$pp[[key]] - m$pm[[key]] - m$mp[[key]] + m$mm[[key]]) / (4 * he * hf)
  }
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (tu in tuples) {
      k2 <- d2expected_key(params, a, b, tu[1L], tu[2L])
      out[[d4expected_key(params, a, b, tu[1L], tu[2L], tu[3L], tu[4L])]] <-
        stencil(tu[3L], tu[4L], k2)
    }
  }
  out
}
