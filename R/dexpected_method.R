#' @include expect_method.R expected_derivatives.R
NULL

# The derivatives of the expected-information matrix under an explicit method.
#
# A family that writes dE[l_ab]/dtheta_c in closed form returns it when `method`
# is NULL. An explicit method gives the matrix a meaning of its own, and the
# derivatives are those of that matrix:
#
# * expect_integrate() and expect_mc() define E_ab(theta) = E_theta[f_ab], the
#   expectation of the integrand f_ab of the identity (l_ab for "hessian",
#   -l_a l_b for "score"). The density moves with theta, so
#   d/dtheta_c E[f] = E[f_c + f l_c], and the operator D_c f = f_c + f l_c is
#   applied once per differentiating parameter.
# * expect_opg() defines the matrix at the data, the per-observation
#   -l_a l_b. The data do not move with theta, so the derivative is the
#   Leibniz derivative of the integrand alone, D_c f = f_c.
#
# Both are polynomials in the observed derivatives l_T of the log-density, so a
# term is a coefficient and a list of index tuples. The expectation of the
# polynomial is taken by the method.

#' Polynomial Algebra of Observed Derivatives
#'
#' @description
#' Helpers for the derivatives of the expected-information matrix under an
#' explicit method. A polynomial is a list of terms, each a list with a numeric
#' `coef` and `f`, a list of sorted integer index tuples, one per factor
#' \eqn{\ell_T}.
#'
#' @details
#' `poly_base()` returns the integrand of the identity for the pair \eqn{(a,b)}:
#' \eqn{\ell_{ab}} for `"hessian"` and \eqn{-\ell_a\ell_b} for `"score"`.
#' `poly_d()` applies the derivative with respect to the parameter `cc`: the
#' Leibniz derivative, which replaces each factor \eqn{\ell_T} by
#' \eqn{\ell_{T\cup c}}, and, when `moved` is `TRUE`, also the term
#' \eqn{f\,\ell_c} that comes from the density. `poly_order()` is the largest
#' derivative order that appears, and `poly_eval()` evaluates the polynomial
#' with a function that returns the observed derivative of a tuple.
#'
#' @param identity `"hessian"` or `"score"`.
#' @param a,b Parameter indices of the pair.
#' @param poly A polynomial.
#' @param cc The index of the differentiating parameter.
#' @param moved Whether the density moves with the parameter.
#' @param getter A function of a sorted index tuple returning the observed
#'   derivative, a numeric vector.
#'
#' @return `poly_base()` and `poly_d()` return a polynomial, `poly_order()` an
#'   integer and `poly_eval()` a numeric vector.
#'
#' @seealso [dexpected_by_method()]
#' @name dexpected_poly
#' @keywords internal
NULL

#' @rdname dexpected_poly
poly_base <- function(identity, a, b) {
  if (identical(identity, "hessian")) {
    list(list(coef = 1, f = list(sort(c(a, b)))))
  } else {
    list(list(coef = -1, f = list(a, b)))
  }
}

#' @rdname dexpected_poly
poly_d <- function(poly, cc, moved) {
  out <- list()
  for (tm in poly) {
    for (i in seq_along(tm$f)) {
      f <- tm$f
      f[[i]] <- sort(c(f[[i]], cc))
      out[[length(out) + 1L]] <- list(coef = tm$coef, f = f)
    }
    if (moved) {
      out[[length(out) + 1L]] <- list(coef = tm$coef, f = c(tm$f, list(cc)))
    }
  }
  out
}

#' @rdname dexpected_poly
poly_order <- function(poly) {
  max(vapply(poly, function(tm) max(lengths(tm$f)), numeric(1)))
}

#' @rdname dexpected_poly
poly_eval <- function(poly, getter) {
  acc <- 0
  for (tm in poly) {
    v <- tm$coef
    for (fac in tm$f) v <- v * getter(fac)
    acc <- acc + v
  }
  acc
}

#' Observed Derivatives by Index Tuple
#'
#' @description
#' Returns a function that gives the observed derivative of the log-density for
#' a sorted index tuple at fixed observations and parameters, computing each
#' order once.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param ys The observations.
#' @param th The parameters, aligned.
#'
#' @return A function of an integer tuple.
#'
#' @seealso [dexpected_by_method()]
#' @keywords internal
observed_getter <- function(distrib, ys, th) {
  params <- distrib@params
  cache <- vector("list", 5L)
  function(tuple) {
    k <- length(tuple)
    if (is.null(cache[[k]])) cache[[k]] <<- observed_deriv(distrib, ys, th, k)
    key <- if (k == 1L) params[tuple] else if (k == 2L) {
      hess_pair_name(params, tuple[1L], tuple[2L])
    } else paste(params[tuple], collapse = "_")
    cache[[k]][[key]]
  }
}

#' The Components of a Derivative of the Expected Information
#'
#' @description
#' Enumerates, for a derivative order, the pairs \eqn{(a,b)} and the multisets
#' of differentiating parameters, with the key under which
#' [distrib_dexpected_hessian()] and its relatives name each component.
#'
#' @param params The parameter names.
#' @param m The derivative order, 0 to 4.
#'
#' @return A list of components, each with `key`, `a`, `b` and `S`.
#'
#' @seealso [dexpected_by_method()]
#' @keywords internal
dexpected_components <- function(params, m) {
  p <- length(params)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  tuples <- if (m == 0L) list(integer(0)) else if (m == 1L) as.list(seq_len(p)) else
    if (m == 2L) lapply(seq_len(nrow(pairs)), function(r) as.integer(pairs[r, ])) else
      deriv_indices(params, m)
  ekey <- function(a, b, S) {
    switch(length(S) + 1L,
      hess_pair_name(params, a, b),
      dexpected_key(params, a, b, S[1L]),
      d2expected_key(params, a, b, S[1L], S[2L]),
      d3expected_key(params, a, b, S[1L], S[2L], S[3L]),
      d4expected_key(params, a, b, S[1L], S[2L], S[3L], S[4L]))
  }
  out <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    for (S in tuples) {
      out[[length(out) + 1L]] <- list(key = ekey(a, b, S), a = a, b = b, S = S)
    }
  }
  out
}

#' Derivatives of the Expected Information Under an Explicit Method
#'
#' @description
#' Returns the derivative of order `m` of the matrix that an
#' [ExpectMethod] defines, with respect to the parameters, on the parameter
#' scale or the link scale. It is the route of
#' [distrib_dexpected_hessian()], [distrib_d2expected_hessian()],
#' [distrib_d3expected_hessian()] and [distrib_d4expected_hessian()] when
#' `method` is an object.
#'
#' @details
#' For [expect_integrate()] and [expect_mc()] the matrix is
#' \eqn{E_{ab}(\theta) = \mathbb{E}_\theta[f_{ab}]} and each derivative moves
#' the density as well as the integrand, so that the integrand of the
#' derivative along \eqn{c} is \eqn{\partial_c f + f\,\ell_c}. For
#' [expect_opg()] the matrix is \eqn{f_{ab}} evaluated at the observations, the
#' data do not move, and the derivative is the Leibniz derivative of
#' \eqn{f_{ab}} alone: \eqn{-(\ell_{ac}\ell_b + \ell_a\ell_{bc})} at order 1 and
#' \eqn{-(\ell_{acd}\ell_b + \ell_{ac}\ell_{bd} + \ell_{ad}\ell_{bc} +
#' \ell_a\ell_{bcd})} at order 2. The integrand \eqn{f_{ab}} is
#' \eqn{\ell_{ab}} for `identity = "hessian"` and \eqn{-\ell_a\ell_b} for
#' `identity = "score"`; [expect_opg()] uses the second.
#'
#' The derivative of order `m` of \eqn{f_{ab}} needs observed derivatives of
#' order up to `m + 1` for the score identity and `m + 2` for the Hessian
#' identity. The package provides five orders, so the Hessian identity at
#' `m = 4` signals an error that names `identity = "score"`.
#'
#' On the link scale the parameter-scale derivatives of every order up to `m`
#' are computed by the same method and carried across by
#' [dexpected_link_higher()].
#'
#' @param distrib A univariate distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned.
#' @param scale `"parameter"` or `"link"`.
#' @param m The derivative order, 1 to 4.
#' @param method An object inheriting from [ExpectMethod].
#'
#' @return A named list keyed as [dexpected_names()], [d2expected_names()],
#'   [d3expected_names()] or [d4expected_names()].
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
dexpected_by_method <- function(distrib, y, theta, scale, m, method) {
  if (S7::S7_inherits(distrib, multivariate_distrib)) {
    stop(sprintf(paste0(
      "'%s' is multivariate, and the derivatives of its expected information ",
      "have no method other than the closed form of the family."),
      distrib@distrib_name), call. = FALSE)
  }
  params <- distrib@params
  n <- length(y)
  opg <- S7::S7_inherits(method, ExpectOPG)
  identity <- if (opg) "score" else method@identity
  need <- m + if (identical(identity, "hessian")) 2L else 1L
  if (need > 5L) {
    stop(sprintf(paste0(
      "The derivative of order %d of the expected information needs ",
      "observed derivatives of order %d for identity = \"hessian\", and ",
      "the package provides five. Use identity = \"score\"."), m, need),
      call. = FALSE)
  }
  if (identical(identity, "hessian")) {
    check_kink_identity(distrib, 2L, if (S7::S7_inherits(method, ExpectMC))
      "expect_mc" else "expect_integrate")
  }
  nat <- vector("list", m + 1L)
  for (j in 0:m) {
    comps <- dexpected_components(params, j)
    polys <- lapply(comps, function(cp) {
      pl <- poly_base(identity, cp$a, cp$b)
      for (cc in cp$S) pl <- poly_d(pl, cc, moved = !opg)
      pl
    })
    vals <- dexpected_evaluate(distrib, y, theta, polys, method)
    names(vals) <- vapply(comps, function(cp) cp$key, character(1))
    canon <- switch(j + 1L, hess_names(params), dexpected_names(params),
                    d2expected_names(params), d3expected_names(params),
                    d4expected_names(params))
    nat[[j + 1L]] <- vals[canon]
  }
  if (identical(scale, "parameter")) return(nat[[m + 1L]])
  h <- inverse_link_derivs(distrib, theta, m + 1L)
  dexpected_link_higher(params, nat, h, m)
}

#' Evaluate Polynomials of Observed Derivatives by a Method
#'
#' @description
#' Computes, for each polynomial, its value at the observations
#' ([expect_opg()]), its expectation by quadrature or summation
#' ([expect_integrate()]) or its Monte Carlo mean ([expect_mc()]).
#'
#' @param distrib A univariate distribution.
#' @param y A numeric vector of observations; only its length is used by the
#'   methods that take an expectation.
#' @param theta A named list of parameters, aligned.
#' @param polys A list of polynomials.
#' @param method An object inheriting from [ExpectMethod].
#'
#' @return A list of numeric vectors, each of length `length(y)`.
#'
#' @seealso [dexpected_by_method()]
#' @keywords internal
dexpected_evaluate <- function(distrib, y, theta, polys, method) {
  n <- length(y)
  if (S7::S7_inherits(method, ExpectOPG)) {
    getter <- observed_getter(distrib, y, theta)
    return(lapply(polys, function(pl) rep_len(poly_eval(pl, getter), n)))
  }
  if (S7::S7_inherits(method, ExpectIntegrate)) {
    old <- set_integration_settings(method)
    on.exit(set_integration_settings(old))
    return(lapply(polys, function(pl) {
      v <- expectation(distrib, function(y, theta) {
        poly_eval(pl, observed_getter(distrib, y, theta))
      }, theta)
      rep(v, length.out = n)
    }))
  }
  p <- distrib@n_params
  th_rows <- transpose_params(expand_params(theta[seq_len(p)]))
  per_row <- lapply(th_rows, function(r) {
    th <- as.list(r)
    ys <- distrib_rng(distrib, method@nsim, th)
    getter <- observed_getter(distrib, ys, th)
    vapply(polys, function(pl) base::mean(poly_eval(pl, getter)), numeric(1))
  })
  m <- do.call(rbind, per_row)
  lapply(seq_along(polys), function(k) rep(unname(m[, k]), length.out = n))
}
