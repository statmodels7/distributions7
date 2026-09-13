#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R
NULL

# The derivative of the expected information in the parameters.
#
# It exists for one consumer: a marginal criterion whose penalized matrix is
# K = -E[l''] + S enters through its determinant, so its gradient asks for
# dK/dbeta, which is dE[l'']/deta contracted with the design. With the OBSERVED
# information that object is -l''', which every family carries; with the
# expected one it is not, because differentiating an expectation moves the
# measure as well as the integrand,
#
#   d/dtheta_c E[l_ab] = E[l_abc] + E[l_ab l_c],
#
# and the second term is a mixed moment no Bartlett identity isolates (the
# third ties the SYMMETRIZED sum, not the single term).
#
# What makes the identity unnecessary for almost every family is that E[l_ab]
# is itself an explicit function of theta wherever the family wrote its
# expected information out: its derivative is then an ordinary derivative of a
# known formula, and ONE central difference of that formula is a single stencil
# on an analytic quantity -- the same licence the skew t has for its degrees of
# freedom, and not the nested differencing the package forbids.
#
# Where the expected information is itself an integral the licence lapses, and
# the reason is COST rather than accuracy. Measured at 100 observations, the
# six families that approximate it cost 1880 to 147300 ms against a median of
# 0.183 ms for the thirty-four that do not -- four orders of magnitude -- so a
# derivative asking for 2p of those calls per evaluation is not a slower route
# but an unusable one.

#' The Derivative of the Expected Information
#'
#' @description
#' \eqn{\partial\,\mathbb{E}[\ell_{ab}]/\partial\theta_c}, one component per
#' pair \eqn{(a,b)} and differentiating parameter \eqn{c}.
#'
#' @details
#' The components are symmetric in \eqn{(a,b)} and NOT in \eqn{c}: writing
#' \deqn{\frac{\partial}{\partial\theta_c}\mathbb{E}[\ell_{ab}]
#'   = \mathbb{E}[\ell_{abc}] + \mathbb{E}[\ell_{ab}\ell_{c}],}
#' the first term is fully symmetric and the second is not, so the result is
#' keyed by [dexpected_names()] rather than by the sorted triples
#' [deriv_names()] uses at order three.
#'
#' **The default method differences the family's own expected
#' information**, one central stencil per parameter, which is a single
#' difference of an analytic quantity wherever that quantity is a written-out
#' formula. It is refused where it is not: see
#' [has_exact_expected_hessian()].
#'
#' On `scale = "link"` the difference is taken along the free scale of the
#' parameter being differentiated, and the expected information is read on the
#' link scale at each of the two points, so the chain rule is never written out
#' here and cannot disagree with the one
#' [distrib_expected_hessian()] already applies.
#'
#' @param distrib A distribution object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, each of length 1 or
#'   `length(y)`.
#' @param scale `"parameter"` for \eqn{\partial/\partial\theta_c},
#'   `"link"` for \eqn{\partial/\partial\eta_c} of the link-scale expected
#'   information.
#' @param approx,nsim Passed to [distrib_expected_hessian()].
#' @param ... Passed to methods.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`dexpected_names(distrib@params)`][dexpected_names].
#'
#' @examples
#' d <- gaussian1_distrib()
#' str(distrib_dexpected_hessian(d, 0, list(mu = 0, sigma = 1)))
#'
#' @seealso [distrib_expected_hessian()],
#'   [dexpected_names()], [has_exact_expected_hessian()]
#'
#' @export
distrib_dexpected_hessian <- S7::new_generic(
  "distrib_dexpected_hessian", "distrib",
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000,
           ...) {
    args <- check_derivative_args(distrib, y, theta)
    y <- args$y
    theta <- args$theta
    S7::S7_dispatch()
  })


#' @title Default Derivative of the Expected Information
#' @name distrib_dexpected_hessian.distrib
#' @description
#' One central difference of [distrib_expected_hessian()] per
#' parameter, refused where that quantity is itself approximated.
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either `"parameter"` or `"link"`.
#' @param approx,nsim Passed through.
#' @param ... Unused.
#' @return A named list keyed as [dexpected_names()].
#' @keywords internal
S7::method(distrib_dexpected_hessian, distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  numerical_dexpected_hessian(distrib, y, theta, match.arg(scale),
                              match.arg(approx), nsim)
}


#' The Names of the Expected Information's Derivative
#'
#' @description
#' One key per pair \eqn{(a,b)} and differentiating parameter \eqn{c}, built by
#' joining [hess_names()] with the parameter differentiated in.
#'
#' @details
#' The keys are BUILT and never parsed, which is the package's rule wherever a
#' component name is a concatenation of parameter names: a parameter whose own
#' name contains an underscore makes the string ambiguous to read back, and
#' [dexpected_key()] exists so that a consumer composes the same
#' string this function enumerates.
#'
#' @param params A character vector of parameter names, in the family's order.
#'
#' @return A character vector, `length(hess_names(params)) *
#'   length(params)` long.
#'
#' @examples
#' dexpected_names(c("mu", "sigma"))
#'
#' @seealso [dexpected_key()], [hess_names()]
#'
#' @export
dexpected_names <- function(params) {
  hn <- hess_names(params)
  as.vector(t(outer(hn, params, paste, sep = "_")))
}


#' The Key of One Component of the Expected Information's Derivative
#'
#' @description
#' The name under which [distrib_dexpected_hessian()] returns
#' \eqn{\partial\,\mathbb{E}[\ell_{ab}]/\partial\theta_c}.
#'
#' @param params A character vector of parameter names, in the family's order.
#' @param a,b Indices into `params`; their order does not matter, the
#'   component being symmetric in them.
#' @param k The index of the parameter differentiated in, which does matter.
#'
#' @return A single string.
#'
#' @examples
#' dexpected_key(c("mu", "sigma"), 1, 2, 2)
#'
#' @seealso [dexpected_names()]
#'
#' @export
dexpected_key <- function(params, a, b, k) {
  nm <- hess_names(params)
  # hess_names() lists the diagonal first and then the off-diagonal pairs in
  # the family's own order, so the pair is looked up both ways round rather
  # than assumed to be sorted -- the trap deriv_index_list() records.
  want <- paste(params[c(a, b)], collapse = "_")
  if (!want %in% nm) want <- paste(params[c(b, a)], collapse = "_")
  if (!want %in% nm) {
    stop(sprintf("No Hessian component for '%s' and '%s'.",
                 params[a], params[b]), call. = FALSE)
  }
  paste0(want, "_", params[k])
}


#' Differencing the Expected Information Once
#'
#' @description
#' The default route behind [distrib_dexpected_hessian()]: a central
#' difference in each parameter of the family's own expected information.
#'
#' @details
#' The step is [fd_steps()]', which shrinks near a finite boundary so
#' that both evaluation points stay strictly inside the parameter's open
#' domain. On the link scale the domain is the whole line and no clamp is
#' needed, so the step is the plain relative one.
#'
#' **It refuses rather than approximating an approximation.** Where the
#' expected information is obtained by quadrature or by simulation, this would
#' be a difference of a difference, which the package forbids everywhere, and
#' it would cost 2p of the dearest call the family has.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either `"parameter"` or `"link"`.
#' @param approx,nsim Passed to [distrib_expected_hessian()].
#' @param h_rel The relative step, a cube root of machine epsilon by default,
#'   the value at which a central difference balances truncation against
#'   rounding.
#'
#' @return A named list keyed as [dexpected_names()].
#'
#' @seealso [distrib_dexpected_hessian()], [fd_steps()]
#'
#' @keywords internal
numerical_dexpected_hessian <- function(distrib, y, theta,
                                        scale = c("parameter", "link"),
                                        approx = "bartlett", nsim = 10000,
                                        h_rel = numericals7::fd_step(1, 1L)) {
  scale <- match.arg(scale)
  if (!has_exact_expected_hessian(distrib)) {
    stop(sprintf(paste0(
      "'%s' approximates its expected information rather than writing it\n",
      "  out, so differencing it would be a difference of a difference. A\n",
      "  family that carries dE[l'']/dtheta must either write its expected\n",
      "  information in closed form or register its own\n",
      "  distrib_dexpected_hessian() method."), distrib@distrib_name),
      call. = FALSE)
  }
  params <- distrib@params
  theta <- align_theta(distrib, theta)
  link <- identical(scale, "link")
  hn <- hess_names(params)
  out <- stats::setNames(vector("list", length(hn) * length(params)),
                         dexpected_names(params))
  # A Monte Carlo expected information is not the same number twice, so the
  # self-consistency reading fd_stable_step() takes would be of that noise
  # rather than of the step: measured on pig1, the quotient moves by 2.75
  # relative between two evaluations at one step, against 17.3 between the step
  # and its half. The step is therefore chosen only where the quantity being
  # differenced is deterministic.
  choose_step <- !identical(approx, "mc")

  for (k in seq_along(params)) {
    p <- params[k]
    E <- function(th) distrib_expected_hessian(distrib, y, th, scale = scale,
                                               approx = approx, nsim = nsim)
    # the point the difference is taken at, and the two points it is taken
    # between, on whichever scale the derivative was asked for
    if (link) {
      lk <- distrib@link_params[[p]]
      base <- linkfunctions7::linkfun(lk, theta[[p]])
      h <- h_rel * pmax(1, abs(base))
      shift <- function(hh) {
        tu <- td <- theta
        tu[[p]] <- linkfunctions7::linkinv(lk, base + hh)
        td[[p]] <- linkfunctions7::linkinv(lk, base - hh)
        list(tu, td)
      }
    } else {
      shift <- function(hh) {
        tu <- td <- theta
        tu[[p]] <- theta[[p]] + hh
        td[[p]] <- theta[[p]] - hh
        list(tu, td)
      }
    }
    quotient <- function(hh) {
      s <- shift(hh)
      a <- E(s[[1L]]); b <- E(s[[2L]])
      stats::setNames(lapply(hn, function(nm) (a[[nm]] - b[[nm]]) / (2 * hh)), hn)
    }
    # The link scale carries no finite bound, so there is nothing to choose
    # there and the magnitude step stands.
    q <- if (link || !choose_step) {
      if (!link) h <- fd_steps(theta[[p]], distrib@params_bounds[[p]], h_rel)
      quotient(h)
    } else {
      fd_stable_step(quotient, theta[[p]], distrib@params_bounds[[p]], h_rel)$value
    }
    for (nm in hn) {
      out[[paste0(nm, "_", p)]] <- q[[nm]]
    }
  }
  out
}


#' The Second Derivative of the Expected Information
#'
#' @description
#' \eqn{\partial^2\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d},
#' one component per pair \eqn{(a,b)} and per pair \eqn{(c,d)}.
#'
#' @details
#' Differentiating the identity of [distrib_dexpected_hessian()] once more
#' moves the measure a second time,
#' \deqn{\frac{\partial^2}{\partial\theta_c\,\partial\theta_d}\mathbb{E}[\ell_{ab}]
#'   = \mathbb{E}[\ell_{abcd}] + \mathbb{E}[\ell_{abd}\ell_{c}]
#'   + \mathbb{E}[\ell_{abc}\ell_{d}] + \mathbb{E}[\ell_{ab}\ell_{cd}]
#'   + \mathbb{E}[\ell_{ab}\ell_{c}\ell_{d}],}
#' and no Bartlett identity isolates those moments. A family supplies the
#' components as ordinary derivatives of its own written-out expected
#' information; the components are symmetric in \eqn{(a,b)} and in \eqn{(c,d)}
#' separately, and are keyed by [d2expected_names()].
#'
#' **There is no numerical default.** A difference of
#' [distrib_dexpected_hessian()], whose own default is already a difference,
#' would be the nested differencing the package forbids, so the base method
#' signals an error and a family that does not register one has no second
#' derivative. The families that do are `gaussian1_distrib()`,
#' `poisson_distrib()`, `gamma1_distrib()`, `negbin2_distrib()` and
#' `beta1_distrib()`, each from a compiled kernel.
#'
#' On `scale = "link"` the expected information is
#' \eqn{F_{ab} = \mathbb{E}[\ell_{ab}]\,h_a' h_b'}, with no term in \eqn{h''}
#' because \eqn{\mathbb{E}[\ell_a] = 0}, and its derivatives follow by
#' Leibniz's rule, written once in [dexpected_link()].
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
#'   [`d2expected_names(distrib@params)`][d2expected_names].
#'
#' @examples
#' d <- gaussian1_distrib()
#' str(distrib_d2expected_hessian(d, 0, list(mu = 0, sigma = 1)))
#'
#' @seealso [distrib_dexpected_hessian()], [d2expected_names()]
#'
#' @export
distrib_d2expected_hessian <- S7::new_generic(
  "distrib_d2expected_hessian", "distrib",
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000,
           ...) {
    args <- check_derivative_args(distrib, y, theta)
    y <- args$y
    theta <- args$theta
    S7::S7_dispatch()
  })


#' @title Default Second Derivative of the Expected Information
#' @name distrib_d2expected_hessian.distrib
#' @description
#' Signals an error: a family without its own method has no second derivative
#' of the expected information, a difference of a difference being refused.
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @return Does not return.
#' @keywords internal
S7::method(distrib_d2expected_hessian, distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  stop(sprintf(paste0(
    "'%s' has no analytic second derivative of its expected information,\n",
    "  and differencing the first would be a difference of a difference."),
    distrib@distrib_name), call. = FALSE)
}


#' The Names of the Expected Information's Second Derivative
#'
#' @description
#' One key per pair \eqn{(a,b)} and pair \eqn{(c,d)}, each pair spelled as
#' [hess_names()] spells it, joined `ab` first.
#'
#' @param params A character vector of parameter names, in the family's order.
#'
#' @return A character vector, `length(hess_names(params))^2` long.
#'
#' @examples
#' d2expected_names(c("mu", "sigma"))
#'
#' @seealso [d2expected_key()], [dexpected_names()]
#'
#' @export
d2expected_names <- function(params) {
  hn <- hess_names(params)
  as.vector(t(outer(hn, hn, paste, sep = "_")))
}


#' The Key of One Component of the Expected Information's Second Derivative
#'
#' @description
#' The name under which [distrib_d2expected_hessian()] returns
#' \eqn{\partial^2\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d}.
#'
#' @param params A character vector of parameter names, in the family's order.
#' @param a,b Indices of the information's pair; their order does not matter.
#' @param c,d Indices of the parameters differentiated in; their order does not
#'   matter either.
#'
#' @return A single string.
#'
#' @examples
#' d2expected_key(c("mu", "sigma"), 1, 1, 2, 2)
#'
#' @seealso [d2expected_names()], [dexpected_key()]
#'
#' @export
d2expected_key <- function(params, a, b, c, d) {
  paste0(hess_pair_name(params, a, b), "_", hess_pair_name(params, c, d))
}


#' The Hessian Name of a Pair of Parameters
#'
#' @description
#' `params[a]_params[b]` in whichever order [hess_names()] lists it.
#'
#' @param params A character vector of parameter names.
#' @param a,b Indices into `params`.
#'
#' @return A single string.
#'
#' @keywords internal
hess_pair_name <- function(params, a, b) {
  nm <- hess_names(params)
  want <- paste(params[c(a, b)], collapse = "_")
  if (!want %in% nm) want <- paste(params[c(b, a)], collapse = "_")
  if (!want %in% nm) {
    stop(sprintf("No Hessian component for '%s' and '%s'.",
                 params[a], params[b]), call. = FALSE)
  }
  want
}


#' The Analytic Derivatives of the Expected Information, on Either Scale
#'
#' @description
#' Reads a family's compiled kernel at order 1 or 2 on the parameter scale and
#' carries the result onto the link scale where it is asked for.
#'
#' @param distrib A distribution object.
#' @param y,theta As the generic takes them.
#' @param scale `"parameter"` or `"link"`.
#' @param order `1L` for [distrib_dexpected_hessian()], `2L` for
#'   [distrib_d2expected_hessian()].
#' @param threads The thread count passed to the kernel and to the family's
#'   expected information.
#' @param kern A function of the order returning the kernel's named list on the
#'   parameter scale.
#'
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#'
#' @seealso [dexpected_link()]
#'
#' @keywords internal
dexpected_analytic <- function(distrib, y, theta, scale, order, threads, kern) {
  params <- distrib@params
  d1 <- kern(1L)[dexpected_names(params)]
  d2 <- if (order == 2L) kern(2L)[d2expected_names(params)] else NULL
  if (identical(scale, "parameter")) return(if (order == 1L) d1 else d2)
  E <- distrib_expected_hessian(distrib, y, theta, scale = "parameter",
                                threads = threads)
  h <- inverse_link_derivs(distrib, theta, order + 1L)
  dexpected_link(params, E, d1, d2, h, order)
}


#' The Derivatives of the Link-Scale Expected Information
#'
#' @description
#' Leibniz's rule on \eqn{F_{ab} = E_{ab}\,h_a' h_b'}, written once for every
#' family.
#'
#' @details
#' Each parameter moves with its own coordinate only, so with
#' \eqn{u_{ab} = h_a' h_b'},
#' \deqn{\partial_c u = [a{=}c]\,h_a'' h_b' + [b{=}c]\,h_a' h_b'',}
#' \deqn{\partial_{cd} u = [a{=}c{=}d]\,h_a''' h_b' + ([a{=}c][b{=}d] +
#'   [b{=}c][a{=}d])\,h_a'' h_b'' + [b{=}c{=}d]\,h_a' h_b''',}
#' and
#' \deqn{\partial_c F = (\partial_c E)\,h_c' u + E\,\partial_c u,}
#' \deqn{\partial_{cd} F = (\partial_{cd} E)\,h_c' h_d' u
#'   + [c{=}d]\,(\partial_c E)\,h_c'' u + (\partial_c E)\,h_c'\,\partial_d u
#'   + (\partial_d E)\,h_d'\,\partial_c u + E\,\partial_{cd} u.}
#'
#' @param params The parameter names.
#' @param E The expected information on the parameter scale.
#' @param d1 Its first derivatives, keyed as [dexpected_names()].
#' @param d2 Its second derivatives, keyed as [d2expected_names()], or `NULL`
#'   at order 1.
#' @param h The inverse link's derivatives, from [inverse_link_derivs()], to
#'   order `order + 1`.
#' @param order `1L` or `2L`.
#'
#' @return A named list on the link scale.
#'
#' @keywords internal
dexpected_link <- function(params, E, d1, d2, h, order) {
  p <- length(params)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  h1 <- function(i) h[[i]][[1L]]
  h2 <- function(i) h[[i]][[2L]]
  h3 <- function(i) h[[i]][[3L]]
  du <- function(a, b, c) {
    (if (a == c) h2(a) * h1(b) else 0) + (if (b == c) h1(a) * h2(b) else 0)
  }
  if (order == 1L) {
    out <- stats::setNames(vector("list", length(d1)), names(d1))
    for (r in seq_len(nrow(pairs))) {
      a <- pairs[r, 1L]; b <- pairs[r, 2L]
      Eab <- E[[hess_pair_name(params, a, b)]]
      u <- h1(a) * h1(b)
      for (c in seq_len(p)) {
        key <- dexpected_key(params, a, b, c)
        out[[key]] <- d1[[key]] * h1(c) * u + Eab * du(a, b, c)
      }
    }
    return(out)
  }
  d2u <- function(a, b, c, d) {
    (if (a == c && a == d) h3(a) * h1(b) else 0) +
      (if (a == c && b == d) h2(a) * h2(b) else 0) +
      (if (b == c && a == d) h2(a) * h2(b) else 0) +
      (if (b == c && b == d) h1(a) * h3(b) else 0)
  }
  out <- stats::setNames(vector("list", length(d2)), names(d2))
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    Eab <- E[[hess_pair_name(params, a, b)]]
    u <- h1(a) * h1(b)
    for (s in seq_len(nrow(pairs))) {
      c <- pairs[s, 1L]; d <- pairs[s, 2L]
      dc <- d1[[dexpected_key(params, a, b, c)]]
      dd <- d1[[dexpected_key(params, a, b, d)]]
      key <- d2expected_key(params, a, b, c, d)
      out[[key]] <- d2[[key]] * h1(c) * h1(d) * u +
        (if (c == d) dc * h2(c) * u else 0) +
        dc * h1(c) * du(a, b, d) + dd * h1(d) * du(a, b, c) +
        Eab * d2u(a, b, c, d)
    }
  }
  out
}
