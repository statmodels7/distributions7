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
#' keyed by \code{\link{dexpected_names}} rather than by the sorted triples
#' \code{\link{deriv_names}} uses at order three.
#'
#' \strong{The default method differences the family's own expected
#' information}, one central stencil per parameter, which is a single
#' difference of an analytic quantity wherever that quantity is a written-out
#' formula. It is refused where it is not: see
#' \code{\link{has_exact_expected_hessian}}.
#'
#' On \code{scale = "link"} the difference is taken along the free scale of the
#' parameter being differentiated, and the expected information is read on the
#' link scale at each of the two points, so the chain rule is never written out
#' here and cannot disagree with the one
#' \code{\link{distrib_expected_hessian}} already applies.
#'
#' @param distrib A distribution object inheriting from \code{distrib}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, each of length 1 or
#'   \code{length(y)}.
#' @param scale \code{"parameter"} for \eqn{\partial/\partial\theta_c},
#'   \code{"link"} for \eqn{\partial/\partial\eta_c} of the link-scale expected
#'   information.
#' @param approx,nsim Passed to \code{\link{distrib_expected_hessian}}.
#' @param ... Passed to methods.
#'
#' @return A named list of numeric vectors, keyed as
#'   \code{\link{dexpected_names}(distrib@params)}.
#'
#' @examples
#' d <- gaussian1_distrib()
#' str(distrib_dexpected_hessian(d, 0, list(mu = 0, sigma = 1)))
#'
#' @seealso \code{\link{distrib_expected_hessian}},
#'   \code{\link{dexpected_names}}, \code{\link{has_exact_expected_hessian}}
#'
#' @export
distrib_dexpected_hessian <- S7::new_generic(
  "distrib_dexpected_hessian", "distrib",
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000,
           ...) {
    args <- check_derivative_args(distrib, y, theta)
    y <- args$y
    theta <- args$theta
    S7::S7_dispatch()
  })


#' @title Default Derivative of the Expected Information
#' @name distrib_dexpected_hessian.distrib
#' @description
#' One central difference of \code{\link{distrib_expected_hessian}} per
#' parameter, refused where that quantity is itself approximated.
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param approx,nsim Passed through.
#' @param ... Unused.
#' @return A named list keyed as \code{\link{dexpected_names}}.
#' @keywords internal
S7::method(distrib_dexpected_hessian, distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  numerical_dexpected_hessian(distrib, y, theta, match.arg(scale),
                              match.arg(approx), nsim)
}


#' The Names of the Expected Information's Derivative
#'
#' @description
#' One key per pair \eqn{(a,b)} and differentiating parameter \eqn{c}, built by
#' joining \code{\link{hess_names}} with the parameter differentiated in.
#'
#' @details
#' The keys are BUILT and never parsed, which is the package's rule wherever a
#' component name is a concatenation of parameter names: a parameter whose own
#' name contains an underscore makes the string ambiguous to read back, and
#' \code{\link{dexpected_key}} exists so that a consumer composes the same
#' string this function enumerates.
#'
#' @param params A character vector of parameter names, in the family's order.
#'
#' @return A character vector, \code{length(hess_names(params)) *
#'   length(params)} long.
#'
#' @examples
#' dexpected_names(c("mu", "sigma"))
#'
#' @seealso \code{\link{dexpected_key}}, \code{\link{hess_names}}
#'
#' @export
dexpected_names <- function(params) {
  hn <- hess_names(params)
  as.vector(t(outer(hn, params, paste, sep = "_")))
}


#' The Key of One Component of the Expected Information's Derivative
#'
#' @description
#' The name under which \code{\link{distrib_dexpected_hessian}} returns
#' \eqn{\partial\,\mathbb{E}[\ell_{ab}]/\partial\theta_c}.
#'
#' @param params A character vector of parameter names, in the family's order.
#' @param a,b Indices into \code{params}; their order does not matter, the
#'   component being symmetric in them.
#' @param k The index of the parameter differentiated in, which does matter.
#'
#' @return A single string.
#'
#' @examples
#' dexpected_key(c("mu", "sigma"), 1, 2, 2)
#'
#' @seealso \code{\link{dexpected_names}}
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
#' The default route behind \code{\link{distrib_dexpected_hessian}}: a central
#' difference in each parameter of the family's own expected information.
#'
#' @details
#' The step is \code{\link{fd_steps}}', which shrinks near a finite boundary so
#' that both evaluation points stay strictly inside the parameter's open
#' domain. On the link scale the domain is the whole line and no clamp is
#' needed, so the step is the plain relative one.
#'
#' \strong{It refuses rather than approximating an approximation.} Where the
#' expected information is obtained by quadrature or by simulation, this would
#' be a difference of a difference, which the package forbids everywhere, and
#' it would cost 2p of the dearest call the family has.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param approx,nsim Passed to \code{\link{distrib_expected_hessian}}.
#' @param h_rel The relative step, a cube root of machine epsilon by default,
#'   which is what a central difference balances.
#'
#' @return A named list keyed as \code{\link{dexpected_names}}.
#'
#' @seealso \code{\link{distrib_dexpected_hessian}}, \code{\link{fd_steps}}
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
  for (k in seq_along(params)) {
    p <- params[k]
    # the point the difference is taken at, and the two points it is taken
    # between, on whichever scale the derivative was asked for
    if (link) {
      lk <- distrib@link_params[[p]]
      base <- linkfunctions7::linkfun(lk, theta[[p]])
      h <- h_rel * pmax(1, abs(base))
      up <- linkfunctions7::linkinv(lk, base + h)
      dn <- linkfunctions7::linkinv(lk, base - h)
    } else {
      h <- fd_steps(theta[[p]], distrib@params_bounds[[p]], h_rel)
      up <- theta[[p]] + h
      dn <- theta[[p]] - h
    }
    tu <- td <- theta
    tu[[p]] <- up
    td[[p]] <- dn
    a <- distrib_expected_hessian(distrib, y, tu, scale = scale,
                                  approx = approx, nsim = nsim)
    b <- distrib_expected_hessian(distrib, y, td, scale = scale,
                                  approx = approx, nsim = nsim)
    for (nm in hn) {
      out[[paste0(nm, "_", p)]] <- (a[[nm]] - b[[nm]]) / (2 * h)
    }
  }
  out
}
