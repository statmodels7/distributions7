#' @include theta2_families.R cross_derivatives_families.R reparametrize.R
NULL

# The mixed grid carried through a reparametrization.
#
# A derivative in the response does not see the parameters, so a
# reparametrization only ever acts on the theta side. At first order that is
# the chain rule mapped_cross_y already takes; at second order it is the
# ordinary two-term expansion,
#
#   d2 f / da db = sum_ij (dth_i/da)(dth_j/db) f_ij + sum_i (d2th_i/da db) f_i
#
# with f_i the parent's first-order mixed component and f_ij its second-order
# one. Both tables are the map's own, which reparam_tables() already keys by
# the indices of the new parameters -- "1" for a first partial and "1,2" for a
# second -- so nothing is differentiated here either.

#' @title A Second-Order Mixed Derivative Through a Map
#'
#' @description
#' Carries the parent's paired components onto a new parametrization by the
#' second-order chain rule. It is the arithmetic behind every
#' [distrib_grad_y_hess()] and [distrib_hess_y_hess()] method of a
#' reparametrized family, and it differentiates nothing itself.
#'
#' @details
#' # The expansion
#'
#' A derivative in the RESPONSE does not see the parameters, so a
#' reparametrization acts on the \eqn{\theta} side alone and the ordinary
#' two-term expansion applies:
#' \deqn{\frac{\partial^2 f}{\partial\alpha_a \partial\alpha_b}
#'   = \sum_{i,j} \frac{\partial\theta_i}{\partial\alpha_a}
#'     \frac{\partial\theta_j}{\partial\alpha_b} f_{ij}
#'   + \sum_i \frac{\partial^2\theta_i}{\partial\alpha_a\partial\alpha_b} f_i,}
#' with \eqn{f_i} the parent's first-order mixed component and \eqn{f_{ij}} its
#' second-order one. Both partial tables are the MAP's, already keyed by
#' [reparam_tables()] under the indices of the new parameters, `"1"` for a
#' first partial and `"1,2"` for a second, so nothing is differentiated here.
#'
#' # What `first` and `second` are
#'
#' They are the same quantity at two orders in \eqn{\theta}, and which
#' quantity depends on the caller: [distrib_cross_y()] with
#' [distrib_grad_y_hess()] for the third-order derivative, and
#' [distrib_cross2_y()] with [distrib_hess_y_hess()] for the fourth.
#'
#' A missing key in a map's table is an exact zero and is skipped, so a map
#' that is affine in one of its coordinates costs nothing for it.
#'
#' @param distrib The distribution in the NEW parametrization. Only its
#'   `params` are read, to key the result.
#' @param parent The parent distribution, whose `params` key the inputs.
#' @param th_par A named list of the parent's parameters, evaluated at the new
#'   ones.
#' @param maps The map's keyed partial tables, from [reparam_tables()]: one
#'   entry per parent parameter, each a list keyed by new-parameter index.
#' @param y A numeric vector of observations. Used only for its length,
#'   through the components.
#' @param first The parent's first-order components, keyed by the parent's
#'   parameters.
#' @param second The parent's second-order components, keyed by the parent's
#'   parameter pairs. Either ordering of a pair's key is accepted.
#'
#' @return A named list with one numeric vector per unordered pair of the NEW
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\theta_i} a parameter
#' of the PARENT and \eqn{\alpha_a} one of the new parametrization, so that
#' \eqn{\theta = \theta(\alpha)} is the map. \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} are the first and second derivatives in the response.
#'
#' @seealso [mapped_cross2_y()] for the first-order counterpart,
#'   [reparam_tables()] for the partial tables, and
#'   [distrib_grad_y_hess.ReparamContinuousDistrib()], its caller.
#'
#' @keywords internal
#'
#' @examples
#' # gaussian2 carries (mu, sigma2) where its parent carries (mu, sigma).
#' d <- gaussian2_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma2 = 1.44)
#'
#' g <- distrib_grad_y_hess(d, y, theta)
#' vapply(g, function(z) z[1], numeric(1))
#'
#' # Against a numerical Hessian of the response gradient in the NEW
#' # parameters, which shares none of the map's arithmetic.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(0.3, 1.44))
mapped_theta2 <- function(distrib, parent, th_par, maps, y, first, second) {
  old <- parent@params
  new_params <- distrib@params
  prs_old <- hess_pairs(old)
  prs_new <- hess_pairs(new_params)
  zero <- 0 * first[[1L]]
  key_old <- function(a, b) {
    k <- paste(old[a], old[b], sep = "_")
    if (k %in% names(second)) k else paste(old[b], old[a], sep = "_")
  }
  stats::setNames(lapply(names(prs_new), function(nm) {
    ij <- match(prs_new[[nm]], new_params)
    a <- ij[1L]
    b <- ij[2L]
    s <- zero
    for (k in seq_along(old)) {
      # the second partial of the map, keyed by the new indices in order
      key <- paste(sort(c(a, b)), collapse = ",")
      v2 <- maps[[k]][[key]]
      if (!is.null(v2)) s <- s + first[[k]] * v2
      da <- maps[[k]][[as.character(a)]]
      if (is.null(da)) next
      for (l in seq_along(old)) {
        db <- maps[[l]][[as.character(b)]]
        if (is.null(db)) next
        s <- s + second[[key_old(k, l)]] * da * db
      }
    }
    s
  }), names(prs_new))
}


#' @title Second-Order Mixed Derivatives of a Reparametrized Distribution
#' @name distrib_grad_y_hess.ReparamContinuousDistrib
#'
#' @description
#' Carries the parent's paired components onto the new parametrization by the
#' second-order chain rule, through [mapped_theta2()]: the parent's
#' second-order components multiplied by two first partials of the map, plus
#' its first-order components multiplied by the map's second partials. The
#' fourth-order method is the same body reading [distrib_cross2_y()] and
#' [distrib_hess_y_hess()] of the parent instead.
#'
#' @details
#' A response derivative does not interact with a reparametrization of
#' \eqn{\theta}, so nothing on the \eqn{y} side is touched and the whole
#' correction is the map's. The partials come from [reparam_tables()], which a
#' family supplies through `map_derivs` or, failing that, through one stencil
#' per partial.
#'
#' @param distrib A `ReparamContinuousDistrib` object, from [reparametrize()].
#' @param y A numeric vector of observations.
#' @param theta A named list of the NEW parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per unordered pair of the new
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\theta_i} a parameter
#' of the PARENT and \eqn{\alpha_a} one of the new parametrization, so that
#' \eqn{\theta = \theta(\alpha)} is the map. \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} are the first and second derivatives in the response.
#'
#' @seealso [mapped_theta2()], which does the work,
#'   [distrib_cross2_y.ReparamContinuousDistrib()] for the first order in
#'   \eqn{\theta}, and [reparametrize()] for the wrapper.
#'
#' @examples
#' # A gaussian obtained in its variance rather than written out.
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma2 = 1.44)
#'
#' g3 <- distrib_grad_y_hess(d, y, theta)
#' vapply(g3, function(z) z[1], numeric(1))
#'
#' # Against the family written out by hand, which shares none of the map's
#' # arithmetic. The map's partials are stencils here, no map_derivs having
#' # been given, so the two agree to about 1e-09.
#' max(abs(unlist(g3) -
#'         unlist(distrib_grad_y_hess(gaussian2_distrib(), y, theta))))
#'
#' # And against a numerical Hessian of the response gradient.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(0.3, 1.44))
#'
#' # The fourth order takes the same route one derivative up in y.
#' g4 <- distrib_hess_y_hess(d, y, theta)
#' max(abs(unlist(g4) -
#'         unlist(distrib_hess_y_hess(gaussian2_distrib(), y, theta))))
#'
#' @keywords internal
S7::method(distrib_grad_y_hess, ReparamContinuousDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    parent <- distrib@parent_distrib
    th <- reparam_theta(distrib, theta)
    maps <- reparam_tables(distrib, theta)
    mapped_theta2(distrib, parent, th, maps, y,
                  distrib_cross_y(parent, y, th),
                  distrib_grad_y_hess(parent, y, th))
  }

#' @rdname distrib_grad_y_hess.ReparamContinuousDistrib
#' @name distrib_hess_y_hess.ReparamContinuousDistrib
#' @keywords internal
S7::method(distrib_hess_y_hess, ReparamContinuousDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    parent <- distrib@parent_distrib
    th <- reparam_theta(distrib, theta)
    maps <- reparam_tables(distrib, theta)
    mapped_theta2(distrib, parent, th, maps, y,
                  distrib_cross2_y(parent, y, th),
                  distrib_hess_y_hess(parent, y, th))
  }

#' @title Second-Response Mixed Derivatives of a Reparametrized Distribution
#' @name distrib_cross2_y.ReparamContinuousDistrib
#'
#' @description
#' Carries the parent's block onto the new parametrization by the FIRST-order
#' chain rule, exactly as [distrib_cross_y()] is carried: a derivative in the
#' response does not interact with a reparametrization of the parameters, so
#' only the map's first partials enter and its second ones never appear.
#'
#' @param distrib A `ReparamContinuousDistrib` object, from [reparametrize()].
#' @param y A numeric vector of observations.
#' @param theta A named list of the NEW parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per new parameter, keyed by
#'   `distrib@params`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\theta_i} a parameter
#' of the PARENT and \eqn{\alpha_a} one of the new parametrization, so that
#' \eqn{\theta = \theta(\alpha)} is the map. \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} are the first and second derivatives in the response.
#'
#' @seealso [mapped_cross2_y()], which does the work,
#'   [distrib_grad_y_hess.ReparamContinuousDistrib()] for the second order in
#'   \eqn{\theta}, and [reparametrize()] for the wrapper.
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma2 = 1.44)
#'
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#'
#' # Against the family written out by hand, and against a numerical
#' # derivative of the response Hessian.
#' max(abs(unlist(distrib_cross2_y(d, y, theta)) -
#'         unlist(distrib_cross2_y(gaussian2_distrib(), y, theta))))
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::grad(f, c(0.3, 1.44))
#'
#' @keywords internal
S7::method(distrib_cross2_y, ReparamContinuousDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    parent <- distrib@parent_distrib
    th <- reparam_theta(distrib, theta)
    mapped_cross2_y(distrib, parent, th, reparam_tables(distrib, theta), y)
  }


#' @title Second-Response Mixed Derivatives Through a Map
#'
#' @description
#' Carries the parent's [distrib_cross2_y()] block onto a new parametrization
#' by the first-order chain rule,
#' \eqn{\sum_i (\partial\theta_i/\partial\alpha_a)\,
#' \partial\ell^{(yy)}/\partial\theta_i}. It is the same expansion
#' [mapped_cross_y()] takes, read one derivative further in the response.
#'
#' @details
#' Only the map's FIRST partials enter, so the second-partial tables
#' [mapped_theta2()] needs are never touched. A missing key is an exact zero
#' and is skipped.
#'
#' @param distrib The distribution in the new parametrization. Only its
#'   `params` are read, to key the result.
#' @param parent The parent distribution, which supplies the block.
#' @param th_par A named list of the parent's parameters, evaluated at the new
#'   ones.
#' @param maps The map's keyed partial tables, from [reparam_tables()].
#' @param y A numeric vector of observations.
#'
#' @return A named list with one numeric vector per new parameter, keyed by
#'   `distrib@params`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\theta_i} a parameter
#' of the PARENT and \eqn{\alpha_a} one of the new parametrization, so that
#' \eqn{\theta = \theta(\alpha)} is the map. \eqn{\ell^{(y)}} and
#' \eqn{\ell^{(yy)}} are the first and second derivatives in the response.
#'
#' @seealso [mapped_theta2()] for the second order in \eqn{\theta},
#'   [distrib_cross2_y.ReparamContinuousDistrib()], its caller, and
#'   [reparam_tables()] for the partial tables.
#'
#' @keywords internal
#'
#' @examples
#' # gaussian2 carries (mu, sigma2) against a parent in (mu, sigma), so the
#' # scale component picks up d(sigma)/d(sigma2) = 1 / (2 sqrt(sigma2)).
#' d <- gaussian2_distrib()
#' y <- c(-0.7, 0.3, 1.4)
#' theta <- list(mu = 0.3, sigma2 = 1.44)
#'
#' vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#'
#' # The parent's own block, carried by that one factor.
#' par_block <- distrib_cross2_y(gaussian1_distrib(), y,
#'                               list(mu = 0.3, sigma = 1.2))
#' c(mu = par_block$mu[1], sigma2 = par_block$sigma[1] / (2 * sqrt(1.44)))
mapped_cross2_y <- function(distrib, parent, th_par, maps, y) {
  cy <- distrib_cross2_y(parent, y, th_par)
  new_params <- distrib@params
  zero <- 0 * cy[[1L]]
  out <- lapply(seq_along(new_params), function(i) {
    s <- zero
    for (k in seq_along(cy)) {
      v <- maps[[k]][[as.character(i)]]
      if (!is.null(v)) s <- s + cy[[k]] * v
    }
    s
  })
  stats::setNames(out, new_params)
}


# --- the Laplace, which was registered at first order and not at second -----

# Both aliases are carried by loc_scale_grad_y_hess's own page, in
# theta2_families.R, where the shared body is documented.
S7::method(distrib_grad_y_hess, LaplaceDistrib) <- loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, LaplaceDistrib) <- loc_scale_hess_y_hess
