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

#' A Second-Order Mixed Derivative Through a Map
#'
#' @description
#' The parent's paired components carried onto a new parametrization by the
#' second-order chain rule.
#'
#' @param distrib The distribution in the new parametrization.
#' @param parent The parent distribution.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param y A numeric vector of observations.
#' @param first The parent's first-order components, keyed by parameter.
#' @param second The parent's second-order components, keyed by pair.
#'
#' @return A named list keyed by the new parameters' pairs.
#'
#' @seealso \code{\link{mapped_cross_y}}, \code{\link{reparam_tables}}
#' @keywords internal
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
#' @description
#' The parent's paired components carried by the second-order chain rule on
#' the map, and its first-order ones by the map's second partials.
#' @param distrib A reparametrized distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the new parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list keyed by parameter pair.
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
#' @description
#' The parent's block carried by the first-order chain rule, exactly as
#' \code{\link{distrib_cross_y}} is: a derivative in the response does not
#' interact with a reparametrization of the parameters.
#' @param distrib A reparametrized distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the new parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_cross2_y, ReparamContinuousDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    parent <- distrib@parent_distrib
    th <- reparam_theta(distrib, theta)
    mapped_cross2_y(distrib, parent, th, reparam_tables(distrib, theta), y)
  }


#' Second-Response Mixed Derivatives Through a Map
#'
#' @description
#' The same first-order chain rule \code{\link{mapped_cross_y}} takes, on the
#' second response derivative.
#'
#' @param distrib The distribution in the new parametrization.
#' @param parent The parent distribution.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param y A numeric vector of observations.
#'
#' @return A named list with one numeric vector per new parameter.
#'
#' @seealso \code{\link{mapped_cross_y}}
#' @keywords internal
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

#' @name distrib_grad_y_hess.LaplaceDistrib
#' @rdname loc_scale_grad_y_hess
#' @keywords internal
S7::method(distrib_grad_y_hess, LaplaceDistrib) <- loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, LaplaceDistrib) <- loc_scale_hess_y_hess
