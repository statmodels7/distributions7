#' @include distrib.R generics.R utility_functions.R
NULL

#' S7 Class for Multivariate Distributions
#'
#' @description
#' A subclass of \code{distrib} for distributions whose observations are
#' vectors. The response is an \eqn{n \times p} matrix, one row per
#' observation, and \code{n_dim} records \eqn{p}.
#'
#' @details
#' The parameters of a multivariate distribution are still \strong{scalars} as
#' far as the rest of the package is concerned. A mean vector contributes
#' \eqn{p} of them and a covariance matrix contributes the free values of the
#' \pkg{covstructs7} structure that parametrises it, so \code{theta} remains
#' the named list of numbers every generic already understands, and the
#' derivative bookkeeping -- \code{\link{deriv_names}}, the Hessian keys, the
#' link scale, \code{\link{fit_distrib}} -- needs no special case. The
#' constraint on the matrix lives inside the structure rather than in a link,
#' which is why the links of a multivariate distribution are all the identity:
#' the free values are already unconstrained.
#'
#' This class is a sibling of \code{\link{continuous_distrib}} rather than a
#' subclass of it, and deliberately. The defaults registered there are
#' one-dimensional -- a cdf by quadrature over an interval, a quantile by root
#' finding, a generator by ratio-of-uniforms on a scalar density -- and none of
#' them means anything in \eqn{p} dimensions. Inheriting them would offer
#' answers that do not exist.
#'
#' Parameters that vary from observation to observation are not supported here.
#' A distribution whose parameters depend on covariates is a model, which is
#' the layer above this one.
#'
#' @inheritParams distrib
#' @param n_dim The dimension \eqn{p} of an observation.
#'
#' @return An object of class \code{multivariate_distrib}.
#'
#' @seealso \code{\link{mvgaussian_distrib}}, \code{\link{n_obs}}
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' S7::S7_inherits(d, multivariate_distrib)
#' c(n_dim = d@n_dim, n_params = d@n_params)
#'
#' @export
multivariate_distrib <- S7::new_class("multivariate_distrib",
  parent = distrib,
  properties = list(n_dim = S7::class_integer),
  validator = function(self) {
    if (length(self@n_dim) != 1L || self@n_dim < 1L) {
      return("@n_dim must be a single positive integer")
    }
    if (!identical(self@dimension, "multivariate")) {
      return("@dimension must be 'multivariate'")
    }
    NULL
  }
)


#' How Many Observations a Response Holds
#'
#' @description
#' The number of observations in \code{y}: its length for a univariate
#' distribution, and the number of rows for a multivariate one.
#'
#' @details
#' Every place that used to write \code{length(y)} goes through this instead,
#' because for a matrix response \code{length(y)} counts entries rather than
#' observations, and the recycling checks built on it would ask for parameters
#' of length \eqn{np}.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y The response.
#'
#' @return A single integer.
#'
#' @seealso \code{\link{multivariate_distrib}}
#'
#' @examples
#' n_obs(gaussian_distrib(), c(1, 2, 3))
#' n_obs(mvgaussian_distrib(2), matrix(0, 5, 2))
#'
#' @export
n_obs <- function(distrib, y) {
  if (S7::S7_inherits(distrib, multivariate_distrib)) {
    if (is.null(dim(y))) {
      return(if (length(y) == 0L) 0L else 1L)
    }
    return(nrow(y))
  }
  length(y)
}


#' Coerce a Multivariate Response to a Matrix
#'
#' @description
#' Puts \code{y} in the \eqn{n \times p} form every multivariate method expects,
#' and refuses a response of the wrong width.
#'
#' @details
#' A plain vector of the right length is read as a single observation, which is
#' what makes \code{distrib_pdf(d, c(0, 0), theta)} mean what a reader expects
#' for a two-dimensional distribution. Anything else must already be a matrix
#' with one column per coordinate: guessing at the orientation of an
#' \eqn{n \times p} matrix would silently transpose a square sample.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param y The response.
#'
#' @return A numeric matrix with \code{distrib@n_dim} columns.
#'
#' @keywords internal
as_mv_matrix <- function(distrib, y) {
  p <- distrib@n_dim
  if (is.null(dim(y))) {
    if (length(y) == 0L) return(matrix(numeric(0), 0L, p))
    if (length(y) != p) {
      stop(sprintf(paste0(
        "'y' is a vector of length %d, which is read as one observation, but\n",
        "  '%s' has dimension %d. Supply a matrix with %d columns for several\n",
        "  observations."
      ), length(y), distrib@distrib_name, p, p), call. = FALSE)
    }
    return(matrix(y, nrow = 1L))
  }
  y <- as.matrix(y)
  if (ncol(y) != p) {
    stop(sprintf(
      "'y' has %d column(s) but '%s' has dimension %d.",
      ncol(y), distrib@distrib_name, p
    ), call. = FALSE)
  }
  if (!is.numeric(y)) stop("'y' must be numeric.", call. = FALSE)
  y
}


#' Require Scalar Parameters
#'
#' @description
#' Refuses a \code{theta} whose components are not single numbers.
#'
#' @details
#' The multivariate families of this package take one parameter value for the
#' whole sample. Vectorised parameters are what a regression supplies, and a
#' distribution that accepted them would be doing the model layer's work with
#' none of its bookkeeping.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#'
#' @return A numeric vector of the parameter values, in declaration order.
#'
#' @keywords internal
mv_flat_theta <- function(distrib, theta) {
  v <- theta[seq_len(distrib@n_params)]
  lens <- lengths(v)
  if (any(lens != 1L)) {
    bad <- distrib@params[which(lens != 1L)[1L]]
    stop(sprintf(paste0(
      "Parameter '%s' has length %d. A multivariate distribution takes one\n",
      "  value per parameter for the whole sample; parameters that vary by\n",
      "  observation belong to a model."
    ), bad, lens[which(lens != 1L)[1L]]), call. = FALSE)
  }
  stats::setNames(as.numeric(unlist(v)), distrib@params)
}


#' Refuse a Quantity That Has No Multivariate Counterpart
#'
#' @description
#' Raises the error a multivariate distribution gives for the one-dimensional
#' quantities: the distribution function, the quantile function and the
#' expectation by quadrature.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param what The name of the quantity.
#' @param why One sentence saying what is missing.
#'
#' @return Never returns; raises an error.
#'
#' @keywords internal
mv_refuse <- function(distrib, what, why) {
  stop(sprintf(
    "%s() is not defined for '%s': %s", what, distrib@distrib_name, why
  ), call. = FALSE)
}


#' @title No Distribution Function in Several Dimensions
#' @name distrib_cdf.multivariate_distrib
#' @description
#' Refused. The distribution function of a multivariate law is an integral over
#' an orthant, which has no closed form for the gaussian and no
#' one-dimensional fallback to stand in for it, and the quadrature registered
#' on \code{\link{continuous_distrib}} integrates over an interval.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param q A numeric matrix of quantiles.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(distrib_cdf, multivariate_distrib) <- function(distrib, q, theta, ...) {
  mv_refuse(
    distrib, "distrib_cdf",
    "the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback."
  )
}

#' @title No Quantile Function in Several Dimensions
#' @name distrib_quantile.multivariate_distrib
#' @description
#' Refused. A quantile is defined by inverting a distribution function on the
#' line; in several dimensions the ordering that would define it does not
#' exist.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(distrib_quantile, multivariate_distrib) <- function(distrib, p, theta, ...) {
  mv_refuse(
    distrib, "distrib_quantile",
    "a quantile inverts an ordering of the line, and there is none in several dimensions."
  )
}


#' @title Expected Information of a Multivariate Distribution
#' @name distrib_expected_hessian.multivariate_distrib
#'
#' @description
#' Fallback for a multivariate family with no closed form: the expectation is
#' taken over draws from the distribution itself.
#'
#' @details
#' The one-dimensional routes do not survive the move to \eqn{p} dimensions.
#' \code{"integrate"} builds its quadrature over an interval and is refused
#' here; \code{"bartlett"} in the univariate package reaches
#' \code{\link{expectation}}, which is that same quadrature. What does
#' generalise is sampling, so both remaining routes draw from the family's own
#' generator and differ in what they average:
#'
#' \code{"mc"} averages the observed Hessian, \eqn{\mathbb{E}[\ell^{(ij)}]}
#' directly. \code{"bartlett"} and \code{"opg"} average the outer product of
#' the score and negate it, which is the second Bartlett identity
#' \eqn{\mathcal{I} = \mathbb{E}[s s^\top]}; it needs no second derivative at
#' all, and is the only route that survives a family whose observed Hessian is
#' degenerate.
#'
#' Both are Monte Carlo, so both carry an error of order
#' \eqn{1/\sqrt{\texttt{nsim}}}, and a fit that uses one is doing Fisher
#' scoring with a noisy information. That is a deliberate choice a caller
#' makes, which is why \code{\link{fit_distrib}} refuses the argument for a
#' family that has an exact expression.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param y An \eqn{n \times p} matrix of observations; only its row count is
#'   used, the expectation being over the distribution rather than the data.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param approx One of \code{"bartlett"} (equivalently \code{"opg"}) or
#'   \code{"mc"}; \code{"integrate"} is refused.
#' @param nsim Monte Carlo sample size.
#' @param ... Unused.
#'
#' @return A named list keyed as \code{\link{hess_names}(distrib@params)}.
#' @keywords internal
S7::method(distrib_expected_hessian, multivariate_distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  approx <- match.arg(approx)
  if (identical(approx, "integrate")) {
    mv_refuse(
      distrib, "distrib_expected_hessian",
      'approx = "integrate" builds a quadrature over an interval, which has no counterpart in several dimensions. Use "mc" or "bartlett".'
    )
  }

  n <- n_obs(distrib, y)
  params <- distrib@params
  nm <- hess_names(params)
  big <- distrib_rng(distrib, nsim, theta)

  vals <- if (identical(approx, "mc")) {
    h <- distrib_hessian(distrib, big, theta)
    vapply(nm, function(k) mean(h[[k]]), numeric(1))
  } else {
    s <- do.call(cbind, distrib_gradient(distrib, big, theta))
    pr <- hess_pairs(params)
    vapply(nm, function(k) {
      ij <- match(pr[[k]], params)
      -mean(s[, ij[1L]] * s[, ij[2L]])
    }, numeric(1))
  }

  stats::setNames(lapply(vals, function(v) rep(v, n)), nm)
}

#' Prefix a Structure's Free Names with the Matrix They Describe
#'
#' @description
#' Returns the structure's free names with \code{"sigma_"} or \code{"omega_"}
#' in front, according to which side of the model the structure parametrises.
#'
#' @details
#' The name of a free value says how the matrix is built, not which matrix it
#' is, so a covariance structure and a precision structure of the same family
#' produce identical names. They are different models --- the inverse of a
#' compound-symmetry matrix is compound symmetry while the inverse of an AR(1)
#' is not AR(1) --- and a printed table that does not distinguish them leaves
#' the reader to guess. The prefix is applied by the distribution rather than
#' by the structure, because the structure does not know which side it has been
#' handed to.
#'
#' @param free_names The structure's free names.
#' @param inverted Whether the structure parametrises the precision.
#'
#' @return A character vector.
#'
#' @keywords internal
mv_prefixed_names <- function(free_names, inverted = FALSE) {
  paste0(if (isTRUE(inverted)) "omega_" else "sigma_", free_names)
}
