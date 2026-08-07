#' @include distrib.R generics.R utility_functions.R cross_derivatives.R moments.R
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
#' \pkg{parameters7} structure that parametrizes it, so \code{theta} remains
#' the named list of numbers every generic already understands, and the
#' derivative bookkeeping -- \code{\link{deriv_names}}, the Hessian keys, the
#' link scale, \code{\link{fit_distrib}} -- needs no special case. The
#' constraint on the matrix lives inside the matrix parameter rather than in a link,
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
#' n_obs(gaussian1_distrib(), c(1, 2, 3))
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
#' and rejects a response of the wrong width.
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
#' Rejects a \code{theta} whose components are not single numbers.
#'
#' @details
#' The multivariate families of this package take one parameter value for the
#' whole sample. Vectorized parameters are what a regression supplies, and a
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
#' Rejected. The distribution function of a multivariate law is an integral over
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
#' Rejected. A quantile is defined by inverting a distribution function on the
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

#' @title Response Derivatives of a Multivariate Distribution
#' @name distrib_grad_y.multivariate_distrib
#' @description
#' Rejected on the base class rather than served numerically: the univariate
#' fallbacks difference along a line, and the derivative of a multivariate
#' log-density in its response is a vector (a matrix at second order) whose
#' shape the base class cannot guess. A family that has the closed form
#' registers it, as the gaussian and the Student t do.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return No return value; called for its error.
#' @keywords internal
S7::method(distrib_grad_y, multivariate_distrib) <- function(distrib, y, theta, ...) {
  mv_refuse(
    distrib, "distrib_grad_y",
    "the univariate numerical fallback differences along a line; register a closed form on the family."
  )
}

#' @rdname distrib_grad_y.multivariate_distrib
#' @name distrib_hess_y.multivariate_distrib
S7::method(distrib_hess_y, multivariate_distrib) <- function(distrib, y, theta, ...) {
  mv_refuse(
    distrib, "distrib_hess_y",
    "the univariate numerical fallback differences along a line; register a closed form on the family."
  )
}

#' @rdname distrib_grad_y.multivariate_distrib
#' @name distrib_cross_y.multivariate_distrib
S7::method(distrib_cross_y, multivariate_distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mv_refuse(
    distrib, "distrib_cross_y",
    "the mixed response-parameter block of a multivariate family is a matrix per observation, and no consumer fixes its shape yet."
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
#' \code{"integrate"} builds its quadrature over an interval and is rejected
#' here; \code{"bartlett"} in the univariate package reaches
#' \code{\link{expectation}}, which is that same quadrature. What does
#' generalize is sampling, so both remaining routes draw from the family's own
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
#' makes, which is why \code{\link{fit_distrib}} rejects the argument for a
#' family that has an exact expression.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param y An \eqn{n \times p} matrix of observations; only its row count is
#'   used, the expectation being over the distribution rather than the data.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param approx One of \code{"bartlett"} (equivalently \code{"opg"}) or
#'   \code{"mc"}; \code{"integrate"} is rejected.
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
#' Returns the matrix parameter's free names with \code{"sigma_"} or \code{"omega_"}
#' in front, according to which side of the model the matrix parameter parametrizes.
#'
#' @details
#' The name of a free value says how the matrix is built, not which matrix it
#' is, so a covariance structure and a precision structure of the same family
#' produce identical names. They are different models --- the inverse of a
#' compound-symmetry matrix is compound symmetry while the inverse of an AR(1)
#' is not AR(1) --- and a printed table that does not distinguish them leaves
#' the reader to guess. The prefix is applied by the distribution rather than
#' by the matrix parameter, because the matrix parameter does not know which side it has been
#' handed to.
#'
#' @param free_names The matrix parameter's free names.
#' @param inverted Whether the matrix parameter parametrizes the precision.
#'
#' @return A character vector.
#'
#' @keywords internal
mv_prefixed_names <- function(free_names, inverted = FALSE) {
  paste0(if (isTRUE(inverted)) "omega_" else "sigma_", free_names)
}


#' The Support Points of a Discrete Multivariate Distribution
#'
#' @description
#' The points a discrete multivariate distribution places mass on, as a matrix
#' with one row per point.
#'
#' @details
#' A univariate discrete distribution needs no such generic: its support is a
#' stretch of the integers and the package walks it. On a vector the support is
#' a set whose shape depends on the family --- the multinomial's is the
#' compositions of its size --- and enumerating it is what lets an expectation
#' be an exact sum and the validator check the total mass by addition rather
#' than by sampling.
#'
#' The base class rejects. A continuous family has no such set, and a discrete
#' one whose support is infinite has no finite matrix to return; either way an
#' answer would be a fiction, and the caller is better told.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters. Families whose support does not
#'   depend on them ignore it.
#' @param ... Passed to methods.
#'
#' @return A matrix with one row per support point and one column per
#'   coordinate.
#'
#' @seealso \code{\link{multinomial_distrib}}, \code{\link[numericals7]{compositions}}
#'
#' @examples
#' d <- multinomial_distrib(3, size = 4)
#' nrow(mv_support(d, list(probs_alr1 = 0, probs_alr2 = 0)))
#'
#' @export
mv_support <- S7::new_generic("mv_support", "distrib",
  function(distrib, theta, ...) S7::S7_dispatch()
)

#' @title No Enumerable Support
#' @name mv_support.multivariate_distrib
#' @description The base-class method, which rejects.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @seealso \code{\link{mv_support}}
#' @keywords internal
S7::method(mv_support, multivariate_distrib) <- function(distrib, theta, ...) {
  stop(sprintf(paste0(
    "'%s' does not enumerate a support. A continuous family has no such set,\n",
    "  and a discrete one whose support is infinite has no finite matrix to\n",
    "  return; a family that does have one registers this generic."
  ), distrib@distrib_name), call. = FALSE)
}


#' A Proposal for Integrating a Multivariate Density
#'
#' @description
#' Draws from a distribution that dominates the family, together with the
#' log-density of those draws, so that an importance-sampling estimate of
#' \eqn{\int f = 1} can be formed.
#'
#' @details
#' The default proposal is a gaussian with the same mean and twice the
#' covariance, which serves any family supported on all of \eqn{\mathbb{R}^p}.
#' The inflation matters: a proposal equal to the distribution itself makes
#' every ratio one and certifies nothing.
#'
#' A family whose support is a lower-dimensional subset of \eqn{\mathbb{R}^p},
#' such as a Dirichlet on the simplex, must register its own method, because a
#' proposal spread over the ambient space places no mass at all on the support
#' and the estimate would be zero. The draws and the log-density must be taken
#' with respect to the same dominating measure the family's density is written
#' against.
#'
#' This is consumed by \code{\link{check_distrib}} and by nothing else. A
#' discrete family does not need it, its normalization being an exact sum over
#' \code{\link{mv_support}}.
#'
#' @param distrib An object inheriting from class
#'   \code{\link{multivariate_distrib}}.
#' @param theta A named list or vector of parameters.
#' @param n The number of draws.
#' @param ... Passed to methods.
#'
#' @return A list with \code{y}, a matrix of \code{n} draws, and \code{logd},
#'   their log-density under the proposal.
#'
#' @seealso \code{\link{check_distrib}}, \code{\link{mv_support}}
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)
#' set.seed(1)
#' str(mv_reference_draw(d, theta, 5))
#'
#' @export
mv_reference_draw <- S7::new_generic("mv_reference_draw", "distrib",
  function(distrib, theta, n, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  }
)

#' @title An Inflated Gaussian Proposal
#' @name mv_reference_draw.multivariate_distrib
#' @description
#' The default: a gaussian with the family's mean and twice its covariance.
#' Requires the covariance to be non-singular, which is what a support of full
#' dimension gives.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#' @param n The number of draws.
#' @param ... Unused.
#' @return A list with the draws \code{y} and their log-density \code{logd}.
#' @seealso \code{\link{mv_reference_draw}}
#' @keywords internal
S7::method(mv_reference_draw, multivariate_distrib) <- function(distrib, theta, n, ...) {
  p <- distrib@n_dim
  mu <- as.numeric(mean(distrib, theta))
  sg <- as.matrix(variance(distrib, theta)) * 2
  l <- t(chol(sg))
  z <- matrix(stats::rnorm(n * p), n, p)
  y <- sweep(z %*% t(l), 2L, mu, "+")
  # The whitened residuals are L^{-1} r, so the system is LOWER triangular and
  # forwardsolve is what solves it; backsolve on the transpose solves L' x = b,
  # which is a different vector with the same shape.
  q <- colSums(forwardsolve(l, t(sweep(y, 2L, mu, "-")))^2)
  list(y = y, logd = -0.5 * (p * log(2 * pi) + 2 * sum(log(diag(l))) + q))
}




#' The Mean Vector and Covariance a Parameter List Describes
#'
#' @description
#' Assembles the mean vector and the covariance matrix of a multivariate
#' distribution from its flat parameter list.
#'
#' @details
#' The parameters of a multivariate distribution are scalars, so that every
#' generic of the package can index them, and these two functions put them back
#' into the shapes a reader thinks in. \code{mv_sigma()} returns the matrix the
#' PARAMETRIZATION carries, whichever side the matrix parameter describes: the
#' covariance for a gaussian, and the scale matrix for a Student t, whose
#' covariance is \eqn{\nu\Sigma/(\nu-2)} and does not exist below two degrees
#' of freedom. The moment is \code{\link{variance}}, and keeping the two apart
#' is what lets a heavy-tailed family be described at all.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list or vector of parameters.
#'
#' Both are generics whose base-class method rejects: not every multivariate
#' family has a location, and one that does not should say so rather than
#' hand back its first p parameters under a name that does not fit them.
#'
#' @return A numeric vector of length \eqn{p} for \code{mv_location()}, and a
#'   \eqn{p \times p} matrix for \code{mv_sigma()}.
#'
#' @seealso \code{\link{mvgaussian_distrib}}
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)
#' mv_location(d, theta)
#' mv_sigma(d, theta)
#'
#' @export
mv_location <- S7::new_generic("mv_location", "distrib", function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' @title No Location Without a Family That Has One
#' @name mv_location.multivariate_distrib
#' @description
#' Rejected. Not every multivariate family has a location: a Dirichlet is
#' described by concentrations and a Wishart by a scale and a count, and
#' handing back the first \eqn{p} parameters under the name of a mean would be
#' a wrong answer in the shape of a right one.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(mv_location, multivariate_distrib) <- function(distrib, theta) {
  mv_refuse(
    distrib, "mv_location",
    "this family has no location parameter. A family that has one registers a method."
  )
}

#' The First p Parameters, Read as a Location
#'
#' @description
#' The helper the elliptical families implement \code{\link{mv_location}} with:
#' the first \eqn{p} entries of the flat parameter vector, labeled by
#' coordinate.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters, already aligned.
#'
#' @return A named numeric vector of length \eqn{p}.
#'
#' @keywords internal
mv_leading_location <- function(distrib, theta) {
  v <- mv_flat_theta(distrib, theta)
  stats::setNames(
    unname(v[seq_len(distrib@n_dim)]), paste0("v", seq_len(distrib@n_dim))
  )
}

#' @rdname mv_location
#' @export
mv_sigma <- S7::new_generic("mv_sigma", "distrib", function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' A Marginal of a Multivariate Distribution
#'
#' @description
#' Returns the distribution of a subset of the coordinates, together with the
#' parameters that describe it.
#'
#' @details
#' A marginal is not available in general: integrating a density over the
#' coordinates one is not interested in has no closed form for most families.
#' It is available for the elliptical ones, where the marginal belongs to the
#' same family with the mean and the matrix subsetted, and those are the ones
#' this generic has methods for. For a family without one the generic signals
#' an error rather than approximating, since a quadrature over the discarded
#' coordinates would be a different object under the same name.
#'
#' This is what makes a picture of a multivariate distribution possible at all:
#' a panel of a pairs plot shows a marginal, so the plot exists exactly when
#' the marginals do.
#'
#' @param distrib An object inheriting from class
#'   \code{\link{multivariate_distrib}}.
#' @param theta A named list or vector of parameters.
#' @param which An integer vector of coordinates to keep.
#' @param ... Passed to methods.
#'
#' @return A list with \code{distrib}, the marginal distribution object, and
#'   \code{theta}, its parameters.
#'
#' @seealso \code{\link{plot.multivariate_distrib}}
#'
#' @examples
#' d <- mvgaussian_distrib(3)
#' theta <- as.list(stats::setNames(c(1, 2, 3, 0, 0, 0, 0.3, 0.2, 0.1), d@params))
#'
#' # the marginal of the first two coordinates is a two-dimensional gaussian
#' m <- mv_marginal(d, theta, c(1, 2))
#' mv_sigma(m$distrib, m$theta)
#'
#' # and it is the corresponding block of the full covariance
#' mv_sigma(d, theta)[1:2, 1:2]
#'
#' @export
mv_marginal <- S7::new_generic("mv_marginal", "distrib",
  function(distrib, theta, which, ...) {
    theta <- align_theta(distrib, theta)
    which <- as.integer(which)
    if (!length(which) || anyNA(which) ||
      any(which < 1L) || any(which > distrib@n_dim) || anyDuplicated(which)) {
      stop(sprintf(
        "'which' must be distinct coordinates in 1:%d.", distrib@n_dim
      ), call. = FALSE)
    }
    S7::S7_dispatch()
  }
)

#' @title No Marginal Without a Closed Form
#' @name mv_marginal.multivariate_distrib
#' @description
#' Rejected. Integrating out the other coordinates has no general closed form,
#' and a numerical marginal would be a different object with the same name.
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#' @param which An integer vector of coordinates.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @keywords internal
S7::method(mv_marginal, multivariate_distrib) <- function(distrib, theta, which, ...) {
  mv_refuse(
    distrib, "mv_marginal",
    "integrating out the other coordinates has no closed form for this family, and a numerical marginal would not be the same object."
  )
}


#' @title No Skewness Without Saying Which One
#' @name skewness.multivariate_distrib
#' @description
#' Rejected. A scalar skewness for a vector response is not one quantity but a
#' choice among several -- Mardia's, Malkovich-Afifi's, or the vector of
#' coordinatewise marginal skewnesses -- and they do not agree. Returning any
#' of them under the bare name would be a wrong answer in the shape of a right
#' one, so the caller names the quantity it wants instead. Note that
#' \code{\link{mv_marginal}} is not a way round this for an elliptical family,
#' whose marginal is a smaller multivariate distribution and rejects in turn;
#' it is for the Dirichlet and the multinomial, whose marginals are univariate.
#' @param x A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @seealso \code{\link{mv_marginal}}
#' @keywords internal
S7::method(skewness, multivariate_distrib) <- function(x, theta, ...) {
  mv_refuse(
    x, "skewness",
    "a vector response has no single skewness: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal skewnesses are different quantities. Ask a univariate family instead."
  )
}

#' @title No Kurtosis Without Saying Which One
#' @name kurtosis.multivariate_distrib
#' @description
#' Rejected, for the reason given at
#' \code{\link[=skewness.multivariate_distrib]{skewness()}}.
#' @param x A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @seealso \code{\link{mv_marginal}}
#' @keywords internal
S7::method(kurtosis, multivariate_distrib) <- function(x, theta, ...) {
  mv_refuse(
    x, "kurtosis",
    "a vector response has no single kurtosis: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal kurtosises are different quantities. Ask a univariate family instead."
  )
}
