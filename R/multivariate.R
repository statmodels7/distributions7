#' @include distrib.R generics.R utility_functions.R cross_derivatives.R moments.R
NULL

#' @title S7 Class for Multivariate Distributions
#'
#' @description
#' The subclass of `distrib` for distributions whose observations are vectors.
#' The response is an \eqn{n \times p} matrix, one row per observation, and
#' `n_dim` records \eqn{p}. It is a SIBLING of [continuous_distrib()] and
#' [discrete_distrib()], not a subclass of either, so it inherits none of their
#' one-dimensional defaults; the distribution function, the quantile function
#' and the response derivatives are refused on this class instead of being
#' approximated.
#'
#' @details
#' # The parameters are still scalars
#'
#' As far as the rest of the package is concerned, a multivariate
#' distribution's parameters are SCALARS. A mean vector contributes \eqn{p} of
#' them and a covariance matrix contributes the free values of the
#' \pkg{parameters7} parametrization that carries it, so `theta` remains the
#' named list of numbers every generic already understands, and the derivative
#' bookkeeping needs no special case: [deriv_names()], the Hessian keys, the
#' link scale and [fit_distrib()] all work unchanged.
#'
#' The constraint on the matrix lives inside the parametrization, not in a
#' link, so every link of a multivariate distribution is the identity: the free
#' values are already unconstrained.
#'
#' # Why it is a sibling of the one-dimensional classes
#'
#' The defaults registered on [continuous_distrib()] are one-dimensional: a
#' distribution function by quadrature over an interval, a quantile by root
#' finding, a generator by ratio-of-uniforms on a scalar density. None of them
#' means anything in \eqn{p} dimensions, and inheriting them would offer
#' answers that do not exist.
#'
#' # No parameter varies by observation
#'
#' A parametrization describes one matrix for the whole sample, and
#' [mv_flat_theta()] rejects a `theta` component longer than one. A
#' distribution whose parameters depend on covariates is a model, which is the
#' layer above this one.
#'
#' @param n_dim The dimension \eqn{p} of one observation. A single positive
#'   integer; the validator rejects anything else, and also rejects a
#'   `dimension` property other than `"multivariate"`.
#' @param ... The properties of the parent `distrib` class, listed under Value.
#'
#' @return An S7 object of class `multivariate_distrib`, inheriting from
#'   `distrib`. Beyond the parent's `distrib_name`, `dimension`, `bounds`,
#'   `params`, `params_interpretation`, `n_params`, `params_bounds`,
#'   `link_params` and `params_smooth`, it carries `n_dim`.
#'
#' @seealso [mvgaussian_distrib()] and [mvstudent_t_distrib()] for the
#'   elliptical families, [dirichlet_distrib()] and [multinomial_distrib()] for
#'   the simplex-valued ones, [n_obs()] for the observation count, and
#'   [mv_summary()] for the quantities a fit reports.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' S7::S7_inherits(d, multivariate_distrib)
#' c(n_dim = d@n_dim, n_params = d@n_params)
#'
#' # It is a sibling of continuous_distrib, not a subclass, so none of the
#' # one-dimensional defaults is inherited.
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # Every link is the identity: the constraint lives in the matrix
#' # parametrization, which needs no link to express it.
#' vapply(d@link_params, function(l) l@link_name, character(1))
#'
#' # And the four shipped families all sit here.
#' vapply(list(mvgaussian_distrib(2), mvstudent_t_distrib(2),
#'             dirichlet_distrib(3), multinomial_distrib(3, size = 5)),
#'        function(x) S7::S7_inherits(x, multivariate_distrib), TRUE)
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


#' @title How Many Observations a Response Holds
#'
#' @description
#' Returns the number of observations in `y`: its length for a univariate
#' distribution, and its number of ROWS for a multivariate one. Every place in
#' the package that would otherwise write `length(y)` goes through this,
#' because for a matrix response `length(y)` counts entries, and a recycling
#' check built on it would ask for parameters of length \eqn{np}.
#'
#' @details
#' A wrong length here is quiet. The multivariate gaussian's expected
#' information once built a zero vector with `length(y)` to stand for a
#' parameter that does not vary; that vector came out \eqn{np} long, recycled
#' against the \eqn{p}-long components, and inflated every diagonal entry of
#' the information by a factor of \eqn{p}, so every standard error of a
#' multivariate fit was \eqn{\sqrt{p}} too small. Nothing failed. Anywhere a
#' matrix response meets code written for a vector, `length()` is a defect and
#' this is the question to ask.
#'
#' @param distrib An object inheriting from `distrib`. Only its class is used,
#'   to decide which reading of `y` applies.
#' @param y The response: a numeric vector for a univariate distribution, or an
#'   \eqn{n \times p} numeric matrix for a multivariate one.
#'
#' @return A single integer.
#'
#' @seealso [multivariate_distrib()] and [as_mv_matrix()], which puts a
#'   response in the shape this counts.
#'
#' @examples
#' # A univariate response is counted by length.
#' n_obs(gaussian1_distrib(), c(1.2, -0.4, 0.8))
#'
#' # A multivariate one by rows, where length() would give the entry count.
#' y <- matrix(0, 5, 2)
#' c(n_obs = n_obs(mvgaussian_distrib(2), y), length = length(y))
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


#' @title Coerce a Multivariate Response to a Matrix
#'
#' @description
#' Puts `y` in the \eqn{n \times p} form every multivariate method expects. A
#' matrix of the right width is returned unchanged; a plain vector of length
#' \eqn{p} is read as a SINGLE observation and returned as a one-row matrix,
#' the reading a caller asking for a density at one point wants. Anything else
#' is an error naming the length it was given and the dimension it should have
#' had.
#'
#' @param distrib A [multivariate_distrib()] object, from which `n_dim` is
#'   read.
#' @param y A numeric matrix with `distrib@n_dim` columns, or a numeric vector
#'   of length `distrib@n_dim`. A vector of any other length is an error, as is
#'   a matrix of the wrong width.
#'
#' @return A numeric matrix with `distrib@n_dim` columns.
#'
#' @seealso [n_obs()] for the row count and [multivariate_distrib()] for the
#'   response convention.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#'
#' # A vector of length p is one observation.
#' dim(distributions7:::as_mv_matrix(d, c(1, -1)))
#'
#' # A matrix of the right width is returned as it stands.
#' dim(distributions7:::as_mv_matrix(d, matrix(0, 4, 2)))
#'
#' # Any other length is an error naming both numbers.
#' try(distributions7:::as_mv_matrix(d, c(1, 2, 3)))
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


#' @title Require Scalar Parameters, and Flatten Them
#'
#' @description
#' Rejects a `theta` whose components are not single numbers, and returns the
#' rest as one named numeric vector in `distrib@params` order. The families
#' here take ONE parameter value for the whole sample: a parametrization
#' describes one matrix, not one per row, and a parameter that varies by
#' observation belongs to a model, which is the layer above this one.
#'
#' @details
#' A univariate family reads a `theta` component of length \eqn{n} as one value
#' per observation, so without this check a caller who wrote that here would
#' get silent recycling against the \eqn{p} columns instead of an error. The
#' message names the offending parameter.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters, already aligned by
#'   [align_theta()]. Every component must have length 1.
#'
#' @return A named numeric vector of length `distrib@n_params`, named and
#'   ordered as `distrib@params`.
#'
#' @seealso [align_theta()], which normalizes the names first, and
#'   [multivariate_distrib()] for the convention.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' distributions7:::mv_flat_theta(d, distributions7:::align_theta(d, theta))
#'
#' # A component longer than one is rejected by name, where a univariate
#' # family would read it as one value per observation.
#' bad <- theta
#' bad$mu1 <- c(0, 1)
#' try(distributions7:::mv_flat_theta(d, bad))
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


#' @title Reject a Quantity That Has No Multivariate Counterpart
#'
#' @description
#' Raises the error a multivariate distribution gives for a one-dimensional
#' quantity, in one wording: the name of the quantity, the family it was asked
#' of, and a sentence saying what is missing. Every refusal on this class goes
#' through it, so the eight of them read alike and none of them silently
#' returns a number of the wrong shape.
#'
#' @param distrib A [multivariate_distrib()] object, whose `distrib_name` is
#'   quoted in the message.
#' @param what The name of the quantity, a single string WITHOUT parentheses:
#'   `"distrib_cdf"`, `"skewness"`. The `()` is appended here, so passing them
#'   gives `distrib_cdf()()`.
#' @param why One sentence saying what is missing, a single string with no
#'   leading capital. It is placed after a colon and carries its own final
#'   period.
#'
#' @return Never returns: it always signals an error, with `call. = FALSE` so
#'   the message is not prefixed by the internal call.
#'
#' @seealso [distrib_cdf.multivariate_distrib()],
#'   [distrib_quantile.multivariate_distrib()],
#'   [skewness.multivariate_distrib()] and
#'   [kurtosis.multivariate_distrib()] for the refusals it produces.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#'
#' # The wording every refusal on this class shares. The name carries no
#' # parentheses; they are appended here.
#' try(distributions7:::mv_refuse(d, "some_quantity",
#'                                "there is no such thing in p dimensions."))
#'
#' # Which is what a caller sees from the generics themselves.
#' try(distrib_cdf(d, rbind(c(0, 0)),
#'                 list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
#'                      sigma_log_L2 = 0, sigma_L2.1 = 0)))
#'
#' @keywords internal
mv_refuse <- function(distrib, what, why) {
  stop(sprintf(
    "%s() is not defined for '%s': %s", what, distrib@distrib_name, why
  ), call. = FALSE)
}


#' @title No Distribution Function in Several Dimensions
#' @name distrib_cdf.multivariate_distrib
#'
#' @description
#' Signals an error. The distribution function of a multivariate law is
#' \eqn{P(Y_1 \le q_1, \ldots, Y_p \le q_p)}, an integral over an orthant.
#' There is no closed form for it in general, and the one-dimensional fallback
#' has no counterpart: [continuous_distrib()]'s default integrates the density
#' over an interval of the line, and an orthant is not one. A numerical orthant
#' probability is a separate piece of work with its own accuracy question, and
#' this package does not attempt it.
#'
#' @details
#' The refusal is on the class, so it covers a family written elsewhere as
#' well as the four that ship. A family that CAN answer registers its own
#' method and the refusal never runs.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param q A numeric matrix of quantiles, one row per point. Not examined: the
#'   error is raised before it is read.
#' @param theta A named list of parameters. Not examined.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Never returns: it always signals an error naming the family.
#'
#' @seealso [distrib_quantile.multivariate_distrib()], refused for a related
#'   reason, [mv_marginal()], which for some families returns a univariate
#'   distribution that does answer, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#'
#' try(distrib_cdf(d, rbind(c(0, 0)), theta))
#'
#' # A gaussian's marginal is a ONE-DIMENSIONAL MvGaussianDistrib, so it
#' # inherits the same refusal.
#' m <- mv_marginal(d, theta, 1)
#' class(m$distrib)[1]
#' try(distrib_cdf(m$distrib, 0, m$theta))
#'
#' # A Dirichlet's marginal is a genuine univariate beta, and answers.
#' b <- mv_marginal(dirichlet_distrib(3),
#'                  list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
#' class(b$distrib)[1]
#' distrib_cdf(b$distrib, 0.4, b$theta)
#'
#' @keywords internal
S7::method(distrib_cdf, multivariate_distrib) <- function(distrib, q, theta, ...) {
  mv_refuse(
    distrib, "distrib_cdf",
    "the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback."
  )
}

#' @title No Quantile Function in Several Dimensions
#' @name distrib_quantile.multivariate_distrib
#'
#' @description
#' Signals an error. A quantile is defined by inverting a distribution function
#' along an ORDERING, and \eqn{\mathbb{R}^p} has none for \eqn{p > 1}. The
#' obstruction is deeper than the one that stops
#' [distrib_cdf.multivariate_distrib()]: there a number exists and is hard to
#' compute, here there is no number to compute. Several quantities go by the
#' name in the literature, and they disagree.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param p A numeric vector of probabilities. Not examined: the error is
#'   raised before it is read.
#' @param theta A named list of parameters. Not examined.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Never returns: it always signals an error naming the family.
#'
#' @seealso [distrib_cdf.multivariate_distrib()] for the related refusal,
#'   [mv_marginal()], which for some families returns a univariate
#'   distribution that does answer, and [distrib_quantile()] for the
#'   generic.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#'
#' try(distrib_quantile(d, 0.5, theta))
#'
#' # A gaussian's marginal is a one-dimensional MvGaussianDistrib and refuses
#' # too; a Dirichlet's is a univariate beta and answers.
#' b <- mv_marginal(dirichlet_distrib(3),
#'                  list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
#' distrib_quantile(b$distrib, c(0.025, 0.5, 0.975), b$theta)
#'
#' @keywords internal
S7::method(distrib_quantile, multivariate_distrib) <- function(distrib, p, theta, ...) {
  mv_refuse(
    distrib, "distrib_quantile",
    "a quantile inverts an ordering of the line, and there is none in several dimensions."
  )
}

#' @title No Numerical Response Gradient in Several Dimensions
#' @name distrib_grad_y.multivariate_distrib
#'
#' @description
#' Signals an error. The univariate numerical fallback differences the
#' log-density along a LINE, which in \eqn{p} dimensions gives a directional
#' derivative: the number it produces is of the wrong SHAPE, not merely
#' inaccurate. The refusal is a design decision. A family that can supply
#' \eqn{\partial\ell/\partial y} registers its own method, and the two
#' elliptical families do.
#'
#' @details
#' [check_distrib()] consults [has_mv_grad_y()] and skips the
#' response-derivative checks where no method is registered, so not
#' registering one costs nothing but the checks it would have earned.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param y An \eqn{n \times p} numeric matrix of observations. Not examined.
#' @param theta A named list of parameters. Not examined.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Never returns: it always signals an error naming the family. A
#'   family that registers a method returns an \eqn{n \times p} matrix.
#'
#' @seealso [distrib_grad_y.MvGaussianDistrib()] and
#'   [distrib_grad_y.MvStudentTDistrib()] for the two families that answer,
#'   [distrib_hess_y.multivariate_distrib()] for the second-order refusal, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' # The Dirichlet registers no response gradient, so the class refuses.
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' try(distrib_grad_y(d, y, theta))
#'
#' # The gaussian registers one and answers with an n by p matrix.
#' g <- mvgaussian_distrib(2)
#' th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'            sigma_L2.1 = 0.5)
#' dim(distrib_grad_y(g, rbind(c(1, -1), c(0, 0)), th))
#'
#' @keywords internal
S7::method(distrib_grad_y, multivariate_distrib) <- function(distrib, y, theta, ...) {
  mv_refuse(
    distrib, "distrib_grad_y",
    "the univariate numerical fallback differences along a line; register a closed form on the family."
  )
}

#' @title No Numerical Response Hessian in Several Dimensions
#' @name distrib_hess_y.multivariate_distrib
#'
#' @description
#' Signals an error, for the reason
#' [distrib_grad_y.multivariate_distrib()] gives: the univariate numerical
#' fallback differences along a line and would produce a scalar where a
#' \eqn{p \times p} matrix is wanted. A family that can supply
#' \eqn{\partial^2\ell/\partial y\,\partial y^\top} registers its own method.
#'
#' @details
#' The two families that answer return DIFFERENT shapes, and a consumer has to
#' allow for both: the gaussian's Hessian does not depend on the observation
#' and comes back as one \eqn{p \times p} matrix, while the Student t's does
#' and comes back as a \eqn{p \times p \times n} array.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param y An \eqn{n \times p} numeric matrix of observations. Not examined.
#' @param theta A named list of parameters. Not examined.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Never returns: it always signals an error naming the family.
#'
#' @seealso [distrib_hess_y.MvGaussianDistrib()] and
#'   [distrib_hess_y.MvStudentTDistrib()] for the two families that answer and
#'   for the shape difference, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' try(distrib_hess_y(d, y, theta))
#'
#' # The two families that answer return different shapes.
#' th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'            sigma_L2.1 = 0.5)
#' yy <- rbind(c(1, -1), c(0, 0))
#' dim(distrib_hess_y(mvgaussian_distrib(2), yy, th))
#' dim(distrib_hess_y(mvstudent_t_distrib(2), yy, c(th, list(nu = 6))))
#'
#' @keywords internal
S7::method(distrib_hess_y, multivariate_distrib) <- function(distrib, y, theta, ...) {
  mv_refuse(
    distrib, "distrib_hess_y",
    "the univariate numerical fallback differences along a line; register a closed form on the family."
  )
}

#' @title No Numerical Mixed Response Derivative in Several Dimensions
#' @name distrib_cross_y.multivariate_distrib
#'
#' @description
#' Signals an error. Unlike the two refusals beside it, the obstruction here is
#' not that a line differences wrongly: it is that no consumer has fixed the
#' shape. The univariate generic returns one number per observation per
#' parameter; the multivariate quantity is one \eqn{n \times p} matrix per
#' parameter, and a family that supplies it registers a method returning that,
#' as the gaussian and the Student t do.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param y An \eqn{n \times p} numeric matrix of observations. Not examined.
#' @param theta A named list of parameters. Not examined.
#' @param scale One of `"parameter"` or `"link"`, handled by the generic before
#'   dispatch. Not examined.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Never returns: it always signals an error naming the family and the
#'   shape a method should produce.
#'
#' @seealso [distrib_cross_y.MvGaussianDistrib()] and
#'   [distrib_cross_y.MvStudentTDistrib()] for the two families that answer,
#'   and [distrib_cross_y()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' try(distrib_cross_y(d, y, theta))
#'
#' # The shape a method returns: one n by p matrix per parameter.
#' g <- mvgaussian_distrib(2)
#' th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'            sigma_L2.1 = 0.5)
#' cy <- distrib_cross_y(g, rbind(c(1, -1), c(0, 0)), th)
#' c(entries = length(cy), rows = nrow(cy[[1]]), cols = ncol(cy[[1]]))
#'
#' @keywords internal
S7::method(distrib_cross_y, multivariate_distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mv_refuse(
    distrib, "distrib_cross_y",
    "register a closed form on the family, returning one n-by-p matrix per parameter as the gaussian and the Student t do."
  )
}


#' @title Expected Information of a Multivariate Distribution, by Sampling
#' @name distrib_expected_hessian.multivariate_distrib
#'
#' @description
#' The fallback for a multivariate family with no closed form: the expectation
#' of the observed Hessian, taken over a sample drawn from the family. All four
#' shipped families register their own closed forms, so this runs only for a
#' family written elsewhere.
#'
#' @details
#' # Why the one-dimensional routes do not survive
#'
#' [continuous_distrib()]'s `"integrate"` splits a one-dimensional integral at
#' quantiles, and there are no quantiles here; `"bartlett"` would need the
#' score's own higher derivatives assembled over a \eqn{p}-dimensional
#' partition sum. What is left is sampling, and both admissible values of
#' `approx` route to it.
#'
#' # The result is a sample, not an integral
#'
#' It carries Monte Carlo error of order `nsim^(-1/2)`, and two calls under
#' different seeds give two answers. Set a seed before calling if the result
#' must be reproducible, and expect a fit that inverts this matrix to move with
#' it. `"bartlett"` averages the outer product of the score and `"mc"` averages
#' the observed Hessian; the two agree in expectation by the second Bartlett
#' identity, and differ by sampling error at any finite `nsim`.
#'
#' @param distrib A [multivariate_distrib()] object with no closed form of its
#'   own.
#' @param y An \eqn{n \times p} numeric matrix of observations; only its row
#'   count is used, the expectation being over the law.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx One of `"bartlett"` (the default, equivalently `"opg"`) or
#'   `"mc"`. `"integrate"` is not available here and falls through to sampling.
#' @param nsim The Monte Carlo sample size. Defaults to `10000`. The error
#'   falls as its square root, so ten times the accuracy costs a hundred times
#'   the draws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed as
#'   [`hess_names(distrib@params)`][hess_names], each vector constant.
#'
#' @seealso [distrib_expected_hessian.MvGaussianDistrib()] and
#'   [distrib_expected_hessian.MvStudentTDistrib()] for the closed forms,
#'   [expected_hessian_exact()], which says which route a family takes, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' # Every shipped family overrides, so reach the fallback directly.
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' base <- S7::method(distrib_expected_hessian, multivariate_distrib)
#'
#' set.seed(1)
#' sampled <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
#'                 approx = "mc", nsim = 4000)
#' exact <- distrib_expected_hessian(d, matrix(0, 3, 2), theta)
#' round(rbind(sampled = vapply(sampled, function(z) z[1], numeric(1)),
#'             closed = vapply(exact, function(z) z[1], numeric(1))), 3)
#'
#' # Two seeds give two answers, this being a sample.
#' set.seed(2)
#' again <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
#'               approx = "mc", nsim = 4000)
#' c(first = sampled$mu1_mu1[1], second = again$mu1_mu1[1])
#'
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

#' @title Prefix a Parametrization's Free Names with the Matrix They Describe
#'
#' @description
#' Returns the matrix parametrization's free names with `"sigma_"` or
#' `"omega_"` in front, so that a printed parameter table says which matrix the
#' coordinate belongs to. A free value's name says how the matrix is BUILT, not
#' which matrix it is, so the same parametrization on the two sides of a model
#' would otherwise give two genuinely different models the same parameter
#' names.
#'
#' @details
#' The prefix is applied by the DISTRIBUTION. The parametrization does not
#' know which side of a model it has been handed to, so it cannot apply one.
#'
#' @param free_names The parametrization's own free names, a character vector.
#' @param inverted Logical of length 1. `TRUE` for a precision, giving
#'   `"omega_"`; `FALSE`, the default, for a covariance or a scale matrix,
#'   giving `"sigma_"`.
#'
#' @return A character vector as long as `free_names`.
#'
#' @seealso [mvgaussian_distrib()], whose two forms this distinguishes, and
#'   [parameters7::log_cholesky()] for a source of free names.
#'
#' @examples
#' distributions7:::mv_prefixed_names(c("log_L1", "log_L2", "L2.1"))
#' distributions7:::mv_prefixed_names(c("log_L1", "L2.1"), inverted = TRUE)
#'
#' # Which is what separates the two parametrizations of one law.
#' mvgaussian_distrib(2)@params
#' mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))@params
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
#' For the multinomial on \eqn{p} coordinates with \eqn{n} trials the support
#' is the weak compositions of \eqn{n} into \eqn{p} parts,
#'
#' \deqn{\mathcal{S} = \Bigl\{y \in \mathbb{N}_0^{p} :
#'   \textstyle\sum_{j=1}^{p} y_j = n\Bigr\},
#'   \qquad \lvert\mathcal{S}\rvert = \binom{n + p - 1}{p - 1},}
#'
#' enumerated by [numericals7::compositions()]. Every expectation is
#' then the finite sum \eqn{\sum_{y \in \mathcal{S}} g(y) f(y; \theta)}.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters. Families whose support does not
#'   depend on them ignore it.
#' @param ... Passed to methods.
#'
#' @return A matrix with one row per support point and one column per
#'   coordinate.
#'
#' @seealso [multinomial_distrib()], [numericals7::compositions()]
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
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @seealso [mv_support()]
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
#' This is consumed by [check_distrib()] and by nothing else. A
#' discrete family does not need it, its normalization being an exact sum over
#' [mv_support()].
#'
#' @param distrib An object inheriting from class
#'   [multivariate_distrib()].
#' @param theta A named list or vector of parameters.
#' @param n The number of draws.
#' @param ... Passed to methods.
#'
#' @return A list with `y`, a matrix of `n` draws, and `logd`,
#'   their log-density under the proposal.
#'
#' @seealso [check_distrib()], [mv_support()]
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
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters.
#' @param n The number of draws.
#' @param ... Unused.
#' @return A list with the draws `y` and their log-density `logd`.
#' @seealso [mv_reference_draw()]
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
#' into the shapes a reader thinks in. `mv_sigma()` returns the matrix the
#' PARAMETRIZATION carries, whichever side the matrix parameter describes: the
#' covariance for a gaussian, and the scale matrix for a Student t, whose
#' covariance is \eqn{\nu\Sigma/(\nu-2)} and does not exist below two degrees
#' of freedom. The moment is [variance()], and keeping the two apart
#' is what lets a heavy-tailed family be described at all.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list or vector of parameters.
#'
#' Both are generics whose base-class method rejects: not every multivariate
#' family has a location, and one that does not should say so rather than
#' hand back its first p parameters under a name that does not fit them.
#'
#' @return A numeric vector of length \eqn{p} for `mv_location()`, and a
#'   \eqn{p \times p} matrix for `mv_sigma()`.
#'
#' @seealso [mvgaussian_distrib()]
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
#' @param distrib A [multivariate_distrib()] object.
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
#' The helper the elliptical families implement [mv_location()] with:
#' the first \eqn{p} entries of the flat parameter vector, labeled by
#' coordinate.
#'
#' @param distrib A [multivariate_distrib()] object.
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
#' For a subset \eqn{A} of the coordinates the marginal density is
#'
#' \deqn{f_A(y_A; \theta)
#'   = \int_{\mathbb{R}^{p - \lvert A \rvert}}
#'     f(y; \theta)\, \mathrm{d}y_{A^{c}},}
#'
#' which for a Gaussian is Gaussian with \eqn{\mu_A} and \eqn{\Sigma_{AA}},
#' and for a Student t is a Student t with \eqn{\mu_A}, \eqn{\Sigma_{AA}}
#' and the **same** \eqn{\nu} (conditioning on the mixing variable of
#' the scale-mixture representation leaves a Gaussian, whose marginal is
#' Gaussian, and the mixture is then taken back). For a Dirichlet the
#' marginal of one coordinate is Beta with the same concentration
#' \eqn{\phi}, and for a multinomial it is binomial.
#'
#' A marginal is not available in general: integrating a density over the
#' coordinates one is not interested in has no closed form for most families.
#' Every family the package ships has one, the elliptical ones by subsetting
#' the mean and the matrix and the two simplex families by the identities
#' above. For a family without one the generic signals an error rather than
#' approximating, since a quadrature over the discarded coordinates would be
#' a different object under the same name.
#'
#' This is what makes a picture of a multivariate distribution possible at all:
#' a panel of a pairs plot shows a marginal, so the plot exists exactly when
#' the marginals do.
#'
#' @param distrib An object inheriting from class
#'   [multivariate_distrib()].
#' @param theta A named list or vector of parameters.
#' @param which An integer vector of coordinates to keep.
#' @param ... Passed to methods.
#'
#' @return A list with `distrib`, the marginal distribution object, and
#'   `theta`, its parameters.
#'
#' @seealso [plot.multivariate_distrib()]
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
#' @param distrib A [multivariate_distrib()] object.
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
#' [mv_marginal()] is not a way round this for an elliptical family,
#' whose marginal is a smaller multivariate distribution and rejects in turn;
#' it is for the Dirichlet and the multinomial, whose marginals are univariate.
#' @param x A [multivariate_distrib()] object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @seealso [mv_marginal()]
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
#' [`skewness()`][skewness.multivariate_distrib].
#' @param x A [multivariate_distrib()] object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return Never returns; raises an error.
#' @seealso [mv_marginal()]
#' @keywords internal
S7::method(kurtosis, multivariate_distrib) <- function(x, theta, ...) {
  mv_refuse(
    x, "kurtosis",
    "a vector response has no single kurtosis: Mardia's, Malkovich-Afifi's and the vector of coordinatewise marginal kurtosises are different quantities. Ask a univariate family instead."
  )
}
