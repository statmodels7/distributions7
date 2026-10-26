#' @include distrib.R generics.R utility_functions.R cross_derivatives.R moments.R
NULL

#' @title S7 Class for Multivariate Distributions
#'
#' @description
#' Defines the subclass of `distrib` for distributions whose observations are
#' vectors.
#' The response is an \eqn{n \times p} matrix, one row per observation, and
#' `n_dim` records \eqn{p}. The class is a sibling of [continuous_distrib()] and
#' [discrete_distrib()] and not a subclass of either, so it inherits none of
#' their one-dimensional defaults. The methods registered on this class signal
#' an error for the distribution function, the quantile function and the
#' response derivatives; the gaussian and the Student t register response
#' derivatives of their own.
#'
#' @details
#' # The parameters are scalars
#'
#' The parameters of a multivariate distribution are scalars, as for every
#' other distribution of the package. A mean vector contributes \eqn{p} of them
#' and a matrix contributes the free values of the \pkg{parameters7}
#' parametrization that carries it, so `theta` is a named list of numbers and
#' [deriv_names()], the Hessian keys, the link scale and [fit_distrib()] apply
#' without a special case.
#'
#' The constraint on a matrix or on a probability vector is carried by its
#' parametrization, so the free values of a matrix and of a simplex take the
#' identity link. A scalar parameter keeps a link of its own: the degrees of
#' freedom of the Student t and the concentration of the Dirichlet take the log
#' link by default.
#'
#' # The one-dimensional defaults
#'
#' The defaults registered on [continuous_distrib()] are one-dimensional: a
#' distribution function by quadrature over an interval, a quantile by root
#' finding and a generator for a scalar density. They have no counterpart in
#' \eqn{p} dimensions, so the class does not inherit them.
#'
#' # Parameters that vary by observation
#'
#' A parametrization describes one matrix for the whole sample, and
#' [mv_flat_theta()] rejects a `theta` component longer than one. A
#' distribution whose parameters depend on covariates is a model, which belongs
#' to the layer above this package.
#'
#' @param n_dim The dimension \eqn{p} of one observation, a single positive
#'   integer. The validator rejects any other value, and also a `dimension`
#'   property other than `"multivariate"`.
#' @inheritParams distrib
#'
#' @return An S7 object of class `multivariate_distrib`, inheriting from
#'   `distrib`. Beyond the parent's `distrib_name`, `dimension`, `bounds`,
#'   `params`, `params_interpretation`, `n_params`, `params_bounds`,
#'   `link_params` and `params_smooth`, it carries `n_dim`.
#'
#' @seealso [mvgaussian1_distrib()] and [mvstudent_t1_distrib()] for the
#'   elliptical families, [dirichlet_distrib()] and [multinomial_distrib()] for
#'   the simplex-valued ones, [n_obs()] for the observation count, and
#'   [mv_summary()] for the quantities that a fit reports.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' S7::S7_inherits(d, multivariate_distrib)
#' c(n_dim = d@n_dim, n_params = d@n_params)
#'
#' # The class is a sibling of continuous_distrib and not a subclass.
#' c(multivariate = S7::S7_inherits(d, multivariate_distrib),
#'   continuous = S7::S7_inherits(d, continuous_distrib))
#'
#' # The free values of the matrix take the identity link; the degrees of
#' # freedom of a Student t take the log link.
#' vapply(d@link_params, function(l) l@link_name, character(1))
#' vapply(mvstudent_t1_distrib(2)@link_params, function(l) l@link_name,
#'        character(1))
#'
#' # The shipped multivariate families all inherit from the class.
#' vapply(list(mvgaussian1_distrib(2), mvstudent_t1_distrib(2),
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


#' @title Number of Observations in a Response
#'
#' @description
#' Returns the number of observations in `y`: its length for a univariate
#' distribution, and its number of rows for a multivariate one. For a
#' multivariate distribution a vector without a `dim` attribute counts as one
#' observation, or none when it is empty, whatever its length. The functions of
#' the package that count observations call this one, because for a matrix
#' response `length(y)` counts the \eqn{np} entries, and a recycling check built
#' on it would require parameters of length \eqn{np}.
#'
#' @param distrib An object inheriting from `distrib`. Only its class is used,
#'   to select the reading of `y` that applies.
#' @param y The response: a numeric vector for a univariate distribution, or an
#'   \eqn{n \times p} numeric matrix for a multivariate one.
#'
#' @return A single integer.
#'
#' @seealso [multivariate_distrib()] and [as_mv_matrix()], which puts a
#'   response in the shape that this function counts.
#'
#' @examples
#' # A univariate response is counted by length.
#' n_obs(gaussian1_distrib(), c(1.2, -0.4, 0.8))
#'
#' # A multivariate response is counted by rows, where length() gives the
#' # number of entries.
#' y <- matrix(0, 5, 2)
#' c(n_obs = n_obs(mvgaussian1_distrib(2), y), length = length(y))
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
#' Puts `y` in the \eqn{n \times p} form that every multivariate method
#' expects. A matrix of the right width is returned unchanged, and a numeric
#' vector of length \eqn{p} is read as a single observation and returned as a
#' one-row matrix, which is the reading of a density at one point. A vector of
#' any other length is an error naming its length and the dimension, and a
#' response that is not numeric is an error as well.
#'
#' @param distrib A [multivariate_distrib()] object, from which `n_dim` is
#'   read.
#' @param y A numeric matrix with `distrib@n_dim` columns, or a numeric vector
#'   of length `distrib@n_dim`. An empty vector gives a matrix with no rows. A
#'   vector of any other length, a matrix of the wrong width and a non-numeric
#'   `y` are errors.
#'
#' @return A numeric matrix with `distrib@n_dim` columns.
#'
#' @seealso [n_obs()] for the row count and [multivariate_distrib()] for the
#'   response convention.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
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
  if (!is.numeric(y)) stop("'y' must be numeric.", call. = FALSE)
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
  y
}


#' @title Require Scalar Parameters, and Flatten Them
#'
#' @description
#' Rejects a `theta` whose components are not single numbers, and returns the
#' rest as one named numeric vector in `distrib@params` order. The multivariate
#' families take one value per parameter for the whole sample: a
#' parametrization describes one matrix, and a parameter that varies by
#' observation belongs to a model, which is the layer above this package.
#'
#' @details
#' A univariate family reads a `theta` component of length \eqn{n} as one value
#' per observation. Without this check the same component given to a
#' multivariate family would be recycled against the \eqn{p} columns without an
#' error. The message names the parameter.
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
#' d <- mvgaussian1_distrib(2)
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


#' @title Signal the Error for a Quantity Without a Multivariate Counterpart
#'
#' @description
#' Signals the error of a multivariate distribution for a quantity that has no
#' multivariate counterpart. The message has one wording and consists of the
#' name of the quantity, the family, and a sentence that states what is
#' missing. Every method of the class that signals such an error calls this
#' function, so the messages share one form and a method does not return a
#' number of the wrong shape.
#'
#' @param distrib A [multivariate_distrib()] object, whose `distrib_name` is
#'   quoted in the message.
#' @param what The name of the quantity, a single string without parentheses,
#'   such as `"distrib_cdf"` or `"skewness"`. The `()` is appended here.
#' @param why One sentence that states what is missing, a single string without
#'   a leading capital. It is placed after a colon and carries its own final
#'   period.
#'
#' @return Does not return: it signals an error, with `call. = FALSE` so that
#'   the message is not prefixed by the internal call.
#'
#' @seealso [distrib_cdf.multivariate_distrib()],
#'   [distrib_quantile.multivariate_distrib()],
#'   [skewness.multivariate_distrib()] and
#'   [kurtosis.multivariate_distrib()] for errors that it produces.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#'
#' # This call shows the wording of these errors. The name carries no
#' # parentheses; they are appended here.
#' try(distributions7:::mv_refuse(d, "some_quantity",
#'                                "there is no such thing in p dimensions."))
#'
#' # A generic signals the error with the same wording.
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


#' @title Distribution Function of a Multivariate Distribution
#' @name distrib_cdf.multivariate_distrib
#'
#' @description
#' Signals an error. The distribution function of a multivariate law is
#' \eqn{P(Y_1 \le q_1, \ldots, Y_p \le q_p)}, an integral over an orthant. It
#' has no closed form in general, and the default of [continuous_distrib()],
#' which integrates the density over an interval of the line, has no
#' counterpart in several dimensions. A numerical orthant probability is a
#' separate computation with an accuracy question of its own, and the package
#' does not provide it.
#'
#' @details
#' The method is registered on the class, so it applies to a family written
#' elsewhere as well as to the shipped ones. A family that can compute the
#' distribution function registers its own method.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param q A numeric matrix of quantiles, one row per point. It is not read,
#'   because the error is signalled first.
#' @param theta A named list of parameters. It is not read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [distrib_quantile.multivariate_distrib()] for the quantile function,
#'   [mv_marginal()], which for the Dirichlet and the multinomial returns a
#'   univariate distribution with a distribution function, and [distrib_cdf()]
#'   for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#'
#' try(distrib_cdf(d, rbind(c(0, 0)), theta))
#'
#' # The marginal of a gaussian is a one-dimensional MvGaussian1Distrib, which
#' # signals the same error.
#' m <- mv_marginal(d, theta, 1)
#' class(m$distrib)[1]
#' try(distrib_cdf(m$distrib, 0, m$theta))
#'
#' # The marginal of a Dirichlet is a univariate beta, which has a distribution
#' # function.
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

#' @title Quantile Function of a Multivariate Distribution
#' @name distrib_quantile.multivariate_distrib
#'
#' @description
#' Signals an error. A quantile inverts a distribution function along an
#' ordering of the line, and \eqn{\mathbb{R}^p} has no ordering compatible with
#' its vector structure for \eqn{p > 1}, so the multivariate quantile is not one
#' object: several quantities go by that name in the literature, and they
#' differ.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param p A numeric vector of probabilities. It is not read, because the
#'   error is signalled first.
#' @param theta A named list of parameters. It is not read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [distrib_cdf.multivariate_distrib()] for the distribution function,
#'   [mv_marginal()], which for the Dirichlet and the multinomial returns a
#'   univariate distribution with a quantile function, and [distrib_quantile()]
#'   for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#'
#' try(distrib_quantile(d, 0.5, theta))
#'
#' # The marginal of a gaussian signals the same error; the marginal of a
#' # Dirichlet is a univariate beta and has a quantile function.
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

#' @title Response Gradient of a Multivariate Family Without a Method
#' @name distrib_grad_y.multivariate_distrib
#'
#' @description
#' Signals an error. The univariate numerical fallback differences the
#' log-density along a line, which in \eqn{p} dimensions gives a directional
#' derivative, a number of the wrong shape. A family that can supply
#' \eqn{\partial\ell/\partial y} registers its own method, as the gaussian and
#' the Student t do.
#'
#' @details
#' [check_distrib()] tests for a method of the family's own with
#' [has_mv_grad_y()] and leaves out the response-derivative check when there is
#' none.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param y An \eqn{n \times p} numeric matrix of observations. It is not read.
#' @param theta A named list of parameters. It is not read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family. A family
#'   that registers a method returns an \eqn{n \times p} matrix.
#'
#' @seealso [distrib_grad_y.MvGaussianDistrib()] and
#'   [distrib_grad_y.MvStudentTDistrib()] for the two families with a method,
#'   [distrib_hess_y.multivariate_distrib()] for the second order, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' # The Dirichlet registers no response gradient, so the method of the class
#' # signals an error.
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' try(distrib_grad_y(d, y, theta))
#'
#' # The gaussian registers one, which returns an n by p matrix.
#' g <- mvgaussian1_distrib(2)
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

#' @title Response Hessian of a Multivariate Family Without a Method
#' @name distrib_hess_y.multivariate_distrib
#'
#' @description
#' Signals an error, for the same reason as
#' [distrib_grad_y.multivariate_distrib()]: the univariate numerical fallback
#' differences along a line and would produce a scalar where a
#' \eqn{p \times p} matrix is required. A family that can
#' supply \eqn{\partial^2\ell/\partial y\,\partial y^\top} registers its own
#' method.
#'
#' @details
#' The two families with a method return different shapes: the gaussian's
#' Hessian does not depend on the observation and is returned as one
#' \eqn{p \times p} matrix, while the Student t's does and is returned as a
#' \eqn{p \times p \times n} array.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param y An \eqn{n \times p} numeric matrix of observations. It is not read.
#' @param theta A named list of parameters. It is not read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [distrib_hess_y.MvGaussianDistrib()] and
#'   [distrib_hess_y.MvStudentTDistrib()] for the two families with a method,
#'   and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' try(distrib_hess_y(d, y, theta))
#'
#' # The two families with a method return different shapes.
#' th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'            sigma_L2.1 = 0.5)
#' yy <- rbind(c(1, -1), c(0, 0))
#' dim(distrib_hess_y(mvgaussian1_distrib(2), yy, th))
#' dim(distrib_hess_y(mvstudent_t1_distrib(2), yy, c(th, list(nu = 6))))
#'
#' @keywords internal
S7::method(distrib_hess_y, multivariate_distrib) <- function(distrib, y, theta, ...) {
  mv_refuse(
    distrib, "distrib_hess_y",
    "the univariate numerical fallback differences along a line; register a closed form on the family."
  )
}

#' @title Mixed Response Derivative of a Multivariate Family Without a Method
#' @name distrib_cross_y.multivariate_distrib
#'
#' @description
#' Signals an error. The univariate generic returns one number per observation
#' and per parameter, while the multivariate quantity is one
#' \eqn{n \times p} matrix per parameter, and the class has no fallback of that
#' shape. A family that supplies it registers a method returning one matrix per
#' parameter, as the gaussian and the Student t do.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param y An \eqn{n \times p} numeric matrix of observations. It is not read.
#' @param theta A named list of parameters. It is not read.
#' @param scale One of `"parameter"` or `"link"`. The generic applies the link
#'   scale to the result of a method. It is not read here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family and the
#'   shape that a method produces.
#'
#' @seealso [distrib_cross_y.MvGaussianDistrib()] and
#'   [distrib_cross_y.MvStudentTDistrib()] for the two families with a method,
#'   and [distrib_cross_y()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' try(distrib_cross_y(d, y, theta))
#'
#' # A method returns one n by p matrix per parameter.
#' g <- mvgaussian1_distrib(2)
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
#' Computes the expected information of a multivariate family without a closed
#' form as an average over a sample drawn from the family. Every shipped
#' multivariate class
#' registers its own closed form, so this method runs only for a family written
#' elsewhere.
#'
#' @details
#' # The routes
#'
#' The `"integrate"` route of [continuous_distrib()] splits a one-dimensional
#' integral at quantiles, and a multivariate family has no quantiles, so
#' `method = expect_integrate()` signals an error. The other methods sample:
#' [expect_opg()] and [expect_mc()] with `identity = "score"` average the
#' negative outer product of the score over simulated draws, and [expect_mc()]
#' with `identity = "hessian"` averages the observed Hessian. The two averages
#' agree in expectation by the second Bartlett identity and differ by sampling
#' error at any finite `nsim`. The sample size is that of [expect_mc()] and
#' 10000 for [expect_opg()]. The routine [mv_expected_sampled()] does the
#' sampling and also serves an explicit `method` for a multivariate family
#' that has a closed form.
#'
#' # Monte Carlo error
#'
#' The result carries Monte Carlo error of order `nsim^(-1/2)`, and two calls
#' under different seeds give different values. A seed set before the call
#' makes the result reproducible, and a fit that inverts this matrix moves with
#' it.
#'
#' @param distrib A [multivariate_distrib()] object without a closed form of
#'   its own.
#' @param y An \eqn{n \times p} numeric matrix of observations; only its row
#'   count is used, because the expectation is over the law.
#' @param theta A named list of parameters, in which each component is a
#'   single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The generic
#'   applies the link scale to the result of the method.
#' @param method `NULL` (the default), which selects [expect_opg()], or an
#'   object inheriting from [ExpectMethod]. [expect_opg()] and [expect_mc()]
#'   sample, and [expect_integrate()] signals an error. See
#'   [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed as
#'   [`hess_names(distrib@params)`][hess_names], in which each vector is
#'   constant.
#'
#' @seealso [distrib_expected_hessian.MvGaussianDistrib()] and
#'   [distrib_expected_hessian.MvStudentTDistrib()] for two closed forms,
#'   [expected_hessian_exact()], which reports whether a family has a closed
#'   form, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' # Every shipped family has its own method, so the fallback is called
#' # directly.
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' base <- S7::method(distrib_expected_hessian, multivariate_distrib)
#'
#' set.seed(1)
#' sampled <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
#'                 method = expect_mc(nsim = 4000))
#' exact <- distrib_expected_hessian(d, matrix(0, 3, 2), theta)
#' round(rbind(sampled = vapply(sampled, function(z) z[1], numeric(1)),
#'             closed = vapply(exact, function(z) z[1], numeric(1))), 3)
#'
#' # Two seeds give two values, because the result is a sample.
#' set.seed(2)
#' again <- base(d, matrix(0, 3, 2), distributions7:::align_theta(d, theta),
#'               method = expect_mc(nsim = 4000))
#' c(first = sampled$mu1_mu1[1], second = again$mu1_mu1[1])
#'
#' @keywords internal
S7::method(distrib_expected_hessian, multivariate_distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    method = NULL, ...) {
  check_method(method)
  mv_expected_sampled(distrib, y, theta, 2L, resolve_method(method, 2L))
}

#' @title Expected Derivatives of a Multivariate Distribution, by Sampling
#' @name mv_expected_sampled
#'
#' @description
#' Computes the expected derivative of an order for a multivariate family as an
#' average over a sample drawn from the family. It serves the fallback of
#' [distrib_expected_hessian.multivariate_distrib()] and an explicit
#' `method` given for a multivariate family that has a closed form.
#'
#' @details
#' A multivariate family has no quantiles, so [expect_integrate()] signals an
#' error. [expect_opg()] and [expect_mc()] with `identity = "score"` average
#' minus the products of the score identity (at order 2, the negative outer
#' product of the score) over `method_nsim(method)` draws, and [expect_mc()]
#' with `identity = "hessian"` averages the observed derivative of the order.
#' The averages agree in expectation and differ by sampling error at any finite
#' `nsim`; the sample size is 10000 for [expect_opg()].
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param y The observations; only the number of rows is used, because the
#'   expectation is over the law.
#' @param theta A named list of parameters, in which each component is a single
#'   number.
#' @param order The derivative order, 2 to 4.
#' @param method An object inheriting from [ExpectMethod].
#'
#' @return A named list of numeric vectors of length \eqn{n}, in which each
#'   vector is constant.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
mv_expected_sampled <- function(distrib, y, theta, order, method) {
  if (S7::S7_inherits(method, ExpectIntegrate)) {
    mv_refuse(
      distrib, "distrib_expected_hessian",
      "expect_integrate() builds a quadrature over an interval, which has no counterpart in several dimensions. Use expect_mc() or expect_opg()."
    )
  }
  n <- n_obs(distrib, y)
  params <- distrib@params
  nms <- if (order == 2L) hess_names(params) else deriv_names(params, order)
  big <- distrib_rng(distrib, method_nsim(method), theta)
  score <- S7::S7_inherits(method, ExpectOPG) ||
    identical(method@identity, "score")
  vals <- if (score && order == 2L) {
    s <- do.call(cbind, distrib_gradient(distrib, big, theta))
    pr <- hess_pairs(params)
    vapply(nms, function(k) {
      ij <- match(pr[[k]], params)
      -mean(s[, ij[1L]] * s[, ij[2L]])
    }, numeric(1))
  } else {
    d <- if (score) score_identity_terms(distrib, big, theta, order)
         else observed_deriv(distrib, big, theta, order)
    vapply(nms, function(k) mean(d[[k]]), numeric(1))
  }
  stats::setNames(lapply(vals, function(v) rep(v, n)), nms)
}

#' @title A Factor of the Covariance for Sampling
#'
#' @description
#' Returns a matrix \eqn{F} with \eqn{F F^\top = \Sigma}, read from the matrix
#' parametrization's own factor instead of a Cholesky decomposition of the
#' assembled matrix. Where the parametrization carries \eqn{\Sigma}, \eqn{F} is
#' its factor \eqn{L}; where it carries the precision \eqn{\Omega = L L^\top},
#' \eqn{F = L^{-\top}}.
#'
#' @details
#' A covariance that is nearly singular can round to a matrix that `chol()`
#' rejects although the parametrization describes a valid one: at a
#' log-Cholesky free value \eqn{\log L_{22} = -20} the entry \eqn{\Sigma_{22} =
#' 1 + e^{-40}} rounds to 1. The factor of the parametrization is exact there.
#'
#' @param s A \pkg{parameters7} matrix parametrization of full rank.
#' @param eta Its free vector.
#' @param inverted Logical of length 1, `TRUE` when `s` carries the precision.
#'
#' @return A \eqn{p \times p} numeric matrix without dimnames.
#'
#' @seealso [distrib_rng.MvGaussianDistrib()] and
#'   [distrib_rng.MvStudentTDistrib()], the two callers, and
#'   [parameters7::param_factor()] for the factor.
#'
#' @examples
#' s <- parameters7::log_cholesky(2)
#' eta <- c(0, -20, 1)
#' f <- distributions7:::mv_sampling_factor(s, eta, FALSE)
#' f
#' # chol() of the assembled matrix rejects it.
#' try(chol(parameters7::param_value(s, eta)))
#'
#' @keywords internal
mv_sampling_factor <- function(s, eta, inverted) {
  f <- unname(unclass(parameters7::param_factor(s, eta)))
  # Omega = L L' gives Sigma = L^-T L^-1, whose factor is L^-T
  if (isTRUE(inverted)) f <- backsolve(t(f), diag(nrow(f)))
  f
}


#' @title The Quadratic Form of a Residual, Through the Chart's Factor
#'
#' @description
#' Computes \eqn{q_i = r_i^\top \Sigma^{-1} r_i} for each row of `r` as the
#' squared norm of a triangular solve or product with the matrix
#' parametrization's own factor: \eqn{\lVert L^{-1} r_i\rVert^2} where the
#' parametrization carries \eqn{\Sigma = L L^\top}, and
#' \eqn{\lVert L^\top r_i\rVert^2} where it carries the precision
#' \eqn{\Omega = L L^\top}.
#'
#' @details
#' Near a singular \eqn{\Sigma} the inverse has large entries of nearly rank
#' one, and \eqn{r^\top \Sigma^{-1} r} formed with the explicit inverse loses
#' the digits that an observation close to the degenerate direction has: at a
#' log-Cholesky free value \eqn{\log L_{22} = -11.5} the relative error was of
#' order \eqn{10^{-7}}. The squared norm adds positive terms only.
#'
#' @param r An \eqn{n \times p} numeric matrix of residuals.
#' @param s A \pkg{parameters7} matrix parametrization of full rank.
#' @param eta Its free vector.
#' @param inverted Logical of length 1, `TRUE` when `s` carries the precision.
#'
#' @return A numeric vector of length \eqn{n}.
#'
#' @seealso [mv_sampling_factor()] for the same factor used to draw, and
#'   [mvt_weights()] and [distrib_pdf.MvGaussianDistrib()] for consumers.
#'
#' @examples
#' s <- parameters7::log_cholesky(2)
#' eta <- c(0.1, -0.2, 0.4)
#' r <- rbind(c(1, -1), c(0.5, 2))
#' distributions7:::mv_quad(r, s, eta, FALSE)
#' rowSums((r %*% solve(parameters7::param_value(s, eta))) * r)
#'
#' @keywords internal
mv_quad <- function(r, s, eta, inverted) {
  if (!nrow(r)) return(numeric(0))
  l <- unname(unclass(parameters7::param_factor(s, eta)))
  u <- if (isTRUE(inverted)) crossprod(l, t(r)) else forwardsolve(l, t(r))
  colSums(u^2)
}


#' @title Prefix Free Names with the Matrix That a Parametrization Describes
#'
#' @description
#' Returns the free names of a matrix parametrization with `"sigma_"` or
#' `"omega_"` in front, so that a printed parameter table names the matrix to
#' which each coordinate belongs. The free name of a value records how the
#' matrix is built and not which matrix it is, so without the prefix the same
#' parametrization of a covariance and of a precision would give two different
#' models the same parameter names.
#'
#' @details
#' The distribution applies the prefix, because the parametrization does not
#' record which side of a model it describes.
#'
#' @param free_names The parametrization's own free names, a character vector.
#' @param inverted Logical of length 1. It is `TRUE` for a precision, which
#'   gives `"omega_"`, and `FALSE` (the default) for a covariance or a scale
#'   matrix, which gives `"sigma_"`.
#'
#' @return A character vector as long as `free_names`.
#'
#' @seealso [mvgaussian1_distrib()], whose two forms this distinguishes, and
#'   [parameters7::log_cholesky()] for a source of free names.
#'
#' @examples
#' distributions7:::mv_prefixed_names(c("log_L1", "log_L2", "L2.1"))
#' distributions7:::mv_prefixed_names(c("log_L1", "L2.1"), inverted = TRUE)
#'
#' # The prefixes separate the two parametrizations of one law.
#' mvgaussian1_distrib(2)@params
#' mvgaussian2_distrib(2, parameters7::log_cholesky(2))@params
#'
#' @keywords internal
mv_prefixed_names <- function(free_names, inverted = FALSE) {
  paste0(if (isTRUE(inverted)) "omega_" else "sigma_", free_names)
}


#' @title The Matrix, Its Inverse and Their Derivatives
#'
#' @description
#' Computes the matrix quantities from which the elliptical families evaluate,
#' for a matrix parametrization that carries either \eqn{\Sigma} or
#' \eqn{\Sigma^{-1}}: the matrix, its inverse, the log-determinant of
#' \eqn{\Sigma}, and the first and second derivative arrays of \eqn{\Sigma} and
#' of \eqn{\Sigma^{-1}}.
#'
#' @details
#' Where the parametrization carries \eqn{\Sigma}, the derivatives of
#' \eqn{\Sigma} are its own arrays and those of \eqn{\Sigma^{-1}} come from
#' [parameters7::param_inv_d1()] and [parameters7::param_inv_d2()], which are
#' exact for [parameters7::log_cholesky()] where the matrix is nearly singular.
#' Where it carries \eqn{\Sigma^{-1}}, the derivatives of the inverse are the
#' parametrization's own arrays, and those of \eqn{\Sigma} follow from the
#' chain rule for an inverse,
#' \deqn{\partial_k\Sigma = -\Sigma A_k \Sigma, \qquad
#'   \partial_{kl}\Sigma = \Sigma\left(A_l\Sigma A_k + A_k\Sigma A_l
#'     - A_{kl}\right)\Sigma,}
#' with \eqn{A_k} and \eqn{A_{kl}} the parametrization's arrays.
#'
#' The derivatives of the log-determinant are returned for \eqn{\Sigma} in
#' either case, as `dlogdet` and `d2logdet`. Where the parametrization carries
#' \eqn{M = \Sigma^{-1}}, \eqn{\log|\Sigma| = -\log|M|}, so the
#' parametrization's own derivatives are negated here, once for both families.
#'
#' @param s A \pkg{parameters7} matrix parametrization.
#' @param eta Its free vector.
#' @param inverted `TRUE` when `s` carries \eqn{\Sigma^{-1}}.
#' @param derivs,derivs2 Logical flags that request the first and the second
#'   derivative arrays. `derivs2 = TRUE` computes the first derivatives as
#'   well.
#'
#' @return A list with `sigma`, `sigma_inv`, `logdet` (of \eqn{\Sigma}) and
#'   `inverted`; with `derivs` or `derivs2` also `a` and `ai`, the first
#'   derivative arrays of \eqn{\Sigma} and of \eqn{\Sigma^{-1}}, and `dlogdet`;
#'   with `derivs2` also `a2`, `ai2` and `d2logdet`.
#'
#' @seealso [mvg_pieces()] and [mvt_pieces()], the two callers.
#'
#' @keywords internal
mv_matrix_pieces <- function(s, eta, inverted, derivs = FALSE,
                             derivs2 = FALSE) {
  m <- parameters7::param_value(s, eta)
  ld <- parameters7::param_logdet(s, eta)

  # the inverse from param_solve() is symmetric only to rounding (1e-16 off
  # the diagonal), and every array built on it inherits the asymmetry
  inv <- parameters7::param_solve(s, eta)
  inv <- (inv + t(inv)) / 2
  if (isTRUE(inverted)) {
    sigma <- inv
    sigma_inv <- m
    logdet <- -ld
  } else {
    sigma <- m
    sigma_inv <- inv
    logdet <- ld
  }
  out <- list(sigma = unname(sigma), sigma_inv = unname(sigma_inv),
              logdet = logdet, inverted = isTRUE(inverted))

  if (derivs || derivs2) {
    d <- parameters7::param_d1(s, eta)
    # THE DERIVATIVES OF THE INVERSE, read from the chart and never formed as
    # -Sigma^-1 A Sigma^-1. Near a singular covariance Sigma^-1 is large and
    # nearly of rank one, and that sandwich cancels: at a log-Cholesky free
    # value log L22 = -11.5 it carried a relative error of 1e-7, which a trace
    # of products of order 1e10 turned into an error of order one in the
    # exact outer gradient. On the precision side the chart's own derivatives
    # are those of the inverse and are exact.
    out$ai <- lapply(if (isTRUE(inverted)) d else
      parameters7::param_inv_d1(s, eta), unname)
    if (isTRUE(inverted)) {
      d <- lapply(d, function(ak) -(sigma %*% ak %*% sigma))
    }
    out$a <- lapply(d, unname)
    dld <- parameters7::param_dlogdet(s, eta)
    out$dlogdet <- if (isTRUE(inverted)) -dld else dld
  }
  if (derivs2) {
    d2ld <- parameters7::param_d2logdet(s, eta)
    out$d2logdet <- if (isTRUE(inverted)) -d2ld else d2ld
    d2 <- parameters7::param_d2(s, eta)
    out$ai2 <- if (isTRUE(inverted)) d2 else parameters7::param_inv_d2(s, eta)
    if (isTRUE(inverted)) {
      idx <- parameters7::param_tuple_indices(s)
      a1 <- parameters7::param_d1(s, eta)
      d2 <- lapply(seq_along(idx), function(i) {
        k <- idx[[i]][1L]
        l <- idx[[i]][2L]
        sigma %*% (a1[[l]] %*% sigma %*% a1[[k]] +
          a1[[k]] %*% sigma %*% a1[[l]] - d2[[i]]) %*% sigma
      })
      names(d2) <- parameters7::param_tuple_names(s)
    }
    out$a2 <- lapply(d2, unname)
  }
  out
}


#' @title The Support Points of a Discrete Multivariate Distribution
#'
#' @description
#' Returns the points on which a discrete multivariate distribution places
#' mass, as a matrix with one row per point. With the support enumerated, every
#' expectation is an exact sum, and [check_distrib()] verifies the total mass
#' by addition instead of by sampling.
#'
#' @details
#' # Why the generic exists
#'
#' A univariate discrete distribution does not need it, since its support is a
#' range of integers, which the package traverses. The support of a vector is a
#' set whose shape depends on the family, so each discrete multivariate family
#' enumerates its own.
#'
#' # The multinomial's support
#'
#' For \eqn{p} coordinates and \eqn{n} trials it is the set of weak
#' compositions of \eqn{n} into \eqn{p} parts,
#' \deqn{\mathcal{S} = \Bigl\{y \in \mathbb{N}_0^{p} :
#'   \textstyle\sum_{j=1}^{p} y_j = n\Bigr\},
#'   \qquad \lvert\mathcal{S}\rvert = \binom{n + p - 1}{p - 1},}
#' enumerated by [numericals7::compositions()]. Every expectation is then the
#' finite sum \eqn{\sum_{y \in \mathcal{S}} g(y) f(y; \theta)}. The count grows
#' quickly: 15 points at \eqn{n = 4, p = 3}, and 5151 at \eqn{n = 100, p = 3}.
#'
#' # The method of the class
#'
#' It signals an error. A continuous family has no such set, and a discrete
#' family whose support is infinite has no finite matrix to return.
#' [has_mv_support()] tests whether a family has a method of its own.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list or vector of parameters. The shipped method does
#'   not read it, because the support of the multinomial does not depend on
#'   them.
#' @param ... Passed to methods. The shipped method does not read it.
#'
#' @return A numeric matrix with one row per support point and one column per
#'   coordinate.
#'
#' @section Notation:
#' \eqn{\mathcal{S}} is the support, \eqn{p} the dimension, \eqn{n} the number
#' of trials of a multinomial and \eqn{f} the mass function.
#'
#' @seealso [multinomial_distrib()], the one family with a method,
#'   [has_mv_support()] for the test, [numericals7::compositions()] for the
#'   enumeration, and [check_distrib()] for the consumer.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 4)
#' theta <- list(probs_alr1 = 0.4, probs_alr2 = -0.3)
#' sup <- mv_support(d, theta)
#'
#' # Each row is one point, and each row sums to the number of trials.
#' dim(sup)
#' unique(rowSums(sup))
#' head(sup, 4)
#'
#' # The count is the number of weak compositions.
#' c(got = nrow(sup), expected = choose(4 + 3 - 1, 3 - 1))
#'
#' # The mass over the support sums to one exactly.
#' sum(distrib_pdf(d, sup, theta))
#'
#' # A continuous family has no such set, and the call signals an error.
#' try(mv_support(mvgaussian1_distrib(2),
#'                list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
#'                     sigma_log_L2 = 0, sigma_L2.1 = 0)))
#'
#' @export
mv_support <- S7::new_generic("mv_support", "distrib",
  function(distrib, theta, ...) S7::S7_dispatch()
)

#' @title Support of a Family Without an Enumerable One
#' @name mv_support.multivariate_distrib
#'
#' @description
#' Signals an error. A continuous family places no mass on a finite set of
#' points, and a discrete family whose support is infinite has no finite matrix
#' to return. A family with a finite support registers its own method; among
#' the shipped families the multinomial does.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param theta A named list of parameters. It is not read, because the error
#'   is signalled first.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [mv_support()] for the generic and the multinomial's method,
#'   [has_mv_support()] for the test that avoids this error, and
#'   [mv_reference_draw()], the route that [check_distrib()] takes for a
#'   continuous family.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0)
#' try(mv_support(d, theta))
#'
#' # check_distrib() applies this test first.
#' c(gaussian = distributions7:::has_mv_support(d),
#'   multinomial = distributions7:::has_mv_support(
#'     multinomial_distrib(3, size = 4)))
#'
#' @keywords internal
S7::method(mv_support, multivariate_distrib) <- function(distrib, theta, ...) {
  stop(sprintf(paste0(
    "'%s' does not enumerate a support. A continuous family has no such set,\n",
    "  and a discrete one whose support is infinite has no finite matrix to\n",
    "  return; a family that does have one registers this generic."
  ), distrib@distrib_name), call. = FALSE)
}


#' @title A Proposal for Integrating a Multivariate Density
#'
#' @description
#' Draws from a distribution that dominates the family, and returns the
#' log-density of the draws under that distribution, so that the
#' importance-sampling estimate of \eqn{\int f = 1} can be formed.
#' [check_distrib()] is its only consumer; a discrete family takes the exact
#' sum over [mv_support()] instead.
#'
#' @details
#' # The default proposal
#'
#' The method of the class draws from a gaussian with the family's mean and
#' twice its covariance, which requires a finite, non-singular covariance. The
#' inflation keeps the estimate informative: for a gaussian family, a proposal
#' equal to the family makes every importance ratio one.
#'
#' # Families with a proposal of their own
#'
#' A Dirichlet lives on the simplex, a set of measure zero in
#' \eqn{\mathbb{R}^p}, where a gaussian spread over the ambient space places no
#' mass; its covariance is singular, and the estimate from the default proposal
#' is far from one. The Dirichlet registers the uniform distribution on the
#' simplex, whose density is the constant \eqn{\Gamma(p)}. A Student t has no
#' finite covariance at \eqn{\nu \le 2}, and its tails are heavier than a
#' gaussian's at every \eqn{\nu}; it registers a Student t with half its degrees
#' of freedom and twice its scale matrix.
#'
#' The draws and their log-density are taken against the dominating measure
#' with respect to which the family's density is written.
#'
#' @param distrib An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, in which each component
#'   is a single number.
#' @param n The number of draws, a single positive whole number.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A named list with `y`, an \eqn{n \times p} numeric matrix of draws,
#'   and `logd`, a numeric vector of length \eqn{n} holding their log-density
#'   under the proposal.
#'
#' @section Notation:
#' \eqn{f} is the family's density, \eqn{p} the dimension and \eqn{\Gamma} the
#' gamma function.
#'
#' @seealso [check_distrib()], the only consumer,
#'   [mv_reference_draw.multivariate_distrib()] for the default,
#'   [mv_reference_draw.MvStudentTDistrib()] for the Student t's, and
#'   [mv_support()], the exact route of a discrete family.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0)
#' set.seed(1)
#' str(mv_reference_draw(d, theta, 5))
#'
#' # This is the estimate for which the proposal is built: the density over
#' # the proposal, averaged.
#' set.seed(2)
#' r <- mv_reference_draw(d, theta, 20000)
#' mean(exp(distrib_pdf(d, r$y, theta, log = TRUE) - r$logd))
#'
#' # A Dirichlet registers its own proposal, the uniform distribution on the
#' # simplex.
#' dd <- dirichlet_distrib(3)
#' thd <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(3)
#' rd <- mv_reference_draw(dd, thd, 5)
#' round(rd$y, 4)
#' unique(round(rd$logd, 10))
#' log(gamma(3))
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
#'
#' @description
#' Draws from the default proposal, a gaussian with the family's mean and twice
#' its covariance, as \eqn{\mu + Lz} with \eqn{LL^\top = 2\Sigma}. For a
#' gaussian family a proposal with the family's own covariance would make every
#' importance ratio one, and the inflation keeps the ratio informative.
#'
#' @details
#' The covariance must be finite and non-singular, which a family supported on
#' all of \eqn{\mathbb{R}^p} with finite second moments has. A family whose
#' support is a lower-dimensional set, or whose covariance is infinite,
#' registers its own method (see [mv_reference_draw()]).
#'
#' The log-density is computed from the Cholesky factor \eqn{L}: the whitened
#' residual \eqn{L^{-1}r} solves a lower triangular system, which
#' [base::forwardsolve()] evaluates.
#'
#' @param distrib A [multivariate_distrib()] object with a mean and a finite
#'   covariance.
#' @param theta A named list of parameters, already aligned by the generic.
#' @param n The number of draws, a single positive whole number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `y`, an \eqn{n \times p} numeric matrix, and
#'   `logd`, a numeric vector of length \eqn{n}.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance, \eqn{L} a lower Cholesky
#' factor of \eqn{2\Sigma} and \eqn{p} the dimension.
#'
#' @seealso [mv_reference_draw()] for the generic and the proposals of the
#'   families that register their own, and [check_distrib()] for the consumer.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#'
#' set.seed(1)
#' r <- mv_reference_draw(d, theta, 4)
#' r$y
#'
#' # logd is the log-density of the inflated gaussian at those draws.
#' S2 <- 2 * variance(d, theta)
#' mu <- c(1, -1)
#' ref <- -0.5 * (2 * log(2 * pi) + log(det(S2)) + mahalanobis(r$y, mu, S2))
#' all.equal(r$logd, as.numeric(ref))
#'
#' # The sample is twice as spread as the family.
#' set.seed(2)
#' round(rbind(proposal = diag(var(mv_reference_draw(d, theta, 20000)$y)),
#'             family = diag(variance(d, theta))), 3)
#'
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




#' @title The Location Vector That a Parameter List Describes
#'
#' @description
#' Returns the location of a multivariate distribution as a numeric vector of
#' length \eqn{p}. The parameters of a multivariate distribution are scalars, so
#' that every generic of the package can index them, and this generic returns
#' the location in the shape of the response.
#'
#' @details
#' The elliptical families return the center of symmetry of the density, which
#' is the mean for a gaussian and, for a Student t, the mean only when
#' \eqn{\nu > 1}; they register [mv_leading_location()]. The Dirichlet returns
#' its mean vector, and the multinomial its mean vector of counts, the
#' probabilities times the number of trials. The method of the class signals an
#' error for a family without a method of its own, because the first \eqn{p}
#' parameters of a family are not in general a location. [base::mean()] is the
#' generic for the moment.
#'
#' @param distrib An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, in which each component
#'   is a single number. The generic aligns it before dispatch.
#'
#' @return A numeric vector of length \eqn{p}. The elliptical families name it
#'   `v1`, ..., `vp` after the coordinates of the response; the Dirichlet and
#'   the multinomial return it without names.
#'
#' @section Notation:
#' \eqn{\mu} is the location and \eqn{p} the dimension of one observation.
#'
#' @seealso [mv_sigma()] for the matrix, [mv_leading_location()] for the
#'   implementation that the elliptical families use, [base::mean()] for the
#'   moment, and [mv_location.multivariate_distrib()] for the method of the
#'   class.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' mv_location(d, theta)
#'
#' # A Student t has a location at every nu, and a mean only above nu = 1.
#' t2 <- mvstudent_t1_distrib(2)
#' th <- c(theta, list(nu = 0.8))
#' rbind(location = mv_location(t2, th), mean = mean(t2, th))
#'
#' # A Dirichlet returns its mean vector, a point of the simplex.
#' mv_location(dirichlet_distrib(3),
#'             list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8))
#'
#' @export
mv_location <- S7::new_generic("mv_location", "distrib", function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' @title Location of a Family Without a Method
#' @name mv_location.multivariate_distrib
#'
#' @description
#' Signals an error. The first \eqn{p} parameters of a family are not a
#' location in general, so the class has no default. The shipped families
#' register their own methods: the elliptical families register
#' [mv_leading_location()], and the Dirichlet and the multinomial register
#' methods that return their mean vector. A family written elsewhere registers
#' its own method.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param theta A named list of parameters. It is not read, because the error is
#'   signalled first.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [mv_location()] for the generic, [mv_leading_location()] for the
#'   implementation that the elliptical families use, and
#'   [mv_location.MvGaussianDistrib()] for a method.
#'
#' @examples
#' # Every shipped family registers a method, so the method of the class is
#' # called directly.
#' base <- S7::method(mv_location, multivariate_distrib)
#' try(base(dirichlet_distrib(3),
#'          list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)))
#'
#' # The Dirichlet's own method returns its mean vector.
#' mv_location(dirichlet_distrib(3),
#'             list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8))
#'
#' @keywords internal
S7::method(mv_location, multivariate_distrib) <- function(distrib, theta) {
  mv_refuse(
    distrib, "mv_location",
    "this family has no location parameter. A family that has one registers a method."
  )
}

#' @title The First p Parameters, Read as a Location
#'
#' @description
#' Reads the first \eqn{p} components of `theta` as the location vector. It is
#' the implementation with which both elliptical families register
#' [mv_location()]. It is correct where the leading \eqn{p} parameters of a
#' family are its location, which is a property of the family and not of the
#' shape of `theta`.
#'
#' @param distrib A [multivariate_distrib()] object whose leading \eqn{p}
#'   parameters are the location.
#' @param theta A named list of parameters, already aligned.
#'
#' @return A numeric vector of length `distrib@n_dim`, named `v1`, ..., `vp`.
#'
#' @seealso [mv_location()] for the generic and
#'   [mv_location.MvGaussianDistrib()] for a registration that uses this.
#'
#' @examples
#' d <- mvgaussian1_distrib(3)
#' theta <- as.list(stats::setNames(c(1, -2, 0.5, rep(0, 6)), d@params))
#' th <- distributions7:::align_theta(d, theta)
#' distributions7:::mv_leading_location(d, th)
#'
#' # mv_location() returns the same vector for this family, because its method
#' # is this function.
#' mv_location(d, theta)
#'
#' @keywords internal
mv_leading_location <- function(distrib, theta) {
  v <- mv_flat_theta(distrib, theta)
  stats::setNames(
    unname(v[seq_len(distrib@n_dim)]), paste0("v", seq_len(distrib@n_dim))
  )
}

#' @title The Matrix That a Parameter List Describes
#'
#' @description
#' Returns the matrix that describes the spread of a multivariate distribution,
#' as a \eqn{p \times p} numeric matrix: the covariance for the gaussian, the
#' Dirichlet and the multinomial, and the scale matrix for the Student t. The
#' Student t's covariance is \eqn{\nu\Sigma/(\nu-2)} for \eqn{\nu > 2} and is
#' infinite for \eqn{\nu \le 2}, while its scale matrix exists at every
#' \eqn{\nu}; [variance()] is the generic for the second moment.
#'
#' @details
#' Where the parametrization carries the precision, as for
#' `mvgaussian2_distrib()` and `mvstudent_t2_distrib()`, the matrix is inverted
#' here, so the result is the covariance or the scale matrix either way.
#'
#' The class has no registered method, so a family without a method of its
#' own fails in dispatch.
#'
#' @param distrib An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, in which each component
#'   is a single number. The generic aligns it before dispatch.
#'
#' @return A \eqn{p \times p} symmetric numeric matrix: positive definite, with
#'   both dimnames `v1`, ..., `vp`, for the gaussian and the Student t; of rank
#'   \eqn{p - 1} and without dimnames for the Dirichlet and the multinomial,
#'   whose coordinates sum to a constant.
#'
#' @section Notation:
#' \eqn{\Sigma} is the matrix that the parametrization carries, \eqn{\nu} the
#' degrees of freedom of a Student t and \eqn{p} the dimension.
#'
#' @seealso [variance()] for the moment, [mv_location()] for the location,
#'   [mv_summary()] for the standard deviations and correlations of a fit, and
#'   [mv_sigma.MvGaussianDistrib()] and [mv_sigma.MvStudentTDistrib()] for the
#'   methods of the gaussian and the Student t.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' mv_sigma(d, theta)
#'
#' # For a gaussian the matrix is the variance; for a Student t it is not.
#' t2 <- mvstudent_t1_distrib(2)
#' th <- c(theta, list(nu = 6))
#' all.equal(mv_sigma(d, theta), variance(d, theta))
#' all.equal(variance(t2, th), (6 / 4) * mv_sigma(t2, th))
#'
#' # At two degrees of freedom the covariance is infinite and the scale matrix
#' # is finite.
#' th2 <- th; th2$nu <- 2
#' c(scale = mv_sigma(t2, th2)[1, 1], covariance = variance(t2, th2)[1, 1])
#'
#' @export
mv_sigma <- S7::new_generic("mv_sigma", "distrib", function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  S7::S7_dispatch()
})

#' @title A Marginal of a Multivariate Distribution
#'
#' @description
#' Returns the distribution of a subset of the coordinates, together with its
#' parameters. For a subset \eqn{A},
#' \deqn{f_A(y_A; \theta) = \int_{\mathbb{R}^{\lvert A^{c}\rvert}}
#'     f(y; \theta)\, \mathrm{d}y_{A^{c}},}
#' which the shipped families have in closed form: any subset for the gaussian
#' and the Student t, and a single coordinate for the Dirichlet and the
#' multinomial. The method of the class signals an error, because a quadrature
#' over the discarded coordinates would give a different object under the same
#' name.
#'
#' @details
#' # What each family returns
#'
#' A gaussian's marginal is a gaussian with \eqn{\mu_A} and \eqn{\Sigma_{AA}}.
#' A Student t's is a Student t with \eqn{\mu_A}, \eqn{\Sigma_{AA}} and the same
#' \eqn{\nu}: conditionally on the mixing variable of the scale-mixture
#' representation the vector is gaussian, whose marginal is gaussian, and
#' mixing again gives a Student t. A Dirichlet's marginal of one coordinate is a
#' beta with the same concentration \eqn{\phi}, and a multinomial's is a
#' binomial. For several coordinates these two signal an error, because the
#' subvector has the family's form only once a coordinate for the remaining
#' mass is added.
#'
#' # The returned object
#'
#' The method returns a new distribution of the reduced dimension, whose
#' parameters are its own. The elliptical families return an object of their own
#' class on an unstructured log-Cholesky matrix, carrying the covariance or the
#' precision as the parent does, so a gaussian's one-coordinate marginal is a
#' one-dimensional `MvGaussian1Distrib` and signals an error for the
#' distribution function. The simplex-valued families return univariate objects,
#' a beta and a binomial, which have every method of a univariate family.
#'
#' # The plot
#'
#' A panel of the pairs plot of [plot.multivariate_distrib()] is a marginal of
#' one or two coordinates, so the plot is available for the gaussian and the
#' Student t, and signals the error above for the Dirichlet and the
#' multinomial.
#'
#' @param distrib An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, in which each component
#'   is a single number. The generic aligns it before dispatch.
#' @param which A vector of distinct whole numbers in `1:distrib@n_dim`, the
#'   coordinates to keep. The generic checks it before dispatch and signals an
#'   error naming the range for an empty, missing, non-integer, duplicated or
#'   out-of-range value. The Dirichlet and the multinomial accept a single
#'   coordinate.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A named list with `distrib`, the marginal distribution object, and
#'   `theta`, its parameters as a named list.
#'
#' @section Notation:
#' \eqn{A} is the retained subset of coordinates, \eqn{A^c} its complement,
#' \eqn{\mu} the location, \eqn{\Sigma} the matrix that the family carries,
#' \eqn{\nu} a Student t's degrees of freedom and \eqn{\phi} a Dirichlet's
#' concentration.
#'
#' @seealso [plot.multivariate_distrib()], whose panels are these marginals,
#'   and [mv_marginal.MvGaussianDistrib()],
#'   [mv_marginal.MvStudentTDistrib()] and
#'   [mv_marginal.multivariate_distrib()] for the methods and the error.
#'
#' @examples
#' d <- mvgaussian1_distrib(3)
#' theta <- as.list(stats::setNames(
#'   c(1, 2, 3, 0, 0, 0, 0.3, 0.2, 0.1), d@params))
#'
#' # The marginal of two coordinates is a two-dimensional gaussian on the
#' # corresponding block of the covariance.
#' m <- mv_marginal(d, theta, c(1, 2))
#' mv_sigma(m$distrib, m$theta)
#' mv_sigma(d, theta)[1:2, 1:2]
#'
#' # A Student t's marginal keeps the same degrees of freedom.
#' t3 <- mvstudent_t1_distrib(3)
#' th <- as.list(stats::setNames(c(unlist(theta), 5), t3@params))
#' c(full = th$nu, marginal = mv_marginal(t3, th, c(1, 3))$theta$nu)
#'
#' # A Dirichlet's is a beta with the same concentration, a univariate object.
#' b <- mv_marginal(dirichlet_distrib(3),
#'                  list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
#' class(b$distrib)[1]
#' b$theta
#'
#' # 'which' is checked before dispatch.
#' try(mv_marginal(d, theta, c(1, 1)))
#' try(mv_marginal(d, theta, 4))
#'
#' @export
mv_marginal <- S7::new_generic("mv_marginal", "distrib",
  function(distrib, theta, which, ...) {
    theta <- align_theta(distrib, theta)
    whole <- is.numeric(which) && all(is.finite(which)) &&
      all(which == round(which))
    which <- if (whole) as.integer(which) else NA_integer_
    if (!length(which) || anyNA(which) ||
      any(which < 1L) || any(which > distrib@n_dim) || anyDuplicated(which)) {
      stop(sprintf(
        "'which' must be distinct coordinates in 1:%d.", distrib@n_dim
      ), call. = FALSE)
    }
    S7::S7_dispatch()
  }
)

#' @title Marginal of a Family Without a Closed Form
#' @name mv_marginal.multivariate_distrib
#'
#' @description
#' Signals an error. Integrating the other coordinates out has no general
#' closed form, and a numerical marginal would be a grid of values and not a
#' distribution object to which the generics apply. A family whose marginals
#' are known registers its own method, as every shipped family does.
#'
#' @param distrib A [multivariate_distrib()] object with no method of its own.
#' @param theta A named list of parameters. It is not read, because the error is
#'   signalled first.
#' @param which An integer vector of coordinates, already validated by the
#'   generic. It is not read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [mv_marginal()] for the generic and what each shipped family
#'   returns, and [plot.multivariate_distrib()], whose panels are marginals.
#'
#' @examples
#' # Every shipped family registers a method, so the method of the class is
#' # called directly.
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0)
#' base <- S7::method(mv_marginal, multivariate_distrib)
#' try(base(d, distributions7:::align_theta(d, theta), 1L))
#'
#' # The family's own method returns the marginal.
#' mv_marginal(d, theta, 1)$distrib@n_dim
#'
#' @keywords internal
S7::method(mv_marginal, multivariate_distrib) <- function(distrib, theta, which, ...) {
  mv_refuse(
    distrib, "mv_marginal",
    "integrating out the other coordinates has no closed form for this family, and a numerical marginal would not be the same object."
  )
}


#' @title Skewness of a Multivariate Distribution
#' @name skewness.multivariate_distrib
#'
#' @description
#' Signals an error. A scalar skewness of a vector response is a choice among
#' several quantities: Mardia's, Malkovich and Afifi's, and the vector of the
#' coordinatewise marginal skewnesses are different numbers, and the method
#' does not return one of them under the bare name.
#'
#' @details
#' [mv_marginal()] gives a coordinatewise skewness where the marginal is a
#' univariate object. A Dirichlet's marginal is a beta and has a skewness; an
#' elliptical family's marginal is a one-dimensional object of its own class
#' and signals the same error. Applying the generic to a marginal specifies
#' which quantity is computed.
#'
#' @param x A [multivariate_distrib()] object.
#' @param theta A named list of parameters. It is not read, because the error is
#'   signalled first.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [kurtosis.multivariate_distrib()] for the fourth moment,
#'   [mv_marginal()] for the coordinatewise route, and [skewness()] for the
#'   generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' try(skewness(d, theta))
#'
#' # A Dirichlet's marginal is a univariate beta and has a skewness. An
#' # elliptical family's marginal is a one-dimensional object of its own class
#' # and signals the same error.
#' m <- mv_marginal(dirichlet_distrib(3),
#'                  list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
#' skewness(m$distrib, m$theta)
#' try(skewness(mv_marginal(d, theta, 1)$distrib,
#'              mv_marginal(d, theta, 1)$theta))
#'
#' @keywords internal
S7::method(skewness, multivariate_distrib) <- function(x, theta, ...) {
  mv_refuse(
    x, "skewness",
    "a vector response has no single skewness: Mardia's, Malkovich and Afifi's and the vector of coordinatewise marginal skewnesses are different quantities. Use a univariate family instead."
  )
}

#' @title Kurtosis of a Multivariate Distribution
#' @name kurtosis.multivariate_distrib
#'
#' @description
#' Signals an error, for the same reason as [skewness.multivariate_distrib()]:
#' a scalar kurtosis of a vector response is a choice among Mardia's,
#' Malkovich and Afifi's and the vector of the coordinatewise marginal
#' kurtoses, which are different numbers. A coordinatewise kurtosis is the
#' kurtosis of a marginal, where the marginal is a univariate object: a
#' Dirichlet's is a beta and has one, while an elliptical family's is a
#' one-dimensional object of its own class and signals the same error.
#'
#' @param x A [multivariate_distrib()] object.
#' @param theta A named list of parameters. Not read: the error is signalled
#'   first.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return Does not return: it signals an error naming the family.
#'
#' @seealso [skewness.multivariate_distrib()] for the third moment,
#'   [mv_marginal()] for the coordinatewise route, and [kurtosis()] for the
#'   generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5, nu = 6)
#' try(kurtosis(d, theta))
#'
#' # A Student t's marginal is a one-dimensional MvStudentT1Distrib, which
#' # signals the same error.
#' try(kurtosis(mv_marginal(d, theta, 1)$distrib,
#'              mv_marginal(d, theta, 1)$theta))
#'
#' # A Dirichlet's marginal is a univariate beta and has a kurtosis.
#' b <- mv_marginal(dirichlet_distrib(3),
#'                  list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8), 1)
#' kurtosis(b$distrib, b$theta)
#'
#' @keywords internal
S7::method(kurtosis, multivariate_distrib) <- function(x, theta, ...) {
  mv_refuse(
    x, "kurtosis",
    "a vector response has no single kurtosis: Mardia's, Malkovich and Afifi's and the vector of coordinatewise marginal kurtoses are different quantities. Use a univariate family instead."
  )
}
