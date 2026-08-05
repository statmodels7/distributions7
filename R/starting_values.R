#' @include distrib.R generics.R multivariate.R mvgaussian_distrib.R mvstudent_t_distrib.R
NULL

#' A Starting Value Drawn from the Data
#'
#' @description
#' Returns a starting value for \code{\link{fit_distrib}}, computed from the
#' response rather than guessed.
#'
#' @details
#' A starting value that ignores the data is a starting value that can be
#' arbitrarily far from the answer, and how far decides whether the fit takes
#' one step or never arrives. A four-dimensional gaussian fitted to the iris
#' measurements is the plain case: started at the origin of the unconstrained
#' scale, which is a unit covariance and a zero mean, Newton with the expected
#' information spends five hundred iterations and stops at a log-likelihood of
#' \eqn{-836}; started at the sample mean and the sample covariance it
#' converges in one iteration to \eqn{-379.9146}, which is the exact maximum.
#' Nothing about the arithmetic changed.
#'
#' The default method returns random parameters, as before, so a distribution
#' that says nothing loses nothing. A family that can do better says so by
#' registering a method: an exact maximum likelihood estimator where one is
#' known, a method-of-moments estimator otherwise, or the estimate of a simpler
#' family the harder one contains.
#'
#' A method may return several starting values, as a list, and
#' \code{fit_distrib()} will try each; the first is the one it prefers.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y The response.
#' @param n_start How many starting values are wanted. A method free to supply
#'   only one may ignore it.
#' @param ... Passed to methods.
#'
#' @return A list of named parameter lists, on the \strong{parameter} scale.
#'
#' @seealso \code{\link{fit_distrib}}, \code{\link{generate_random_theta}}
#'
#' @examples
#' set.seed(1)
#' d <- mvgaussian_distrib(2)
#' y <- distrib_rng(d, 200, list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0,
#'                               sigma_log_L2 = 0, sigma_L2.1 = 0.5))
#'
#' # the gaussian knows its own maximum likelihood estimate, so the fit starts
#' # there and has nothing left to do
#' start <- distrib_start(d, y)[[1]]
#' mv_location(d, start)
#' colMeans(y)
#'
#' @export
distrib_start <- S7::new_generic("distrib_start", "distrib",
  function(distrib, y, n_start = 5L, ...) {
    S7::S7_dispatch()
  }
)

#' @title Random Starting Values
#' @name distrib_start.distrib
#' @description
#' The default: \code{n_start} draws from \code{\link{generate_random_theta}},
#' which uses the parameter domains and not the data. A family with a better
#' idea registers its own method.
#' @param distrib A \code{\link{distrib}} object.
#' @param y The response, unused here.
#' @param n_start How many to draw.
#' @param ... Unused.
#' @return A list of named parameter lists.
#' @keywords internal
S7::method(distrib_start, distrib) <- function(distrib, y, n_start = 5L, ...) {
  lapply(seq_len(max(1L, n_start)), function(i) generate_random_theta(distrib))
}


#' The Moment Estimates a Multivariate Family Starts From
#'
#' @description
#' Returns the sample mean and the sample covariance of a multivariate
#' response, with the covariance made safely positive definite.
#'
#' @details
#' A sample covariance is singular when there are fewer observations than
#' coordinates, and nearly singular when two coordinates almost coincide.
#' Either way a structure cannot be inverted onto it, so the eigenvalues are
#' floored at a small multiple of the largest before the matrix is handed on.
#' The floor moves a starting value, which is allowed to be approximate; it
#' would not be allowed anywhere the answer is reported.
#'
#' @param y The response, an \eqn{n \times p} matrix.
#' @param p The dimension.
#'
#' @return A list with \code{mu} and \code{sigma}.
#'
#' @keywords internal
mv_moment_start <- function(y, p) {
  y <- as.matrix(y)
  mu <- colMeans(y)
  n <- nrow(y)
  s <- if (n > 1L) crossprod(sweep(y, 2L, mu)) / n else diag(p)
  ev <- eigen(s, symmetric = TRUE)
  floor_at <- 1e-6 * max(ev$values, 1)
  if (min(ev$values) < floor_at) {
    s <- ev$vectors %*% diag(pmax(ev$values, floor_at), p) %*% t(ev$vectors)
    s <- (s + t(s)) / 2
  }
  list(mu = unname(mu), sigma = unname(s))
}


#' Project a Matrix onto What a Structure Can Represent
#'
#' @description
#' Returns the free vector of the matrix parameter whose matrix is closest to the one
#' supplied, or the matrix parameter's own inverse map when it has one.
#'
#' @details
#' \code{\link[parameters7]{param_free}} is exact or refused: a structure that
#' cannot represent the matrix says so rather than returning something
#' plausible. That is the right contract for reporting an estimate and the
#' wrong one for choosing where to begin, so a refusal here falls back to a
#' short numerical search over the free values, which is allowed to be
#' approximate because a starting value is.
#'
#' @param s A \pkg{parameters7} structure.
#' @param m The matrix to represent.
#'
#' @return A numeric vector of length \code{s@n_free}.
#'
#' @keywords internal
param_free_or_fit <- function(s, m) {
  eta <- tryCatch(parameters7::param_free(s, m), error = function(e) NULL)
  if (!is.null(eta) && all(is.finite(eta))) return(unname(eta))

  # Least squares on the entries, from the zero vector. This is a starting
  # value for a starting value, and its accuracy does not matter.
  obj <- function(v) {
    mm <- tryCatch(parameters7::param_value(s, v), error = function(e) NULL)
    if (is.null(mm) || anyNA(mm)) return(1e10)
    sum((mm - m)^2)
  }
  r <- tryCatch(
    stats::optim(rep(0, s@n_free), obj, method = "BFGS",
                 control = list(maxit = 200))$par,
    error = function(e) rep(0, s@n_free)
  )
  unname(r)
}


#' @title The Maximum Likelihood Estimate as a Starting Value
#' @name distrib_start.MvGaussianDistrib
#' @description
#' The sample mean and the sample covariance, which for an unstructured
#' covariance are the maximum likelihood estimate itself, so the fit begins at
#' the answer and confirms it in one step. For a structured covariance they are
#' the closest thing the matrix parameter can represent, which is a good deal nearer
#' than the origin.
#' @details
#' When the matrix parameter parametrizes the precision the sample covariance is
#' inverted first, since that is the matrix the matrix parameter has to represent.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y The response.
#' @param n_start Unused; one starting value is enough when it is the estimate.
#' @param ... Unused.
#' @return A list with one named parameter list.
#' @keywords internal
S7::method(distrib_start, MvGaussianDistrib) <- function(distrib, y, n_start = 5L, ...) {
  p <- distrib@n_dim
  m <- mv_moment_start(y, p)
  target <- if (isTRUE(distrib@inverted)) solve(m$sigma) else m$sigma
  eta <- param_free_or_fit(distrib@param, target)
  list(as.list(stats::setNames(c(m$mu, eta), distrib@params)))
}


#' @title The Gaussian Estimate as a Starting Value for a t
#' @name distrib_start.MvStudentTDistrib
#' @description
#' The sample mean and the sample covariance, with the degrees of freedom set
#' where the family is heavy tailed but its second moment exists. The gaussian
#' estimate is the limit of this family as \eqn{\nu} grows, so it is the right
#' place to start looking for a finite one.
#' @details
#' The scale matrix is the covariance divided by \eqn{\nu/(\nu-2)}, and that
#' factor is applied, since a starting value that confused the two would begin
#' with a scale a third too large.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param y The response.
#' @param n_start Unused.
#' @param ... Unused.
#' @return A list with one named parameter list.
#' @keywords internal
S7::method(distrib_start, MvStudentTDistrib) <- function(distrib, y, n_start = 5L, ...) {
  p <- distrib@n_dim
  m <- mv_moment_start(y, p)
  nu0 <- 8
  eta <- param_free_or_fit(distrib@param, m$sigma * (nu0 - 2) / nu0)
  list(as.list(stats::setNames(c(m$mu, eta, nu0), distrib@params)))
}
