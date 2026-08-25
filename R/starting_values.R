#' @include distrib.R generics.R multivariate.R mvgaussian_distrib.R mvstudent_t_distrib.R pig1_distrib.R pig2_distrib.R
NULL

#' @title A Starting Value Drawn From the Data
#'
#' @description
#' Returns the starting values [fit_distrib()] begins from, computed from the
#' response wherever the family knows how. The result is a list of named
#' parameter lists on the **parameter** scale; the fit tries them in order and
#' stops at the first that converges, so the first element is the one the
#' family prefers.
#'
#' @details
#' A starting value that ignores the data can be arbitrarily far from the
#' answer, and how far decides whether the fit takes one step or never
#' arrives. A four-dimensional Gaussian fitted to the iris measurements is the
#' plain case: started at the origin of the unconstrained scale, which is a
#' unit covariance and a zero mean, Newton with the expected information spends
#' five hundred iterations and stops at a log-likelihood of \eqn{-836};
#' started at the sample mean and the sample covariance it converges in one
#' iteration to \eqn{-379.9146}, the exact maximum. Nothing about the
#' arithmetic changed.
#'
#' The base method draws each parameter from its own domain and never reads
#' `y`, so a family that registers nothing loses nothing. A family with a
#' better estimator registers a method: an exact maximum likelihood estimator
#' where one is known, a method-of-moments estimator otherwise, or the estimate
#' of a simpler family the harder one contains. The two univariate methods
#' route through [start_from_moments()], which uses the family's own moment
#' inversion where [moment_estimates()] has one and reads
#' `params_interpretation` where it does not.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y The response: a numeric vector, or an \eqn{n \times p} matrix for a
#'   multivariate family. A method may ignore it, and the base one does.
#' @param n_start How many starting values are wanted, a single positive
#'   integer. Defaults to 5. A method that returns its own estimate returns one
#'   and ignores this, which the two multivariate methods and the two
#'   Poisson-inverse Gaussian ones do.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A list of named parameter lists on the parameter scale, each
#'   complete: one component per entry of `distrib@params`, named and ordered
#'   as `distrib@params`. Every value is strictly inside its parameter's
#'   bounds, which the validator treats as open.
#'
#' @examples
#' set.seed(1)
#' d <- mvgaussian_distrib(2)
#' y <- distrib_rng(d, 200, list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0,
#'                               sigma_log_L2 = 0, sigma_L2.1 = 0.5))
#'
#' # An unstructured covariance has a closed-form estimate, so the fit starts
#' # at the answer and has nothing left to do.
#' start <- distrib_start(d, y)[[1]]
#' rbind(start = mv_location(d, start), sample = colMeans(y))
#'
#' # A univariate family asked for several starts gets one from the data and
#' # the rest at random, so a caller asking for more still explores.
#' set.seed(2)
#' s <- distrib_start(gaussian1_distrib(), rnorm(200, 900, 170), n_start = 3)
#' length(s)
#' vapply(s, function(th) unlist(th), numeric(2))
#'
#' @seealso [fit_distrib()], which calls this;
#'   [start_from_moments()] for the univariate route;
#'   [moment_estimates()] for the family-by-family inversions;
#'   [generate_random_theta()] for the draw the base method makes.
#' @export
distrib_start <- S7::new_generic("distrib_start", "distrib",
  function(distrib, y, n_start = 5L, ...) {
    S7::S7_dispatch()
  }
)

#' @title Random Starting Values
#' @name distrib_start.distrib
#'
#' @description
#' The base method: `n_start` independent draws from
#' [generate_random_theta()], which samples each parameter inside its own
#' domain and never looks at the response. It is what a family gets when it
#' registers no method of its own, and it is adequate wherever the parameters
#' are of order one whatever the data, which is true of a shape, a dispersion
#' or a probability and false of a location or a scale.
#'
#' @param distrib A [distrib()] object.
#' @param y The response. Unused here, and accepted only because the generic
#'   passes it.
#' @param n_start How many starting values to draw, a single positive integer.
#'   Defaults to 5. A value below 1 is raised to 1.
#' @param ... Unused.
#'
#' @return A list of `n_start` named parameter lists on the parameter scale.
#'
#' @seealso [distrib_start()] for the generic;
#'   [start_from_moments()], which the univariate classes register instead;
#'   [generate_random_theta()] for the draw itself.
#' @keywords internal
S7::method(distrib_start, distrib) <- function(distrib, y, n_start = 5L, ...) {
  lapply(seq_len(max(1L, n_start)), function(i) generate_random_theta(distrib))
}


#' @title The Moment Estimates a Multivariate Family Starts From
#'
#' @description
#' Returns the sample mean and the sample covariance of a multivariate
#' response, with the covariance made safely positive definite. The covariance
#' divides by \eqn{n} rather than \eqn{n - 1}, which is the maximum likelihood
#' estimator and therefore the point a Gaussian fit is looking for.
#'
#' @details
#' A sample covariance is singular when there are fewer observations than
#' coordinates, and nearly singular when two coordinates almost coincide.
#' Either way a `parameters7` structure cannot be inverted onto it, so the
#' eigenvalues are floored at a small multiple of the largest before the matrix
#' is handed on. The floor moves a starting value, which is allowed to be
#' approximate; it would not be allowed anywhere the answer is reported.
#'
#' A response with one row has no covariance at all, and the identity is
#' returned in its place.
#'
#' @param y The response, an \eqn{n \times p} matrix, or anything
#'   [base::as.matrix()] turns into one.
#' @param p The dimension, a single positive integer. Used only for the
#'   one-row fallback.
#'
#' @return A list with two components: `mu`, a numeric vector of length
#'   \eqn{p}, and `sigma`, a \eqn{p \times p} positive definite matrix.
#'
#' @seealso [distrib_start.MvGaussianDistrib()] and
#'   [distrib_start.MvStudentTDistrib()], the two callers;
#'   [param_free_or_fit()], which carries the matrix onto a structure.
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


#' @title A Free Vector Representing a Matrix, Exactly or Approximately
#'
#' @description
#' Returns the free vector of a `parameters7` structure whose matrix is `m`.
#' `parameters7::param_free()` answers exactly where the structure can
#' represent the matrix, which for an unstructured covariance is always. Where
#' it cannot, that call signals an error and this falls back to least squares
#' on the entries, started from the zero vector. A compound-symmetric
#' structure asked for an arbitrary covariance is the ordinary case.
#'
#' The fallback is a starting value for a starting value and its accuracy does
#' not matter. It is capped at 200 BFGS iterations and returns the zero vector
#' if even that fails, so the caller always receives a usable vector.
#'
#' @param s A \pkg{parameters7} structure, supplying `n_free`, `param_free()`
#'   and `param_value()`.
#' @param m The matrix to represent, of the structure's own dimension.
#'
#' @return An unnamed numeric vector of length `s@n_free`.
#'
#' @seealso [mv_moment_start()], which supplies `m`;
#'   [parameters7::param_free()] for the exact route.
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
#'
#' @description
#' Returns the sample mean and the sample covariance, carried onto the
#' distribution's matrix parametrization. For an unstructured covariance those
#' **are** the maximum likelihood estimate, so the fit begins at the answer
#' and confirms it in one step; for a structured one they are the closest
#' matrix the parametrization can represent, which is a great deal nearer than
#' the origin.
#'
#' One starting value is returned whatever `n_start` asks for. There is nothing
#' to explore when the first point is the estimate.
#'
#' @details
#' Where the object parametrizes the **precision**, the sample covariance is
#' inverted first: the matrix the structure has to represent is
#' \eqn{\hat\Sigma^{-1}}, not \eqn{\hat\Sigma}.
#'
#' The covariance divides by \eqn{n}, which is the maximum likelihood
#' estimator, and its eigenvalues are floored by [mv_moment_start()] so that a
#' singular sample covariance still produces a usable point. Carrying it onto
#' the structure goes through [param_free_or_fit()], which is exact where
#' `parameters7::param_free()` succeeds and a least-squares fit where it does
#' not.
#'
#' @param distrib An [MvGaussianDistrib()] object.
#' @param y The response, an \eqn{n \times p} matrix.
#' @param n_start Ignored: one starting value is returned, and it is the
#'   estimate.
#' @param ... Unused.
#'
#' @return A list of length 1 holding one named parameter list: the \eqn{p}
#'   location components followed by the structure's free values, named and
#'   ordered as `distrib@params`.
#'
#' @seealso [distrib_start()] for the generic;
#'   [distrib_start.MvStudentTDistrib()], which starts from this;
#'   [mv_moment_start()] and [param_free_or_fit()] for the two steps.
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
#'
#' @description
#' Returns the sample mean and the sample covariance with the degrees of
#' freedom set at 8, which is heavy-tailed enough to be worth fitting and light
#' enough for the second moment to exist. The Gaussian estimate is this
#' family's limit as \eqn{\nu} grows, so it is where a finite \eqn{\nu} is
#' looked for from.
#'
#' One starting value is returned whatever `n_start` asks for.
#'
#' @details
#' The scale matrix is not the covariance. For a \eqn{t} with \eqn{\nu}
#' degrees of freedom \eqn{\mathrm{Var}(Y) = \nu\Sigma/(\nu-2)}, so the sample
#' covariance is multiplied by \eqn{(\nu_0-2)/\nu_0 = 0.75} before it is
#' carried onto the structure. A starting value that confused the two would
#' begin with a scale a third too large.
#'
#' @param distrib An [MvStudentTDistrib()] object.
#' @param y The response, an \eqn{n \times p} matrix.
#' @param n_start Ignored: one starting value is returned.
#' @param ... Unused.
#'
#' @return A list of length 1 holding one named parameter list: the \eqn{p}
#'   location components, the structure's free values, then `nu` at 8, named
#'   and ordered as `distrib@params`.
#'
#' @seealso [distrib_start()] for the generic;
#'   [distrib_start.MvGaussianDistrib()], the limit this starts from;
#'   [mv_sigma()] for the scale matrix and [variance()] for the covariance.
#' @keywords internal
S7::method(distrib_start, MvStudentTDistrib) <- function(distrib, y, n_start = 5L, ...) {
  p <- distrib@n_dim
  m <- mv_moment_start(y, p)
  nu0 <- 8
  eta <- param_free_or_fit(distrib@param, m$sigma * (nu0 - 2) / nu0)
  list(as.list(stats::setNames(c(m$mu, eta, nu0), distrib@params)))
}


#' @title Poisson-Inverse Gaussian Starting Values
#' @name distrib_start.Pig1Distrib
#'
#' @description
#' Returns the method-of-moments estimate. The family has mean \eqn{\mu} and
#' variance \eqn{\mu + \sigma\mu^2}, so setting both equal to the sample gives
#' \eqn{\hat\mu = \bar y} and \eqn{\hat\sigma = (s^2 - \bar y)/\bar y^2}
#' directly, with no root to find.
#'
#' Both are floored: \eqn{\mu} just above zero, and \eqn{\sigma} at
#' \eqn{10^{-3}} when the sample is **underdispersed** and the inversion
#' returns a negative number. A Poisson sample is the case that produces it,
#' and \eqn{10^{-3}} is where this family is nearly Poisson, which is the right
#' place to start from there.
#'
#' @param distrib A `Pig1Distrib` object.
#' @param y A numeric vector of counts.
#' @param n_start Ignored: one moment start is returned.
#' @param ... Unused.
#'
#' @return A list of length 1 holding one named parameter list with components
#'   `mu` and `sigma`.
#'
#' @examples
#' set.seed(4)
#' d <- pig1_distrib()
#' y <- distrib_rng(d, 5000, list(mu = 4, sigma = 0.5))
#'
#' # The inversion recovers the parameters it was drawn from.
#' unlist(distrib_start(d, y)[[1]])
#'
#' # It is exactly the sample moments, read through mean and variance.
#' c(mu = mean(y), sigma = (var(y) - mean(y)) / mean(y)^2)
#'
#' # An underdispersed sample would give a negative sigma, so it is floored
#' # where the family is nearly Poisson.
#' set.seed(5)
#' unlist(distrib_start(d, rpois(2000, 4))[[1]])
#'
#' @seealso [pig1_distrib()] for the family;
#'   [distrib_start.Pig2Distrib()], the same estimate on the orthogonal chart;
#'   [distrib_start()] for the generic.
S7::method(distrib_start, Pig1Distrib) <- function(distrib, y, n_start = 5L, ...) {
  mu <- max(mean(y), 1e-8)
  list(list(mu = mu, sigma = max((stats::var(y) - mu) / mu^2, 1e-3)))
}

#' @title Orthogonal Poisson-Inverse Gaussian Starting Values
#' @name distrib_start.Pig2Distrib
#'
#' @description
#' Returns the moment estimate of [`pig1()`][distrib_start.Pig1Distrib] carried
#' onto this chart. The two families are the same law, so the estimate is the
#' same estimate: \eqn{\hat\mu = \bar y} and
#' \eqn{\hat\sigma = (s^2 - \bar y)/\bar y^2} as before, then
#' \deqn{\alpha = \frac{\sqrt{1 + 2\sigma\mu}}{\sigma}.}
#'
#' The same two floors apply, so an underdispersed sample gives a very large
#' \eqn{\alpha}, which is where this parametrization puts the Poisson limit.
#'
#' @param distrib A `Pig2Distrib` object.
#' @param y A numeric vector of counts.
#' @param n_start Ignored: one moment start is returned.
#' @param ... Unused.
#'
#' @return A list of length 1 holding one named parameter list with components
#'   `mu` and `alpha`.
#'
#' @examples
#' set.seed(4)
#' d <- pig2_distrib()
#' alpha0 <- sqrt(1 + 2 * 0.5 * 4) / 0.5
#' y <- distrib_rng(d, 5000, list(mu = 4, alpha = alpha0))
#'
#' # The inversion recovers what the sample was drawn from.
#' rbind(start = unlist(distrib_start(d, y)[[1]]),
#'       truth = c(mu = 4, alpha = alpha0))
#'
#' # It is pig1's estimate carried through the map.
#' s <- (var(y) - mean(y)) / mean(y)^2
#' c(alpha = sqrt(1 + 2 * s * mean(y)) / s)
#'
#' @seealso [pig2_distrib()] for the family;
#'   [distrib_start.Pig1Distrib()] for the estimate this maps;
#'   [distrib_start()] for the generic.
S7::method(distrib_start, Pig2Distrib) <- function(distrib, y, n_start = 5L, ...) {
  mu <- max(mean(y), 1e-8)
  sg <- max((stats::var(y) - mu) / mu^2, 1e-3)
  list(list(mu = mu, alpha = sqrt(1 + 2 * sg * mu) / sg))
}


#' @title Where a Univariate Family Starts, From the Data
#'
#' @description
#' Returns `n_start` starting values whose first is computed from the response
#' and whose rest are random draws. It is the method both
#' `continuous_distrib` and `discrete_distrib` register, so it is what every
#' univariate family in the package uses.
#'
#' The data-based value comes from [moment_estimates()] where the family has an
#' entry there, which 37 of the 42 univariate families do. The other five fall
#' back to reading `params_interpretation`: a parameter meaning a location is
#' started at the sample median, one meaning a spread at the sample standard
#' deviation or its square, one meaning degrees of freedom at what the sample
#' kurtosis implies, and a shape, a dispersion or a probability keeps its
#' draw, those being of order one whatever the data.
#'
#' @details
#' # What this replaced, and why it was necessary
#' The base method draws each parameter from its own domain and never reads
#' `y`. That is adequate while the response is of order one and fails
#' completely when it is not: on a response of mean 919 and standard deviation
#' 169 the draws are of order one, the first Newton step is taken where the
#' residuals are hundreds of standard deviations wide, and the scale runs to
#' the largest representable double. Measured on a Gaussian, [fit_distrib()]
#' recovers \eqn{N(5, 2)} and \eqn{N(50, 20)} and fails on \eqn{N(500, 200)}:
#' the defect was a threshold in the scale of the data, not in the family.
#'
#' # Degrees of freedom from the kurtosis
#' A \eqn{t} of \eqn{\nu} degrees has excess kurtosis \eqn{6/(\nu-4)}, so a
#' sample kurtosis above 0.05 inverts to \eqn{\nu = 6/\hat\gamma_2 + 4}, capped
#' at 100; a sample no heavier-tailed than a Gaussian starts at 30. Starting
#' large matters. A `student_t2`, whose scale parameter is the standard
#' deviation, has a degenerate ridge at its lower bound of \eqn{\nu = 2} where
#' the scale runs to infinity as \eqn{\nu} falls, and a run started small
#' slides down it: measured on 610 abdominal circumferences, whose excess
#' kurtosis is \eqn{-1.09}, the random draws put \eqn{\nu} between 2.8 and 8
#' and the fit reached the boundary at a log-likelihood of \eqn{-3688.28},
#' while any start of 2.5 or more with the location and scale at their sample
#' values reaches \eqn{-3600.71}, the Gaussian limit and the true maximum.
#'
#' # The clamp, and why only the first start
#' Every value is moved strictly inside its parameter's bounds before it is
#' returned, because a sample median can land exactly on the boundary of a
#' support and the validator treats bounds as open. Only the first start is
#' data-based; the rest stay random, so a caller asking for several still
#' explores, and [fit_distrib()] reaches them when the first fails.
#'
#' @param distrib A univariate distribution, supplying `params`,
#'   `params_interpretation` and `params_bounds`.
#' @param y The response, a numeric vector. Non-finite entries are dropped, and
#'   a sample with fewer than two usable values gets the random draws alone.
#' @param n_start How many starting values are wanted, a single positive
#'   integer. Defaults to 5, and a value below 1 is raised to 1.
#' @param ... Unused.
#'
#' @return A list of `n_start` named parameter lists on the parameter scale,
#'   the first from the data where one could be computed.
#'
#' @seealso [moment_estimates()] for the family-by-family inversions;
#'   [distrib_start()] for the generic;
#'   [clamp_to_bounds()] for the boundary rule;
#'   [fit_distrib()], which consumes the list.
#' @keywords internal
start_from_moments <- function(distrib, y, n_start = 5L, ...) {
  out <- lapply(seq_len(max(1L, n_start)),
                function(i) generate_random_theta(distrib))
  yy <- as.numeric(y)
  yy <- yy[is.finite(yy)]
  if (length(yy) < 2L) return(out)

  # The family's own moment estimate where it has one, which is an estimate
  # rather than a reading of what the parameters are called. What follows is
  # the fallback for the families that do not.
  me <- tryCatch(moment_estimates(distrib, yy), error = function(e) NULL)
  if (!is.null(me)) {
    out[[1L]] <- me
    return(out)
  }

  loc <- stats::median(yy)
  sdev <- stats::sd(yy)
  if (!is.finite(sdev) || sdev <= 0) sdev <- 1
  pos <- yy[yy > 0]
  lmean <- if (length(pos) > 1L) mean(log(pos)) else NA_real_
  lvar <- if (length(pos) > 1L) stats::var(log(pos)) else NA_real_

  # Degrees of freedom from the sample kurtosis, a t of nu degrees having
  # excess kurtosis 6/(nu - 4). Where the sample is no heavier-tailed than a
  # gaussian the answer is "large", and starting large matters: a
  # student_t2, whose scale parameter is the STANDARD DEVIATION, has a
  # degenerate ridge at its lower bound of 2 where the scale runs to
  # infinity as nu falls, and a run started at a small nu slides down it.
  # Measured on 610 abdominal circumferences, whose excess kurtosis is
  # -1.09: from the random draws, nu between 2.8 and 8, the fit reached the
  # boundary at a log-likelihood of -3688.28; from any start of 2.5 or more
  # with the location and scale at their sample values it reaches -3600.71,
  # the gaussian limit, which is the true maximum.
  z <- (yy - mean(yy)) / sdev
  kurt <- mean(z^4) - 3
  ndf <- if (is.finite(kurt) && kurt > 0.05) min(6 / kurt + 4, 100) else 30

  params <- distrib@params
  ip <- tolower(distrib@params_interpretation)
  bnds <- distrib@params_bounds
  th <- out[[1L]]
  for (j in seq_along(params)) {
    p <- params[j]
    what <- if (j <= length(ip)) ip[[j]] else ""
    v <- switch(
      what,
      "mean" = loc,
      "location" = loc,
      "mean direction" = loc,
      "standard deviation" = sdev,
      "scale" = sdev,
      "variance" = sdev^2,
      "mean (log scale)" = lmean,
      "variance (log scale)" = lvar,
      "degrees of freedom" = ndf,
      NA_real_)
    if (!is.finite(v)) next
    b <- bnds[[p]]
    if (!is.null(b) && length(b) == 2L) {
      eps <- .Machine$double.eps
      if (is.finite(b[1L]) && v <= b[1L]) {
        v <- if (b[1L] == 0) max(sdev, 1) * eps else b[1L] + abs(b[1L]) * eps
      }
      if (is.finite(b[2L]) && v >= b[2L]) v <- b[2L] - abs(b[2L]) * eps
    }
    th[[p]] <- v
  }
  out[[1L]] <- th
  out
}

#' @rdname start_from_moments
#' @name distrib_start.continuous_distrib
#' @keywords internal
S7::method(distrib_start, continuous_distrib) <- start_from_moments

#' @rdname start_from_moments
#' @name distrib_start.discrete_distrib
#' @keywords internal
S7::method(distrib_start, discrete_distrib) <- start_from_moments


#' @title Method-of-Moments Estimates, Family by Family
#'
#' @description
#' Returns the parameters a family's own first two moments imply for a sample,
#' in closed form where the inversion has one, and `NULL` where this family has
#' no entry. **37 of the 42 univariate families have one.**
#'
#' A starting value should be an estimate. For most families the moment
#' estimate is one line: the sample mean and variance are set equal to the
#' family's own and the pair is inverted. [fit_distrib()] then refines it by
#' maximum likelihood, and \pkg{statmodels7} takes the result as the intercept
#' of each equation.
#'
#' @details
#' # How the inversions are written and checked
#' Each is written against the family's own [mean()] and [variance()], not from
#' memory, and the tests check it the same way: a family's moment estimate
#' applied to a large sample drawn from a known parameter must return that
#' parameter.
#'
#' # Where a moment does not exist, and where the inversion is not closed
#' A family with no moments takes the robust analogue. The Cauchy uses the
#' median and half the interquartile range, which are its location and its
#' scale. Where the inversion needs a root, the standard approximation
#' stands in: the Weibull's shape from the coefficient of variation, the von
#' Mises's concentration from the mean resultant length. A starting value is
#' allowed to be approximate.
#'
#' Five families have no entry and fall through to
#' [start_from_moments()]'s reading of `params_interpretation`: `gengamma1`,
#' `gengamma2`, `skewt`, `pseudohuber` and `enet`. Their inversions are neither
#' closed nor standard: two gamma-function ratios in two shapes for the
#' generalized gamma, a moment that does not identify the pair for the skew
#' \eqn{t}, no elementary moments at all for the elastic net.
#'
#' # A fixed constant in the name
#' A family carrying one announces it there, as in `"beta-binomial [size=10]"`,
#' so the bracketed part is stripped before the name is used as a lookup key.
#' Without that no beta-binomial of any size would match its own entry, and
#' silently: a missing key is a legal fallback.
#'
#' @param distrib A univariate distribution, read for its `distrib_name` and
#'   its parameter names.
#' @param y The response, a numeric vector, already filtered to finite values
#'   by the caller.
#'
#' @return A named list of parameters on the parameter scale, one component per
#'   entry of `distrib@params`, or `NULL` where this family has no entry.
#'
#' @seealso [start_from_moments()], the only caller and the fallback;
#'   [distrib_start()] for the generic; [mean()] and [variance()], which the
#'   inversions are written against.
#' @keywords internal
moment_estimates <- function(distrib, y) {
  # A family carrying a fixed constant announces it in its name, as in
  # "beta-binomial [size=10]", so the bracketed part is dropped before the
  # name is used as a key: the moment estimate depends on the size through
  # distrib@size and not through the spelling.
  nm <- gsub("[^a-z0-9]", "", tolower(sub("\\[.*$", "", distrib@distrib_name)))
  yy <- as.numeric(y)
  yy <- yy[is.finite(yy)]
  if (length(yy) < 2L) return(NULL)
  m <- mean(yy)
  v <- stats::var(yy)
  if (!is.finite(v) || v <= 0) v <- max(1e-8, abs(m))
  s <- sqrt(v)
  # kurtosis identifies a t's degrees of freedom, excess = 6/(nu - 4)
  z <- (yy - m) / s
  kurt <- mean(z^4) - 3
  nu_k <- if (is.finite(kurt) && kurt > 0.05) min(6 / kurt + 4, 100) else 30
  # skewness identifies the shape of the skew normal in either chart; the
  # family cannot reach |gamma1| beyond 0.9953, so the sample value is held
  # inside that
  skew <- mean(z^3)
  skew <- if (is.finite(skew)) min(max(skew, -0.9), 0.9) else 0
  pos <- yy[yy > 0]
  ok_pos <- length(pos) > 1L

  out <- switch(
    nm,
    gaussian1 = list(mu = m, sigma = s),
    gaussian2 = list(mu = m, sigma2 = v),
    gaussian3 = list(mu = m, tau = 1 / v),
    # no moments exist; the location and scale are what these estimate
    cauchy = list(mu = stats::median(yy),
                  sigma = max(stats::IQR(yy) / 2, 1e-8)),
    laplace = list(mu = stats::median(yy), sigma = sqrt(v / 2)),
    laplace2 = list(mu = stats::median(yy), lambda = sqrt(2 / v)),
    logistic = list(mu = m, sigma = sqrt(3 * v) / pi),
    # var = sigma^2 nu / (nu - 2), the scale being smaller than the sd
    studentt1 = list(mu = m, sigma = s * sqrt(max(nu_k - 2, 1e-3) / nu_k),
                     nu = nu_k),
    # parametrized BY the standard deviation, so sigma is s itself
    studentt2 = list(mu = m, sigma = s, nu = nu_k),
    # var = phi mu^2
    gamma1 = list(mu = m, phi = v / m^2),
    gamma2 = list(mu = m, sigma2 = v),
    exponential = list(mu = m),
    chisq = list(mu = m),
    poisson = list(mu = m),
    geometric = list(mu = m),
    # var = mu (1 + theta)
    negbin1 = list(mu = m, theta = max(v / m - 1, 1e-3)),
    # var = mu + mu^2 / theta
    negbin2 = list(mu = m, theta = if (v > m) m^2 / (v - m) else 1e3),
    # var = mu (1 - mu) / (1 + phi)
    beta1 = list(mu = min(max(m, 1e-6), 1 - 1e-6),
                 phi = max(m * (1 - m) / v - 1, 1e-3)),
    beta2 = {
      k <- max(m * (1 - m) / v - 1, 1e-3)
      list(alpha = max(m * k, 1e-3), beta = max((1 - m) * k, 1e-3))
    },
    bernoulli = list(mu = min(max(m, 1e-6), 1 - 1e-6)),
    binomial = list(mu = min(max(m, 1e-6), 1 - 1e-6)),
    # shape from the coefficient of variation, the standard approximation
    weibull1 = {
      sh <- max((s / m)^(-1.086), 1e-2)
      list(mu = m / gamma(1 + 1 / sh), sigma = sh)
    },
    # the same shape, and the mean is the mean
    weibull3 = list(mean = m, sigma = max((s / m)^(-1.086), 1e-2)),
    # The centred parametrization IS the first three moments.
    skewnormal2 = list(mu = m, sigma = s, gamma1 = skew),
    # and the direct one carries them across the map of reparam_maps.R:
    # with r the real cube root of 2 gamma1/(4 - pi) and b^2 = 2/pi,
    # xi = m - s r, omega = s sqrt(1 + r^2), alpha = r/sqrt(b^2 + (b^2-1) r^2)
    skewnormal1 = {
      r <- sign(skew) * abs(2 * skew / (4 - pi))^(1 / 3)
      b2 <- 2 / pi
      list(mu = m - s * r, sigma = s * sqrt(1 + r^2),
           alpha = r / sqrt(b2 + (b2 - 1) * r^2))
    },
    # var = n mu (1-mu) [1 + (n-1) rho] with rho = sigma/(1+sigma), so the
    # intra-class correlation is read off the variance and inverted
    betabinomial = ,
    betabinom2 = {
      nsz <- distrib@size
      mu <- min(max(m / nsz, 1e-6), 1 - 1e-6)
      rho <- if (nsz > 1) (v / (nsz * mu * (1 - mu)) - 1) / (nsz - 1) else 0
      rho <- min(max(rho, 1e-8), 1 - 1e-8)
      sg <- rho / (1 - rho)
      if (identical(nm, "betabinomial")) list(mu = mu, sigma = sg)
      else list(alpha = mu / sg, beta = (1 - mu) / sg)
    },
    lognormal1 = if (ok_pos)
      list(mu = mean(log(pos)), sigma2 = max(stats::var(log(pos)), 1e-8))
      else NULL,
    lognormal2 = list(mean = m, var = v),
    # var = phi mu^3
    invgauss1 = list(mu = m, phi = v / m^3),
    # var = mu^3 / lambda
    invgauss2 = list(mu = m, lambda = m^3 / v),
    # var = pi^2 sigma^2 / 6, mean = mu + gamma sigma
    gumbel = {
      sg <- sqrt(6 * v) / pi
      list(mu = m - 0.5772156649015329 * sg, sigma = sg)
    },
    # mean = sigma/(1 - xi), var = sigma^2 / ((1-xi)^2 (1-2xi))
    generalizedpareto = {
      xi <- (1 - m^2 / v) / 2
      xi <- min(max(xi, -0.45), 0.45)
      list(sigma = max(m * (1 - xi), 1e-8), xi = xi)
    },
    # var = mu (1 + sigma mu)
    poissoninversegaussian = list(mu = m,
                                  sigma = max((v / m - 1) / m, 1e-6)),
    # the same estimate through the orthogonal chart, inverting
    # sigma = (mu + sqrt(mu^2 + alpha^2))/alpha^2 for alpha
    poissoninversegaussianorthogonal = {
      sg <- max((v / m - 1) / m, 1e-6)
      list(mu = m, alpha = sqrt(2 * m * sg + 1) / sg)
    },
    # circular: the resultant length gives the concentration
    vonmises1 = ,
    vonmises = {
      R <- sqrt(mean(cos(yy))^2 + mean(sin(yy))^2)
      k <- if (R < 0.53) 2 * R + R^3 + 5 * R^5 / 6 else
        if (R < 0.85) -0.4 + 1.39 * R + 0.43 / (1 - R) else
          1 / (R^3 - 4 * R^2 + 3 * R)
      list(mu = atan2(mean(sin(yy)), mean(cos(yy))),
           kappa = min(max(k, 1e-3), 1e3))
    },
    vonmises2 = {
      R <- sqrt(mean(cos(yy))^2 + mean(sin(yy))^2)
      list(mu = atan2(mean(sin(yy)), mean(cos(yy))),
           rho = min(max(R, 1e-6), 1 - 1e-6))
    },
    NULL)
  if (is.null(out)) return(NULL)
  if (!all(vapply(out, function(z) length(z) == 1L && is.finite(z),
                  logical(1)))) return(NULL)
  # the family's own order, and strictly inside its bounds
  out <- out[distrib@params]
  if (anyNA(names(out))) return(NULL)
  clamp_to_bounds(out, distrib)
}

#' @title Move a Starting Value Strictly Inside the Parameter Bounds
#'
#' @description
#' Returns `theta` with every component moved strictly inside its own bounds.
#' A moment estimate can land exactly on a boundary. A sample of non-negative
#' counts gives a median of zero and a sample of proportions a maximum of one,
#' and `align_theta()` treats `params_bounds` as **open**, so such a value is
#' rejected by the first generic that sees it.
#'
#' A value at a finite non-zero bound is moved in by 16 machine epsilons
#' relative to the bound. A value at a bound of exactly zero cannot be, there
#' being no relative scale there, so it is replaced by \eqn{10^{-8}} times the
#' larger of the value's own magnitude and 1. Both are far below anything a
#' first optimizer step resolves.
#'
#' @param theta A named list of parameters on the parameter scale. Components
#'   with no entry in `params_bounds`, and bounds that are not a pair, are left
#'   alone.
#' @param distrib The distribution whose `params_bounds` apply.
#'
#' @return `theta`, each component strictly inside its own bounds.
#'
#' @seealso [start_from_moments()] and [moment_estimates()], the callers;
#'   [align_theta()], which enforces the open bounds this exists to satisfy.
#' @keywords internal
clamp_to_bounds <- function(theta, distrib) {
  bnds <- distrib@params_bounds
  eps <- .Machine$double.eps
  for (p in names(theta)) {
    b <- bnds[[p]]
    if (is.null(b) || length(b) != 2L) next
    v <- theta[[p]]
    if (is.finite(b[1L]) && v <= b[1L]) {
      v <- if (b[1L] == 0) max(abs(v), 1) * 1e-8 else b[1L] + abs(b[1L]) * 16 * eps
    }
    if (is.finite(b[2L]) && v >= b[2L]) v <- b[2L] - abs(b[2L]) * 16 * eps
    theta[[p]] <- v
  }
  theta
}
