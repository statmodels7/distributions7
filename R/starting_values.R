#' @include distrib.R generics.R multivariate.R mvgaussian_distrib.R mvstudent_t_distrib.R pig1_distrib.R pig2_distrib.R
NULL

#' A Starting Value Drawn from the Data
#'
#' @description
#' Returns a starting value for [fit_distrib()], computed from the
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
#' that registers nothing loses nothing. A family with a better estimator
#' registers a method: an exact maximum likelihood estimator where one is
#' known, a method-of-moments estimator otherwise, or the estimate of a simpler
#' family the harder one contains.
#'
#' A method may return several starting values, as a list, and
#' `fit_distrib()` will try each; the first is the one it prefers.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y The response.
#' @param n_start How many starting values are wanted. A method free to supply
#'   only one may ignore it.
#' @param ... Passed to methods.
#'
#' @return A list of named parameter lists, on the **parameter** scale.
#'
#' @seealso [fit_distrib()], [generate_random_theta()]
#'
#' @examples
#' set.seed(1)
#' d <- mvgaussian_distrib(2)
#' y <- distrib_rng(d, 200, list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0,
#'                               sigma_log_L2 = 0, sigma_L2.1 = 0.5))
#'
#' # the gaussian has a closed-form maximum likelihood estimate, so the fit
#' # starts there and has nothing left to do
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
#' The default: `n_start` draws from [generate_random_theta()],
#' which uses the parameter domains and not the data. A family with a better
#' idea registers its own method.
#' @param distrib A [distrib()] object.
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
#' @return A list with `mu` and `sigma`.
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
#' [parameters7::param_free()] is exact or rejected: a parameter that
#' cannot represent the matrix signals an error rather than returning
#' something plausible. That is the right contract for reporting an estimate and the
#' wrong one for choosing where to begin, so a rejection here falls back to a
#' short numerical search over the free values, which is allowed to be
#' approximate because a starting value is.
#'
#' @param s A \pkg{parameters7} structure.
#' @param m The matrix to represent.
#'
#' @return A numeric vector of length `s@n_free`.
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
#' @param distrib A [MvGaussianDistrib()] object.
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
#' @param distrib A [MvStudentTDistrib()] object.
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


#' @title Poisson-Inverse Gaussian Starting Values
#' @name distrib_start.Pig1Distrib
#' @description Method of moments: the sample mean for \eqn{\mu}, and
#' \eqn{(s^2 - \bar y)/\bar y^2} for \eqn{\sigma}, floored just above zero
#' when the sample is underdispersed.
#' @param distrib A `Pig1Distrib` object.
#' @param y A numeric vector of observations.
#' @param n_start Ignored; one moment start is returned.
#' @param ... Unused.
#' @return A list with one named parameter list.
#' @seealso [pig1_distrib()]
S7::method(distrib_start, Pig1Distrib) <- function(distrib, y, n_start = 5L, ...) {
  mu <- max(mean(y), 1e-8)
  list(list(mu = mu, sigma = max((stats::var(y) - mu) / mu^2, 1e-3)))
}

#' @title Orthogonal Poisson-Inverse Gaussian Starting Values
#' @name distrib_start.Pig2Distrib
#' @description The moment start of
#' [`pig1()`][distrib_start.Pig1Distrib], mapped onto
#' \eqn{\alpha = \sqrt{1 + 2\sigma\mu}/\sigma}.
#' @param distrib A `Pig2Distrib` object.
#' @param y A numeric vector of observations.
#' @param n_start Ignored; one moment start is returned.
#' @param ... Unused.
#' @return A list with one named parameter list.
#' @seealso [pig2_distrib()]
S7::method(distrib_start, Pig2Distrib) <- function(distrib, y, n_start = 5L, ...) {
  mu <- max(mean(y), 1e-8)
  sg <- max((stats::var(y) - mu) / mu^2, 1e-3)
  list(list(mu = mu, alpha = sqrt(1 + 2 * sg * mu) / sg))
}


#' Where a Univariate Family Starts, From the Data
#'
#' @description
#' A starting value whose location and spread are the response's, for any
#' univariate family that declares what its parameters mean.
#'
#' @details
#' The base method draws each parameter from its own domain and never looks
#' at `y`, which is fine while the response is of order one and fails
#' completely when it is not: on a response of mean 919 and standard
#' deviation 169 the draws are of order one, the first Newton step is taken
#' from a point where the residuals are hundreds of standard deviations
#' wide, and the scale runs to the largest representable double. Measured
#' on a gaussian, `fit_distrib()` recovers N(5, 2) and N(50, 20) and
#' fails on N(500, 200): the defect is a threshold in the scale of the data,
#' not in the family.
#'
#' What makes a general fix possible is that every shipped family already
#' declares `params_interpretation`. A parameter that means a location
#' is started at the sample median, one that means a spread at the sample
#' standard deviation or its square, and one whose meaning is a shape, a
#' dispersion or a probability is left to the draw, those being of order one
#' whatever the data. A family declaring nothing recognizable loses nothing:
#' it keeps the draw it had.
#'
#' The result is CLAMPED strictly inside each parameter's bounds, because a
#' sample median can land exactly on the boundary of a support and the
#' validator treats bounds as open. Only the first start is data-based; the
#' rest stay random, so a caller asking for several still explores.
#'
#' @param distrib A univariate distribution.
#' @param y The response.
#' @param n_start How many starting values.
#' @param ... Unused.
#'
#' @return A list of named lists, one per start.
#'
#' @seealso [distrib_start()], [fit_distrib()]
#'
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


#' Method-of-Moments Estimates, Family by Family
#'
#' @description
#' The parameters a family's own first two moments imply for a sample, in
#' closed form where the inversion has one.
#'
#' @details
#' A starting value should be an estimate, not a guess, and for most
#' families the moment estimate is one line: the mean and the variance of
#' the sample are set equal to the family's own and the pair is inverted.
#' [fit_distrib()] then refines it by maximum likelihood, and
#' \pkg{statmodels7} takes the result as the intercept of each equation.
#'
#' The inversions are written against each family's `mean()` and
#' `variance()` rather than from memory, and the tests check them that
#' way: a family's moment estimate applied to a large sample from a known
#' parameter must return that parameter.
#'
#' Where the family has no moments the robust analogue is used instead --
#' the Cauchy takes the median and half the interquartile range, which is
#' what its location and scale are. Where the inversion needs a root the
#' standard approximation is used, the Weibull's shape from the coefficient
#' of variation and the von Mises's concentration from the mean resultant
#' length; a starting value is allowed to be approximate. Families whose
#' inversion is neither closed nor standard -- the generalized gamma, the
#' skew normal in its direct parametrization, the skew t, the elastic net --
#' are not listed and fall back to
#' [start_from_moments()]'s reading of
#' `params_interpretation`.
#'
#' @param distrib A univariate distribution.
#' @param y The response.
#'
#' @return A named list of parameters, or `NULL` where this family has
#'   no entry.
#'
#' @seealso [distrib_start()], [start_from_moments()]
#'
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

#' Move a Starting Value Strictly Inside the Parameter Bounds
#'
#' @param theta A named list.
#' @param distrib The distribution whose bounds apply.
#'
#' @return `theta`, each element strictly inside its own bounds.
#'
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
