#' @include distrib.R generics.R multivariate.R moments.R
NULL

#' @title S7 Class for the Multinomial Distribution
#' @name MultinomialDistrib
#'
#' @description A subclass of \code{multivariate_distrib} representing the
#'   multinomial distribution, the first multivariate family here whose
#'   support is a set of points rather than a region.
#' @inheritParams multivariate_distrib
#' @param size The number of trials, a constant of the distribution.
#' @param param The \pkg{parameters7} simplex carrying the probabilities.
#' @return An object of class \code{MultinomialDistrib}.
#' @seealso \code{\link{multinomial_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_expected_hessian.MultinomialDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.MultinomialDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.MultinomialDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.MultinomialDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_rng.MultinomialDistrib]{distrib_rng()}},
#'   \code{\link[=mv_location.MultinomialDistrib]{mv_location()}},
#'   \code{\link[=mv_marginal.MultinomialDistrib]{mv_marginal()}},
#'   \code{\link[=mv_sigma.MultinomialDistrib]{mv_sigma()}},
#'   \code{\link[=mv_support.MultinomialDistrib]{mv_support()}}
#'
#' Everything else is inherited from \code{\link{multivariate_distrib}}.
MultinomialDistrib <- S7::new_class("MultinomialDistrib",
  parent = multivariate_distrib,
  properties = list(size = S7::class_numeric, param = S7::class_any)
)

#' The Pieces a Multinomial Derivative Needs
#'
#' @description
#' The probability vector and the simplex's first two derivatives, computed
#' once per call.
#'
#' @param distrib A \code{\link{MultinomialDistrib}} object.
#' @param theta A named list of parameters.
#'
#' @return A list with \code{prob}, \code{A}, \code{B}, \code{idx} and
#'   \code{k}.
#'
#' @seealso \code{\link{multinomial_distrib}}
#'
#' @keywords internal
mn_parts <- function(distrib, theta) {
  s <- distrib@param
  eta <- mv_flat_theta(distrib, theta)
  list(prob = as.numeric(parameters7::param_value(s, eta)),
       A = do.call(cbind, parameters7::param_d1(s, eta)),
       B = parameters7::param_d2(s, eta),
       idx = parameters7::param_tuple_indices(s, 2L),
       k = s@n_free)
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Multinomial Probability Mass Function
#' @name distrib_pdf.MultinomialDistrib
#' @description
#' \deqn{P(Y = y) = \dfrac{n!}{\prod_j y_j!}\prod_j p_j^{y_j}}
#' for a row of non-negative integers summing to \eqn{n}.
#' @param distrib A \code{MultinomialDistrib} object.
#' @param y A matrix with one row per observation, each summing to the size.
#' @param theta A named list of parameters.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector, one entry per row of \code{y}.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(distrib_pdf, MultinomialDistrib) <- function(distrib, y, theta, log = FALSE) {
  p <- mn_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  n <- distrib@size
  out <- lgamma(n + 1) - rowSums(lgamma(y + 1)) +
    as.numeric(y %*% log(p$prob))
  bad <- apply(y, 1L, function(r) {
    any(r < 0) || any(r != round(r)) || abs(sum(r) - n) > 1e-8
  })
  out[bad] <- -Inf
  if (log) out else exp(out)
}

#' @title Multinomial Random Generation
#' @name distrib_rng.MultinomialDistrib
#' @description \code{\link[stats]{rmultinom}}, transposed so that one row is
#'   one observation.
#' @param distrib A \code{MultinomialDistrib} object.
#' @param n The number of draws.
#' @param theta A named list of parameters.
#' @return A matrix with \code{n} rows.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(distrib_rng, MultinomialDistrib) <- function(distrib, n, theta) {
  p <- mn_parts(distrib, theta)
  t(stats::rmultinom(n, size = distrib@size, prob = p$prob))
}

#' @title The Support Points of a Multinomial
#' @name mv_support.MultinomialDistrib
#' @description
#' Every vector of non-negative integers summing to the size, enumerated as a
#' matrix with one row per point.
#' @param distrib A \code{MultinomialDistrib} object.
#' @param theta Ignored; the support does not depend on the parameters.
#' @param ... Unused.
#' @return A matrix with \code{choose(n + p - 1, p - 1)} rows.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(mv_support, MultinomialDistrib) <- function(distrib, theta, ...) {
  compositions(distrib@size, distrib@n_dim)
}

#' @title Multinomial Analytical Gradient
#' @name distrib_gradient.MultinomialDistrib
#' @description
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k}
#'       = \sum_j \dfrac{y_j}{p_j}A_{jk}, \qquad
#'       A = \dfrac{\partial p}{\partial\eta}}
#' @param distrib A \code{MultinomialDistrib} object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list, one component per parameter.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(distrib_gradient, MultinomialDistrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"), ...) {
  p <- mn_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  yp <- sweep(y, 2L, p$prob, "/")
  stats::setNames(
    lapply(seq_len(p$k), function(k) as.numeric(yp %*% p$A[, k])),
    distrib@params
  )
}

#' @title Multinomial Analytical Observed Hessian
#' @name distrib_hessian.MultinomialDistrib
#' @description
#' \deqn{\ell^{(\eta_k\eta_l)} = \sum_j\left(\dfrac{y_j}{p_j}B_{j,kl}
#'       - \dfrac{y_j}{p_j^2}A_{jk}A_{jl}\right)}
#' @param distrib A \code{MultinomialDistrib} object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(distrib_hessian, MultinomialDistrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  p <- mn_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  yp <- sweep(y, 2L, p$prob, "/")
  yp2 <- sweep(y, 2L, p$prob^2, "/")
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    bi <- which(vapply(p$idx, function(t) {
      length(t) == 2L && identical(sort(as.integer(t)),
                                   sort(as.integer(c(a, b))))
    }, logical(1)))[1L]
    as.numeric(yp %*% p$B[[bi]]) -
      as.numeric(yp2 %*% (p$A[, a] * p$A[, b]))
  }), nm)
}

#' @title Multinomial Analytical Expected Hessian
#' @name distrib_expected_hessian.MultinomialDistrib
#' @description
#' Closed form. \eqn{\mathbb{E}[y_j] = n p_j}, so the first term becomes
#' \eqn{n\sum_j B_{j,kl}}, which vanishes because the probabilities sum to a
#' constant, and
#' \deqn{\mathbb{E}[\ell^{(\eta_k\eta_l)}]
#'       = -n\sum_j \dfrac{A_{jk}A_{jl}}{p_j}}
#' @param distrib A \code{MultinomialDistrib} object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(distrib_expected_hessian, MultinomialDistrib) <- function(distrib, y, theta,
                                                                      scale = c("parameter", "link"),
                                                                      approx = c("bartlett", "integrate", "mc", "opg"),
                                                                      nsim = 10000, ...) {
  p <- mn_parts(distrib, theta)
  n <- n_obs(distrib, y)
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    rep(-distrib@size * sum(p$A[, a] * p$A[, b] / p$prob), n)
  }), nm)
}

#' @title Multinomial Mean Vector
#' @name mv_location.MultinomialDistrib
#' @description \eqn{n p}, the size times the probability vector.
#' @param distrib A \code{MultinomialDistrib} object.
#' @param theta A named list of parameters.
#' @return A numeric vector summing to the size.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(mv_location, MultinomialDistrib) <- function(distrib, theta) {
  distrib@size * mn_parts(distrib, theta)$prob
}

#' @title Multinomial Covariance Matrix
#' @name mv_sigma.MultinomialDistrib
#' @description
#' \deqn{\operatorname{Cov}(Y_i, Y_j) = n(\delta_{ij}p_i - p_i p_j)}
#' singular by construction, the coordinates summing to the size.
#' @param distrib A \code{MultinomialDistrib} object.
#' @param theta A named list of parameters.
#' @return A symmetric matrix.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(mv_sigma, MultinomialDistrib) <- function(distrib, theta) {
  p <- mn_parts(distrib, theta)$prob
  distrib@size * (diag(p, nrow = length(p)) - tcrossprod(p))
}

#' @title Mean of a Multinomial
#' @name mean.MultinomialDistrib
#' @description The vector \eqn{np}, which \code{\link{mv_location}} returns.
#' @param x A \code{MultinomialDistrib} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector summing to the size.
#' @seealso \code{\link{multinomial_distrib}}
#' @keywords internal
S7::method(mean, MultinomialDistrib) <- function(x, theta, ...) {
  mv_location(x, theta)
}

#' @title Variance of a Multinomial
#' @name variance.MultinomialDistrib
#' @description The covariance matrix \code{\link{mv_sigma}} carries, singular
#'   because the counts sum to the size.
#' @param x A \code{MultinomialDistrib} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A symmetric \eqn{p \times p} matrix.
#' @seealso \code{\link{multinomial_distrib}}
#' @keywords internal
S7::method(variance, MultinomialDistrib) <- function(x, theta, ...) {
  mv_sigma(x, theta)
}

#' @title Multinomial Marginal
#' @name mv_marginal.MultinomialDistrib
#' @description
#' A single coordinate is \eqn{\mathrm{Binomial}(n, p_j)}, the other outcomes
#' collapsing into a single failure.
#' @param distrib A \code{MultinomialDistrib} object.
#' @param theta A named list of parameters.
#' @param which The coordinate wanted.
#' @param ... Unused.
#' @return A list with the marginal \code{distrib} and its \code{theta}.
#' @seealso \code{\link{multinomial_distrib}}
S7::method(mv_marginal, MultinomialDistrib) <- function(distrib, theta, which, ...) {
  if (length(which) != 1L) {
    stop(paste0(
      "A multinomial marginal is returned one coordinate at a time. Several\n",
      "  coordinates are again multinomial, but only after the remaining\n",
      "  outcomes are collapsed into a category of their own, which is a\n",
      "  different object from a sub-vector."
    ), call. = FALSE)
  }
  list(distrib = binomial_distrib(size = distrib@size),
       theta = list(mu = mn_parts(distrib, theta)$prob[which]))
}

# --- CONSTRUCTOR WRAPPER ---

#' Multinomial Distribution Object
#'
#' @description
#' Creates a distribution object for the multinomial distribution,
#' parametrised by a probability vector on the simplex.
#'
#' @param n_dim The number of categories \eqn{p}.
#' @param size The number of trials \eqn{n}. A constant of the distribution
#'   rather than a parameter, as for \code{\link{binomial_distrib}}.
#' @param probs A \pkg{parameters7} \code{\link[parameters7]{simplex}} of the
#'   same dimension. Defaults to \code{parameters7::simplex(n_dim)}.
#'
#' @details
#' The first family here that is multivariate and \strong{discrete}, so its
#' support is a finite set of points --- every vector of non-negative integers
#' summing to \eqn{n} --- rather than a region. \code{\link{mv_support}}
#' enumerates them, which is what lets an expectation be an exact sum and the
#' validator check the total mass by addition rather than by sampling.
#'
#' The probabilities are carried by a \pkg{parameters7} simplex and flattened
#' into scalars with identity links, exactly as a covariance is for the
#' multivariate gaussian: the constraint that they be positive and sum to one
#' lives in the parameter, where a scalar link could not express it.
#'
#' \strong{Probability mass function:}
#' \deqn{P(Y=y) = \dfrac{n!}{\prod_j y_j!}\prod_j p_j^{y_j}}
#'
#' \strong{Score and information.} With \eqn{A = \partial p/\partial\eta},
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k} = \sum_j \dfrac{y_j}{p_j}A_{jk},
#'       \qquad
#'       \mathbb{E}[\ell^{(\eta_k\eta_l)}]
#'         = -n\sum_j \dfrac{A_{jk}A_{jl}}{p_j}}
#' The expected form is closed because \eqn{\mathbb{E}[y_j] = np_j} turns the
#' second-derivative term into \eqn{n\sum_j B_{j,kl}}, which vanishes: the
#' probabilities sum to one, so every derivative of their sum is zero.
#'
#' \strong{Moments:} mean \eqn{np} and
#' \eqn{\operatorname{Cov}(Y_i,Y_j) = n(\delta_{ij}p_i - p_ip_j)}, singular by
#' construction.
#'
#' \strong{The marginals are binomial}, coordinate \eqn{j} being
#' \eqn{\mathrm{Binomial}(n, p_j)} with the other categories collapsed into a
#' single failure.
#'
#' @return An S7 object of class \code{MultinomialDistrib}.
#'
#' @seealso \code{\link{binomial_distrib}}, \code{\link{dirichlet_distrib}}
#'   for the family that is conjugate to it,
#'   \code{\link[parameters7]{simplex}}
#'
#' @importFrom linkfunctions7 identity_link
#' @importFrom stats rmultinom
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' d@params
#'
#' theta <- as.list(stats::setNames(c(0.3, -0.2), d@params))
#' mv_location(d, theta)
#'
#' # the support is a finite set of points, so the mass sums exactly
#' supp <- mv_support(d, theta)
#' c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, theta)))
#'
#' @export
multinomial_distrib <- function(n_dim, size,
                                probs = parameters7::simplex(n_dim)) {
  p <- as.integer(n_dim)
  if (length(p) != 1L || is.na(p) || p < 2L) {
    stop("'n_dim' must be a single integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
      size < 1 || size != round(size)) {
    stop("'size' must be a single positive integer.", call. = FALSE)
  }
  if (!S7::S7_inherits(probs, parameters7::parameter)) {
    stop("'probs' must be a parameters7 parameter.", call. = FALSE)
  }
  if (length(parameters7::param_value(probs, numeric(probs@n_free))) != p) {
    stop("The probability parameter does not have the stated dimension.",
         call. = FALSE)
  }

  params <- paste0("probs_", probs@free_names)
  MultinomialDistrib(
    distrib_name = sprintf("multinomial [%dd, size=%g, probs=%s]", p, size,
                           probs@param_name),
    dimension = "multivariate", n_dim = p, bounds = c(0, size),
    params = params,
    params_interpretation = stats::setNames(
      rep("probability", probs@n_free), params
    ),
    n_params = probs@n_free,
    params_bounds = stats::setNames(
      rep(list(c(-Inf, Inf)), probs@n_free), params
    ),
    link_params = stats::setNames(
      rep(list(linkfunctions7::identity_link()), probs@n_free), params
    ),
    size = as.numeric(size),
    param = probs
  )
}
