#' @include distrib.R generics.R multivariate.R moments.R
NULL

#' @title S7 Class for the Dirichlet Distribution
#' @name DirichletDistrib
#'
#' @description A subclass of \code{multivariate_distrib} representing the
#'   Dirichlet distribution on the simplex, written in a mean vector and a
#'   concentration.
#' @inheritParams multivariate_distrib
#' @param param The \pkg{parameters7} simplex carrying the mean.
#' @return An object of class \code{DirichletDistrib}.
#' @seealso \code{\link{dirichlet_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_expected_hessian.DirichletDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.DirichletDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.DirichletDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.DirichletDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_rng.DirichletDistrib]{distrib_rng()}},
#'   \code{\link[=mv_location.DirichletDistrib]{mv_location()}},
#'   \code{\link[=mv_marginal.DirichletDistrib]{mv_marginal()}},
#'   \code{\link[=mv_sigma.DirichletDistrib]{mv_sigma()}}
#'
#' Everything else is inherited from \code{\link{multivariate_distrib}}, which
#' rejects the distribution function and the quantile.
DirichletDistrib <- S7::new_class("DirichletDistrib",
  parent = multivariate_distrib,
  properties = list(param = S7::class_any)
)

#' The Pieces a Dirichlet Derivative Needs
#'
#' @description
#' The mean vector, the concentration, the shapes \eqn{\alpha = \phi\mu} and
#' the simplex's first two derivatives, computed once per call.
#'
#' @details
#' Two identities keep every formula short, both
#' following from \eqn{\sum_j \mu_j = 1} differentiated once and twice: the
#' columns of \eqn{A = \partial\mu/\partial\eta} sum to zero, and so does
#' every second-derivative vector. They are what make the expected information
#' closed form, since the terms carrying the data drop out under expectation.
#'
#' @param distrib A \code{\link{DirichletDistrib}} object.
#' @param theta A named list of parameters.
#'
#' @return A list with \code{mu}, \code{phi}, \code{alpha}, \code{A} (a
#'   \eqn{p \times (p-1)} matrix), \code{B} (the second derivatives, keyed by
#'   tuple) and \code{idx} (their index tuples).
#'
#' @seealso \code{\link{dirichlet_distrib}}
#'
#' @keywords internal
dir_parts <- function(distrib, theta) {
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  k <- s@n_free
  eta <- v[seq_len(k)]
  phi <- v[[k + 1L]]
  mu <- as.numeric(parameters7::param_value(s, eta))
  list(mu = mu, phi = phi, alpha = phi * mu,
       A = do.call(cbind, parameters7::param_d1(s, eta)),
       B = parameters7::param_d2(s, eta),
       idx = parameters7::param_tuple_indices(s, 2L),
       k = k)
}

#' The Position of a Second-Derivative Block
#'
#' @description
#' Locates the \eqn{(k, l)} entry in the list of second derivatives a
#' \pkg{parameters7} parameter returns, whose elements are keyed by unordered
#' tuple rather than by position.
#'
#' @param idx The tuple index list, as returned by
#'   \code{parameters7::param_tuple_indices()}.
#' @param k,l The two free-value positions.
#'
#' @return An integer position into the second-derivative list.
#'
#' @seealso \code{\link{dirichlet_distrib}}
#'
#' @keywords internal
dir_b_index <- function(idx, k, l) {
  which(vapply(idx, function(t) {
    length(t) == 2L && identical(sort(as.integer(t)), sort(as.integer(c(k, l))))
  }, logical(1)))[1L]
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Dirichlet Density
#' @name distrib_pdf.DirichletDistrib
#' @description
#' \deqn{f(y) = \dfrac{\Gamma(\phi)}{\prod_j \Gamma(\alpha_j)}
#'       \prod_j y_j^{\alpha_j - 1}, \qquad \alpha = \phi\mu}
#' evaluated on the open simplex, one row of \code{y} per observation.
#' @param distrib A \code{DirichletDistrib} object.
#' @param y A matrix with one row per observation, each summing to one.
#' @param theta A named list of parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector, one entry per row of \code{y}.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(distrib_pdf, DirichletDistrib) <- function(distrib, y, theta, log = FALSE) {
  p <- dir_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  const <- lgamma(p$phi) - sum(lgamma(p$alpha))
  # The support is tested before the logarithm rather than after: log() of a
  # negative coordinate warns, and a density evaluated off its support is zero
  # rather than a numerical complaint.
  ok <- apply(y, 1L, function(r) all(r > 0) && abs(sum(r) - 1) <= 1e-8)
  out <- rep(-Inf, nrow(y))
  if (any(ok)) {
    out[ok] <- const + as.numeric(log(y[ok, , drop = FALSE]) %*% (p$alpha - 1))
  }
  if (log) out else exp(out)
}

#' @title Dirichlet Random Generation
#' @name distrib_rng.DirichletDistrib
#' @description
#' Independent Gamma draws with the shapes \eqn{\alpha_j}, normalized by their
#' sum, which is the representation the family is defined by.
#' @param distrib A \code{DirichletDistrib} object.
#' @param n The number of draws.
#' @param theta A named list of parameters.
#' @return A matrix with \code{n} rows, each summing to one.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(distrib_rng, DirichletDistrib) <- function(distrib, n, theta) {
  p <- dir_parts(distrib, theta)
  d <- length(p$alpha)
  g <- matrix(stats::rgamma(n * d, rep(p$alpha, each = n)), nrow = n)
  g / rowSums(g)
}

#' @title Dirichlet Analytical Gradient
#' @name distrib_gradient.DirichletDistrib
#' @description
#' With \eqn{g_j = \log y_j - \psi(\alpha_j)} and
#' \eqn{A = \partial\mu/\partial\eta},
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k} = \phi \sum_j g_j A_{jk},
#'       \qquad
#'       \dfrac{\partial\ell}{\partial\phi} = \psi(\phi) + \sum_j g_j \mu_j}
#' @param distrib A \code{DirichletDistrib} object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list, one component per parameter.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(distrib_gradient, DirichletDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  p <- dir_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  g <- sweep(t(log(y)), 1L, digamma(p$alpha), "-")
  out <- vector("list", distrib@n_params)
  for (k in seq_len(p$k)) out[[k]] <- p$phi * as.numeric(crossprod(p$A[, k], g))
  out[[p$k + 1L]] <- digamma(p$phi) + as.numeric(crossprod(p$mu, g))
  stats::setNames(out, distrib@params)
}

#' @title Dirichlet Analytical Observed Hessian
#' @name distrib_hessian.DirichletDistrib
#' @description
#' The second derivatives of the same expressions, with
#' \eqn{t_j = \psi'(\alpha_j)}:
#' \deqn{\ell^{(\eta_k\eta_l)} = \phi\sum_j\left(-t_j\phi A_{jk}A_{jl}
#'         + g_j B_{j,kl}\right), \qquad
#'       \ell^{(\phi\phi)} = \psi'(\phi) - \sum_j t_j \mu_j^2}
#' the last free of the data, the family being an exponential family in
#' \eqn{\log y}.
#' @param distrib A \code{DirichletDistrib} object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(distrib_hessian, DirichletDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  p <- dir_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  n <- nrow(y)
  g <- sweep(t(log(y)), 1L, digamma(p$alpha), "-")
  tt <- trigamma(p$alpha)
  kk <- p$k
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  ix <- match(distrib@params, distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    if (a <= kk && b <= kk) {
      bi <- dir_b_index(p$idx, a, b)
      base <- -p$phi * p$phi * sum(tt * p$A[, a] * p$A[, b])
      rep(base, n) + p$phi * as.numeric(crossprod(p$B[[bi]], g))
    } else if (a > kk && b > kk) {
      rep(trigamma(p$phi) - sum(tt * p$mu^2), n)
    } else {
      j <- if (a > kk) b else a
      base <- -p$phi * sum(tt * p$mu * p$A[, j])
      rep(base, n) + as.numeric(crossprod(p$A[, j], g))
    }
  }), nm)
}

#' @title Dirichlet Analytical Expected Hessian
#' @name distrib_expected_hessian.DirichletDistrib
#' @description
#' Closed form. \eqn{\mathbb{E}[\log y_j] = \psi(\alpha_j) - \psi(\phi)}, so
#' \eqn{\mathbb{E}[g_j] = -\psi(\phi)} for every \eqn{j}; the columns of
#' \eqn{A} and every second-derivative vector of the simplex sum to zero, so
#' each data-carrying term drops out and
#' \deqn{\mathbb{E}[\ell^{(\eta_k\eta_l)}] = -\phi^2\sum_j t_j A_{jk}A_{jl},
#'       \qquad
#'       \mathbb{E}[\ell^{(\eta_k\phi)}] = -\phi\sum_j t_j \mu_j A_{jk}}
#' @param distrib A \code{DirichletDistrib} object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(distrib_expected_hessian, DirichletDistrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  p <- dir_parts(distrib, theta)
  n <- n_obs(distrib, y)
  tt <- trigamma(p$alpha)
  kk <- p$k
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    val <- if (a <= kk && b <= kk) {
      -p$phi * p$phi * sum(tt * p$A[, a] * p$A[, b])
    } else if (a > kk && b > kk) {
      trigamma(p$phi) - sum(tt * p$mu^2)
    } else {
      j <- if (a > kk) b else a
      -p$phi * sum(tt * p$mu * p$A[, j])
    }
    rep(val, n)
  }), nm)
}

#' @title Dirichlet Mean Vector
#' @name mv_location.DirichletDistrib
#' @description The simplex the mean parameter carries, which is the mean.
#' @param distrib A \code{DirichletDistrib} object.
#' @param theta A named list of parameters.
#' @return A numeric vector summing to one.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(mv_location, DirichletDistrib) <- function(distrib, theta) {
  dir_parts(distrib, theta)$mu
}

#' @title Dirichlet Covariance Matrix
#' @name mv_sigma.DirichletDistrib
#' @description
#' \deqn{\operatorname{Cov}(Y_i, Y_j)
#'       = \dfrac{\delta_{ij}\mu_i - \mu_i\mu_j}{\phi + 1}}
#' singular by construction, the coordinates summing to a constant.
#' @param distrib A \code{DirichletDistrib} object.
#' @param theta A named list of parameters.
#' @return A symmetric matrix.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(mv_sigma, DirichletDistrib) <- function(distrib, theta) {
  p <- dir_parts(distrib, theta)
  (diag(p$mu, nrow = length(p$mu)) - tcrossprod(p$mu)) / (p$phi + 1)
}

#' @title Mean of a Dirichlet
#' @name mean.DirichletDistrib
#' @description The mean vector, which is a parameter of the family.
#' @param x A \code{DirichletDistrib} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector summing to one.
#' @seealso \code{\link{dirichlet_distrib}}
#' @keywords internal
S7::method(mean, DirichletDistrib) <- function(x, theta, ...) {
  mv_location(x, theta)
}

#' @title Variance of a Dirichlet
#' @name variance.DirichletDistrib
#' @description The covariance matrix \code{\link{mv_sigma}} carries, singular
#'   because the coordinates sum to one.
#' @param x A \code{DirichletDistrib} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A symmetric \eqn{p \times p} matrix.
#' @seealso \code{\link{dirichlet_distrib}}
#' @keywords internal
S7::method(variance, DirichletDistrib) <- function(x, theta, ...) {
  mv_sigma(x, theta)
}

#' @title A Uniform Proposal on the Simplex
#' @name mv_reference_draw.DirichletDistrib
#' @description
#' Draws from the uniform distribution on the simplex, which is the Dirichlet
#' with every shape equal to one and has the constant density \eqn{\Gamma(p)}
#' with respect to the same dominating measure the family's density is written
#' against. The base class's gaussian proposal lives in \eqn{\mathbb{R}^p} and
#' would place no mass at all on the simplex.
#' @param distrib A \code{DirichletDistrib} object.
#' @param theta A named list of parameters.
#' @param n The number of draws.
#' @param ... Unused.
#' @return A list with the draws \code{y} and their log-density \code{logd}.
#' @seealso \code{\link{dirichlet_distrib}}
#' @keywords internal
S7::method(mv_reference_draw, DirichletDistrib) <- function(distrib, theta, n, ...) {
  p <- distrib@n_dim
  g <- matrix(stats::rexp(n * p), nrow = n)
  list(y = g / rowSums(g), logd = rep(lgamma(p), n))
}

#' @title Dirichlet Marginal
#' @name mv_marginal.DirichletDistrib
#' @description
#' A single coordinate is \eqn{\mathrm{Beta}(\alpha_j, \phi-\alpha_j)}, which
#' in this package's mean-and-precision parametrization of the Beta is simply
#' mean \eqn{\mu_j} and precision \eqn{\phi}: the concentration is shared by
#' every marginal. A group of coordinates is again Dirichlet, but only after
#' the remaining mass is collapsed into one of its own, so that case is
#' rejected.
#' @param distrib A \code{DirichletDistrib} object.
#' @param theta A named list of parameters.
#' @param which The coordinate wanted.
#' @param ... Unused.
#' @return A list with the marginal \code{distrib} and its \code{theta}.
#' @seealso \code{\link{dirichlet_distrib}}
S7::method(mv_marginal, DirichletDistrib) <- function(distrib, theta, which, ...) {
  if (length(which) != 1L) {
    stop(paste0(
      "A Dirichlet marginal is returned one coordinate at a time. Several\n",
      "  coordinates are again Dirichlet, but only after the remaining mass is\n",
      "  collapsed into a coordinate of its own, which is a different object\n",
      "  from a sub-vector and would be misleading to return under this name."
    ), call. = FALSE)
  }
  p <- dir_parts(distrib, theta)
  # The Beta here is written in a mean and a precision, exactly as this family
  # is, so the marginal needs no reparametrization at all: coordinate j has
  # shapes (alpha_j, phi - alpha_j), hence mean mu_j and precision phi. The
  # CONCENTRATION IS SHARED by every marginal, as the multivariate t's degrees
  # of freedom are.
  list(distrib = beta1_distrib(),
       theta = list(mu = p$mu[which], phi = p$phi))
}

# --- CONSTRUCTOR WRAPPER ---

#' Dirichlet Distribution Object
#'
#' @description
#' Creates a distribution object for the Dirichlet distribution on the
#' simplex, parametrized by a mean vector and a concentration \eqn{\phi}.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{\Gamma(\phi)}{\prod_{j=1}^{p}\Gamma(\phi\mu_j)}\prod_{j=1}^{p} y_j^{\phi\mu_j - 1}}
#' on the simplex \eqn{y_j > 0}, \eqn{\sum_j y_j = 1}, with
#'
#' \deqn{\mathbb{E}[Y_j] = \mu_j, \qquad \operatorname{Var}(Y_j) = \frac{\mu_j(1-\mu_j)}{\phi + 1}.}
#'
#' The mean lies on a \pkg{parameters7} simplex and \eqn{\phi} is the
#' concentration, shared by every coordinate; the marginals are Beta
#' with that same \eqn{\phi}.
#'
#' @param n_dim The number of coordinates \eqn{p}.
#' @param mean A \pkg{parameters7} \code{\link[parameters7]{simplex}} of the
#'   same dimension, carrying the mean. Defaults to
#'   \code{parameters7::simplex(n_dim)}.
#' @param link_phi A link function object for \eqn{\phi}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The first family here that is multivariate and \strong{not elliptical}, and
#' therefore the second real test of that layer: there is no location and
#' scale to separate, the support is a simplex rather than a Euclidean space,
#' and the covariance is singular by construction because the coordinates sum
#' to one.
#'
#' The parametrization follows the same design as the multivariate gaussian's.
#' The constrained object --- here a point of the simplex rather than a
#' positive definite matrix --- is carried by a \pkg{parameters7} parameter and
#' \strong{flattened into scalars} with identity links, so every generic of the
#' package indexes it as it always did. The shapes are
#' \eqn{\alpha = \phi\mu}.
#'
#' \strong{Density:}
#' \deqn{f(y) = \dfrac{\Gamma(\phi)}{\prod_j\Gamma(\alpha_j)}
#'       \prod_j y_j^{\alpha_j-1}}
#'
#' \strong{Moments:} mean \eqn{\mu} and
#' \eqn{\operatorname{Cov}(Y_i,Y_j) = (\delta_{ij}\mu_i - \mu_i\mu_j)/(\phi+1)},
#' so \eqn{\phi} is a precision: the larger it is, the tighter the draws about
#' the mean.
#'
#' \strong{The expected information is closed form}, which two identities make
#' possible. Differentiating \eqn{\sum_j\mu_j = 1} once and twice shows that
#' the columns of \eqn{A = \partial\mu/\partial\eta} sum to zero and so does
#' every second-derivative vector of the simplex; and
#' \eqn{\mathbb{E}[\log y_j] = \psi(\alpha_j) - \psi(\phi)} makes
#' \eqn{\mathbb{E}[g_j] = -\psi(\phi)} the same for every \eqn{j}. Every term
#' carrying the data is therefore a constant times one of those zero sums, and
#' drops out.
#'
#' \strong{The marginals are Beta}, coordinate \eqn{j} being
#' \eqn{\mathrm{Beta}(\alpha_j, \phi-\alpha_j)}, so
#' \code{\link{mv_marginal}} returns an object rather than signaling an error --- which
#' is what makes this family a useful test of that generic rather than another
#' rejection. Several coordinates together are again Dirichlet, but only after
#' the remaining mass is collapsed into a coordinate of its own, so that case
#' is rejected rather than returned under a name that would mislead.
#'
#' The distribution function and the quantile are rejected by
#' \code{\link{multivariate_distrib}}, as for every family of that class.
#'
#' @return An S7 object of class \code{DirichletDistrib}.
#'
#' @seealso \code{\link{beta1_distrib}} for the two-coordinate case seen on the
#'   line, \code{\link{mvgaussian_distrib}},
#'   \code{\link[parameters7]{simplex}}
#'
#' @importFrom linkfunctions7 log_link identity_link
#' @importFrom stats rgamma
#' @examples
#' d <- dirichlet_distrib(3)
#' d@params
#'
#' theta <- as.list(stats::setNames(c(0.3, -0.2, log(12)), d@params))
#' mv_location(d, theta)
#' round(mv_sigma(d, theta), 5)
#'
#' # the marginals are Beta, so a panel of a pairs plot is a real object
#' mv_marginal(d, theta, which = 1)$theta
#'
#' @export
dirichlet_distrib <- function(n_dim, mean = parameters7::simplex(n_dim),
                              link_phi = log_link()) {
  p <- as.integer(n_dim)
  if (length(p) != 1L || is.na(p) || p < 2L) {
    stop("'n_dim' must be a single integer of at least 2.", call. = FALSE)
  }
  if (!S7::S7_inherits(mean, parameters7::parameter)) {
    stop("'mean' must be a parameters7 parameter.", call. = FALSE)
  }
  if (length(parameters7::param_value(mean, numeric(mean@n_free))) != p) {
    stop(sprintf(
      "The mean parameter produces %d coordinates but the distribution has %d.",
      length(parameters7::param_value(mean, numeric(mean@n_free))), p
    ), call. = FALSE)
  }

  params <- c(paste0("mean_", mean@free_names), "phi")
  n_par <- length(params)

  DirichletDistrib(
    distrib_name = sprintf("dirichlet [%dd, mean=%s]", p, mean@param_name),
    dimension = "multivariate", n_dim = p, bounds = c(0, 1),
    params = params,
    params_interpretation = stats::setNames(
      c(rep("mean", mean@n_free), "concentration"), params
    ),
    n_params = n_par,
    params_bounds = stats::setNames(
      c(rep(list(c(-Inf, Inf)), mean@n_free), list(c(0, Inf))), params
    ),
    link_params = stats::setNames(
      c(rep(list(linkfunctions7::identity_link()), mean@n_free),
        list(link_phi)), params
    ),
    param = mean
  )
}

