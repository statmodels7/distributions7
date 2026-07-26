#' @include distrib.R generics.R utility_functions.R numerical_functions.R zero_inflated.R zero_adjusted.R

# ===========================================================================
# Truncation.
#
# Unlike the zero wrappers, truncation adds **no parameter**: the truncation
# points are known constants, like a binomial's `size`. What it adds instead is
# a normalising constant that depends on theta,
#
#     Z(theta) = P(lower <= Y <= upper) = F(upper; theta) - F(lower^-; theta),
#
# and therefore contributes to every derivative of the log-likelihood:
#
#     l_T(y) = l(y) - log Z.
#
# Differentiating log Z is the whole content of this file. The key step is that
# its derivatives are expectations **under the truncated distribution itself**:
#
#     d_i Z / Z     = E_T[s_i]                       =: m_i
#     d_ij Z / Z    = E_T[H_ij + s_i s_j]            =: M_ij
#
# because d_i Z = \int_T d_i f = \int_T f s_i, and likewise at second order with
# d_ij f / f = H_ij + s_i s_j. Hence
#
#     d_i    l_T = s_i(y) - m_i
#     d_ij   l_T = H_ij(y) - (M_ij - m_i m_j)
#     E[d_ij l_T] = -( E_T[s_i s_j] - m_i m_j ) = -Cov_T(s_i, s_j).
#
# The last line is the second Bartlett identity for the truncated model, which
# is a useful consistency check rather than a separate derivation: the expected
# Hessian is minus the covariance of the parent's score under truncation.
# ===========================================================================

#' @title S7 Class for Truncated Continuous Distributions
#' @name TruncatedContinuousDistrib
#'
#' @description
#' A subclass of \code{continuous_distrib} representing a continuous distribution
#' restricted to \eqn{[\ell, u]} and renormalised.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{continuous_distrib} object.
#' @param lower,upper The truncation points.
#' @seealso \code{\link{truncated}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_atoms.TruncatedContinuousDistrib]{distrib_atoms()}},
#'   \code{\link[=distrib_cdf.TruncatedContinuousDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.TruncatedContinuousDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.TruncatedContinuousDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.TruncatedContinuousDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.TruncatedContinuousDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.TruncatedContinuousDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.TruncatedContinuousDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.TruncatedContinuousDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.TruncatedContinuousDistrib]{distrib_rng()}},
#'   \code{\link[=expectation.TruncatedContinuousDistrib]{expectation()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
TruncatedContinuousDistrib <- S7::new_class("TruncatedContinuousDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib,
    lower = S7::class_numeric,
    upper = S7::class_numeric
  )
)

#' @title S7 Class for Truncated Discrete Distributions
#' @name TruncatedDiscreteDistrib
#'
#' @description
#' A subclass of \code{discrete_distrib} representing a discrete distribution
#' restricted to the lattice points in \eqn{[\ell, u]} and renormalised. Both
#' endpoints are \strong{included}.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{discrete_distrib} object.
#' @param lower,upper The truncation points, included in the support.
#' @seealso \code{\link{truncated}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.TruncatedDiscreteDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.TruncatedDiscreteDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.TruncatedDiscreteDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.TruncatedDiscreteDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.TruncatedDiscreteDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.TruncatedDiscreteDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.TruncatedDiscreteDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{discrete_distrib}}.
TruncatedDiscreteDistrib <- S7::new_class("TruncatedDiscreteDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib,
    lower = S7::class_numeric,
    upper = S7::class_numeric
  )
)

# --------------------------------------------------------------------------
# Internals
# --------------------------------------------------------------------------

is_truncated <- function(distrib) {
  S7::S7_inherits(distrib, TruncatedContinuousDistrib) ||
    S7::S7_inherits(distrib, TruncatedDiscreteDistrib)
}

# P(Y = x) under the parent: the pmf for a lattice distribution, the atom's
# probability for a mixed one, and zero for an ordinary continuous distribution.
#
# This is the one quantity that separates the two truncation classes, and getting
# it wrong for a *mixed* parent is subtle: the cdf of zero_adjusted(gamma) already
# includes the point mass at zero, so F(0) is not F(0^-) even though the parent is
# a continuous_distrib. Truncating it above, with the atom retained, then loses
# exactly that mass from the normalising constant -- and the resulting density
# integrates to something other than one while every formula still looks right.
parent_mass_at <- function(distrib, x, theta) {
  parent <- distrib@parent_distrib
  if (S7::S7_inherits(distrib, TruncatedDiscreteDistrib)) {
    return(distrib_pdf(parent, x, theta))
  }
  at <- distrib_atoms(parent, theta)
  if (length(at$y) && any(at$y == x)) sum(at$p[at$y == x]) else 0
}

# F(lower^-) and Z = F(upper) - F(lower^-), vectorized in theta.
#
# The lower endpoint is *included* in the truncated support, so any mass sitting
# exactly on it has to be added back: F(lower^-) = F(lower) - P(Y = lower).
trunc_constants <- function(distrib, theta) {
  parent <- distrib@parent_distrib
  lo <- distrib@lower
  up <- distrib@upper

  Fl <- if (is.infinite(lo)) 0 else {
    distrib_cdf(parent, lo, theta) - parent_mass_at(distrib, lo, theta)
  }
  Fu <- if (is.infinite(up)) 1 else distrib_cdf(parent, up, theta)

  Z <- Fu - Fl
  if (any(!is.finite(Z)) || any(Z <= 0)) {
    stop(sprintf(paste0(
      "The truncation interval [%s, %s] carries no probability under these ",
      "parameter values (computed mass %s). A truncated distribution is not ",
      "defined there."
    ), format(lo), format(up), format(min(Z))), call. = FALSE)
  }
  list(Fl = Fl, Z = Z)
}

# TRUE where y lies in the truncated support.
trunc_inside <- function(distrib, y) {
  y >= distrib@lower & y <= distrib@upper
}

# E_T[s_i] for every parameter: the m_i of the header comment.
trunc_score_mean <- function(distrib, theta) {
  parent <- distrib@parent_distrib
  out <- lapply(distrib@params, function(p) {
    expectation(distrib, function(y, theta) distrib_gradient(parent, y, theta)[[p]], theta)
  })
  stats::setNames(out, distrib@params)
}

# E_T[s_i s_j] for every Hessian component.
trunc_score_prod_mean <- function(distrib, theta) {
  parent <- distrib@parent_distrib
  pairs <- hess_pairs(distrib@params)
  out <- lapply(pairs, function(pr) {
    expectation(distrib, function(y, theta) {
      g <- distrib_gradient(parent, y, theta)
      g[[pr[1]]] * g[[pr[2]]]
    }, theta)
  })
  stats::setNames(out, names(pairs))
}

# E_T[H_ij] for every Hessian component.
trunc_hess_mean <- function(distrib, theta) {
  parent <- distrib@parent_distrib
  nms <- hess_names(distrib@params)
  out <- lapply(nms, function(nm) {
    expectation(distrib, function(y, theta) distrib_hessian(parent, y, theta)[[nm]], theta)
  })
  stats::setNames(out, nms)
}

# --------------------------------------------------------------------------
# Shared method bodies. Truncation treats the two kinds of parent identically
# once trunc_constants() has resolved the one place they differ, so the bodies
# are written once and registered on both classes.
# --------------------------------------------------------------------------

trunc_pdf <- function(distrib, y, theta, log = FALSE) {
  Z <- trunc_constants(distrib, theta)$Z
  ld <- distrib_pdf(distrib@parent_distrib, y, theta, log = TRUE) - log(Z)
  outside <- !trunc_inside(distrib, y)
  if (any(outside)) ld[rep_len(outside, length(ld))] <- -Inf
  if (log) ld else exp(ld)
}

trunc_cdf <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  cs <- trunc_constants(distrib, theta)
  res <- (distrib_cdf(distrib@parent_distrib, q, theta) - cs$Fl) / cs$Z
  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

trunc_quantile <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  p <- pmin(pmax(p, 0), 1)

  # Expand p and theta to a common length before combining them: Fl and Z follow
  # the length of theta, p its own, and multiplying the two directly would recycle
  # silently when they disagree.
  all_params <- expand_params(c(list(.p = p), theta))
  p <- all_params$.p
  th <- all_params[distrib@params]

  cs <- trunc_constants(distrib, th)
  # Inverse transform on the parent: F_T(q) = p  <=>  F(q) = F(lower^-) + p Z.
  # The generalized inverse of a lattice cdf satisfies the same relation, so the
  # discrete case needs no separate treatment.
  distrib_quantile(distrib@parent_distrib, pmin(pmax(cs$Fl + p * cs$Z, 0), 1), th)
}

trunc_rng <- function(distrib, n, theta) {
  trunc_quantile(distrib, stats::runif(n), theta)
}

trunc_gradient <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  m <- trunc_score_mean(distrib, theta)
  g <- distrib_gradient(distrib@parent_distrib, y, theta)
  stats::setNames(lapply(distrib@params, function(p) g[[p]] - m[[p]]), distrib@params)
}

trunc_hessian <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  n <- length(y)
  m <- trunc_score_mean(distrib, theta)
  EH <- trunc_hess_mean(distrib, theta)
  ES <- trunc_score_prod_mean(distrib, theta)
  h <- distrib_hessian(distrib@parent_distrib, y, theta)
  pairs <- hess_pairs(distrib@params)

  res <- lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    # M_ij = E_T[H_ij + s_i s_j];  d_ij log Z = M_ij - m_i m_j
    M <- EH[[nm]] + ES[[nm]]
    h[[nm]] - M + m[[pr[1]]] * m[[pr[2]]]
  })
  expand_params(stats::setNames(res, names(pairs))[hess_names(distrib@params)], n)
}

trunc_expected_hessian <- function(distrib, y, theta, scale = c("parameter", "link"),
                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                   nsim = 10000, ...) {
  n <- length(y)
  m <- trunc_score_mean(distrib, theta)
  ES <- trunc_score_prod_mean(distrib, theta)
  pairs <- hess_pairs(distrib@params)

  # E[H_T] = -Cov_T(s_i, s_j): the parent's expected Hessian cancels exactly
  # against the term it contributes to the second derivative of log Z.
  res <- lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    -(ES[[nm]] - m[[pr[1]]] * m[[pr[2]]])
  })
  expand_params(stats::setNames(res, names(pairs))[hess_names(distrib@params)], n)
}

# --------------------------------------------------------------------------
# Registration
# --------------------------------------------------------------------------

#' @title Truncated Probability Density Function
#' @name distrib_pdf.TruncatedContinuousDistrib
#' @description
#' \deqn{f_T(y) = \dfrac{f(y;\theta)}{Z(\theta)}\ \ \ (\ell \le y \le u), \qquad 0 \text{ otherwise}}
#' with \eqn{Z(\theta) = F(u;\theta) - F(\ell;\theta)}.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @seealso \code{\link{truncated}}
S7::method(distrib_pdf, TruncatedContinuousDistrib) <- trunc_pdf

#' @title Truncated Probability Mass Function
#' @name distrib_pdf.TruncatedDiscreteDistrib
#' @description
#' \deqn{P_T(Y = y) = \dfrac{f(y;\theta)}{Z(\theta)}\ \ \ (\ell \le y \le u), \qquad 0 \text{ otherwise}}
#' with \eqn{Z(\theta) = F(u;\theta) - F(\ell;\theta) + f(\ell;\theta)}, the mass
#' at the lower endpoint being added back because the endpoint is included.
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @seealso \code{\link{truncated}}
S7::method(distrib_pdf, TruncatedDiscreteDistrib) <- trunc_pdf

#' @title Truncated Cumulative Distribution Function (Continuous)
#' @name distrib_cdf.TruncatedContinuousDistrib
#' @description \deqn{F_T(q) = \dfrac{F(q;\theta) - F(\ell^-;\theta)}{Z(\theta)}} clamped to \eqn{[0,1]}.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @seealso \code{\link{truncated}}
S7::method(distrib_cdf, TruncatedContinuousDistrib) <- trunc_cdf

#' @title Truncated Cumulative Distribution Function (Discrete)
#' @name distrib_cdf.TruncatedDiscreteDistrib
#' @description \deqn{F_T(q) = \dfrac{F(q;\theta) - F(\ell^-;\theta)}{Z(\theta)}} clamped to \eqn{[0,1]}.
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @seealso \code{\link{truncated}}
S7::method(distrib_cdf, TruncatedDiscreteDistrib) <- trunc_cdf

#' @title Truncated Quantile Function (Continuous)
#' @name distrib_quantile.TruncatedContinuousDistrib
#' @description \deqn{Q_T(p) = Q\!\left(F(\ell^-;\theta) + p\,Z(\theta)\right)}
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @seealso \code{\link{truncated}}
S7::method(distrib_quantile, TruncatedContinuousDistrib) <- trunc_quantile

#' @title Truncated Quantile Function (Discrete)
#' @name distrib_quantile.TruncatedDiscreteDistrib
#' @description \deqn{Q_T(p) = Q\!\left(F(\ell^-;\theta) + p\,Z(\theta)\right)}
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @seealso \code{\link{truncated}}
S7::method(distrib_quantile, TruncatedDiscreteDistrib) <- trunc_quantile

#' @title Truncated Random Number Generator (Continuous)
#' @name distrib_rng.TruncatedContinuousDistrib
#' @description
#' Inverse transform sampling on the parent: \eqn{Y = Q(F(\ell^-) + U Z)} with
#' \eqn{U \sim \mathrm{Unif}(0,1)}. Exact, and unlike rejection sampling it always
#' terminates in one pass however small \eqn{Z} is.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A named list of the parent's parameters.
#' @seealso \code{\link{truncated}}
S7::method(distrib_rng, TruncatedContinuousDistrib) <- trunc_rng

#' @title Truncated Random Number Generator (Discrete)
#' @name distrib_rng.TruncatedDiscreteDistrib
#' @description Inverse transform sampling on the parent, exact for a lattice cdf.
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A named list of the parent's parameters.
#' @seealso \code{\link{truncated}}
S7::method(distrib_rng, TruncatedDiscreteDistrib) <- trunc_rng

#' @title Truncated Analytical Gradient (Continuous)
#' @name distrib_gradient.TruncatedContinuousDistrib
#' @description
#' \deqn{\dfrac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad m_i = \mathbb{E}_T[s_i]}
#' the parent's score recentred at its mean over the truncated support.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{truncated}}
S7::method(distrib_gradient, TruncatedContinuousDistrib) <- trunc_gradient

#' @title Truncated Analytical Gradient (Discrete)
#' @name distrib_gradient.TruncatedDiscreteDistrib
#' @description
#' \deqn{\dfrac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad m_i = \mathbb{E}_T[s_i]}
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{truncated}}
S7::method(distrib_gradient, TruncatedDiscreteDistrib) <- trunc_gradient

#' @title Truncated Analytical Observed Hessian (Continuous)
#' @name distrib_hessian.TruncatedContinuousDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j} = H_{ij}(y) - M_{ij} + m_i m_j}
#' with \eqn{M_{ij} = \mathbb{E}_T[H_{ij} + s_i s_j]}.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{truncated}}
S7::method(distrib_hessian, TruncatedContinuousDistrib) <- trunc_hessian

#' @title Truncated Analytical Observed Hessian (Discrete)
#' @name distrib_hessian.TruncatedDiscreteDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j} = H_{ij}(y) - M_{ij} + m_i m_j}
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{truncated}}
S7::method(distrib_hessian, TruncatedDiscreteDistrib) <- trunc_hessian

#' @title Truncated Analytical Expected Hessian (Continuous)
#' @name distrib_expected_hessian.TruncatedContinuousDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j}\right]
#'   = -\operatorname{Cov}_T(s_i, s_j)}
#' the covariance of the parent's score under the truncated distribution.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{truncated}}
S7::method(distrib_expected_hessian, TruncatedContinuousDistrib) <- trunc_expected_hessian

#' @title Truncated Analytical Expected Hessian (Discrete)
#' @name distrib_expected_hessian.TruncatedDiscreteDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j}\right]
#'   = -\operatorname{Cov}_T(s_i, s_j)}
#' @param distrib A \code{TruncatedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{truncated}}
S7::method(distrib_expected_hessian, TruncatedDiscreteDistrib) <- trunc_expected_hessian

#' @title Atoms of a Truncated Continuous Distribution
#' @name distrib_atoms.TruncatedContinuousDistrib
#' @description
#' Truncation preserves the parent's atoms that survive it, rescaled by
#' \eqn{1/Z}. This matters only when the parent is itself mixed, as
#' \code{\link{zero_adjusted}()} of a continuous distribution is.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param theta A named list of the parent's parameters.
#' @return A list with components \code{y} and \code{p}.
#' @seealso \code{\link{truncated}}, \code{\link{distrib_atoms}}
S7::method(distrib_atoms, TruncatedContinuousDistrib) <- function(distrib, theta) {
  at <- distrib_atoms(distrib@parent_distrib, theta)
  if (!length(at$y)) return(at)
  keep <- trunc_inside(distrib, at$y)
  Z <- trunc_constants(distrib, theta)$Z
  list(y = at$y[keep], p = at$p[keep] / Z[1])
}

#' @title Expectation for Truncated Continuous Distributions
#' @name expectation.TruncatedContinuousDistrib
#' @description
#' The inherited continuous method integrates the density over
#' \eqn{[\ell, u]}, which is correct unless the parent carries point masses ---
#' as it does when it is a \code{\link{zero_adjusted}()} continuous distribution.
#' Those masses are added explicitly, exactly as in
#' \code{\link[=expectation.ZeroAdjustedContinuousDistrib]{the untruncated case}}.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param f A function \code{f(y, theta, ...)}.
#' @param theta A named list of the parent's parameters.
#' @param ... Additional arguments passed to \code{f}.
#' @keywords internal
S7::method(expectation, TruncatedContinuousDistrib) <- function(distrib, f, theta, ...) {
  cont <- S7::method(expectation, continuous_distrib)(distrib, f, theta, ...)
  at <- distrib_atoms(distrib, theta)
  if (!length(at$y)) return(cont)
  for (k in seq_along(at$y)) cont <- cont + at$p[k] * f(at$y[k], theta, ...)
  cont
}

#' @title Truncated Continuous Response Gradient
#' @name distrib_grad_y.TruncatedContinuousDistrib
#' @description
#' \eqn{Z} does not depend on \eqn{y}, so inside \eqn{[\ell, u]} the response
#' derivative is the parent's. Outside, the log-density is \eqn{-\infty} and no
#' derivative exists, so \code{NaN} is returned.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector.
#' @seealso \code{\link{truncated}}
S7::method(distrib_grad_y, TruncatedContinuousDistrib) <- function(distrib, y, theta) {
  trunc_y_deriv(distrib, y, theta, distrib_grad_y)
}

#' @title Truncated Continuous Response Hessian
#' @name distrib_hess_y.TruncatedContinuousDistrib
#' @description As \code{\link[=distrib_grad_y.TruncatedContinuousDistrib]{distrib_grad_y()}}, at second order.
#' @param distrib A \code{TruncatedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector.
#' @seealso \code{\link{truncated}}
S7::method(distrib_hess_y, TruncatedContinuousDistrib) <- function(distrib, y, theta) {
  trunc_y_deriv(distrib, y, theta, distrib_hess_y)
}

# Internal: a response derivative of the parent, restricted to the truncated
# support. Defined after the methods so that its comment block cannot capture
# their documentation.
trunc_y_deriv <- function(distrib, y, theta, fun) {
  out <- rep(NaN, length(y))
  ok <- trunc_inside(distrib, y)
  if (any(ok)) {
    th_sub <- lapply(theta, function(x) if (length(x) > 1) x[ok] else x)
    out[ok] <- fun(distrib@parent_distrib, y[ok], th_sub)
  }
  out
}

# --- CONSTRUCTOR WRAPPER ---

# Internal: everything the constructor refuses, and why.
check_truncation_points <- function(distrib, lower, upper, is_disc) {
  b <- distrib@bounds
  nm <- distrib@distrib_name

  for (arg in c("lower", "upper")) {
    v <- if (arg == "lower") lower else upper
    if (is.null(v)) next
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be a single non-missing number.", arg), call. = FALSE)
    }
    if (is_disc && is.finite(v) && v != round(v)) {
      stop(sprintf(paste0(
        "'%s' = %s is not a lattice point. A discrete distribution is supported on\n",
        "  the integers, so a non-integer truncation point is ambiguous: use %s or %s."
      ), arg, format(v), format(floor(v)), format(ceiling(v))), call. = FALSE)
    }
  }

  if (!is.null(lower) && !is.null(upper) && lower >= upper) {
    stop(sprintf("'lower' (%s) must be strictly less than 'upper' (%s).",
                 format(lower), format(upper)), call. = FALSE)
  }

  # A truncation point outside the support removes nothing. This is the
  # truncated(gamma_distrib(), lower = -2) case: the result would be the parent,
  # and a user who wrote that meant something else.
  if (!is.null(lower) && lower <= b[1]) {
    stop(sprintf(paste0(
      "truncated() was given lower = %s, which is at or below the lower bound of the\n",
      "  support of '%s' (%s). Truncating there removes no probability mass and the\n",
      "  result would be the parent distribution. Choose a point strictly inside the\n",
      "  support, or omit 'lower'."
    ), format(lower), nm, format(b[1])), call. = FALSE)
  }
  if (!is.null(upper) && upper >= b[2]) {
    stop(sprintf(paste0(
      "truncated() was given upper = %s, which is at or above the upper bound of the\n",
      "  support of '%s' (%s). Truncating there removes no probability mass and the\n",
      "  result would be the parent distribution. Choose a point strictly inside the\n",
      "  support, or omit 'upper'."
    ), format(upper), nm, format(b[2])), call. = FALSE)
  }
  if (!is.null(lower) && lower >= b[2]) {
    stop(sprintf("lower = %s lies at or above the support of '%s', which ends at %s: the truncated support would be empty.",
                 format(lower), nm, format(b[2])), call. = FALSE)
  }
  if (!is.null(upper) && upper <= b[1]) {
    stop(sprintf("upper = %s lies at or below the support of '%s', which starts at %s: the truncated support would be empty.",
                 format(upper), nm, format(b[1])), call. = FALSE)
  }
  invisible(NULL)
}

#' Truncated Distribution Object
#'
#' @description
#' Restricts an existing distribution to \eqn{[\ell, u]} and renormalises it, so
#' that all the probability mass the parent placed outside the interval is
#' redistributed inside it. Either endpoint may be omitted, giving one-sided
#' truncation; at least one must be given.
#'
#' Works on discrete and continuous parents alike. **Both endpoints are
#' included**, which for a discrete parent is the difference between
#' \code{truncated(poisson_distrib(), lower = 1)} --- the zero-truncated Poisson,
#' supported on \eqn{\{1, 2, \dots\}} --- and truncating above 1.
#'
#' @param distrib An object inheriting from \code{discrete_distrib} or
#'   \code{continuous_distrib}.
#' @param lower,upper The truncation points. Each may be \code{NULL} (no
#'   truncation on that side), and at least one must be supplied. For a discrete
#'   parent both must be whole numbers.
#'
#' @details
#' Write \eqn{Z(\theta) = P(\ell \le Y \le u)} for the retained mass. Then
#' \deqn{f_T(y;\theta) = \frac{f(y;\theta)}{Z(\theta)}\ \ (\ell \le y \le u),
#' \qquad F_T(q) = \frac{F(q;\theta) - F(\ell^-;\theta)}{Z(\theta)},
#' \qquad Q_T(p) = Q\!\left(F(\ell^-;\theta) + pZ(\theta)\right),}
#' with \eqn{F(\ell^-) = F(\ell)} for a continuous parent and
#' \eqn{F(\ell) - f(\ell)} for a discrete one, since the lower endpoint is kept.
#'
#' \strong{Truncation adds no parameter.} The endpoints are known constants, like
#' a binomial's \code{size}, so the truncated distribution has exactly the
#' parent's parameters, domains and links. What it does add is a
#' \eqn{\theta}-dependent normalising constant, and that constant contributes to
#' every derivative of the log-likelihood \eqn{\ell_T = \ell - \log Z}. Writing
#' \deqn{m_i = \mathbb{E}_T[s_i], \qquad M_{ij} = \mathbb{E}_T[H_{ij} + s_is_j],}
#' where the expectations are taken under the \emph{truncated} distribution,
#' \deqn{\frac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad
#' \frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j} = H_{ij}(y) - M_{ij} + m_im_j,}
#' \deqn{\mathbb{E}\left[\frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j}\right]
#' = -\operatorname{Cov}_T(s_i, s_j).}
#' The score is simply the parent's score \emph{recentred} at its truncated mean,
#' and the information is the covariance of that score --- which is the second
#' Bartlett identity for the truncated model, and is used as a consistency check
#' rather than derived separately.
#'
#' \strong{What this costs.} \eqn{m_i} and \eqn{M_{ij}} have no closed form for a
#' general parent, and are obtained by quadrature (continuous) or summation
#' (discrete) through \code{\link{expectation}}. Derivatives of a truncated
#' distribution are therefore substantially more expensive than the parent's, and
#' third and fourth derivatives fall back to finite differences of the analytical
#' Hessian.
#'
#' \strong{What the constructor refuses.}
#' \itemize{
#'   \item Both endpoints \code{NULL}: nothing to do, and silently returning the
#'     parent would hide the mistake.
#'   \item \code{lower >= upper}.
#'   \item A truncation point that removes no mass, such as
#'     \code{truncated(gamma_distrib(), lower = -2)}: the Gamma is supported on
#'     \eqn{(0,\infty)}, so the result would be the Gamma itself.
#'   \item A non-integer endpoint for a discrete parent, which is ambiguous.
#'   \item A discrete truncation leaving too few support points to identify the
#'     parameters: \eqn{k} points carry \eqn{k-1} free probabilities, so
#'     \code{n_params + 1} points are needed.
#'   \item A parent that models a probability of zero ---
#'     \code{\link{zero_inflated}()} or \code{\link{zero_adjusted}()} --- when the
#'     truncation removes \eqn{0} from the support. Truncating zero away cancels
#'     that parameter out of the likelihood entirely, leaving an identically zero
#'     score. Truncating elsewhere, as in
#'     \code{truncated(zero_adjusted(gamma_distrib()), upper = 5)}, is fine and the
#'     point mass is carried through \code{\link{distrib_atoms}}.
#' }
#' Truncating an already truncated distribution is allowed and is collapsed into a
#' single object over the intersection of the two intervals, rather than nested.
#'
#' @return An S7 object of class \code{TruncatedDiscreteDistrib} or
#'   \code{TruncatedContinuousDistrib}.
#'
#' @examples
#' # The zero-truncated Poisson
#' ztp <- truncated(poisson_distrib(), lower = 1)
#' distrib_pdf(ztp, 0:4, list(mu = 2))
#' dpois(1:4, 2) / (1 - dpois(0, 2))
#'
#' # A Gaussian restricted to an interval
#' tn <- truncated(gaussian_distrib(), lower = -1, upper = 2)
#' mean(tn, list(mu = 0, sigma = 1))
#'
#' # A truncation point that removes nothing is refused
#' try(truncated(gamma_distrib(), lower = -2))
#'
#' @seealso \code{\link{zero_inflated}}, \code{\link{zero_adjusted}},
#'   \code{\link{transformation}}, \code{\link{check_distrib}}
#' @export
truncated <- function(distrib, lower = NULL, upper = NULL) {
  is_disc <- S7::S7_inherits(distrib, discrete_distrib)
  if (!is_disc && !S7::S7_inherits(distrib, continuous_distrib)) {
    stop("Input must inherit from 'discrete_distrib' or 'continuous_distrib'.", call. = FALSE)
  }
  if (is.null(lower) && is.null(upper)) {
    stop(paste0(
      "truncated() needs at least one of 'lower' and 'upper'. With neither, the\n",
      "  result would be the parent distribution unchanged; returning it silently\n",
      "  would hide a missing argument rather than report it."
    ), call. = FALSE)
  }

  # Collapse a nested truncation into a single object over the intersection.
  # Nesting is harmless mathematically but pays the quadrature cost twice.
  if (is_truncated(distrib)) {
    lower <- if (is.null(lower)) distrib@lower else max(lower, distrib@lower)
    upper <- if (is.null(upper)) distrib@upper else min(upper, distrib@upper)
    distrib <- distrib@parent_distrib
    is_disc <- S7::S7_inherits(distrib, discrete_distrib)
    if (is.infinite(lower) && lower < 0) lower <- NULL
    if (is.infinite(upper) && upper > 0) upper <- NULL
  }

  check_truncation_points(distrib, lower, upper, is_disc)

  b <- distrib@bounds
  lo <- if (is.null(lower)) b[1] else lower
  up <- if (is.null(upper)) b[2] else upper

  # A zero wrapper's extra parameter lives entirely in the mass at zero. Truncate
  # zero away and it leaves the likelihood: the score is identically zero and the
  # parameter wanders on a flat ridge. Same defect as stacking two zero wrappers.
  if (is_zero_wrapper(distrib) && !(lo <= 0 && up >= 0)) {
    stop(sprintf(paste0(
      "truncated() cannot remove 0 from the support of '%s', which models the\n",
      "  probability of a zero. The truncation constant cancels that parameter out of\n",
      "  the likelihood, leaving an identically zero score. Truncate elsewhere, or\n",
      "  truncate the underlying distribution before wrapping it."
    ), distrib@distrib_name), call. = FALSE)
  }

  if (is_disc) {
    k <- if (all(is.finite(c(lo, up)))) up - lo + 1 else Inf
    if (k < distrib@n_params + 1) {
      stop(sprintf(paste0(
        "Truncating '%s' to [%s, %s] leaves %g support point(s), carrying %g free\n",
        "  probabilities, while the distribution has %d parameter(s). They would not be\n",
        "  identified. At least %d support points are required."
      ), distrib@distrib_name, format(lo), format(up), k, k - 1,
      distrib@n_params, distrib@n_params + 1), call. = FALSE)
    }
  }

  side <- if (is.null(upper)) sprintf("lower=%s", format(lo)) else
    if (is.null(lower)) sprintf("upper=%s", format(up)) else
      sprintf("%s, %s", format(lo), format(up))

  common <- list(
    parent_distrib = distrib,
    lower = lo,
    upper = up,
    distrib_name = sprintf("truncated %s [%s]", distrib@distrib_name, side),
    dimension = distrib@dimension,
    bounds = c(lo, up),
    params = distrib@params,
    params_interpretation = distrib@params_interpretation,
    n_params = distrib@n_params,
    params_bounds = distrib@params_bounds,
    link_params = distrib@link_params,
    params_smooth = param_smoothness(distrib)
  )

  if (is_disc) {
    do.call(TruncatedDiscreteDistrib, common)
  } else {
    do.call(TruncatedContinuousDistrib, common)
  }
}
