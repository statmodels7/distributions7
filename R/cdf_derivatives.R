#' @include distrib.R generics.R utility_functions.R numerical_functions.R
#' @include gaussian1_distrib.R logistic_distrib.R cauchy_distrib.R laplace_distrib.R
#' @include lognormal1_distrib.R invgauss1_distrib.R poisson_distrib.R binomial_distrib.R
#' @include bernoulli_distrib.R negbin2_distrib.R student_t1_distrib.R pseudohuber_distrib.R
NULL

# ===========================================================================
# Derivatives of the distribution function with respect to the parameters.
#
# The identity that governs everything here is one exchange of derivative and
# integral. The region of integration does not depend on theta, so
#
#     d^I F(q) / F(q) = E[ d^I f / f  |  Y <= q ],
#
# a partial expectation of exactly the quantity the Bartlett lemma expands: at
# first order the conditional mean of the score, at second order the conditional
# mean of l^(ij) + l^(i) l^(j).
#
# Two consequences shape the implementation. For a *discrete* distribution the
# expectation is a finite sum, so the identity is not an approximation and the
# default method uses it directly. For a *continuous* one it is an integral over
# a semi-infinite region, and evaluating it by quadrature is both slower and less
# accurate than differencing the cdf, which for every distribution in the
# catalogue is an analytic function; the default there differences the cdf, and
# distributions with a closed form register it.
# ===========================================================================

#' Put CDF Derivatives on the Requested Tail and Scale
#'
#' @description
#' Converts derivatives of \eqn{F} into derivatives of whichever tail was asked
#' for, on the natural or the log scale.
#'
#' @details
#' Every route to a cdf derivative in this file produces derivatives of \eqn{F}
#' itself; the \code{lower.tail} and \code{log} arguments are handled once, here,
#' rather than in each of them. Switching tail flips the sign, since
#' \eqn{S = 1 - F}, and switching to the log scale divides by the probability,
#' which for a second derivative brings in the familiar
#' \eqn{d^2 \log P = d^2 P / P - (dP/P)(dP/P)}.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param Fq The cdf evaluated at the quantile.
#' @param dF1 A named list of first derivatives of \eqn{F}.
#' @param dF2 An optional named list of second derivatives of \eqn{F}.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of derivative component vectors.
#'
#' @seealso \code{\link{distrib_grad_cdf}}, \code{\link{distrib_hess_cdf}}
#' @keywords internal
cdf_tail_scale <- function(distrib, Fq, dF1, dF2 = NULL, lower.tail, log) {
  params <- distrib@params
  P <- if (lower.tail) Fq else 1 - Fq
  sgn <- if (lower.tail) 1 else -1

  d1 <- lapply(dF1, function(v) sgn * v)
  names(d1) <- params
  if (is.null(dF2)) {
    if (!log) return(d1)
    return(stats::setNames(lapply(params, function(p) d1[[p]] / P), params))
  }

  d2 <- lapply(dF2, function(v) sgn * v)
  names(d2) <- names(dF2)
  if (!log) return(d2)

  pairs <- hess_pairs(params)
  stats::setNames(lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    d2[[nm]] / P - (d1[[pr[1]]] / P) * (d1[[pr[2]]] / P)
  }), names(pairs))
}

# --- the exact route, for discrete distributions ----------------------------
#
# d^I F(q) = sum over the support up to q of  f(y) * (d^I f / f)(y).
# The sum is finite whenever the support has a finite lower bound, which every
# discrete distribution in the package has and the class requires.

#' CDF Derivatives of a Discrete Distribution
#'
#' @description
#' Evaluates \eqn{d^I F(q)} for a discrete distribution as the finite sum
#' \eqn{\sum_{y \le q} f(y) \, (d^I f / f)(y)}.
#'
#' @details
#' This is the governing identity
#' \eqn{d^I F(q) / F(q) = \mathbb{E}[d^I f / f \mid Y \le q]}
#' written out. For a discrete family the conditional expectation is a finite sum
#' whenever the support has a finite lower bound -- which the discrete class
#' requires -- so the identity is exact rather than an approximation, and it is
#' used directly.
#'
#' Note for anyone writing a test: checking this against the partial-expectation
#' sum proves nothing, because it is the same sum computed twice. A discrete
#' implementation has to be checked against finite differences of the cdf.
#'
#' @param distrib An object inheriting from class \code{"discrete_distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 or 2.
#'
#' @return A named list of derivative component vectors of \eqn{F}.
#'
#' @seealso \code{\link{cdf_tail_scale}}
#' @keywords internal
discrete_cdf_deriv <- function(distrib, q, theta, order) {
  params <- distrib@params
  lo <- distrib@bounds[1]
  n <- length(q)
  nms <- if (order == 1L) params else hess_names(params)
  pairs <- if (order == 2L) hess_pairs(params) else NULL

  out <- lapply(nms, function(nm) numeric(n))
  names(out) <- nms

  for (k in seq_len(n)) {
    th_k <- lapply(theta[seq_along(params)],
                   function(v) if (length(v) > 1) v[k] else v)
    if (q[k] < lo) next
    grid <- seq(lo, floor(q[k] + 1e-9))
    fy <- distrib_pdf(distrib, grid, th_k)
    g <- distrib_gradient(distrib, grid, th_k)
    if (order == 1L) {
      for (nm in nms) out[[nm]][k] <- sum(fy * g[[nm]])
    } else {
      h <- distrib_hessian(distrib, grid, th_k)
      for (nm in nms) {
        pr <- pairs[[nm]]
        out[[nm]][k] <- sum(fy * (h[[nm]] + g[[pr[1]]] * g[[pr[2]]]))
      }
    }
  }
  out
}

# --- the numerical route, for continuous distributions ---------------------

#' Numerical Derivatives of the Distribution Function
#'
#' @description
#' Central finite differences of \code{\link{distrib_cdf}} with respect to each
#' parameter. These power the default \code{\link{distrib_grad_cdf}} and
#' \code{\link{distrib_hess_cdf}} methods for continuous distributions that do
#' not supply a closed form.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order Either 1 or 2.
#' @param h_rel Numeric. Relative finite-difference step.
#' @param which Character vector of parameter names to differentiate, or
#'   \code{NULL} (default) for all of them. Used at first order by families that
#'   have a closed form for some parameters and not others, so that only the
#'   remaining ones cost a pair of cdf evaluations.
#' @return A named list of derivative components of \eqn{F}, not of its logarithm.
#' @seealso \code{\link{distrib_grad_cdf}}
#' @examples
#' numerical_cdf_deriv(gaussian1_distrib(), 1, list(mu = 0, sigma = 1), order = 1)
#'
#' @export
numerical_cdf_deriv <- function(distrib, q, theta, order = 1L,
                                h_rel = .Machine$double.eps^(1 / (order + 2)),
                                which = NULL) {
  params <- distrib@params
  p <- length(params)
  keep <- if (is.null(which)) params else which
  bump <- function(j, s, hj) {
    t2 <- theta
    t2[[j]] <- theta[[j]] + s * hj
    t2
  }
  # one step per observation, not per parameter: theta may be vectorised, and a
  # step chosen from its first element would be wrong everywhere else.
  hs <- lapply(seq_len(p), function(j) h_rel * pmax(1, abs(theta[[j]])))

  if (order == 1L) {
    js <- match(keep, params)
    out <- lapply(js, function(j) {
      (distrib_cdf(distrib, q, bump(j, 1, hs[[j]])) -
         distrib_cdf(distrib, q, bump(j, -1, hs[[j]]))) / (2 * hs[[j]])
    })
    return(stats::setNames(out, keep))
  }

  pairs <- hess_pairs(params)
  stats::setNames(lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    i <- match(pr[1], params); j <- match(pr[2], params)
    if (i == j) {
      (distrib_cdf(distrib, q, bump(i, 1, hs[[i]])) -
         2 * distrib_cdf(distrib, q, theta) +
         distrib_cdf(distrib, q, bump(i, -1, hs[[i]]))) / hs[[i]]^2
    } else {
      shift <- function(a, b) {
        t2 <- theta
        t2[[i]] <- theta[[i]] + a * hs[[i]]
        t2[[j]] <- theta[[j]] + b * hs[[j]]
        t2
      }
      (distrib_cdf(distrib, q, shift(1, 1)) - distrib_cdf(distrib, q, shift(1, -1)) -
         distrib_cdf(distrib, q, shift(-1, 1)) + distrib_cdf(distrib, q, shift(-1, -1))) /
        (4 * hs[[i]] * hs[[j]])
    }
  }), names(pairs))
}


# --- default methods -------------------------------------------------------

#' @title Default Log-CDF Gradient for Continuous Distributions
#' @name distrib_grad_cdf.continuous_distrib
#' @description Fallback: finite differences of the distribution function.
#' @param distrib A \code{continuous_distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @keywords internal
S7::method(distrib_grad_cdf, continuous_distrib) <- function(distrib, q, theta,
                                                             lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 numerical_cdf_deriv(distrib, q, theta, 1L),
                 NULL, lower.tail, log)
}

#' @title Default Log-CDF Hessian for Continuous Distributions
#' @name distrib_hess_cdf.continuous_distrib
#' @description Fallback: finite differences of the distribution function.
#' @param distrib A \code{continuous_distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @keywords internal
S7::method(distrib_hess_cdf, continuous_distrib) <- function(distrib, q, theta,
                                                             lower.tail = TRUE, log = TRUE) {
  # the first-order part comes from distrib_grad_cdf, so a closed form registered
  # for this family is used here too rather than differenced a second time
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 distrib_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = FALSE),
                 numerical_cdf_deriv(distrib, q, theta, 2L),
                 lower.tail, log)
}

#' @title Log-CDF Gradient for Discrete Distributions
#' @name distrib_grad_cdf.discrete_distrib
#' @description
#' Exact: the partial expectation of the score is a finite sum over the support
#' up to \eqn{q}, so nothing is differenced.
#' @param distrib A \code{discrete_distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @keywords internal
S7::method(distrib_grad_cdf, discrete_distrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 discrete_cdf_deriv(distrib, q, theta, 1L),
                 NULL, lower.tail, log)
}

#' @title Log-CDF Hessian for Discrete Distributions
#' @name distrib_hess_cdf.discrete_distrib
#' @description Exact, by the same finite sum as the gradient.
#' @param distrib A \code{discrete_distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @keywords internal
S7::method(distrib_hess_cdf, discrete_distrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 distrib_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = FALSE),
                 discrete_cdf_deriv(distrib, q, theta, 2L),
                 lower.tail, log)
}


# --- closed forms for location-scale families ------------------------------
#
# If F(y) = G(z) with z = (y - mu)/sigma, then differentiating the composition
# gives everything in terms of the density and its derivative in y:
#
#   dF/dmu     = -f
#   dF/dsigma  = -z f
#   d2F/dmu2   = f * l_y                        l_y = d log f / dy
#   d2F/dmu dsigma = f * (z l_y + 1/sigma)
#   d2F/dsigma2    = f * (z^2 l_y + 2z/sigma)
#
# so a family only has to say that it is location-scale. This covers the
# censored-regression workhorses -- Gaussian, logistic, Cauchy, Laplace -- and
# needs nothing from the family beyond its density and response derivative,
# both of which every distribution already provides.

#' CDF Derivatives of a Location-Scale Family
#'
#' @description
#' Closed-form derivatives of \eqn{F} for a family that is location-scale in its
#' first two parameters.
#'
#' @details
#' With \eqn{z = (q - \mu)/\sigma} and \eqn{\ell_y = \partial \log f/\partial y},
#' \deqn{\partial F/\partial \mu = -f, \qquad \partial F/\partial \sigma = -z f}
#' \deqn{\partial^2 F/\partial \mu^2 = f \ell_y, \qquad
#'       \partial^2 F/\partial \mu \partial \sigma = f (z \ell_y + 1/\sigma), \qquad
#'       \partial^2 F/\partial \sigma^2 = f (z^2 \ell_y + 2z/\sigma).}
#'
#' A family therefore only has to declare that it is location-scale; nothing is
#' needed beyond its density and its response derivative, both of which every
#' distribution already provides. This covers the censored-regression workhorses
#' -- Gaussian, logistic, Cauchy, Laplace.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, location first and scale second.
#' @param order The derivative order, 1 or 2.
#'
#' @return A named list of derivative component vectors of \eqn{F}.
#'
#' @seealso \code{\link{loc_scale_grad_cdf}}, \code{\link{loc_scale_hess_cdf}}
#' @keywords internal
loc_scale_cdf_deriv <- function(distrib, q, theta, order) {
  mu <- theta[[1]]
  s <- theta[[2]]
  z <- (q - mu) / s
  f <- distrib_pdf(distrib, q, theta)
  if (order == 1L) {
    return(stats::setNames(list(-f, -z * f), distrib@params))
  }
  ly <- distrib_grad_y(distrib, q, theta)
  vals <- list(f * ly, f * (z^2 * ly + 2 * z / s), f * (z * ly + 1 / s))
  stats::setNames(vals, hess_names(distrib@params))
}

#' Location-Scale CDF Gradient
#'
#' @description
#' The \code{\link{distrib_grad_cdf}} body shared by the location-scale families:
#' \code{\link{loc_scale_cdf_deriv}} at order 1, put on the requested tail and
#' scale.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of gradient component vectors.
#'
#' @keywords internal
loc_scale_grad_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 loc_scale_cdf_deriv(distrib, q, theta, 1L), NULL, lower.tail, log)
}

#' Location-Scale CDF Hessian
#'
#' @description
#' The \code{\link{distrib_hess_cdf}} body shared by the location-scale families.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of Hessian component vectors.
#'
#' @keywords internal
loc_scale_hess_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 loc_scale_cdf_deriv(distrib, q, theta, 1L),
                 loc_scale_cdf_deriv(distrib, q, theta, 2L), lower.tail, log)
}

#' @title Gaussian Log-CDF Derivatives
#' @name distrib_grad_cdf.Gaussian1Distrib
#' @description
#' Closed form, from the location-scale structure:
#' \eqn{\partial F/\partial\mu = -f(y)} and
#' \eqn{\partial F/\partial\sigma = -z f(y)} with \eqn{z = (y-\mu)/\sigma}.
#' @param distrib A \code{Gaussian1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{gaussian1_distrib}}
S7::method(distrib_grad_cdf, Gaussian1Distrib) <- loc_scale_grad_cdf

#' @title Gaussian Log-CDF Second Derivatives
#' @name distrib_hess_cdf.Gaussian1Distrib
#' @description Closed form, from the same location-scale structure.
#' @param distrib A \code{Gaussian1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{gaussian1_distrib}}
S7::method(distrib_hess_cdf, Gaussian1Distrib) <- loc_scale_hess_cdf

#' @title Logistic Log-CDF Derivatives
#' @name distrib_grad_cdf.LogisticDistrib
#' @description Closed form, from the location-scale structure.
#' @param distrib A \code{LogisticDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_grad_cdf, LogisticDistrib) <- loc_scale_grad_cdf

#' @title Logistic Log-CDF Second Derivatives
#' @name distrib_hess_cdf.LogisticDistrib
#' @description Closed form, from the location-scale structure.
#' @param distrib A \code{LogisticDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{logistic_distrib}}
S7::method(distrib_hess_cdf, LogisticDistrib) <- loc_scale_hess_cdf

#' @title Cauchy Log-CDF Derivatives
#' @name distrib_grad_cdf.CauchyDistrib
#' @description Closed form, from the location-scale structure.
#' @param distrib A \code{CauchyDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_grad_cdf, CauchyDistrib) <- loc_scale_grad_cdf

#' @title Cauchy Log-CDF Second Derivatives
#' @name distrib_hess_cdf.CauchyDistrib
#' @description Closed form, from the location-scale structure.
#' @param distrib A \code{CauchyDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{cauchy_distrib}}
S7::method(distrib_hess_cdf, CauchyDistrib) <- loc_scale_hess_cdf

#' @title Laplace Log-CDF Derivatives
#' @name distrib_grad_cdf.LaplaceDistrib
#' @description
#' Closed form, from the location-scale structure. Note that the second
#' derivatives inherit the kink at \eqn{y = \mu}, where \eqn{\partial\ell/\partial y}
#' does not exist; \code{\link{param_smoothness}} records this.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{b}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_grad_cdf, LaplaceDistrib) <- loc_scale_grad_cdf

#' @title Laplace Log-CDF Second Derivatives
#' @name distrib_hess_cdf.LaplaceDistrib
#' @description
#' Closed form, exact away from \eqn{y = \mu}. At the kink the second derivative
#' of \eqn{F} in \eqn{\mu} genuinely does not exist --- it jumps between
#' \eqn{\pm 1/(2b^{2})} --- so the value returned there is the one-sided limit
#' the sign convention picks out, and is reported rather than smoothed.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{b}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_hess_cdf, LaplaceDistrib) <- loc_scale_hess_cdf


# --- lognormal -------------------------------------------------------------
#
# F(y) = Phi(z) with z = (log y - mu)/sigma, so the family is location-scale on
# the log scale and the same two derivatives appear, expressed through the
# density of Y itself: phi(z) = y sigma f(y), whence
#
#   dF/dmu     = -y f(y)
#   dF/dsigma2 = -y f(y) z / (2 sigma)
#
# the second carrying the extra 1/(2 sigma) because the package parametrises by
# the variance on the log scale rather than the standard deviation.

#' @title Lognormal Log-CDF Gradient
#' @name distrib_grad_cdf.Lognormal1Distrib
#' @description
#' Closed form. On the log scale the lognormal is a location-scale family, so
#' \eqn{\partial F/\partial\mu = -y f(y)} and
#' \eqn{\partial F/\partial\sigma^{2} = -y f(y) z/(2\sigma)} with
#' \eqn{z = (\log y - \mu)/\sigma}.
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{lognormal1_distrib}}
S7::method(distrib_grad_cdf, Lognormal1Distrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  s <- sqrt(theta[[2]])
  z <- (base::log(q) - theta[[1]]) / s
  f <- distrib_pdf(distrib, q, theta)
  d1 <- list(mu = -q * f, sigma2 = -q * f * z / (2 * s))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

# --- inverse Gaussian ------------------------------------------------------
#
# Unusually for a positive family, the distribution function is elementary:
#
#   F(y) = Phi(a) + exp(2/(phi mu)) Phi(b),
#   a = (y/mu - 1)/sqrt(phi y),   b = -(y/mu + 1)/sqrt(phi y),
#
# so it can simply be differentiated. The exponential overflows for small
# phi*mu, so the product exp(2/(phi mu)) Phi(b) is formed on the log scale --
# the factor is huge exactly where Phi(b) is tiny.

#' @title Inverse Gaussian Log-CDF Gradient
#' @name distrib_grad_cdf.InvGauss1Distrib
#' @description
#' Closed form, obtained by differentiating the elementary representation
#' \eqn{F(y) = \Phi(a) + e^{2/(\phi\mu)}\Phi(b)}. The exponential factor is
#' combined with \eqn{\Phi(b)} on the log scale, since it overflows for small
#' \eqn{\phi\mu} exactly where \eqn{\Phi(b)} underflows.
#' @param distrib An \code{InvGauss1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{phi}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{invgauss1_distrib}}
S7::method(distrib_grad_cdf, InvGauss1Distrib) <- function(distrib, q, theta,
                                                          lower.tail = TRUE, log = TRUE) {
  mu <- theta[[1]]
  ph <- theta[[2]]
  rt <- sqrt(ph * q)
  a <- (q / mu - 1) / rt
  b <- -(q / mu + 1) / rt
  lg <- 2 / (ph * mu)
  tail_b <- exp(lg + stats::pnorm(b, log.p = TRUE))   # exp(lg) * Phi(b), stably
  eb <- exp(lg + stats::dnorm(b, log = TRUE))         # exp(lg) * phi(b), stably

  d_mu <- stats::dnorm(a) * (-q / (mu^2 * rt)) + tail_b * (-2 / (ph * mu^2)) +
    eb * (q / (mu^2 * rt))
  d_ph <- stats::dnorm(a) * (-a / (2 * ph)) + tail_b * (-2 / (ph^2 * mu)) +
    eb * (-b / (2 * ph))

  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 list(mu = d_mu, phi = d_ph), NULL, lower.tail, log)
}

# --- discrete families ------------------------------------------------------
#
# For the Poisson the sum defining F telescopes:
#   d/dmu sum_{j<=k} e^-mu mu^j/j!  =  sum_j [-f(j) + f(j-1)]  =  -f(k),
# so the sensitivity of the cdf to the mean is minus the mass at the last point
# retained. The binomial has the companion identity
#   d/dp P(X <= k) = -n dbinom(k, n-1, p),
# and the negative binomial, in the (mu, theta) parametrisation,
#   dF(k)/dmu = -f(k) (k + theta)/(theta + mu),
# which tends to the Poisson result as theta grows. The theta direction is a
# derivative of the incomplete beta in its parameter and has no elementary form,
# so it keeps the exact summation.

#' @title Poisson Log-CDF Gradient
#' @name distrib_grad_cdf.PoissonDistrib
#' @description
#' Closed form: the sum defining \eqn{F} telescopes, leaving
#' \eqn{\partial F(k)/\partial\mu = -f(k)}.
#' @param distrib A \code{PoissonDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list with one element.
#' @seealso \code{\link{poisson_distrib}}
S7::method(distrib_grad_cdf, PoissonDistrib) <- function(distrib, q, theta,
                                                         lower.tail = TRUE, log = TRUE) {
  d1 <- list(mu = -distrib_pdf(distrib, floor(q), theta) * (q >= 0))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Binomial Log-CDF Gradient
#' @name distrib_grad_cdf.BinomialDistrib
#' @description
#' Closed form: \eqn{\partial F(k)/\partial\mu = -n\,\mathrm{dbinom}(k; n-1, \mu)}.
#' @param distrib A \code{BinomialDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list with one element.
#' @seealso \code{\link{binomial_distrib}}
S7::method(distrib_grad_cdf, BinomialDistrib) <- function(distrib, q, theta,
                                                          lower.tail = TRUE, log = TRUE) {
  n <- distrib@size
  d1 <- list(mu = -n * stats::dbinom(floor(q), n - 1, theta[[1]]) * (q >= 0))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Bernoulli Log-CDF Gradient
#' @name distrib_grad_cdf.BernoulliDistrib
#' @description
#' Closed form, the binomial identity at \eqn{n = 1}: the derivative is
#' \eqn{-1} at \eqn{k = 0} and zero at \eqn{k = 1}.
#' @param distrib A \code{BernoulliDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list with one element.
#' @seealso \code{\link{bernoulli_distrib}}
S7::method(distrib_grad_cdf, BernoulliDistrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  d1 <- list(mu = -stats::dbinom(floor(q), 0, theta[[1]]) * (q >= 0))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Negative Binomial Log-CDF Gradient
#' @name distrib_grad_cdf.NegBin2Distrib
#' @description
#' Closed form in \eqn{\mu}, \eqn{\partial F(k)/\partial\mu = -f(k)(k+\theta)/(\theta+\mu)},
#' which reduces to the Poisson identity as \eqn{\theta\to\infty}. The \eqn{\theta}
#' direction is a derivative of the incomplete beta in its parameter, has no
#' elementary form, and keeps the exact summation.
#' @param distrib A \code{NegBin2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{theta}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{negbin2_distrib}}
S7::method(distrib_grad_cdf, NegBin2Distrib) <- function(distrib, q, theta,
                                                        lower.tail = TRUE, log = TRUE) {
  k <- floor(q)
  f <- distrib_pdf(distrib, k, theta)
  d1 <- list(mu = -f * (k + theta[[2]]) / (theta[[2]] + theta[[1]]) * (q >= 0),
             theta = discrete_cdf_deriv(distrib, q, theta, 1L)[["theta"]])
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

# --- families that are location-scale in only some of their parameters -----
#
# The Student t and the pseudo-Huber are location-scale in (mu, sigma) with a
# further shape parameter nu. The two location-scale derivatives are closed
# form; nu is a derivative of a hypergeometric-type integral and is differenced.

#' CDF Gradient When Only Some Parameters Are Location-Scale
#'
#' @description
#' For a family that is location-scale in \eqn{(\mu, \sigma)} but carries a
#' further shape parameter: the two location-scale directions in closed form, the
#' shape direction by finite differences.
#'
#' @details
#' The Student t and the pseudo-Huber are the two. Their shape direction is a
#' derivative of a hypergeometric-type integral with no elementary form, so it is
#' differenced; the other two are exact.
#'
#' For the pseudo-Huber this is an \strong{accuracy} gain rather than merely a
#' speed one. Its cdf is itself a quadrature, so differencing it is good to only
#' about \code{1e-6}, whereas \eqn{\partial F/\partial \mu = -f} is exact.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of gradient component vectors.
#'
#' @keywords internal
partial_loc_scale_grad_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
  params <- distrib@params
  z <- (q - theta[[1]]) / theta[[2]]
  f <- distrib_pdf(distrib, q, theta)
  rest <- params[-(1:2)]
  d1 <- c(list(-f, -z * f), numerical_cdf_deriv(distrib, q, theta, 1L, which = rest))
  names(d1) <- params
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Student t Log-CDF Gradient
#' @name distrib_grad_cdf.StudentT1Distrib
#' @description
#' Closed form in the location and scale, \eqn{-f(y)} and \eqn{-z f(y)}; the
#' degrees of freedom are differenced, having no elementary form.
#' @param distrib A \code{StudentT1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{sigma} and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{student_t1_distrib}}
S7::method(distrib_grad_cdf, StudentT1Distrib) <- partial_loc_scale_grad_cdf

#' @title Pseudo-Huber Log-CDF Gradient
#' @name distrib_grad_cdf.PseudoHuberDistrib
#' @description
#' Closed form in the location and scale; the shape \eqn{\nu} is differenced.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{sigma} and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_grad_cdf, PseudoHuberDistrib) <- partial_loc_scale_grad_cdf
