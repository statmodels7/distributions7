#' @include distrib.R generics.R utility_functions.R numerical_functions.R zero_inflated.R zero_adjusted.R
NULL

# ===========================================================================
# Truncation.
#
# Unlike the zero wrappers, truncation adds **no parameter**: the truncation
# points are known constants, like a binomial's `size`. What it adds instead is
# a normalizing constant that depends on theta,
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
#' Represents a continuous parent restricted to \eqn{[L, U]} and renormalized
#' by the retained mass \eqn{Z(\theta)}. Construct one with [truncated()], which
#' validates the endpoints, collapses a nested truncation and copies the
#' parent's parameter metadata; calling the class directly does none of that.
#'
#' @details
#' # What truncation adds, and what it does not
#'
#' It adds NO parameter. The endpoints are known constants, like a binomial's
#' `size`, so the truncated object carries exactly the parent's `params`,
#' `params_bounds` and `link_params`. What it adds is the
#' \eqn{\theta}-dependent normalizing constant \eqn{Z}, and every derivative of
#' \eqn{\ell_T = \ell - \log Z} carries its contribution.
#'
#' The support does not move with \eqn{\theta}. That is the condition under
#' which \eqn{Z} may be differentiated under the integral sign, and it is what
#' keeps truncation at fixed points a regular problem.
#'
#' # A mixed parent
#'
#' The parent may itself carry point masses, as [zero_adjusted()] of a
#' continuous distribution does. Those atoms survive truncation where they lie
#' inside the interval, rescaled by \eqn{1/Z}, and
#' [distrib_atoms.TruncatedContinuousDistrib()] reports them. Two methods exist
#' only for that case: [expectation.TruncatedContinuousDistrib()] adds the
#' masses to the integral, and [parent_mass_at()] asks the parent for the mass
#' on a single point instead of assuming a continuous parent has none.
#'
#' @inheritParams distrib
#' @param parent_distrib The wrapped `continuous_distrib` object. Its
#'   parameters become the truncated object's, unchanged.
#' @param lower,upper The truncation points, `L` and `U`. Either may be
#'   infinite, giving one-sided truncation, and both are included in the
#'   support.
#'
#' @return An S7 object of class `TruncatedContinuousDistrib`, inheriting from
#'   `continuous_distrib`.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [truncated()] for the constructor, [TruncatedDiscreteDistrib] for
#'   the discrete branch, and [trunc_constants()] for \eqn{Z}.
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_atoms()`][distrib_atoms.TruncatedContinuousDistrib],
#'   [`distrib_cdf()`][distrib_cdf.TruncatedContinuousDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.TruncatedContinuousDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.TruncatedContinuousDistrib],
#'   [`distrib_gradient()`][distrib_gradient.TruncatedContinuousDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.TruncatedContinuousDistrib],
#'   [`distrib_hessian()`][distrib_hessian.TruncatedContinuousDistrib],
#'   [`distrib_pdf()`][distrib_pdf.TruncatedContinuousDistrib],
#'   [`distrib_quantile()`][distrib_quantile.TruncatedContinuousDistrib],
#'   [`distrib_rng()`][distrib_rng.TruncatedContinuousDistrib],
#'   [`expectation()`][expectation.TruncatedContinuousDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' class(tn)[1]
#' c(name = tn@distrib_name, params = paste(tn@params, collapse = ", "))
#' tn@bounds
#'
#' # The parameters are the parent's, unchanged: truncation adds none.
#' identical(tn@params, gaussian1_distrib()@params)
#' identical(tn@params_bounds, gaussian1_distrib()@params_bounds)
#'
#' # The density is the parent's divided by the retained mass.
#' Z <- pnorm(2, 0.3, 1.2) - pnorm(-1, 0.3, 1.2)
#' c(truncated = distrib_pdf(tn, 0, theta),
#'   parent_over_Z = dnorm(0, 0.3, 1.2) / Z)
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
#' Represents a discrete parent restricted to the support points in
#' \eqn{[L, U]} and renormalized by the retained mass \eqn{Z(\theta)}. Both
#' endpoints are INCLUDED, so `truncated(poisson_distrib(), lower = 1)` is the
#' zero-truncated Poisson on \eqn{\{1, 2, \dots\}} and keeps the mass at one.
#' Construct one with [truncated()]; calling the class directly skips every
#' validation the constructor performs.
#'
#' @details
#' The lower endpoint being included is the one place the two truncation
#' classes differ. The tail below the interval is \eqn{F(L^-) = F(L) - f(L)},
#' so the mass sitting exactly on \eqn{L} has to be added back;
#' [trunc_constants()] does that through [parent_mass_at()], and
#' [trunc_mass_derivs()] applies the matching correction to the derivatives of
#' \eqn{Z}.
#'
#' Truncation adds no parameter, so the object carries exactly the parent's
#' `params`, `params_bounds` and `link_params`. What it can remove is
#' IDENTIFIABILITY: \eqn{k} retained support points carry \eqn{k-1} free
#' probabilities, so [truncated()] rejects an interval leaving fewer than
#' `n_params + 1` of them.
#'
#' @inheritParams distrib
#' @param parent_distrib The wrapped `discrete_distrib` object. Its parameters
#'   become the truncated object's, unchanged.
#' @param lower,upper The truncation points, `L` and `U`, both included in the
#'   support. Either may be infinite, and a finite one must be a whole number.
#'
#' @return An S7 object of class `TruncatedDiscreteDistrib`, inheriting from
#'   `discrete_distrib`.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [truncated()] for the constructor, [TruncatedContinuousDistrib] for
#'   the continuous branch, and [trunc_constants()] for \eqn{Z}.
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.TruncatedDiscreteDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.TruncatedDiscreteDistrib],
#'   [`distrib_gradient()`][distrib_gradient.TruncatedDiscreteDistrib],
#'   [`distrib_hessian()`][distrib_hessian.TruncatedDiscreteDistrib],
#'   [`distrib_pdf()`][distrib_pdf.TruncatedDiscreteDistrib],
#'   [`distrib_quantile()`][distrib_quantile.TruncatedDiscreteDistrib],
#'   [`distrib_rng()`][distrib_rng.TruncatedDiscreteDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
#'
#' @examples
#' ztp <- truncated(poisson_distrib(), lower = 1)
#' class(ztp)[1]
#' c(name = ztp@distrib_name, bounds = paste(ztp@bounds, collapse = ", "))
#'
#' # The lower endpoint is IN the support, so this is the zero-truncated
#' # Poisson: mass at 0 is gone, mass at 1 is not.
#' distrib_pdf(ztp, 0:3, list(mu = 2))
#' dpois(1:3, 2) / (1 - dpois(0, 2))
#'
#' # An interval leaving too few support points is rejected by the
#' # constructor, the parameter no longer being identified.
#' try(truncated(bernoulli_distrib(), lower = 1))
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

#' @title Is This Distribution Already Truncated?
#'
#' @description
#' Reports whether an object belongs to either truncated class. [truncated()]
#' asks before wrapping, and collapses a nested truncation into a single
#' object over the intersection of the two intervals. Nesting would be correct
#' but would pay the quadrature cost of [trunc_score_mean()] twice for a law
#' one truncation already describes.
#'
#' @param distrib An object inheriting from class `distrib`. Any other input
#'   returns `FALSE`; nothing raises.
#'
#' @return A single logical.
#'
#' @seealso [truncated()], which uses it, and [TruncatedContinuousDistrib] and
#'   [TruncatedDiscreteDistrib], the two classes it recognizes.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' c(truncated = distributions7:::is_truncated(tn),
#'   parent = distributions7:::is_truncated(gaussian1_distrib()))
#'
#' # Nesting collapses: one object, over the intersection.
#' t2 <- truncated(truncated(gaussian1_distrib(), lower = -1), upper = 2)
#' c(lower = t2@lower, upper = t2@upper)
#' distributions7:::is_truncated(t2@parent_distrib)
is_truncated <- function(distrib) {
  S7::S7_inherits(distrib, TruncatedContinuousDistrib) ||
    S7::S7_inherits(distrib, TruncatedDiscreteDistrib)
}

# P(Y = x) under the parent: the pmf for a discrete distribution, the atom's
# probability for a mixed one, and zero for an ordinary continuous distribution.
#
# This is the one quantity that separates the two truncation classes, and getting
# it wrong for a *mixed* parent is subtle: the cdf of zero_adjusted(gamma) already
# includes the point mass at zero, so F(0) is not F(0^-) even though the parent is
# a continuous_distrib. Truncating it above, with the atom retained, then loses
# exactly that mass from the normalizing constant -- and the resulting density
# integrates to something other than one while every formula still looks right.
#' @title Probability the Parent Puts on a Single Point
#'
#' @description
#' Returns \eqn{P(Y = x)} under the parent: the mass function for a discrete
#' parent, the atom's probability for a mixed one, and zero for an ordinary
#' continuous parent. [trunc_constants()] needs it because the lower endpoint
#' is included in the truncated support, so the tail below the interval is
#' \eqn{F(L^-) = F(L) - P(Y = L)}.
#'
#' @details
#' The question is asked of the OBJECT, through [distrib_atoms()], never of its
#' class. Branching on `discrete_distrib` looks right and is wrong: the
#' distribution function of `zero_adjusted(gamma2_distrib())` already contains
#' the point mass at zero, so \eqn{F(0) \ne F(0^-)} although the object is a
#' `continuous_distrib`. Truncating such a parent from above, with the atom
#' retained, would then drop exactly that mass out of \eqn{Z} while every
#' formula still read correctly, and the density would integrate to something
#' other than one.
#'
#' @param distrib A truncated distribution object, of either class. The branch
#'   is taken on it, and the mass is asked of `distrib@parent_distrib`.
#' @param x A single number, the point to evaluate at.
#' @param theta A named list of the parent's parameters.
#'
#' @return A numeric vector of probabilities, of length one unless a parameter
#'   in `theta` varies by observation.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [distrib_atoms()] for the declaration it reads, and
#'   [trunc_constants()] for its one caller.
#'
#' @keywords internal
#'
#' @examples
#' # A discrete parent: the mass function.
#' ztp <- truncated(poisson_distrib(), lower = 1)
#' c(reported = distributions7:::parent_mass_at(ztp, 1, list(mu = 2)),
#'   dpois = dpois(1, 2))
#'
#' # A mixed parent: the atom, although the object is a continuous_distrib.
#' tz <- truncated(zero_adjusted(gamma2_distrib()), upper = 5)
#' th <- list(mu = 2, sigma2 = 1, za = 0.3)
#' distributions7:::parent_mass_at(tz, 0, th)
#'
#' # An ordinary continuous parent: zero, at every point.
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#' distributions7:::parent_mass_at(tn, 0, theta)
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
#' @title The Truncation Constant and Lower Tail
#'
#' @description
#' Returns the tail below the interval, \eqn{F(L^-)}, and the retained mass
#' \eqn{Z = F(U) - F(L^-)}. Every method of a truncated distribution divides
#' by \eqn{Z}, and [trunc_cdf()] and [trunc_quantile()] shift by \eqn{F(L^-)}
#' as well.
#'
#' @details
#' Both endpoints are INCLUDED in the truncated support, so any mass sitting
#' exactly on the lower one is added back:
#' \eqn{F(L^-) = F(L) - P(Y = L)}. That correction belongs to the ATOM case,
#' not to the discrete case, so it goes through [parent_mass_at()]. An
#' infinite endpoint contributes \eqn{0} or \eqn{1} without a call on the
#' parent.
#'
#' Both quantities are vectorized in \eqn{\theta}, so a parameter varying by
#' observation gives one constant per observation.
#'
#' An interval carrying no probability under the given parameters raises an
#' error naming the interval and the computed mass, the truncated law not
#' being defined there. That is a likely place for a search to wander to, so
#' the message says which endpoints produced it.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters. A component may vary
#'   by observation.
#'
#' @return A named list with `Fl`, the tail below the interval, and `Z`, the
#'   retained mass; each a numeric vector following the length of `theta`.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [parent_mass_at()] for the endpoint correction,
#'   [trunc_mass_derivs()] for the derivatives of \eqn{Z}, and [truncated()].
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' cs <- distributions7:::trunc_constants(tn, theta)
#' unlist(cs)
#' c(Z = cs$Z, direct = pnorm(2, 0.3, 1.2) - pnorm(-1, 0.3, 1.2))
#'
#' # A discrete parent: the mass at the lower endpoint is added back, so the
#' # zero-truncated Poisson keeps everything except dpois(0, mu).
#' ztp <- truncated(poisson_distrib(), lower = 1)
#' cz <- distributions7:::trunc_constants(ztp, list(mu = 2))
#' c(Z = cz$Z, direct = 1 - dpois(0, 2), Fl = cz$Fl)
#'
#' # Vectorized in theta.
#' distributions7:::trunc_constants(tn, list(mu = c(0, 0.5), sigma = 1.2))$Z
#'
#' # An interval carrying no mass is reported, not returned.
#' far <- truncated(gaussian1_distrib(), lower = 100, upper = 200)
#' try(distributions7:::trunc_constants(far, theta))
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

#' @title Which Observations Lie in the Truncated Support
#'
#' @description
#' Tests `y >= L & y <= U`, with both endpoints included. [trunc_pdf()] uses it
#' to set the log-density to \eqn{-\infty} outside the interval,
#' [trunc_y_deriv()] to return `NaN` there, and
#' [distrib_atoms.TruncatedContinuousDistrib()] to keep the parent's atoms that
#' survive.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param y A numeric vector of observations.
#'
#' @return A logical vector as long as `y`.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_pdf()] and [trunc_y_deriv()], its two callers.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' # Both endpoints count as inside.
#' distributions7:::trunc_inside(tn, c(-2, -1, 0, 2, 3))
#'
#' # Which is what makes the density positive at them and zero beyond.
#' distrib_pdf(tn, c(-2, -1, 2, 3), theta)
trunc_inside <- function(distrib, y) {
  y >= distrib@lower & y <= distrib@upper
}

# Derivatives of the retained mass, from the parent's cdf derivatives.
#
#   d^B Z = d^B F(U) - d^B F(L^-)
#
# which turns the truncated expectations E_T[d^B f / f] = d^B Z / Z into two
# calls on the parent rather than one quadrature per component. The lower
# endpoint keeps the same correction as the mass itself: F(L^-) = F(L) - P(Y = L),
# so its derivatives lose the derivatives of that mass, which for a discrete distribution
# parent are f(L) l^(i)(L) and f(L)(l^(ij) + l^(i)l^(j)) at L.
#
# Returns NULL when the route should not be taken, and the caller falls back to
# quadrature. Two situations call for that.
#
# The first is correctness: a mixed parent with an atom sitting exactly on the
# lower endpoint, whose mass derivative is not the parent's pmf.
#
# The second is accuracy, and it is the reason this is a choice rather than a
# replacement. When the parent has no closed-form cdf derivative the route above
# differences its cdf, carrying roughly 1e-8 of relative error into the Hessian,
# where the quadrature it replaces carried 1e-10. That is invisible in the
# Hessian itself but not downstream: numerical_deriv4() differentiates the
# analytical Hessian, so a noisier Hessian degrades the *reference* the
# fourth-order check compares against, and the check fails on code that is
# right. The route is therefore taken only where it is at least as accurate as
# what it replaces -- a parent with a genuine closed form, or a discrete parent,
# whose cdf derivatives are an exact finite sum.
#' @title Can the Parent Supply Exact CDF Derivatives?
#'
#' @description
#' Reports whether the parent has a genuine closed-form derivative of its
#' distribution function at the given order, or is a discrete family, whose cdf
#' derivatives are an exact finite sum. [trunc_mass_derivs()] takes the cheap
#' route only where this is `TRUE`, and falls back to quadrature otherwise.
#'
#' @details
#' # Why the cheap route is gated
#'
#' Reading \eqn{d^B Z} off the parent's cdf derivatives replaces one quadrature
#' per component with two calls on the parent: measured on a truncated
#' gaussian's Hessian, 1.4 ms against 4.9 ms. Where the parent has no closed
#' form, though, that route differences its distribution function and carries
#' roughly \eqn{10^{-8}} of relative error into the Hessian, where the
#' quadrature it replaced carried \eqn{10^{-10}}.
#'
#' That is invisible in the Hessian itself, and visible downstream:
#' [numerical_deriv4()] differentiates the analytical Hessian, so a noisier
#' Hessian degrades the REFERENCE the fourth-order check compares against, and
#' the check then fails on code that is right.
#'
#' # How a genuine method is recognized
#'
#' `attr(m, "signature")[[1]]` is the class the method was registered on, so an
#' inherited fallback answers with a base class. `identical()` on the method
#' object cannot be used for this, S7 wrapping it. The third- and
#' fourth-order defaults sit on `distrib` where the first two sit on
#' `continuous_distrib`, so both base classes are excluded or a stencil would
#' be read as a closed form.
#'
#' The question is asked of the OWNING CLASS, so a family that registers a
#' method combining closed components with a stencil in one direction answers
#' `TRUE`. That is the intended reading: what the gate protects against is the
#' generic differencing the cdf itself.
#'
#' @param parent The parent distribution, before wrapping.
#' @param order The derivative order, an integer from 1 to 4.
#'
#' @return A single logical.
#'
#' @seealso [trunc_mass_derivs()], its one caller, and [distrib_grad_cdf()] for
#'   the generics it asks about.
#'
#' @keywords internal
#'
#' @examples
#' # A discrete family answers TRUE at every order: its cdf derivatives are an
#' # exact sum.
#' vapply(1:4, function(k)
#'   distributions7:::has_exact_cdf_deriv(poisson_distrib(), k), logical(1))
#'
#' # The gaussian writes all four out.
#' vapply(1:4, function(k)
#'   distributions7:::has_exact_cdf_deriv(gaussian1_distrib(), k), logical(1))
#'
#' # The gamma does not, the derivative of an incomplete gamma in its shape
#' # having no elementary form, so truncation falls back to quadrature there.
#' vapply(1:4, function(k)
#'   distributions7:::has_exact_cdf_deriv(gamma2_distrib(), k), logical(1))
has_exact_cdf_deriv <- function(parent, order) {
  if (S7::S7_inherits(parent, discrete_distrib)) return(TRUE)
  gen <- switch(order, distrib_grad_cdf, distrib_hess_cdf,
                distrib_deriv3_cdf, distrib_deriv4_cdf)
  m <- tryCatch(S7::method(gen, S7::S7_class(parent)), error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  reg <- tryCatch(attr(m, "signature")[[1]], error = function(e) NULL)
  if (is.null(reg)) return(FALSE)
  # identical() on a method object does not answer "is this the fallback?" --
  # S7 wraps it -- but the class it was registered on does. The defaults of
  # the third and fourth orders sit on `distrib` rather than on
  # `continuous_distrib`, so both have to be excluded or a stencil would be
  # mistaken for a closed form.
  !is_class(reg, continuous_distrib) && !is_class(reg, distrib)
}

#' @title Derivatives of the Truncation Constant via the Parent's CDF
#'
#' @description
#' Computes \eqn{d^B Z = d^B F(U) - d^B F(L^-)} from the parent's cdf
#' derivatives, or returns `NULL` when that route is not available. Dividing
#' the result by \eqn{Z} gives the truncated expectations
#' [trunc_score_mean()] and [trunc_M()] need, replacing one quadrature per
#' component with two calls on the parent.
#'
#' @details
#' # The discrete correction
#'
#' The lower endpoint is included in the truncated support, so what leaves
#' \eqn{Z} is \eqn{F(L)} minus the mass at \eqn{L}. Its derivatives lose that
#' mass's derivatives with it, \eqn{d^I F(L^-) = d^I F(L) - f(L) B_I}, where
#' \eqn{B_I} is the complete Bell polynomial in the parent's log-mass
#' derivatives that [bell_f_ratio()] assembles.
#'
#' # When it declines
#'
#' `NULL` is returned, and the caller falls back to quadrature, in two
#' situations. The first is accuracy, and is decided by
#' [has_exact_cdf_deriv()]. The second is correctness: a mixed parent with an
#' atom sitting exactly on the lower endpoint, whose mass derivative is not the
#' parent's own mass function.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters.
#' @param order The derivative order, 1 or 2.
#'
#' @return A named list of derivative components of \eqn{Z}, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names], or `NULL` when the
#'   route is declined.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [has_exact_cdf_deriv()] for the gate, [trunc_score_mean()] and
#'   [trunc_M()] for the two callers, and [trunc_constants()] for \eqn{Z}.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' Z <- distributions7:::trunc_constants(tn, theta)$Z
#' dZ <- distributions7:::trunc_mass_derivs(tn, theta, 1L)
#' unlist(dZ) / Z
#'
#' # The same numbers the quadrature route returns, to machine precision.
#' unlist(distributions7:::trunc_score_mean_quad(tn, theta))
#'
#' # A gamma parent has no closed-form cdf derivative, so the route declines
#' # and the caller integrates instead.
#' tg <- truncated(gamma2_distrib(), lower = 0.5, upper = 5)
#' is.null(distributions7:::trunc_mass_derivs(tg, list(mu = 2, sigma2 = 1), 1L))
trunc_mass_derivs <- function(distrib, theta, order) {
  if (!has_exact_cdf_deriv(distrib@parent_distrib, order)) return(NULL)
  parent <- distrib@parent_distrib
  params <- distrib@params
  lo <- distrib@lower
  up <- distrib@upper
  is_disc <- S7::S7_inherits(distrib, TruncatedDiscreteDistrib)
  nms <- deriv_names(params, order)
  idx <- deriv_indices(params, order)

  if (!is_disc && is.finite(lo)) {
    at <- distrib_atoms(parent, theta)
    if (length(at$y) && any(at$y == lo)) return(NULL)
  }

  gen <- switch(order, distrib_grad_cdf, distrib_hess_cdf,
                distrib_deriv3_cdf, distrib_deriv4_cdf)
  zero <- stats::setNames(lapply(nms, function(nm) 0), nms)
  dU <- if (is.infinite(up)) zero else gen(parent, up, theta, log = FALSE)
  dL <- if (is.infinite(lo)) zero else {
    d <- gen(parent, lo, theta, log = FALSE)
    if (is_disc) {
      # the lower point is IN the truncated support, so what leaves Z is
      # F(lo) minus the mass at lo: d^I F(lo^-) = d^I F(lo) - f(lo) B_I,
      # the same complete Bell polynomial the sum itself carries
      f0 <- distrib_pdf(parent, lo, theta)
      ell <- parent_ell(parent, lo, theta, order, params)
      for (m in seq_along(nms)) {
        d[[nms[m]]] <- d[[nms[m]]] - f0 * bell_f_ratio(params[idx[[m]]], ell)
      }
    }
    d
  }
  stats::setNames(lapply(nms, function(nm) dU[[nm]] - dL[[nm]]), nms)
}

#' @title Mean of the Parent's Score Under the Truncated Law
#'
#' @description
#' Computes \eqn{m_i = \mathbb{E}_T[s_i]}, the quantity that recenters the
#' parent's score: the truncated score is \eqn{s_i(y) - m_i}. It is
#' \eqn{d_i Z / Z}, because \eqn{d_i Z = \int_T f s_i}, so it comes from
#' [trunc_mass_derivs()] where those are exact and from
#' [trunc_score_mean_quad()] otherwise.
#'
#' @details
#' \eqn{m_i} is a function of \eqn{\theta} and of the endpoints alone, not of
#' the data, so one evaluation serves a whole sample. It is also what makes
#' derivatives of a truncated distribution dearer than the parent's: the
#' quadrature route costs one integral per parameter.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters.
#'
#' @return A named list of numeric vectors, one component per parameter, in
#'   `distrib@params` order.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_gradient()], which subtracts it, [trunc_M()] for the
#'   second-order counterpart, and [trunc_score_mean_quad()] for the fallback.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' m <- distributions7:::trunc_score_mean(tn, theta)
#' unlist(m)
#'
#' # The truncated score is the parent's, recentered by exactly this.
#' y <- c(-0.5, 0, 1)
#' g <- distrib_gradient(tn, y, theta)
#' gp <- distrib_gradient(gaussian1_distrib(), y, theta)
#' max(abs(unlist(g) - (unlist(gp) - rep(unlist(m), each = length(y)))))
#'
#' # It is also what makes the truncated score have mean zero, which the
#' # untruncated score does not have over the interval.
#' set.seed(1)
#' ys <- distrib_rng(tn, 50000, theta)
#' round(vapply(distrib_gradient(tn, ys, theta), mean, numeric(1)), 3)
trunc_score_mean <- function(distrib, theta) {
  dZ <- trunc_mass_derivs(distrib, theta, 1L)
  if (!is.null(dZ)) {
    Z <- trunc_constants(distrib, theta)$Z
    return(stats::setNames(lapply(distrib@params, function(p) dZ[[p]] / Z), distrib@params))
  }
  trunc_score_mean_quad(distrib, theta)
}

#' @title Second-Order Truncated Moment of the Parent's Derivatives
#'
#' @description
#' Computes \eqn{M_{ij} = \mathbb{E}_T[H_{ij} + s_i s_j]}, which is
#' \eqn{d_{ij} Z / Z}, since \eqn{d_{ij} f / f = H_{ij} + s_i s_j}. It enters
#' the truncated Hessian as
#' \eqn{d_{ij}\ell_T = H_{ij}(y) - M_{ij} + m_i m_j}, the subtracted part being
#' the second derivative of \eqn{\log Z}.
#'
#' @details
#' Where the parent has exact cdf derivatives the whole quantity comes from
#' [trunc_mass_derivs()] at order two. Otherwise it is assembled from two
#' quadratures, [trunc_hess_mean()] and [trunc_score_prod_mean()], which is
#' also why [trunc_expected_hessian()] cannot use the cheap route: it needs the
#' two pieces separately and \eqn{d_{ij} Z} only gives their sum.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters.
#'
#' @return A named list of numeric vectors, one component per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_hessian()], which subtracts it, [trunc_score_mean()] for the
#'   first order, and [trunc_hess_mean()] and [trunc_score_prod_mean()] for the
#'   two halves of the fallback.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' M <- distributions7:::trunc_M(tn, theta)
#' unlist(M)
#'
#' # It is the sum of the two quadratures, whichever route was taken.
#' EH <- distributions7:::trunc_hess_mean(tn, theta)
#' ES <- distributions7:::trunc_score_prod_mean(tn, theta)
#' max(abs(unlist(M) - (unlist(EH) + unlist(ES))))
trunc_M <- function(distrib, theta) {
  dZ <- trunc_mass_derivs(distrib, theta, 2L)
  if (!is.null(dZ)) {
    Z <- trunc_constants(distrib, theta)$Z
    nms <- hess_names(distrib@params)
    return(stats::setNames(lapply(nms, function(nm) dZ[[nm]] / Z), nms))
  }
  EH <- trunc_hess_mean(distrib, theta)
  ES <- trunc_score_prod_mean(distrib, theta)
  nms <- hess_names(distrib@params)
  stats::setNames(lapply(nms, function(nm) EH[[nm]] + ES[[nm]]), nms)
}

#' @title Truncated Score Mean by Quadrature
#'
#' @description
#' Computes \eqn{m_i = \mathbb{E}_T[s_i]} by integrating the parent's score
#' against the truncated law, one integral per parameter. It is the route
#' [trunc_score_mean()] falls back on for a parent whose cdf derivatives are
#' not exact, and is the only route for a gamma or a beta, whose distribution
#' function has no elementary derivative in its shape.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters.
#'
#' @return A named list of numeric vectors, one component per parameter, in
#'   `distrib@params` order.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_score_mean()], which prefers the cdf route where it is
#'   exact, and [expectation()], which performs the integration.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' # The two routes agree; only the cost differs.
#' unlist(distributions7:::trunc_score_mean_quad(tn, theta))
#' unlist(distributions7:::trunc_score_mean(tn, theta))
#'
#' # For a gamma parent this is the only route available.
#' tg <- truncated(gamma2_distrib(), lower = 0.5, upper = 5)
#' unlist(distributions7:::trunc_score_mean_quad(tg, list(mu = 2, sigma2 = 1)))
trunc_score_mean_quad <- function(distrib, theta) {
  parent <- distrib@parent_distrib
  out <- lapply(distrib@params, function(p) {
    expectation(distrib, function(y, theta) distrib_gradient(parent, y, theta)[[p]], theta)
  })
  stats::setNames(out, distrib@params)
}

#' @title Truncated Mean of Products of Scores
#'
#' @description
#' Computes \eqn{\mathbb{E}_T[s_i s_j]} for every unordered pair of parameters,
#' by quadrature. Two consumers need it: [trunc_M()] adds it to
#' [trunc_hess_mean()] where the cdf route is unavailable, and
#' [trunc_expected_hessian()] uses it whatever the parent is, the expected
#' Hessian being \eqn{-\mathrm{Cov}_T(s_i, s_j)}.
#'
#' @details
#' This is the quantity that keeps a truncated expected information expensive
#' even for a parent with closed-form cdf derivatives. Those give
#' \eqn{d_{ij} Z}, which is \eqn{\mathbb{E}_T[H_{ij}] + \mathbb{E}_T[s_i s_j]},
#' and no rearrangement separates the two terms; the covariance needs the
#' second alone.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters.
#'
#' @return A named list of numeric vectors, one component per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_expected_hessian()] and [trunc_M()], its two consumers, and
#'   [trunc_hess_mean()] for the other half of the second-order moment.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' ES <- distributions7:::trunc_score_prod_mean(tn, theta)
#' unlist(ES)
#'
#' # Against the same expectation taken by simulation.
#' set.seed(1)
#' ys <- distrib_rng(tn, 100000, theta)
#' g <- distrib_gradient(gaussian1_distrib(), ys, theta)
#' round(c(mu_mu = mean(g$mu^2), sigma_sigma = mean(g$sigma^2),
#'         mu_sigma = mean(g$mu * g$sigma)), 3)
#' round(unlist(ES), 3)
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

#' @title Truncated Mean of the Parent's Hessian
#'
#' @description
#' Computes \eqn{\mathbb{E}_T[H_{ij}]} for every unordered pair of parameters,
#' by quadrature. [trunc_M()] adds it to [trunc_score_prod_mean()] where the
#' parent's cdf derivatives are not exact, that sum being
#' \eqn{d_{ij} Z / Z}.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param theta A named list of the parent's parameters.
#'
#' @return A named list of numeric vectors, one component per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_M()], its one caller, and [trunc_score_prod_mean()] for the
#'   other half of the sum.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' EH <- distributions7:::trunc_hess_mean(tn, theta)
#' unlist(EH)
#'
#' # Against the same expectation taken by simulation.
#' set.seed(1)
#' ys <- distrib_rng(tn, 100000, theta)
#' H <- distrib_hessian(gaussian1_distrib(), ys, theta)
#' round(vapply(H, mean, numeric(1)), 3)
#' round(unlist(EH), 3)
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

#' @title Density of a Truncated Distribution
#'
#' @description
#' Evaluates the parent's density divided by the retained mass \eqn{Z}, and
#' zero outside \eqn{[L, U]}. This is one of the shared method bodies:
#' truncation treats the two kinds of parent identically once
#' [trunc_constants()] has resolved the one place they differ, so the body is
#' written once and registered on both classes.
#'
#' @details
#' The division is done on the LOG scale and the outside is set to
#' \eqn{-\infty} there, so a point far in a tail keeps its precision instead of
#' underflowing to zero before the division.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param y A numeric vector of observations. A point outside \eqn{[L, U]}
#'   gives `0`, or `-Inf` on the log scale.
#' @param theta A named list of the parent's parameters.
#' @param log Logical, default `FALSE`. When `TRUE` the log-density is
#'   returned.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector, of the length `y` and `theta` recycle to.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [distrib_pdf.TruncatedContinuousDistrib()] and
#'   [distrib_pdf.TruncatedDiscreteDistrib()], the two registrations, and
#'   [trunc_constants()] for \eqn{Z}.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' distrib_pdf(tn, c(-2, 0, 1, 3), theta)
#'
#' # The parent's density divided by the retained mass.
#' Z <- distributions7:::trunc_constants(tn, theta)$Z
#' dnorm(c(0, 1), 0.3, 1.2) / Z
#'
#' # Outside the interval, zero, and -Inf on the log scale.
#' distrib_pdf(tn, c(-2, 3), theta, log = TRUE)
#'
#' # It integrates to one over the interval, which the parent's does not.
#' integrate(function(y) distrib_pdf(tn, y, theta), -1, 2)$value
trunc_pdf <- function(distrib, y, theta, log = FALSE, ...) {
  Z <- trunc_constants(distrib, theta)$Z
  ld <- distrib_pdf(distrib@parent_distrib, y, theta, log = TRUE) - log(Z)
  outside <- !trunc_inside(distrib, y)
  if (any(outside)) ld[rep_len(outside, length(ld))] <- -Inf
  if (log) ld else exp(ld)
}

#' @title Distribution Function of a Truncated Distribution
#'
#' @description
#' Evaluates \eqn{F_T(q) = (F(q) - F(L^-))/Z}, clamped to \eqn{[0, 1]}. One of
#' the shared method bodies, registered on both truncated classes.
#'
#' @details
#' The clamp is there to make the endpoints exact. It returns `0` at and below
#' \eqn{L} and `1` at and above \eqn{U}, where the unclamped ratio would give a
#' small negative number or one plus a rounding.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical, default `TRUE`. When `FALSE` the upper tail
#'   \eqn{1 - F_T(q)} is returned.
#' @param log.p Logical, default `FALSE`. When `TRUE` the probability is
#'   returned on the log scale.
#'
#' @return A numeric vector of probabilities, of the length `q` and `theta`
#'   recycle to.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [distrib_cdf.TruncatedContinuousDistrib()] and
#'   [distrib_cdf.TruncatedDiscreteDistrib()], the two registrations, and
#'   [trunc_quantile()] for the inverse.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' # Exactly 0 and 1 at the endpoints.
#' distrib_cdf(tn, c(-1, 0.3, 2), theta)
#'
#' # Against the ratio written out.
#' cs <- distributions7:::trunc_constants(tn, theta)
#' (pnorm(0.3, 0.3, 1.2) - cs$Fl) / cs$Z
#'
#' # The two tails sum to one, on either scale.
#' distrib_cdf(tn, 0.3, theta) + distrib_cdf(tn, 0.3, theta, lower.tail = FALSE)
#' exp(distrib_cdf(tn, 0.3, theta, log.p = TRUE))
trunc_cdf <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  cs <- trunc_constants(distrib, theta)
  res <- (distrib_cdf(distrib@parent_distrib, q, theta) - cs$Fl) / cs$Z
  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Quantile Function of a Truncated Distribution
#'
#' @description
#' Inverts [trunc_cdf()] through the PARENT's quantile function, no
#' root-finding of its own being needed: \eqn{F_T(q) = p} holds exactly when
#' \eqn{F(q) = F(L^-) + pZ}. The generalized inverse of a discrete
#' distribution function satisfies the same relation, so the discrete case
#' needs no separate treatment and the body is registered on both classes.
#'
#' @details
#' `p` and `theta` are expanded to a common length before being combined.
#' \eqn{F(L^-)} and \eqn{Z} follow the length of `theta` and `p` follows its
#' own, so multiplying them directly would recycle silently where the two
#' disagree.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param p A numeric vector of probabilities, clamped to \eqn{[0, 1]}.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical, default `TRUE`. When `FALSE`, `p` is read as an
#'   upper-tail probability.
#' @param log.p Logical, default `FALSE`. When `TRUE`, `p` is given on the log
#'   scale.
#'
#' @return A numeric vector of quantiles, of the length `p` and `theta` recycle
#'   to.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [distrib_quantile.TruncatedContinuousDistrib()] and
#'   [distrib_quantile.TruncatedDiscreteDistrib()], the two registrations,
#'   [trunc_cdf()] for the inverse, and [trunc_rng()], which draws through it.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' distrib_quantile(tn, c(0.1, 0.5, 0.9), theta)
#'
#' # It inverts the distribution function exactly.
#' p <- c(0.1, 0.5, 0.9)
#' max(abs(distrib_cdf(tn, distrib_quantile(tn, p, theta), theta) - p))
#'
#' # And every quantile lies inside the interval.
#' range(distrib_quantile(tn, seq(0, 1, by = 0.25), theta))
#'
#' # A discrete parent needs no separate treatment.
#' ztp <- truncated(poisson_distrib(), lower = 1)
#' distrib_quantile(ztp, c(0.1, 0.5, 0.9), list(mu = 2))
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
  # The generalized inverse of a discrete cdf satisfies the same relation, so the
  # discrete case needs no separate treatment.
  distrib_quantile(distrib@parent_distrib, pmin(pmax(cs$Fl + p * cs$Z, 0), 1), th)
}

#' @title Random Generation From a Truncated Distribution
#'
#' @description
#' Draws by inverse transform through [trunc_quantile()],
#' \eqn{Y = Q(F(L^-) + V Z)} with \eqn{V} uniform on \eqn{(0, 1)}. The draw is
#' exact and terminates in one pass however small \eqn{Z} is, where rejection
#' sampling from the parent would need \eqn{1/Z} draws on average. One of the
#' shared bodies, registered on both classes.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param n The number of draws, a single positive integer.
#' @param theta A named list of the parent's parameters. A component varying by
#'   observation must have length `n`.
#'
#' @return A numeric vector of `n` draws, every one inside \eqn{[L, U]}.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [distrib_rng.TruncatedContinuousDistrib()] and
#'   [distrib_rng.TruncatedDiscreteDistrib()], the two registrations, and
#'   [trunc_quantile()], which does the work.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' set.seed(1)
#' y <- distrib_rng(tn, 5000, theta)
#' c(inside = all(y >= -1 & y <= 2), min = min(y), max = max(y))
#'
#' # The sample mean tracks the truncated mean, not the parent's.
#' c(sampled = mean(y), truncated = mean(tn, theta), parent = 0.3)
#'
#' # An interval carrying little mass costs no more than any other.
#' far <- truncated(gaussian1_distrib(), lower = 4, upper = 5)
#' set.seed(2)
#' range(distrib_rng(far, 1000, theta))
trunc_rng <- function(distrib, n, theta) {
  trunc_quantile(distrib, stats::runif(n), theta)
}

#' @title Score of a Truncated Distribution
#'
#' @description
#' Computes the parent's score recentered by its truncated mean,
#' \deqn{\frac{\partial \ell_T}{\partial \theta_i} = s_i(y) - m_i, \qquad
#'   m_i = \mathbb{E}_T[s_i],}
#' the subtraction being the derivative of \eqn{-\log Z}. One of the shared
#' bodies, registered on both truncated classes.
#'
#' @details
#' Truncation adds no parameter, the endpoints being constants, so the returned
#' list has exactly the parent's components. What it adds is the
#' \eqn{\theta}-dependent normalizing constant, and the recentering is that
#' constant's whole contribution at first order. The support does not depend on
#' \eqn{\theta}, which licenses differentiating \eqn{Z} under the integral sign
#' and keeps truncation at fixed points a regular problem.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, applied by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one component per parameter, in
#'   `distrib@params` order.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_score_mean()] for \eqn{m_i}, [trunc_hessian()] for the
#'   second order, and [distrib_gradient.TruncatedContinuousDistrib()] and
#'   [distrib_gradient.TruncatedDiscreteDistrib()], the two registrations.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' y <- c(-0.5, 0, 1)
#' g <- distrib_gradient(tn, y, theta)
#' vapply(g, sum, numeric(1))
#'
#' # Against a numerical derivative of the truncated log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- tn@params
#'   sum(distrib_pdf(tn, y, t2, log = TRUE))
#' }
#' max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#'
#' # It is the parent's score, shifted by one number per parameter.
#' gp <- distrib_gradient(gaussian1_distrib(), y, theta)
#' round(unlist(g) - unlist(gp), 8)
trunc_gradient <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  m <- trunc_score_mean(distrib, theta)
  g <- distrib_gradient(distrib@parent_distrib, y, theta)
  stats::setNames(lapply(distrib@params, function(p) g[[p]] - m[[p]]), distrib@params)
}

#' @title Observed Hessian of a Truncated Distribution
#'
#' @description
#' Computes
#' \deqn{\frac{\partial^2 \ell_T}{\partial\theta_i \partial\theta_j}
#'   = H_{ij}(y) - M_{ij} + m_i m_j,}
#' the subtracted part \eqn{M_{ij} - m_i m_j} being the second derivative of
#' \eqn{\log Z}. One of the shared bodies, registered on both truncated
#' classes.
#'
#' @details
#' Only \eqn{H_{ij}(y)} varies with the data; the correction is a function of
#' \eqn{\theta} and the endpoints alone, so it is computed once and recycled to
#' the length of `y`.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, applied by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length `length(y)`, one component
#'   per unordered pair of parameters, keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_M()] and [trunc_score_mean()] for the two corrections,
#'   [trunc_expected_hessian()] for the expectation, and
#'   [distrib_hessian.TruncatedContinuousDistrib()] and
#'   [distrib_hessian.TruncatedDiscreteDistrib()], the two registrations.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' y <- c(-0.5, 0, 1)
#' H <- distrib_hessian(tn, y, theta)
#' vapply(H, sum, numeric(1))
#'
#' # Against a numerical Hessian of the truncated log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- tn@params
#'   sum(distrib_pdf(tn, y, t2, log = TRUE))
#' }
#' Hn <- numDeriv::hessian(ll, unlist(theta))
#' ref <- vapply(distributions7:::hess_pairs(tn@params),
#'               function(q) Hn[match(q[1], tn@params), match(q[2], tn@params)],
#'               numeric(1))
#' max(abs(vapply(H, sum, numeric(1)) - ref))
#'
#' # The correction is one number per component, not one per observation.
#' Hp <- distrib_hessian(gaussian1_distrib(), y, theta)
#' round(unlist(H) - unlist(Hp), 8)
trunc_hessian <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  n <- length(y)
  m <- trunc_score_mean(distrib, theta)
  M <- trunc_M(distrib, theta)
  h <- distrib_hessian(distrib@parent_distrib, y, theta)
  pairs <- hess_pairs(distrib@params)

  res <- lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    # d_ij log Z = M_ij - m_i m_j
    h[[nm]] - M[[nm]] + m[[pr[1]]] * m[[pr[2]]]
  })
  expand_params(stats::setNames(res, names(pairs))[hess_names(distrib@params)], n)
}

#' @title Expected Hessian of a Truncated Distribution
#'
#' @description
#' Computes
#' \deqn{\mathbb{E}\left[\frac{\partial^2 \ell_T}{\partial\theta_i
#'   \partial\theta_j}\right] = -\mathrm{Cov}_T(s_i, s_j),}
#' the second Bartlett identity under the truncated law: the parent's expected
#' Hessian cancels exactly against the term it contributes to the second
#' derivative of \eqn{\log Z}. One of the shared bodies, registered on both
#' truncated classes.
#'
#' @details
#' This needs one quadrature per component even where the parent has exact cdf
#' derivatives. Those give \eqn{d_{ij} Z}, which is
#' \eqn{\mathbb{E}_T[H_{ij}] + \mathbb{E}_T[s_i s_j]}, and the covariance needs
#' the second term on its own.
#'
#' No component depends on the data, so every returned vector is constant.
#' `approx` and `nsim` are accepted so that the signature matches the
#' generic's, and neither is read: the identity above is exact, and what it
#' rests on is computed by quadrature whatever `approx` says.
#'
#' @param distrib A truncated distribution object, of either class.
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, applied by the
#'   generic before dispatch.
#' @param approx Ignored. Present so that the signature matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length `length(y)`, one component
#'   per unordered pair of parameters, keyed as
#'   [`hess_names(distrib@params)`][hess_names], each vector constant.
#'
#' @section Notation:
#' \eqn{L} and \eqn{U} are the truncation endpoints, both included in the
#' support; \eqn{Z(\theta) = P(L \le Y \le U)} is the retained mass; \eqn{f}
#' and \eqn{F} are the parent's density and distribution function; \eqn{s_i}
#' and \eqn{H_{ij}} are the parent's score and observed Hessian; and
#' \eqn{\mathbb{E}_T} is expectation under the truncated law.
#'
#' @seealso [trunc_score_prod_mean()] and [trunc_score_mean()] for the two
#'   pieces of the covariance, [trunc_hessian()] for the observed matrix, and
#'   [distrib_expected_hessian.TruncatedContinuousDistrib()] and
#'   [distrib_expected_hessian.TruncatedDiscreteDistrib()], the two
#'   registrations.
#'
#' @keywords internal
#'
#' @examples
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' EH <- distrib_expected_hessian(tn, 1, theta)
#' unlist(EH)
#'
#' # Minus the covariance of the parent's score under the truncated law.
#' m <- distributions7:::trunc_score_mean(tn, theta)
#' ES <- distributions7:::trunc_score_prod_mean(tn, theta)
#' c(mu_mu = -(ES$mu_mu - m$mu^2),
#'   sigma_sigma = -(ES$sigma_sigma - m$sigma^2),
#'   mu_sigma = -(ES$mu_sigma - m$mu * m$sigma))
#'
#' # And what averaging the observed Hessian over the truncated law gives.
#' set.seed(1)
#' ys <- distrib_rng(tn, 200000, theta)
#' round(vapply(distrib_hessian(tn, ys, theta), mean, numeric(1)), 3)
#' round(vapply(EH, function(z) z[1], numeric(1)), 3)
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
#' \deqn{f_T(y) = \dfrac{f(y;\theta)}{Z(\theta)}\ \ \ (L \le y \le U), \qquad 0 \text{ otherwise}}
#' with \eqn{Z(\theta) = F(U;\theta) - F(L;\theta)}.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [truncated()]
S7::method(distrib_pdf, TruncatedContinuousDistrib) <- trunc_pdf

#' @title Truncated Probability Mass Function
#' @name distrib_pdf.TruncatedDiscreteDistrib
#' @description
#' \deqn{P_T(Y = y) = \dfrac{f(y;\theta)}{Z(\theta)}\ \ \ (L \le y \le U), \qquad 0 \text{ otherwise}}
#' with \eqn{Z(\theta) = F(U;\theta) - F(L;\theta) + f(L;\theta)}, the mass
#' at the lower endpoint being added back because the endpoint is included.
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param log Logical; if `TRUE`, returns the log-probability.
#' @return A numeric vector of density values.
#' @seealso [truncated()]
S7::method(distrib_pdf, TruncatedDiscreteDistrib) <- trunc_pdf

#' @title Truncated Cumulative Distribution Function (Continuous)
#' @name distrib_cdf.TruncatedContinuousDistrib
#' @description \deqn{F_T(q) = \dfrac{F(q;\theta) - F(L^-;\theta)}{Z(\theta)}} clamped to \eqn{[0,1]}.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned as logs.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [truncated()]
S7::method(distrib_cdf, TruncatedContinuousDistrib) <- trunc_cdf

#' @title Truncated Cumulative Distribution Function (Discrete)
#' @name distrib_cdf.TruncatedDiscreteDistrib
#' @description \deqn{F_T(q) = \dfrac{F(q;\theta) - F(L^-;\theta)}{Z(\theta)}} clamped to \eqn{[0,1]}.
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned as logs.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [truncated()]
S7::method(distrib_cdf, TruncatedDiscreteDistrib) <- trunc_cdf

#' @title Truncated Quantile Function (Continuous)
#' @name distrib_quantile.TruncatedContinuousDistrib
#' @description \deqn{Q_T(p) = Q\!\left(F(L^-;\theta) + p\,Z(\theta)\right)}
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if `TRUE`, probabilities are given as logs.
#' @return A numeric vector of quantiles.
#' @seealso [truncated()]
S7::method(distrib_quantile, TruncatedContinuousDistrib) <- trunc_quantile

#' @title Truncated Quantile Function (Discrete)
#' @name distrib_quantile.TruncatedDiscreteDistrib
#' @description \deqn{Q_T(p) = Q\!\left(F(L^-;\theta) + p\,Z(\theta)\right)}
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if `TRUE`, probabilities are given as logs.
#' @return A numeric vector of quantiles.
#' @seealso [truncated()]
S7::method(distrib_quantile, TruncatedDiscreteDistrib) <- trunc_quantile

#' @title Truncated Random Number Generator (Continuous)
#' @name distrib_rng.TruncatedContinuousDistrib
#' @description
#' Inverse transform sampling on the parent: \eqn{Y = Q(F(L^-) + V Z)} with
#' \eqn{V \sim \mathrm{Unif}(0,1)}. Exact, and unlike rejection sampling it always
#' terminates in one pass however small \eqn{Z} is.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector of random draws.
#' @seealso [truncated()]
S7::method(distrib_rng, TruncatedContinuousDistrib) <- trunc_rng

#' @title Truncated Random Number Generator (Discrete)
#' @name distrib_rng.TruncatedDiscreteDistrib
#' @description Inverse transform sampling on the parent, exact for a discrete cdf.
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector of random draws.
#' @seealso [truncated()]
S7::method(distrib_rng, TruncatedDiscreteDistrib) <- trunc_rng

#' @title Truncated Analytical Gradient (Continuous)
#' @name distrib_gradient.TruncatedContinuousDistrib
#' @description
#' \deqn{\dfrac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad m_i = \mathbb{E}_T[s_i]}
#' the parent's score recentered at its mean over the truncated support.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of first derivatives.
#' @seealso [truncated()]
S7::method(distrib_gradient, TruncatedContinuousDistrib) <- trunc_gradient

#' @title Truncated Analytical Gradient (Discrete)
#' @name distrib_gradient.TruncatedDiscreteDistrib
#' @description
#' \deqn{\dfrac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad m_i = \mathbb{E}_T[s_i]}
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of first derivatives.
#' @seealso [truncated()]
S7::method(distrib_gradient, TruncatedDiscreteDistrib) <- trunc_gradient

#' @title Truncated Analytical Observed Hessian (Continuous)
#' @name distrib_hessian.TruncatedContinuousDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j} = H_{ij}(y) - M_{ij} + m_i m_j}
#' with \eqn{M_{ij} = \mathbb{E}_T[H_{ij} + s_i s_j]}.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of second derivatives.
#' @seealso [truncated()]
S7::method(distrib_hessian, TruncatedContinuousDistrib) <- trunc_hessian

#' @title Truncated Analytical Observed Hessian (Discrete)
#' @name distrib_hessian.TruncatedDiscreteDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j} = H_{ij}(y) - M_{ij} + m_i m_j}
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of second derivatives.
#' @seealso [truncated()]
S7::method(distrib_hessian, TruncatedDiscreteDistrib) <- trunc_hessian

#' @title Truncated Analytical Expected Hessian (Continuous)
#' @name distrib_expected_hessian.TruncatedContinuousDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j}\right]
#'   = -\operatorname{Cov}_T(s_i, s_j)}
#' the covariance of the parent's score under the truncated distribution.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [truncated()]
S7::method(distrib_expected_hessian, TruncatedContinuousDistrib) <- trunc_expected_hessian

#' @title Truncated Analytical Expected Hessian (Discrete)
#' @name distrib_expected_hessian.TruncatedDiscreteDistrib
#' @description
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j}\right]
#'   = -\operatorname{Cov}_T(s_i, s_j)}
#' @param distrib A `TruncatedDiscreteDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [truncated()]
S7::method(distrib_expected_hessian, TruncatedDiscreteDistrib) <- trunc_expected_hessian

#' @title Atoms of a Truncated Continuous Distribution
#' @name distrib_atoms.TruncatedContinuousDistrib
#' @description
#' Truncation preserves the parent's atoms that survive it, rescaled by
#' \eqn{1/Z}. This matters only when the parent is itself mixed, as
#' [zero_adjusted()] of a continuous distribution is.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param theta A named list of the parent's parameters.
#' @return A list with components `y` and `p`.
#' @seealso [truncated()], [distrib_atoms()]
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
#' \eqn{[L, U]}, which is correct unless the parent carries point masses ---
#' as it does when it is a [zero_adjusted()] continuous distribution.
#' Those masses are added explicitly, exactly as in
#' [`the untruncated case()`][expectation.ZeroAdjustedContinuousDistrib].
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param f A function `f(y, theta, ...)`.
#' @param theta A named list of the parent's parameters.
#' @param ... Additional arguments passed to `f`.
#' @keywords internal
#' @return A numeric scalar, the expectation of `f` under the distribution.
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
#' \eqn{Z} does not depend on \eqn{y}, so inside \eqn{[L, U]} the response
#' derivative is the parent's. Outside, the log-density is \eqn{-\infty} and no
#' derivative exists, so `NaN` is returned.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector.
#' @seealso [truncated()]
S7::method(distrib_grad_y, TruncatedContinuousDistrib) <- function(distrib, y, theta) {
  trunc_y_deriv(distrib, y, theta, distrib_grad_y)
}

#' @title Truncated Continuous Response Hessian
#' @name distrib_hess_y.TruncatedContinuousDistrib
#' @description As [`distrib_grad_y()`][distrib_grad_y.TruncatedContinuousDistrib], at second order.
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector.
#' @seealso [truncated()]
S7::method(distrib_hess_y, TruncatedContinuousDistrib) <- function(distrib, y, theta) {
  trunc_y_deriv(distrib, y, theta, distrib_hess_y)
}

# Internal: a response derivative of the parent, restricted to the truncated
# support. Defined after the methods so that its comment block cannot capture
# their documentation.
#' Response Derivative of a Truncated Distribution
#'
#' @description
#' Evaluates a response derivative of the parent inside the interval.
#'
#' @details
#' The normalizing constant does not depend on \eqn{y}, so inside the support the
#' parent's derivative is exact and nothing needs correcting.
#'
#' @param distrib A truncated distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param fun The parent's response-derivative function.
#'
#' @return A numeric vector.
#'
#' @keywords internal
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
#' Validate the Truncation Endpoints
#'
#' @description
#' Checks that the interval is well formed and that truncating this parent there
#' leaves an identifiable model.
#'
#' @details
#' One case is rejected outright: truncating zero away from a zero wrapper. The
#' \eqn{(1-\zeta)} factor then cancels between the numerator and \eqn{Z}, so
#' \eqn{\zeta} leaves the likelihood entirely and its score is identically zero
#' -- the same defect as stacking the two zero wrappers, arriving by a different
#' route. Truncating a zero wrapper anywhere else is fine, and the atom is
#' carried through.
#'
#' @param distrib The parent distribution.
#' @param lower The lower endpoint, or `NULL`.
#' @param upper The upper endpoint, or `NULL`.
#'
#' @return Invisibly `NULL`; raises an error on a bad interval.
#'
#' @seealso [truncated()]
#' @keywords internal
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
        "'%s' = %s is not a point of the support. A discrete distribution is supported on\n",
        "  the integers, so a non-integer truncation point is ambiguous: use %s or %s."
      ), arg, format(v), format(floor(v)), format(ceiling(v))), call. = FALSE)
    }
  }

  if (!is.null(lower) && !is.null(upper) && lower >= upper) {
    stop(sprintf("'lower' (%s) must be strictly less than 'upper' (%s).",
                 format(lower), format(upper)), call. = FALSE)
  }

  # A truncation point outside the support removes nothing. This is the
  # truncated(gamma2_distrib(), lower = -2) case: the result would be the parent,
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
#' Restricts an existing distribution to \eqn{[L, U]} and renormalizes it, so
#' that all the probability mass the parent placed outside the interval is
#' redistributed inside it. Either endpoint may be omitted, giving one-sided
#' truncation; at least one must be given.
#'
#' Works on discrete and continuous parents alike. **Both endpoints are
#' included**, which for a discrete parent is the difference between
#' `truncated(poisson_distrib(), lower = 1)` --- the zero-truncated Poisson,
#' supported on \eqn{\{1, 2, \dots\}} --- and truncating above 1.
#'
#' @param distrib An object inheriting from `discrete_distrib` or
#'   `continuous_distrib`.
#' @param lower,upper The truncation points. Each may be `NULL` (no
#'   truncation on that side), and at least one must be supplied. For a discrete
#'   parent both must be whole numbers.
#'
#' @details
#' Write \eqn{Z(\theta) = P(L \le Y \le U)} for the retained mass. Then
#' \deqn{f_T(y;\theta) = \frac{f(y;\theta)}{Z(\theta)}\ \ (L \le y \le U),
#' \qquad F_T(q) = \frac{F(q;\theta) - F(L^-;\theta)}{Z(\theta)},
#' \qquad Q_T(p) = Q\!\left(F(L^-;\theta) + pZ(\theta)\right),}
#' with \eqn{F(L^-) = F(L)} for a continuous parent and
#' \eqn{F(L) - f(L)} for a discrete one, since the lower endpoint is kept.
#'
#' **Truncation adds no parameter.** The endpoints are known constants, like
#' a binomial's `size`, so the truncated distribution has exactly the
#' parent's parameters, domains and links. What it does add is a
#' \eqn{\theta}-dependent normalizing constant, and that constant contributes to
#' every derivative of the log-likelihood \eqn{\ell_T = \ell - \log Z}. Writing
#' \deqn{m_i = \mathbb{E}_T[s_i], \qquad M_{ij} = \mathbb{E}_T[H_{ij} + s_is_j],}
#' where the expectations are taken under the *truncated* distribution,
#' \deqn{\frac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad
#' \frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j} = H_{ij}(y) - M_{ij} + m_im_j,}
#' \deqn{\mathbb{E}\left[\frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j}\right]
#' = -\operatorname{Cov}_T(s_i, s_j).}
#' The score is simply the parent's score *recentered* at its truncated mean,
#' and the information is the covariance of that score --- which is the second
#' Bartlett identity for the truncated model, and is used as a consistency check
#' rather than derived separately.
#'
#' **What this costs.** \eqn{m_i} and \eqn{M_{ij}} have no closed form for a
#' general parent, and are obtained by quadrature (continuous) or summation
#' (discrete) through [expectation()]. Derivatives of a truncated
#' distribution are therefore substantially more expensive than the parent's, and
#' third and fourth derivatives fall back to finite differences of the analytical
#' Hessian.
#'
#' **What the constructor rejects.**
#'
#' - Both endpoints `NULL`: nothing to do, and silently returning the
#'   parent would hide the mistake.
#' - `lower >= upper`.
#' - A truncation point that removes no mass, such as
#'   `truncated(gamma2_distrib(), lower = -2)`: the Gamma is supported on
#'   \eqn{(0,\infty)}, so the result would be the Gamma itself.
#' - A non-integer endpoint for a discrete parent, which is ambiguous.
#' - A discrete truncation leaving too few support points to identify the
#'   parameters: \eqn{k} points carry \eqn{k-1} free probabilities, so
#'   `n_params + 1` points are needed.
#' - A parent that models a probability of zero ---
#'   [zero_inflated()] or [zero_adjusted()] --- when the
#'   truncation removes \eqn{0} from the support. Truncating zero away cancels
#'   that parameter out of the likelihood entirely, leaving an identically zero
#'   score. Truncating elsewhere, as in
#'   `truncated(zero_adjusted(gamma2_distrib()), upper = 5)`, is fine and the
#'   point mass is carried through [distrib_atoms()].
#'
#' Truncating an already truncated distribution is allowed and is collapsed into a
#' single object over the intersection of the two intervals, rather than nested.
#'
#' @return An S7 object of class `TruncatedDiscreteDistrib` or
#'   `TruncatedContinuousDistrib`.
#'
#' @examples
#' # The zero-truncated Poisson
#' ztp <- truncated(poisson_distrib(), lower = 1)
#' distrib_pdf(ztp, 0:4, list(mu = 2))
#' dpois(1:4, 2) / (1 - dpois(0, 2))
#'
#' # A Gaussian restricted to an interval
#' tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
#' mean(tn, list(mu = 0, sigma = 1))
#'
#' # A truncation point that removes nothing is rejected
#' try(truncated(gamma2_distrib(), lower = -2))
#'
#' @seealso [zero_inflated()], [zero_adjusted()],
#'   [transformation()], [check_distrib()]
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
