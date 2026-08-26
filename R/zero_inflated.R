#' @include distrib.R generics.R utility_functions.R
NULL

#' @title S7 Class for Zero-Inflated Distributions
#' @name ZeroInflatedDistrib
#'
#' @description
#' The S7 class of the zero-inflated version of a discrete distribution, with
#' mass function
#' \deqn{P(Y = y) = \zeta\,\mathbb{I}(y = 0) + (1 - \zeta)\, f(y; \theta).}
#' It inherits from `discrete_distrib` and carries the parent's parameters
#' followed by `zi`, which is the probability \eqn{\zeta} of a structural zero
#' and rides a link of its own.
#'
#' Build one with [zero_inflated()], which checks that the parent is discrete,
#' is not already a zero wrapper, and has enough support points for the extra
#' parameter to be identified. This page documents the raw S7 constructor,
#' which checks none of that.
#'
#' @param parent_distrib The wrapped `discrete_distrib` object.
#' @inheritParams distrib
#'
#' @return An S7 object of class `ZeroInflatedDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. It carries `parent_distrib` beside
#'   the parent's `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params` and
#'   `params_smooth`. For an object built by [zero_inflated()] the parameters
#'   are the parent's followed by `zi`, whose bound is \eqn{(0, 1)} and whose
#'   interpretation is `"prob. of structural zero"`.
#'
#' @seealso [zero_inflated()] to build one, [ZeroAdjustedDiscreteDistrib] for
#'   the wrapper that REPLACES the mass at zero where this one adds to it, and
#'   [zero_adjusted()], the constructor to reach for with a continuous
#'   parent.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.ZeroInflatedDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ZeroInflatedDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ZeroInflatedDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ZeroInflatedDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ZeroInflatedDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ZeroInflatedDistrib],
#'   [`distrib_rng()`][distrib_rng.ZeroInflatedDistrib]
#'
#' The third and fourth derivatives come from the shared wrapper machinery in
#' `wrapper_derivatives.R`; everything else is inherited from
#' [discrete_distrib()].
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # The parent's parameters, then zi last with its own link and bound.
#' d@params
#' d@params_bounds$zi
#' vapply(d@link_params, function(l) l@link_name, character(1))
#' d@params_interpretation
#'
#' # Inflation can only ADD zeros: the mass at zero exceeds the parent's.
#' c(inflated = distrib_pdf(d, 0, theta), parent = dpois(0, 3))
ZeroInflatedDistrib <- S7::new_class("ZeroInflatedDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib
  )
)

#' @title Split a Wrapper's Parameter From Its Parent's
#'
#' @description
#' Separates a full `theta` into the parent distribution's parameters and the
#' single mixing parameter the wrapper adds. Both zero wrappers append their
#' parameter LAST, so the split is positional and needs no name matching, which
#' keeps it correct for a parent whose own parameter happens to be called `zi`
#' or `pi`.
#'
#' @param distrib A `ZeroInflatedDistrib` or a zero-adjusted object, whose
#'   `parent_distrib` supplies the names of the first block.
#' @param theta A named list of parameters, already aligned, the parent's
#'   followed by the wrapper's.
#'
#' @return A named list with `orig`, the parent's parameters as a named list,
#'   and `mix`, the wrapper's parameter as a numeric vector.
#'
#' @seealso [zero_inflated()] and [zero_adjusted()], whose methods all begin
#'   with this.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#' p <- distributions7:::split_mix_theta(
#'   d, distributions7:::align_theta(d, theta))
#' str(p)
#'
#' # The first block is exactly what the parent expects.
#' all.equal(distrib_pdf(d@parent_distrib, 0:3, p$orig), dpois(0:3, 3))
#'
#' @keywords internal
split_mix_theta <- function(distrib, theta) {
  n <- distrib@n_params
  list(orig = theta[seq_len(n - 1L)], mix = theta[[n]])
}

# ---------------------------------------------------------------------------
# Validation shared by zero_inflated() and zero_adjusted().
#
# Both wrappers add one parameter to a distribution and both are easy to apply
# where the result is not a model at all. Neither failure is visible at run
# time: the pmf still sums to one, check_distrib() still passes, and the fit
# still converges -- to an arbitrary point of a flat ridge. The constructor is
# the only place where they can be caught, so they are caught there.
# ---------------------------------------------------------------------------

#' @title Number of Points in a Discrete Support
#'
#' @description
#' Returns how many points the distribution's support contains, read off its
#' `bounds`, and `Inf` when the support is unbounded above. It is what
#' [check_support_is_rich_enough()] counts against the parameters a zero
#' wrapper would spend.
#'
#' @param distrib A `distrib` object, whose `bounds` are read.
#'
#' @return A single number: the count of integer points between the bounds
#'   inclusive, or `Inf`.
#'
#' @seealso [check_support_is_rich_enough()], the consumer, and
#'   [zero_inflated()] for the counting rule.
#'
#' @examples
#' # A count family has an unbounded support.
#' distributions7:::n_support_points(poisson_distrib())
#'
#' # A binomial has size + 1 points, and a Bernoulli has two.
#' c(binomial5 = distributions7:::n_support_points(binomial_distrib(size = 5)),
#'   bernoulli = distributions7:::n_support_points(bernoulli_distrib()))
#'
#' @keywords internal
n_support_points <- function(distrib) {
  b <- distrib@bounds
  if (!all(is.finite(b))) Inf else b[2] - b[1] + 1
}

#' @title Does This Distribution Already Model a Probability of Zero
#'
#' @description
#' Answers whether the object is already one of the two zero wrappers, which is
#' how each of them refuses to be applied on top of the other. Stacking them
#' leaves only their combination identified, so the check is what turns a model
#' with a flat direction into an error at construction.
#'
#' @param distrib A `distrib` object.
#'
#' @return `TRUE` for a zero-inflated or a zero-adjusted object of either kind,
#'   `FALSE` otherwise.
#'
#' @seealso [check_not_stacked()], which reports the refusal, and
#'   [zero_inflated()] and [zero_adjusted()] for what stacking would cost.
#'
#' @examples
#' c(plain = distributions7:::is_zero_wrapper(poisson_distrib()),
#'   inflated = distributions7:::is_zero_wrapper(
#'     zero_inflated(poisson_distrib())),
#'   adjusted = distributions7:::is_zero_wrapper(
#'     zero_adjusted(poisson_distrib())))
#'
#' @keywords internal
is_zero_wrapper <- function(distrib) {
  S7::S7_inherits(distrib, ZeroInflatedDistrib) ||
    S7::S7_inherits(distrib, ZeroAdjustedDiscreteDistrib) ||
    S7::S7_inherits(distrib, ZeroAdjustedContinuousDistrib)
}

#' @title Reject the Composition of Two Zero Wrappers
#'
#' @description
#' Signals an error when the parent already models the probability of a zero,
#' which is how both wrappers refuse to stack. The models this rejects are
#' well defined and inestimable, so the constructor is the only place the
#' problem can be caught: nothing at run time reports it.
#'
#' @details
#' Two zero parameters cannot both be identified. Zero-adjusting a
#' zero-inflated parent truncates the zero away, which cancels the
#' \eqn{(1-\zeta)} factor between the numerator and the normalizing constant,
#' so \eqn{\zeta} leaves the likelihood entirely and its score is IDENTICALLY
#' zero. The other order keeps only the combination \eqn{\zeta + (1-\zeta)\pi}.
#' Either way an optimizer wanders along a flat ridge, the mass function sums
#' to one throughout, and [check_distrib()] passes.
#'
#' @param distrib The candidate parent, a `distrib` object.
#' @param fun The name of the calling constructor, a single string, quoted in
#'   the message.
#' @param param The name of the parameter it would add, a single string. Not
#'   currently placed in the message; it is carried for the caller's own
#'   record.
#'
#' @return `NULL`, invisibly, when the parent is not a zero wrapper. Otherwise
#'   it signals an error naming the parent.
#'
#' @seealso [is_zero_wrapper()] for the test, and [zero_inflated()] and
#'   [zero_adjusted()], the two callers.
#'
#' @examples
#' # A plain parent passes silently.
#' distributions7:::check_not_stacked(poisson_distrib(),
#'                                    "zero_inflated()", "zi")
#'
#' # A wrapped one does not, and the two orders are refused alike.
#' try(zero_inflated(zero_adjusted(poisson_distrib())))
#' try(zero_adjusted(zero_inflated(poisson_distrib())))
#'
#' @keywords internal
check_not_stacked <- function(distrib, fun, param) {
  if (is_zero_wrapper(distrib)) {
    stop(sprintf(
      paste0(
        "%s() cannot wrap '%s', which already models the probability of a zero.\n",
        "  Stacking the two leaves only their combination identified: the second\n",
        "  parameter has an identically zero score, and any optimizer will wander\n",
        "  along that ridge. Apply exactly one zero wrapper to a plain distribution."
      ),
      fun, distrib@distrib_name
    ), call. = FALSE)
  }
  if (param %in% distrib@params) {
    stop(sprintf(
      "The parent distribution already has a parameter named '%s'.", param
    ), call. = FALSE)
  }
  invisible(NULL)
}

#' @title Reject a Model With More Parameters Than the Support Can Distinguish
#'
#' @description
#' Signals an error when the parent's support has too few points for the zero
#' wrapper's extra parameter to be identified. A discrete distribution on
#' \eqn{k} points has \eqn{k - 1} free probabilities, and a wrapper spends
#' `n_params + 1` of them, so \eqn{k \ge n_{params} + 2} is necessary. The rule
#' is the SAME for both wrappers.
#'
#' @details
#' What it rules out is exactly the Bernoulli and `binomial_distrib(size = 1)`.
#' Zero-adjusting a Bernoulli leaves the truncated part on the single point
#' \eqn{\{1\}} and `mu` disappears: the mass function is literally the same at
#' `mu = 0.2` and at `mu = 0.9`. None of this is visible at run time, the mass
#' summing to one and [check_distrib()] passing, so the constructor is the only
#' place it can be caught.
#'
#' A large support is NECESSARY without being sufficient. With a mean large
#' enough that the parent puts almost no mass at zero, the extra parameter is
#' weakly identified whatever the support size, which is a question about the
#' data and one this check cannot ask.
#'
#' @param distrib The candidate parent, a `distrib` object.
#' @param fun The name of the calling constructor, a single string, quoted in
#'   the message.
#'
#' @return `NULL`, invisibly, when the support is large enough. Otherwise it
#'   signals an error giving the point count, the parameter count and the
#'   count required.
#'
#' @seealso [n_support_points()] for the count, and [zero_inflated()] and
#'   [zero_adjusted()], the two callers.
#'
#' @examples
#' # A count family has room.
#' distributions7:::check_support_is_rich_enough(poisson_distrib(),
#'                                               "zero_inflated()")
#'
#' # A Bernoulli does not, and neither does a one-trial binomial.
#' try(zero_inflated(bernoulli_distrib()))
#' try(zero_inflated(binomial_distrib(size = 1)))
#'
#' # Five trials is enough: six points against two parameters.
#' zero_inflated(binomial_distrib(size = 5))@params
#'
#' @keywords internal
check_support_is_rich_enough <- function(distrib, fun) {
  k <- n_support_points(distrib)
  needed <- distrib@n_params + 2
  if (k < needed) {
    plural <- function(n, one, many) if (n == 1) one else many
    stop(sprintf(
      paste0(
        "%s() cannot wrap '%s': its support has %g points, so the family has %g free\n",
        "  %s, while the wrapped distribution would have %d %s. They are not\n",
        "  identified -- different parameter values give exactly the same distribution.\n",
        "  A support of at least %d points is required."
      ),
      fun, distrib@distrib_name, k, k - 1,
      plural(k - 1, "probability", "probabilities"),
      distrib@n_params + 1,
      plural(distrib@n_params + 1, "parameter", "parameters"),
      needed
    ), call. = FALSE)
  }
  invisible(NULL)
}

#' @title Zero-Inflated Probability Mass Function
#' @name distrib_pdf.ZeroInflatedDistrib
#'
#' @description
#' Computes the zero-inflated mass function
#' \deqn{P(Y = y) = \zeta\,\mathbb{I}(y = 0) + (1-\zeta)\, f(y; \theta),}
#' with \eqn{f} the parent's mass function. Inflation can only ADD zeros: the
#' mass at zero is \eqn{\zeta + (1-\zeta)f(0)}, which exceeds \eqn{f(0)} for
#' every \eqn{\zeta > 0}, and no single observed zero can be attributed to one
#' mechanism or the other.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param y A numeric vector of observations. A point off the parent's support
#'   returns whatever the parent returns there, scaled by \eqn{1-\zeta}.
#' @param theta A named list with the parent's parameters followed by `zi`,
#'   each a numeric vector of length 1 or of the length of `y`. `zi` must lie
#'   strictly inside \eqn{(0, 1)}.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned. The
#'   logarithm is taken of the mixture, not inside the parent, so it underflows
#'   where the mixture does. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_cdf.ZeroInflatedDistrib()] for the distribution function,
#'   [distrib_gradient.ZeroInflatedDistrib()] for the score,
#'   [distrib_pdf.ZeroAdjustedDiscreteDistrib()] for the wrapper that replaces
#'   the mass at zero, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#'
#' distrib_pdf(d, 0:5, theta)
#'
#' # Which is the mixture written out.
#' y <- 0:5
#' all.equal(distrib_pdf(d, y, theta), 0.25 * (y == 0) + 0.75 * dpois(y, 3))
#'
#' # The mass at zero exceeds the parent's, inflation adding to it.
#' c(inflated = distrib_pdf(d, 0, theta), parent = dpois(0, 3))
#'
#' # And the whole mass function still sums to one.
#' sum(distrib_pdf(d, 0:200, theta))
S7::method(distrib_pdf, ZeroInflatedDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  pars <- split_mix_theta(distrib, theta)
  zi <- pars$mix

  res <- (1 - zi) * distrib_pdf(distrib@parent_distrib, y, pars$orig, log = FALSE) +
    zi * (y == 0)

  if (log) log(res) else res
}

#' @title Zero-Inflated Cumulative Distribution Function
#' @name distrib_cdf.ZeroInflatedDistrib
#'
#' @description
#' Computes
#' \deqn{F_{ZI}(q) = (1-\zeta) F(q; \theta) + \zeta\,\mathbb{I}(q \ge 0)}
#' from the parent's own distribution function, exactly and with no summation.
#' The result is clamped to \eqn{[0, 1]} against rounding.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param q A numeric vector of quantiles. Values below zero give `0`.
#' @param theta A named list with the parent's parameters followed by `zi`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   formed as \eqn{1 - F} and so subject to that subtraction's cancellation
#'   far into the upper tail.
#' @param log.p Logical of length 1. When `TRUE` the logarithm is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities, in \eqn{[0, 1]}.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_pdf.ZeroInflatedDistrib()] for the mass function,
#'   [distrib_quantile.ZeroInflatedDistrib()], which inverts this, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#'
#' q <- c(0, 2, 5)
#' distrib_cdf(d, q, theta)
#'
#' # Which is the parent's, shrunk and shifted.
#' all.equal(distrib_cdf(d, q, theta), 0.75 * ppois(q, 3) + 0.25)
#'
#' # It agrees with the mass function summed, as it must on a lattice.
#' c(cdf = distrib_cdf(d, 5, theta), summed = sum(distrib_pdf(d, 0:5, theta)))
#'
#' # Both tails and the logarithm.
#' distrib_cdf(d, 2, theta, lower.tail = FALSE)
#' distrib_cdf(d, 2, theta, log.p = TRUE)
S7::method(distrib_cdf, ZeroInflatedDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  pars <- split_mix_theta(distrib, theta)
  zi <- pars$mix

  res <- (1 - zi) * distrib_cdf(distrib@parent_distrib, q, pars$orig) + zi * (q >= 0)
  res <- pmin(pmax(res, 0), 1)

  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Zero-Inflated Quantile Function
#' @name distrib_quantile.ZeroInflatedDistrib
#'
#' @description
#' Inverts the mixture distribution function. The quantile is `0` for
#' \eqn{p \le \zeta + (1-\zeta)F(0; \theta)}, which is the whole of the
#' inflated mass at zero, and otherwise the parent's quantile at the rescaled
#' probability \eqn{(p - \zeta)/(1 - \zeta)}.
#'
#' @details
#' The result is a LATTICE quantile and overshoots, as any discrete quantile
#' does: `distrib_cdf(d, distrib_quantile(d, p, theta), theta)` returns the
#' smallest attainable probability at or above `p`, not `p` itself. The
#' inflated atom makes the first step larger than the parent's, so a
#' probability below \eqn{L_0} maps to zero however small it is.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param p A numeric vector of probabilities, clamped to \eqn{[0, 1]} after
#'   the `log.p` and `lower.tail` transformations are applied.
#' @param theta A named list with the parent's parameters followed by `zi`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles, of the recycled length of `p` and
#'   `theta`.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_cdf.ZeroInflatedDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#'
#' distrib_quantile(d, c(0.1, 0.3, 0.9), theta)
#'
#' # Everything below the inflated mass at zero maps to zero.
#' c(mass_at_zero = distrib_pdf(d, 0, theta),
#'   q_below = distrib_quantile(d, 0.2, theta))
#'
#' # The round trip overshoots, as on any lattice: asked for 0.3, the
#' # attainable probability at the returned point is higher.
#' p <- c(0.1, 0.3, 0.9)
#' rbind(asked = p,
#'       reached = distrib_cdf(d, distrib_quantile(d, p, theta), theta))
S7::method(distrib_quantile, ZeroInflatedDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  p <- pmin(pmax(p, 0), 1)

  # Expand everything to a common length for consistent subsetting
  all_params <- expand_params(c(list(.p = p), theta))
  p <- all_params$.p
  th_orig <- all_params[names(pars$orig)]
  zi <- all_params[[distrib@params[distrib@n_params]]]

  F0 <- distrib_cdf(parent, 0, th_orig)
  prob_at_zero <- zi + (1 - zi) * F0

  q_vals <- numeric(length(p))
  idx <- (p > prob_at_zero)

  if (any(idx)) {
    p_trans <- (p[idx] - zi[idx]) / (1 - zi[idx])
    p_trans <- pmin(pmax(p_trans, 0), 1)
    th_sub <- lapply(th_orig, function(x) if (length(x) > 1) x[idx] else x)
    q_vals[idx] <- distrib_quantile(parent, p_trans, th_sub)
  }

  q_vals
}

#' @title Zero-Inflated Random Number Generator
#' @name distrib_rng.ZeroInflatedDistrib
#'
#' @description
#' Draws `n` values from the parent, then replaces a Bernoulli(\eqn{\zeta})
#' fraction of them with structural zeros. It consumes the parent's draws
#' followed by `n` uniforms, so the stream is reproducible under
#' [base::set.seed()] but is not the parent's stream with a filter: the
#' uniforms come after.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param n The number of draws, a single non-negative whole number.
#' @param theta A named list with the parent's parameters followed by `zi`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length `n`.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_pdf.ZeroInflatedDistrib()] for the mass function these
#'   are drawn from, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#'
#' set.seed(1)
#' distrib_rng(d, 10, theta)
#'
#' # A large sample reproduces the inflated mass at zero, which is well above
#' # the parent's.
#' set.seed(1)
#' big <- distrib_rng(d, 20000, theta)
#' c(sampled = mean(big == 0), exact = distrib_pdf(d, 0, theta),
#'   parent = dpois(0, 3))
S7::method(distrib_rng, ZeroInflatedDistrib) <- function(distrib, n, theta, ...) {
  pars <- split_mix_theta(distrib, theta)
  y <- distrib_rng(distrib@parent_distrib, n, pars$orig)
  y[stats::runif(n) < pars$mix] <- 0
  y
}

#' @title Zero-Inflated Score
#' @name distrib_gradient.ZeroInflatedDistrib
#'
#' @description
#' Computes the first derivatives of the zero-inflated log-mass in closed form.
#' For the parent's parameters the score is the parent's own, weighted by the
#' posterior probability that an observed zero came from the parent,
#' \deqn{\frac{\partial \ell}{\partial \theta_i} = w\, s_i(y), \qquad
#'   w = \begin{cases} \dfrac{(1-\zeta)f(0)}{L_0} & y = 0 \\ 1 & y > 0,
#'   \end{cases}}
#' and for the inflation parameter
#' \deqn{\frac{\partial \ell}{\partial \zeta}
#'   = \mathbb{I}(y = 0)\frac{1 - f(0)}{L_0}
#'   - \mathbb{I}(y > 0)\frac{1}{1 - \zeta}.}
#'
#' @details
#' A positive observation carries no information about \eqn{\zeta} beyond the
#' \eqn{(1-\zeta)} factor, which is why its \eqn{\zeta} component is the same
#' number at every such point. A zero carries all of it, and its weight
#' \eqn{w} is what shares the credit between the two mechanisms.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `zi`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. `zi` rides a logit by default, so the two differ.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, in
#'   `distrib@params` order.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score, \eqn{H} its observed
#' Hessian and \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_hessian.ZeroInflatedDistrib()] for the second order,
#'   [distrib_expected_hessian.ZeroInflatedDistrib()] for the information, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#' set.seed(2)
#' y <- distrib_rng(d, 200, theta)
#'
#' g <- distrib_gradient(d, y, theta)
#' vapply(g, sum, numeric(1))
#'
#' # Against a numerical derivative of the log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#'
#' # Every positive observation gives the same zi component, -1 / (1 - zi).
#' c(unique(round(g$zi[y > 0], 12)), theory = -1 / 0.75)
#'
#' # And the parent's component is zero at a zero only in the limit: it is
#' # the parent's score there, weighted down by w.
#' c(weighted = g$mu[y == 0][1], unweighted = distrib_gradient(
#'     poisson_distrib(), 0, list(mu = 3))$mu)
S7::method(distrib_gradient, ZeroInflatedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  zi <- pars$mix

  f0 <- distrib_pdf(parent, 0, pars$orig)
  l0 <- zi + (1 - zi) * f0

  grad_orig <- distrib_gradient(parent, y, pars$orig)
  w <- ifelse(y == 0, ((1 - zi) * f0) / l0, 1)

  res <- lapply(grad_orig, function(g) w * g)
  res[[distrib@params[distrib@n_params]]] <- ifelse(y == 0, (1 - f0) / l0, -1 / (1 - zi))
  res
}

#' @title Zero-Inflated Observed Hessian
#' @name distrib_hessian.ZeroInflatedDistrib
#'
#' @description
#' Computes the second derivatives of the zero-inflated log-mass in closed
#' form. At a POSITIVE observation the parent block is the parent's own
#' Hessian, the mixed block is zero and the \eqn{\zeta} block is
#' \eqn{-1/(1-\zeta)^2}. At a ZERO every block picks up the mixture's
#' correction:
#' \deqn{\ell^{(\theta_i\theta_j)} = w\,H_{ij}(0) + w(1-w)\,s_i(0)s_j(0),
#'   \qquad
#'   \ell^{(\theta_i\zeta)} = -\frac{f(0)\,s_i(0)}{L_0^2},
#'   \qquad
#'   \ell^{(\zeta\zeta)} = -\frac{(1 - f(0))^2}{L_0^2},}
#' with \eqn{w = (1-\zeta)f(0)/L_0} the posterior weight of the parent
#' component.
#'
#' @details
#' The parent block is the ordinary two-component mixture curvature: the
#' weighted Hessian plus the variance of the score across the two components,
#' \eqn{w(1-w)} being that variance's weight when one component's score is
#' zero.
#'
#' The mixed block collapses. Written out it is
#' \eqn{-f'(0)/L_0 - (1-\zeta)f'(0)(1 - f(0))/L_0^2}, and the bracket is
#' exactly one, so the whole thing is \eqn{-f'(0)/L_0^2}.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `zi`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score, \eqn{H} its observed
#' Hessian and \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_gradient.ZeroInflatedDistrib()] for the first order,
#'   [distrib_expected_hessian.ZeroInflatedDistrib()] for the expectation, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#' set.seed(2)
#' y <- distrib_rng(d, 200, theta)
#'
#' H <- distrib_hessian(d, y, theta)
#' vapply(H, sum, numeric(1))
#'
#' # Against a numerical Hessian of the log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' Hn <- numDeriv::hessian(ll, unlist(theta))
#' ref <- vapply(distributions7:::hess_pairs(d@params),
#'               function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
#'               numeric(1))
#' max(abs(vapply(H, sum, numeric(1)) - ref))
#'
#' # At a positive observation the mixed block is exactly zero and the zi
#' # block is a constant; at a zero neither is.
#' rbind(positive = c(mu_zi = H$mu_zi[y > 0][1], zi_zi = H$zi_zi[y > 0][1]),
#'       at_zero = c(mu_zi = H$mu_zi[y == 0][1], zi_zi = H$zi_zi[y == 0][1]))
#' -1 / 0.75^2
S7::method(distrib_hessian, ZeroInflatedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  zi <- pars$mix
  zi_name <- distrib@params[distrib@n_params]
  n <- length(y)

  f0 <- distrib_pdf(parent, 0, pars$orig)
  l0 <- zi + (1 - zi) * f0
  w0 <- ((1 - zi) * f0) / l0

  grad_0 <- distrib_gradient(parent, 0, pars$orig)
  hess_0_obs <- distrib_hessian(parent, 0, pars$orig)
  h_orig <- distrib_hessian(parent, y, pars$orig)
  pairs <- hess_pairs(names(pars$orig))

  res <- list()

  # Block theta-theta: at y=0, w*H(0) + w(1-w)*S(0)S(0)'
  for (nm in names(pairs)) {
    p1 <- pairs[[nm]][1]
    p2 <- pairs[[nm]][2]
    val_y0 <- w0 * hess_0_obs[[nm]] + w0 * (1 - w0) * grad_0[[p1]] * grad_0[[p2]]
    res[[nm]] <- ifelse(y == 0, val_y0, h_orig[[nm]])
  }

  # Block theta-zi: -f(0) S(0) / L0^2 at y=0, 0 otherwise
  factor_mix <- ifelse(y == 0, -f0 / (l0^2), 0)
  for (p in names(pars$orig)) {
    res[[paste0(p, "_", zi_name)]] <- factor_mix * grad_0[[p]]
  }

  # Block zi-zi
  val_zi_0 <- -((1 - f0)^2) / (l0^2)
  val_zi_pos <- -1 / ((1 - zi)^2)
  res[[paste0(zi_name, "_", zi_name)]] <- ifelse(y == 0, val_zi_0, val_zi_pos)

  expand_params(res[hess_names(distrib@params)], n)
}

#' @title Zero-Inflated Expected Information
#' @name distrib_expected_hessian.ZeroInflatedDistrib
#'
#' @description
#' Computes the expectation of the observed Hessian in closed form, by
#' splitting the expectation over the two events \eqn{y = 0} and \eqn{y > 0}.
#' The zero contributes with probability \eqn{L_0} and carries the mixture's
#' corrections; the positive part contributes with probability \eqn{1 - L_0}
#' and carries the parent's own expected Hessian, less what the zero would have
#' contributed to it. No component depends on the data, so every returned
#' vector is constant and `y` is read for its length alone.
#'
#' @details
#' This is EXACT. `approx` and `nsim` are accepted so that the signature
#' matches the generic's, and neither is read. The route
#' works because the parent's expected Hessian is an expectation over its whole
#' support, from which the single point at zero can be subtracted in closed
#' form.
#'
#' @param distrib A `ZeroInflatedDistrib` object, from [zero_inflated()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with the parent's parameters followed by `zi`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored: the expectation is exact. Present so that the
#'   signature matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length `length(y)`, keyed as
#'   [`hess_names(distrib@params)`][hess_names], each vector constant.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\zeta} the probability of a
#' structural zero, \eqn{L_0 = \zeta + (1-\zeta)f(0)} the inflated mass at
#' zero, \eqn{w = (1-\zeta)f(0)/L_0} the posterior probability that an observed
#' zero came from the parent, \eqn{s} the parent's score, \eqn{H} its observed
#' Hessian and \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_hessian.ZeroInflatedDistrib()] for the observed matrix,
#'   [fit_distrib()], whose Fisher scoring inverts this, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.25)
#' set.seed(2)
#' y <- distrib_rng(d, 200, theta)
#'
#' EH <- distrib_expected_hessian(d, y, theta)
#' vapply(EH, function(z) z[1], numeric(1))
#'
#' # Closed form, so two calls agree to the bit and nothing is sampled.
#' identical(EH, distrib_expected_hessian(d, y, theta))
#'
#' # It is what summing the observed Hessian against the mass function gives,
#' # over the support taken far enough out.
#' sup <- 0:400
#' m <- distrib_pdf(d, sup, theta)
#' Hs <- distrib_hessian(d, sup, theta)
#' rbind(summed = vapply(Hs, function(z) sum(z * m), numeric(1)),
#'       closed = vapply(EH, function(z) z[1], numeric(1)))
S7::method(distrib_expected_hessian, ZeroInflatedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  zi <- pars$mix
  zi_name <- distrib@params[distrib@n_params]
  n <- length(y)

  f0 <- distrib_pdf(parent, 0, pars$orig)
  l0 <- zi + (1 - zi) * f0
  w0 <- ((1 - zi) * f0) / l0

  grad_0 <- distrib_gradient(parent, 0, pars$orig)
  hess_0_obs <- distrib_hessian(parent, 0, pars$orig)
  h_orig_exp <- distrib_expected_hessian(parent, y, pars$orig)
  pairs <- hess_pairs(names(pars$orig))

  res <- list()

  # Block theta-theta:
  #   L0 * [w0 H(0) + w0(1-w0) S(0)S(0)'] + (1-zi) * (E[H] - f(0)H(0))
  for (nm in names(pairs)) {
    p1 <- pairs[[nm]][1]
    p2 <- pairs[[nm]][2]
    h_zi_0 <- w0 * hess_0_obs[[nm]] + w0 * (1 - w0) * grad_0[[p1]] * grad_0[[p2]]
    contrib_pos <- h_orig_exp[[nm]] - hess_0_obs[[nm]] * f0
    res[[nm]] <- l0 * h_zi_0 + (1 - zi) * contrib_pos
  }

  # Block theta-zi: -f(0) S(0) / L0
  for (p in names(pars$orig)) {
    res[[paste0(p, "_", zi_name)]] <- rep(-(f0 / l0) * grad_0[[p]], length.out = n)
  }

  # Block zi-zi: negative Fisher information of the mixture weight
  score_zi_0 <- (1 - f0) / l0
  score_zi_pos <- -1 / (1 - zi)
  res[[paste0(zi_name, "_", zi_name)]] <- -(l0 * score_zi_0^2 + (1 - l0) * score_zi_pos^2)

  expand_params(res[hess_names(distrib@params)], n)
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Zero-Inflated Distribution Object
#'
#' @description
#' Wraps a DISCRETE distribution into a mixture that keeps the parent intact
#' and adds a second source of zeros, with probability \eqn{\zeta} carried by
#' the new parameter `zi`:
#' \deqn{P(Y = y; \theta, \zeta) = \begin{cases}
#'   \zeta + (1 - \zeta)f(0; \theta) & y = 0 \\
#'   (1 - \zeta)f(y; \theta) & y > 0. \end{cases}}
#' It is the wrapper to reach for when the data carry MORE zeros than the
#' parent can produce and a zero may plausibly have come either from the count
#' process or from a mechanism that switches it off.
#'
#' @details
#' # Zero-inflation against zero-adjustment
#'
#' The two wrappers differ in what they do to the mass the parent already
#' places at zero. Zero-inflation ADDS to it, so
#' \eqn{P(Y = 0) = \zeta + (1-\zeta)f(0) > f(0)}: the model can produce more
#' zeros than the parent and never fewer, and the observed zeros are a mixture
#' of structural and sampling ones that no single observation can be assigned
#' to. Zero-adjustment REPLACES it: the parent is truncated away from zero and
#' the mass there becomes a free parameter, which may sit above or below
#' \eqn{f(0)}.
#'
#' A hurdle model therefore also handles UNDER-dispersed zeros, and its
#' likelihood factorizes into a binary part and a positive-count part that can
#' be read separately. Zero-inflation keeps the parent's interpretation, so
#' \eqn{\theta} still describes the count process the non-structural
#' observations come from, while the hurdle re-interprets \eqn{\theta} as the
#' parameters of a truncated law.
#'
#' # What the parent must be
#'
#' Zero-inflation adds mass to a zero that already carries some, so it needs a
#' discrete distribution with 0 in its support. A continuous one has
#' \eqn{P(Y = 0) = 0} and nothing to inflate; placing a point mass at zero
#' beside a density is zero-ADJUSTMENT, and [zero_adjusted()] handles it. The
#' constructor also fails where the result would not be identified:
#'
#' - the parent already models a probability of zero, since only the total mass
#'   at zero would then be identified;
#' - the support is too small for one more parameter. A distribution on \eqn{k}
#'   points has \eqn{k-1} free probabilities, so at least `n_params + 2` support
#'   points are needed. This rules out the Bernoulli and
#'   `binomial_distrib(size = 1)`, where `mu` and `zi` between them describe a
#'   single free cell.
#'
#' A large support is NECESSARY without being sufficient. With \eqn{\mu} large
#' enough that \eqn{f(0)} underflows, or \eqn{\zeta} close to zero, the ridge
#' reappears in the data instead of in the model, and [fit_distrib()] reports
#' the standard errors that reveal it.
#'
#' # What the result supports
#'
#' The whole `distrib` contract: the mass function, the distribution function,
#' the quantile function, the generator, the analytic score, the observed and
#' the expected Hessian, and the third and fourth derivatives from the shared
#' wrapper machinery. The moments come from [moment()] numerically.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\theta} its parameters,
#' \eqn{\zeta} the probability of a structural zero and \eqn{k} the number of
#' support points.
#'
#' @param distrib An object inheriting from `discrete_distrib` whose support
#'   includes 0, such as [poisson_distrib()] or [negbin2_distrib()]. A
#'   continuous distribution, an already-wrapped one, and a support of fewer
#'   than `n_params + 2` points are each rejected with an error saying which
#'   condition failed.
#' @param link_zi The link carrying \eqn{\zeta} to the unconstrained scale, a
#'   `linkfunctions7::link` object. Defaults to
#'   [linkfunctions7::logit_link()], which keeps it strictly inside
#'   \eqn{(0, 1)} at every point of the free scale.
#'
#' @return An S7 object of class [ZeroInflatedDistrib], inheriting from
#'   `discrete_distrib`. Its `params` are the parent's followed by `zi`;
#'   `n_params` is the parent's plus one; `params_bounds$zi` is
#'   \eqn{(0, 1)}; `link_params$zi` is `link_zi` and the rest are the parent's;
#'   `params_interpretation` gains `"prob. of structural zero"`; and
#'   `distrib_name` is `"zero-inflated "` followed by the parent's.
#'
#' @seealso [zero_adjusted()] for the hurdle counterpart,
#'   [ZeroInflatedDistrib] for the class, [truncated()], which the hurdle uses
#'   internally, and [check_distrib()] to validate the result.
#'
#' @examples
#' zip <- zero_inflated(poisson_distrib())
#' theta <- list(mu = 3, zi = 0.2)
#' zip@params
#'
#' distrib_pdf(zip, 0:5, theta)
#'
#' # More mass at zero than the Poisson alone can put there, and the rest of
#' # the mass function shrunk by 1 - zi.
#' c(inflated = distrib_pdf(zip, 0, theta), poisson = dpois(0, 3))
#' all.equal(distrib_pdf(zip, 1:5, theta), 0.8 * dpois(1:5, 3))
#'
#' # It is still a distribution: the mass sums to one.
#' sum(distrib_pdf(zip, 0:200, theta))
#'
#' # A fit recovers both parameters, the extra zeros identifying zi.
#' set.seed(1)
#' y <- distrib_rng(zip, 2000, theta)
#' round(coef(fit_distrib(zip, y)), 3)
#'
#' # A Bernoulli has no room for a second parameter, and a wrapper cannot be
#' # stacked on another.
#' try(zero_inflated(bernoulli_distrib()))
#' try(zero_inflated(zero_adjusted(poisson_distrib())))
#'
#' # Nor can a continuous parent be inflated: there is no mass at zero to add
#' # to, and the message names the wrapper that does apply.
#' try(zero_inflated(gaussian1_distrib()))
#'
#' @importFrom linkfunctions7 logit_link
#' @export
zero_inflated <- function(distrib, link_zi = logit_link()) {
  if (S7::S7_inherits(distrib, continuous_distrib)) {
    stop(paste0(
      "zero_inflated() requires a discrete distribution: a continuous one has\n",
      "  P(Y = 0) = 0, so there is no mass at zero to inflate. To place a point mass\n",
      "  at zero alongside a density, use zero_adjusted()."
    ), call. = FALSE)
  }
  if (!S7::S7_inherits(distrib, discrete_distrib)) {
    stop("zero_inflated() requires a discrete distribution ('discrete_distrib').", call. = FALSE)
  }
  check_not_stacked(distrib, "zero_inflated", "zi")
  if (distrib@bounds[1] > 0) {
    stop(sprintf(paste0(
      "zero_inflated() requires 0 in the support of '%s', which starts at %g.\n",
      "  With P(Y = 0) = 0 there is nothing to inflate: the mixture would put all of\n",
      "  its zeros in the added component, which is a zero-adjusted model."
    ), distrib@distrib_name, distrib@bounds[1]), call. = FALSE)
  }
  check_support_is_rich_enough(distrib, "zero_inflated")

  ZeroInflatedDistrib(
    parent_distrib = distrib,
    distrib_name = paste0("zero-inflated ", distrib@distrib_name),
    dimension = distrib@dimension,
    bounds = c(min(0, distrib@bounds[1]), distrib@bounds[2]),

    params = c(distrib@params, "zi"),
    params_interpretation = c(distrib@params_interpretation, zi = "prob. of structural zero"),
    n_params = distrib@n_params + 1,

    params_bounds = c(distrib@params_bounds, list(zi = c(0, 1))),
    link_params = c(distrib@link_params, list(zi = link_zi)),
    params_smooth = c(param_smoothness(distrib), zi = TRUE)
  )
}
