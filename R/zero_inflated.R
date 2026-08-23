#' @include distrib.R generics.R utility_functions.R
NULL

#' @title S7 Class for Zero-Inflated Distributions
#' @name ZeroInflatedDistrib
#'
#' @description
#' A subclass of `discrete_distrib` representing the zero-inflated version of a
#' wrapped discrete distribution: a mixture of a point mass at zero (with probability
#' \eqn{\zeta}) and the original count distribution.
#' @inheritParams distrib
#' @param parent_distrib The wrapped `discrete_distrib` object.
#' @return An object of class `ZeroInflatedDistrib`.
#' @seealso [zero_inflated()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.ZeroInflatedDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ZeroInflatedDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ZeroInflatedDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ZeroInflatedDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ZeroInflatedDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ZeroInflatedDistrib],
#'   [`distrib_rng()`][distrib_rng.ZeroInflatedDistrib]
#'
#' Everything else is inherited from [discrete_distrib()].
ZeroInflatedDistrib <- S7::new_class("ZeroInflatedDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib
  )
)

#' Split a Wrapper's Parameters From Its Parent's
#'
#' @description
#' Separates the full `theta` into the parent distribution's parameters and
#' the single mixture parameter the wrapper adds.
#'
#' @details
#' Both zero wrappers append their parameter last, so the split is positional and
#' does not depend on what that parameter is called -- `zi` for inflation,
#' `pi` for adjustment.
#'
#' @param distrib A zero-inflated or zero-adjusted distribution object.
#' @param theta A named list of parameters, already aligned.
#'
#' @return A list with `orig`, the parent's parameters, and `mix`,
#'   the wrapper's own.
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

#' Number of Points in a Discrete Support
#'
#' @description
#' How many points the distribution's support contains, `Inf` when it is
#' unbounded.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#'
#' @return A single number.
#'
#' @seealso [check_support_is_rich_enough()]
#' @keywords internal
n_support_points <- function(distrib) {
  b <- distrib@bounds
  if (!all(is.finite(b))) Inf else b[2] - b[1] + 1
}

#' Does This Distribution Already Model a Probability of Zero?
#'
#' @description
#' `TRUE` for a distribution produced by [zero_inflated()] or
#' [zero_adjusted()], in either of its two forms.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#'
#' @return A single logical.
#'
#' @seealso [check_not_stacked()]
#' @keywords internal
is_zero_wrapper <- function(distrib) {
  S7::S7_inherits(distrib, ZeroInflatedDistrib) ||
    S7::S7_inherits(distrib, ZeroAdjustedDiscreteDistrib) ||
    S7::S7_inherits(distrib, ZeroAdjustedContinuousDistrib)
}

#' Reject the Composition of Two Zero Wrappers
#'
#' @description
#' Rejects an attempt to wrap a distribution that already models the probability
#' of a zero, and rejects a parameter name the parent has already used.
#'
#' @details
#' Two zero parameters cannot both be identified. Zero-truncating a distribution
#' that already has one removes it from the likelihood entirely -- the factor
#' cancels between the numerator and the truncation constant, leaving its score
#' identically zero -- and mixing a further point mass in only ever shifts the
#' total mass at zero, which one parameter already describes.
#'
#' The distributions this rejects are well-defined but not estimable, and
#' nothing detects that at run time: the pmf sums to one,
#' [check_distrib()] passes, and a fit converges to an arbitrary
#' point of a flat ridge. The constructor is therefore the only place the
#' condition can be enforced.
#'
#' @param distrib The parent distribution being wrapped.
#' @param fun The calling constructor's name, used in the message.
#' @param param The name of the parameter the wrapper wants to add.
#'
#' @return Invisibly `NULL`; raises an error if either condition fails.
#'
#' @seealso [zero_inflated()], [zero_adjusted()]
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

#' Reject a Model With More Parameters Than the Support Can Distinguish
#'
#' @description
#' Enforces the counting rule that makes a zero wrapper identifiable.
#'
#' @details
#' A discrete distribution on \eqn{k} points has \eqn{k-1} free probabilities,
#' and either wrapper spends `n_params + 1` of them, so \eqn{k \ge}
#' `n_params + 2` is necessary. The bound is the same for inflation and for
#' adjustment.
#'
#' What it rules out is exactly the Bernoulli, and
#' `binomial_distrib(size = 1)` with it. Zero-inflating a Bernoulli gives
#' two parameters for the one free cell of \eqn{\{0, 1\}}; zero-adjusting it
#' leaves the truncated part concentrated on \eqn{\{1\}} with no free parameter
#' at all, so \eqn{\mu} disappears from the likelihood and the pmf is literally
#' the same for \eqn{\mu = 0.2} and \eqn{\mu = 0.9}.
#'
#' @param distrib The parent distribution being wrapped.
#' @param fun The calling constructor's name, used in the message.
#'
#' @return Invisibly `NULL`; raises an error when the support is too small.
#'
#' @seealso [n_support_points()], [check_not_stacked()]
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
#' @description
#' \deqn{P(Y=y) = \zeta\,\mathbb{I}(y=0) + (1-\zeta) f(y; \theta)}
#' where \eqn{f} is the PMF of the wrapped distribution.
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @param log Logical; if `TRUE`, returns the log-probability.
#' @return A numeric vector of density values.
#' @seealso [zero_inflated()]
S7::method(distrib_pdf, ZeroInflatedDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  pars <- split_mix_theta(distrib, theta)
  zi <- pars$mix

  res <- (1 - zi) * distrib_pdf(distrib@parent_distrib, y, pars$orig, log = FALSE) +
    zi * (y == 0)

  if (log) log(res) else res
}

#' @title Zero-Inflated Cumulative Distribution Function
#' @name distrib_cdf.ZeroInflatedDistrib
#' @description
#' \deqn{F_{ZI}(q) = (1-\zeta) F(q; \theta) + \zeta\,\mathbb{I}(q \ge 0)}
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned as logs.
#' @return A numeric vector of cumulative probabilities.
#' @seealso [zero_inflated()]
S7::method(distrib_cdf, ZeroInflatedDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  pars <- split_mix_theta(distrib, theta)
  zi <- pars$mix

  res <- (1 - zi) * distrib_cdf(distrib@parent_distrib, q, pars$orig) + zi * (q >= 0)
  res <- pmin(pmax(res, 0), 1)

  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Zero-Inflated Quantile Function
#' @name distrib_quantile.ZeroInflatedDistrib
#' @description
#' Inverts the mixture CDF: the quantile is 0 for \eqn{p \le \zeta + (1-\zeta)F(0;\theta)},
#' otherwise \eqn{Q\left(\dfrac{p-\zeta}{1-\zeta}; \theta\right)}.
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if `TRUE`, probabilities are given as logs.
#' @return A numeric vector of quantiles.
#' @seealso [zero_inflated()]
S7::method(distrib_quantile, ZeroInflatedDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
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
#' @description
#' Draws from the wrapped distribution, then replaces a Bernoulli(\eqn{\zeta}) fraction
#' of the draws with structural zeros.
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @return A numeric vector of random draws.
#' @seealso [zero_inflated()]
S7::method(distrib_rng, ZeroInflatedDistrib) <- function(distrib, n, theta) {
  pars <- split_mix_theta(distrib, theta)
  y <- distrib_rng(distrib@parent_distrib, n, pars$orig)
  y[stats::runif(n) < pars$mix] <- 0
  y
}

#' @title Zero-Inflated Analytical Gradient
#' @name distrib_gradient.ZeroInflatedDistrib
#' @description
#' Score function of the zero-inflated model. For the parent's parameters the score is
#' the parent's score weighted by \eqn{w = (1-\zeta)f(0)/L_0} at \eqn{y=0} (and 1 otherwise),
#' where \eqn{L_0 = \zeta + (1-\zeta)f(0)}. For \eqn{\zeta}:
#' \deqn{\dfrac{\partial \ell}{\partial \zeta} = \mathbb{I}(y=0)\dfrac{1-f(0)}{L_0} - \mathbb{I}(y>0)\dfrac{1}{1-\zeta}}
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @return A list containing the vectors of first derivatives.
#' @seealso [zero_inflated()]
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

#' @title Zero-Inflated Analytical Observed Hessian
#' @name distrib_hessian.ZeroInflatedDistrib
#' @description
#' Observed Hessian of the zero-inflated model, combining the parent's observed Hessian
#' with rank-one corrections at \eqn{y=0} (see the old package documentation for the full
#' derivation).
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @return A list containing the vectors of second derivatives.
#' @seealso [zero_inflated()]
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

#' @title Zero-Inflated Analytical Expected Hessian
#' @name distrib_expected_hessian.ZeroInflatedDistrib
#' @description
#' Expected Hessian (negative Fisher information) of the zero-inflated model, derived by
#' decomposing the expectation over \eqn{y=0} and \eqn{y>0} and using the parent's
#' expected Hessian.
#' @param distrib A `ZeroInflatedDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by `zi`.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso [zero_inflated()]
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

#' Zero-Inflated Distribution Object (Discrete)
#'
#' @description
#' Creates a zero-inflated version of an existing **discrete** distribution: a
#' mixture that keeps the parent intact and adds a second source of zeros, with
#' probability \eqn{\zeta} (parameter `zi`).
#'
#' Zero-inflation is the right wrapper when the data contain *more* zeros than
#' the parent can produce, and a zero can plausibly have come either from the count
#' process or from a separate mechanism that switches it off. If instead the zeros
#' come from one identifiable mechanism and the positive values from another, the
#' appropriate model is the hurdle, [zero_adjusted()].
#'
#' @param distrib An object inheriting from `discrete_distrib` whose support
#'   includes 0, e.g. [poisson_distrib()] or [negbin2_distrib()].
#' @param link_zi A link function object for the zero-inflation probability \eqn{\zeta}.
#'   Defaults to [linkfunctions7::logit_link()].
#'
#' @details
#' \deqn{
#' P(Y=y; \theta, \zeta) =
#' \begin{cases}
#' \zeta + (1 - \zeta)f(0; \theta) & y = 0 \\
#' (1 - \zeta)f(y; \theta) & y > 0
#' \end{cases}
#' }
#'
#' **Zero-inflation versus zero-adjustment.** The two wrappers differ in what
#' they do to the mass the parent already places at zero. Zero-inflation
#' *adds* to it, so \eqn{P(Y = 0) = \zeta + (1-\zeta)f(0) > f(0)}: the model can
#' only ever produce more zeros than the parent, never fewer, and the observed zeros
#' are a mixture of structural and sampling ones that no single observation can be
#' assigned to. Zero-adjustment ([zero_adjusted()]) *replaces* it: the
#' parent is truncated away from zero and the mass at zero becomes a free parameter,
#' which can be above or below \eqn{f(0)}. A hurdle model therefore also handles
#' *under*-dispersed zeros, and its likelihood factorizes into a binary part and
#' a positive-count part that can be read separately. Zero-inflation keeps the
#' parent's interpretation --- \eqn{\theta} still describes the count process the
#' non-structural observations come from --- while the hurdle re-interprets
#' \eqn{\theta} as the parameters of a truncated law.
#'
#' **What the parent must be.** Zero-inflation adds mass to a zero that already
#' carries some, so it requires a discrete distribution with \eqn{0} in its support.
#' A continuous distribution has \eqn{P(Y = 0) = 0} and nothing to inflate; putting a
#' point mass at zero next to a density is zero-*adjustment*, and
#' [zero_adjusted()] handles it. Constructing the object also fails when the
#' result would not be identified:
#'
#' - the parent already models a probability of zero (a wrapper cannot be
#'   stacked on another wrapper: only the total mass at zero would be identified);
#' - the support is too small for one more parameter --- a distribution on
#'   \eqn{k} points has \eqn{k-1} free probabilities, so at least
#'   `n_params + 2` support points are needed. This rules out the Bernoulli
#'   and `binomial_distrib(size = 1)`, where `mu` and `zi` between
#'   them describe a single free cell.
#'
#' A large support is necessary but not sufficient: with \eqn{\mu} large enough that
#' \eqn{f(0)} underflows, or \eqn{\zeta} close to 0, the ridge reappears in the data
#' rather than in the model. [fit_distrib()] reports the standard errors
#' that reveal it.
#'
#' The resulting object supports the full `distrib` API: pdf, cdf, quantile, rng,
#' analytical gradient, observed and expected Hessian (all derived from the parent's),
#' plus numerical moments via [moment()].
#'
#' @return An S7 object of class `ZeroInflatedDistrib` (inheriting from `discrete_distrib`).
#'
#' @examples
#' zip <- zero_inflated(poisson_distrib())
#' distrib_pdf(zip, 0:5, list(mu = 3, zi = 0.2))
#'
#' # More mass at zero than the Poisson alone can put there
#' distrib_pdf(zip, 0, list(mu = 3, zi = 0.2)) > dpois(0, 3)
#'
#' # A Bernoulli has no room for a second parameter
#' try(zero_inflated(bernoulli_distrib()))
#'
#' @seealso [zero_adjusted()] for the hurdle counterpart,
#'   [check_distrib()] to validate the result.
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
