#' @include distrib.R generics.R utility_functions.R numerical_functions.R zero_inflated.R
NULL

#' @title S7 Class for Zero-Adjusted Discrete (Hurdle) Distributions
#' @name ZeroAdjustedDiscreteDistrib
#'
#' @description
#' The S7 class of the hurdle version of a discrete distribution, with mass
#' function
#' \deqn{P(Y = y) = \begin{cases} \pi & y = 0 \\
#'   (1-\pi)\dfrac{f(y; \theta)}{1 - f(0; \theta)} & y > 0. \end{cases}}
#' The mass at zero is REPLACED by the free parameter \eqn{\pi}, carried by
#' `za`, and the parent is truncated away from zero. Unlike
#' [ZeroInflatedDistrib], the model can produce FEWER zeros than the parent as
#' well as more.
#'
#' @details
#' The likelihood factorizes into a binary part in \eqn{\pi} and a
#' positive-count part in \eqn{\theta}. Every mixed block of the Hessian is
#' exactly zero for that reason, and the two halves can be read separately. It
#' also re-interprets \eqn{\theta}: they are now the parameters of a TRUNCATED
#' law, not of the count process the observations come from.
#'
#' Build one with [zero_adjusted()], which checks that the parent is not
#' already a zero wrapper and has enough support points. This page documents
#' the raw S7 constructor, which checks neither.
#'
#' @inheritParams distrib
#' @param parent_distrib The wrapped `discrete_distrib` object.
#'
#' @return An S7 object of class `ZeroAdjustedDiscreteDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. It carries `parent_distrib` beside
#'   the parent's properties. For an object built by [zero_adjusted()] the
#'   parameters are the parent's followed by `za`, whose bound is
#'   \eqn{(0, 1)} and whose interpretation is `"prob. of zero"`.
#'
#' @seealso [zero_adjusted()] to build one, [ZeroInflatedDistrib] for the
#'   wrapper that ADDS to the mass at zero, and
#'   [ZeroAdjustedContinuousDistrib] for the continuous branch, which produces
#'   a mixed distribution: a density with an atom on it.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.ZeroAdjustedDiscreteDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ZeroAdjustedDiscreteDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ZeroAdjustedDiscreteDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ZeroAdjustedDiscreteDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ZeroAdjustedDiscreteDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ZeroAdjustedDiscreteDistrib],
#'   [`distrib_rng()`][distrib_rng.ZeroAdjustedDiscreteDistrib]
#'
#' The third and fourth derivatives come from the shared wrapper machinery in
#' `wrapper_derivatives.R`; everything else is inherited from
#' [discrete_distrib()].
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\pi} the probability of a zero
#' and \eqn{\ell} the log-mass of one observation. The positive part is the
#' parent truncated away from zero, with mass
#' \eqn{f(y)/\{1 - f(0)\}} at \eqn{y > 0}.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#' S7::S7_inherits(d, discrete_distrib)
#' d@params
#'
#' # The mass at zero IS the parameter, and can be set below the parent's,
#' # which zero-inflation cannot do.
#' c(adjusted = distrib_pdf(d, 0, theta), parent = dpois(0, 3))
#' distrib_pdf(d, 0, list(mu = 3, za = 0.01))
#'
#' # The positive part is the parent renormalized away from zero.
#' all.equal(distrib_pdf(d, 1:4, theta),
#'           0.6 * dpois(1:4, 3) / (1 - dpois(0, 3)))
ZeroAdjustedDiscreteDistrib <- S7::new_class("ZeroAdjustedDiscreteDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib
  )
)

#' @title S7 Class for Zero-Adjusted Continuous Distributions
#' @name ZeroAdjustedContinuousDistrib
#'
#' @description
#' The S7 class of a continuous distribution with a point mass at zero:
#' \deqn{P(Y = 0) = \pi, \qquad f_Y(y) = (1-\pi) f(y; \theta) \; (y \ne 0).}
#' The result is a MIXED distribution, a density plus an atom, and it declares
#' that atom through [distrib_atoms()]. That declaration is how
#' [check_distrib()] and [expectation()] learn to treat it as one.
#'
#' @details
#' No truncation is needed here. A continuous parent has \eqn{P(Y = 0) = 0}, so
#' there is no mass to remove before placing the atom, and the density is
#' simply scaled by \eqn{1-\pi}. That is the whole difference from
#' [ZeroAdjustedDiscreteDistrib], whose parent must be truncated away from the
#' point it already occupies.
#'
#' The likelihood factorizes completely: the mixed blocks of the Hessian are
#' exactly zero, \eqn{\pi} is estimated by the proportion of zeros and
#' \eqn{\theta} by the parent's own fit to the non-zero observations.
#'
#' Build one with [zero_adjusted()]. This page documents the raw S7
#' constructor, which validates nothing.
#'
#' @inheritParams distrib
#' @param parent_distrib The wrapped `continuous_distrib` object.
#'
#' @return An S7 object of class `ZeroAdjustedContinuousDistrib`, inheriting
#'   from `continuous_distrib` and from `distrib`. It carries `parent_distrib`
#'   beside the parent's properties, with `za` added last.
#'
#' @seealso [zero_adjusted()] to build one, [ZeroAdjustedDiscreteDistrib] for
#'   the discrete branch, [distrib_atoms.ZeroAdjustedContinuousDistrib()] for
#'   the atom it declares, and [folded()], which REJECTS a parent of this class
#'   for exactly that atom.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_atoms()`][distrib_atoms.ZeroAdjustedContinuousDistrib],
#'   [`distrib_cdf()`][distrib_cdf.ZeroAdjustedContinuousDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.ZeroAdjustedContinuousDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.ZeroAdjustedContinuousDistrib],
#'   [`distrib_gradient()`][distrib_gradient.ZeroAdjustedContinuousDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.ZeroAdjustedContinuousDistrib],
#'   [`distrib_hessian()`][distrib_hessian.ZeroAdjustedContinuousDistrib],
#'   [`distrib_pdf()`][distrib_pdf.ZeroAdjustedContinuousDistrib],
#'   [`distrib_quantile()`][distrib_quantile.ZeroAdjustedContinuousDistrib],
#'   [`distrib_rng()`][distrib_rng.ZeroAdjustedContinuousDistrib],
#'   [`expectation()`][expectation.ZeroAdjustedContinuousDistrib]
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{\pi} the probability of a zero and
#' \eqn{\ell} the log-density of one observation.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#' d@params
#'
#' # It is a MIXED distribution: an atom at zero and a density elsewhere.
#' distrib_atoms(d, theta)
#' c(at_zero = distrib_pdf(d, 0, theta),
#'   elsewhere = distrib_pdf(d, 2, theta),
#'   scaled_parent = 0.7 * dnorm(2, 1, 2))
#'
#' # The density part integrates to 1 - pi, the atom carrying the rest.
#' integrate(function(z) ifelse(z == 0, 0, distrib_pdf(d, z, theta)),
#'           -Inf, Inf)$value
ZeroAdjustedContinuousDistrib <- S7::new_class("ZeroAdjustedContinuousDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib
  )
)

# ============================== DISCRETE (HURDLE) ==============================

#' @title Zero-Adjusted Discrete Probability Mass Function
#' @name distrib_pdf.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Computes the hurdle mass function
#' \deqn{P(Y = 0) = \pi, \qquad
#'   P(Y = y) = (1-\pi)\frac{f(y;\theta)}{1 - f(0;\theta)} \quad (y > 0).}
#' The mass at zero is the parameter itself, so it may sit ABOVE or BELOW the
#' parent's \eqn{f(0)}; the positive part is the parent renormalized away from
#' zero. Everything is computed on the log scale and exponentiated at the end,
#' with the normalizing constant taken through [base::log1p()] so that a parent
#' with almost no mass at zero does not lose digits to `log(1 - f0)`.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `za`,
#'   each a numeric vector of length 1 or of the length of `y`. `za` must lie
#'   strictly inside \eqn{(0, 1)}.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned, which
#'   is the quantity actually computed. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{F} its distribution function,
#' \eqn{\pi} the probability of a zero, \eqn{s} the parent's score and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_cdf.ZeroAdjustedDiscreteDistrib()] for the distribution
#'   function, [distrib_pdf.ZeroInflatedDistrib()] for the wrapper that adds
#'   to the mass at zero where this one replaces it, and [distrib_pdf()] for
#'   the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#'
#' distrib_pdf(d, 0:5, theta)
#'
#' # The mass at zero is the parameter, and the rest is the parent
#' # renormalized away from it.
#' c(at_zero = distrib_pdf(d, 0, theta), parameter = 0.4)
#' all.equal(distrib_pdf(d, 1:5, theta),
#'           0.6 * dpois(1:5, 3) / (1 - dpois(0, 3)))
#'
#' # It can produce FEWER zeros than the parent, which inflation cannot.
#' c(hurdle = distrib_pdf(d, 0, list(mu = 3, za = 0.01)),
#'   parent = dpois(0, 3))
#'
#' # And it is a distribution: the mass sums to one.
#' sum(distrib_pdf(d, 0:200, theta))
S7::method(distrib_pdf, ZeroAdjustedDiscreteDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  za <- pars$mix

  f0 <- distrib_pdf(parent, 0, pars$orig)
  log_res_pos <- log(1 - za) + distrib_pdf(parent, y, pars$orig, log = TRUE) - log1p(-f0)

  log_res <- ifelse(y == 0, rep(log(za), length.out = length(y)), log_res_pos)

  if (log) log_res else exp(log_res)
}

#' @title Zero-Adjusted Discrete Cumulative Distribution Function
#' @name distrib_cdf.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Computes the hurdle distribution function
#' \deqn{F_{ZA}(q) = \pi + (1-\pi)\frac{F(q;\theta) - f(0;\theta)}{1 - f(0;\theta)}
#'   \quad (q \ge 0),}
#' the parent's own distribution function shifted and rescaled by the
#' truncation. It is `0` below zero and is clamped to \eqn{[0, 1]} against
#' rounding.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param q A numeric vector of quantiles. Values below zero give `0`.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   formed as \eqn{1 - F}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm is returned.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities, in \eqn{[0, 1]}.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{F} its distribution function,
#' \eqn{\pi} the probability of a zero, \eqn{s} the parent's score and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_pdf.ZeroAdjustedDiscreteDistrib()] for the mass function,
#'   [distrib_quantile.ZeroAdjustedDiscreteDistrib()], which inverts this, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#'
#' distrib_cdf(d, c(0, 2, 5), theta)
#'
#' # At zero it is the parameter itself, the whole mass there.
#' c(cdf_at_zero = distrib_cdf(d, 0, theta), parameter = 0.4)
#'
#' # It agrees with the mass function summed, as it must on a lattice.
#' c(cdf = distrib_cdf(d, 5, theta), summed = sum(distrib_pdf(d, 0:5, theta)))
#'
#' # Both tails and the logarithm.
#' distrib_cdf(d, 2, theta, lower.tail = FALSE)
#' distrib_cdf(d, 2, theta, log.p = TRUE)
S7::method(distrib_cdf, ZeroAdjustedDiscreteDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  za <- pars$mix

  f0 <- distrib_pdf(parent, 0, pars$orig)
  F_orig <- distrib_cdf(parent, q, pars$orig)

  cdf_trunc <- pmax(0, (F_orig - f0) / (1 - f0))
  res <- za + (1 - za) * cdf_trunc
  res[q < 0] <- 0
  res <- pmin(res, 1)

  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Zero-Adjusted Discrete Quantile Function
#' @name distrib_quantile.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Inverts the hurdle distribution function. The quantile is `0` for
#' \eqn{p \le \pi}, the whole mass at zero, and otherwise the PARENT's quantile
#' at \eqn{u\{1 - f(0)\} + f(0)} with \eqn{u = (p - \pi)/(1 - \pi)}: the
#' probability is rescaled out of the hurdle and back onto the parent's own
#' scale before the parent is asked.
#'
#' @details
#' The result is a LATTICE quantile and overshoots, as any discrete quantile
#' does. Every probability at or below \eqn{\pi} maps to zero, so the first
#' step is as large as the mass at zero.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param p A numeric vector of probabilities, clamped to \eqn{[0, 1]} after
#'   the `log.p` and `lower.tail` transformations are applied.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles, of the recycled length of `p` and
#'   `theta`. Every value is either `0` or a positive support point of the
#'   parent.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{F} its distribution function,
#' \eqn{\pi} the probability of a zero, \eqn{s} the parent's score and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_cdf.ZeroAdjustedDiscreteDistrib()], which this inverts,
#'   and [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#'
#' distrib_quantile(d, c(0.2, 0.5, 0.95), theta)
#'
#' # Everything at or below the mass at zero maps to zero.
#' distrib_quantile(d, c(0.1, 0.4), theta)
#'
#' # The round trip overshoots, as on any lattice.
#' p <- c(0.5, 0.95)
#' rbind(asked = p,
#'       reached = distrib_cdf(d, distrib_quantile(d, p, theta), theta))
S7::method(distrib_quantile, ZeroAdjustedDiscreteDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  p <- pmin(pmax(p, 0), 1)

  all_params <- expand_params(c(list(.p = p), theta))
  p <- all_params$.p
  th_orig <- all_params[names(pars$orig)]
  za <- all_params[[distrib@params[distrib@n_params]]]

  f0 <- distrib_pdf(parent, 0, th_orig)
  q_vals <- numeric(length(p))
  idx <- (p > za)

  if (any(idx)) {
    u <- (p[idx] - za[idx]) / (1 - za[idx])
    f0_sub <- if (length(f0) > 1) f0[idx] else f0
    target <- pmin(u * (1 - f0_sub) + f0_sub, 1 - 1e-10)
    th_sub <- lapply(th_orig, function(x) if (length(x) > 1) x[idx] else x)
    q_vals[idx] <- distrib_quantile(parent, target, th_sub)
  }

  q_vals
}

#' @title Zero-Adjusted Discrete Random Number Generator
#' @name distrib_rng.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Draws zeros with probability \eqn{\pi} and otherwise samples from the
#' ZERO-TRUNCATED parent, so no draw from the positive part is ever zero. The
#' truncated draw is taken by inverting the parent's distribution function at a
#' uniform rescaled onto \eqn{(f(0), 1)}, which needs no rejection loop and so
#' terminates in bounded time however small the positive mass is.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param n The number of draws, a single non-negative whole number.
#' @param theta A named list with the parent's parameters followed by `za`.
#'
#' @return A numeric vector of length `n`, every value a non-negative support
#'   point of the parent.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{F} its distribution function,
#' \eqn{\pi} the probability of a zero, \eqn{s} the parent's score and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_pdf.ZeroAdjustedDiscreteDistrib()] for the mass function
#'   these are drawn from, [truncated()] for the truncation the positive part
#'   uses, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#'
#' set.seed(1)
#' distrib_rng(d, 10, theta)
#'
#' # A large sample reproduces the mass at zero, and no positive draw is zero.
#' set.seed(1)
#' big <- distrib_rng(d, 20000, theta)
#' c(sampled_zero = mean(big == 0), parameter = 0.4,
#'   smallest_positive = min(big[big > 0]))
S7::method(distrib_rng, ZeroAdjustedDiscreteDistrib) <- function(distrib, n, theta) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  za <- pars$mix

  is_zero <- stats::runif(n) < za
  y <- numeric(n)

  if (any(!is_zero)) {
    n_pos <- sum(!is_zero)
    f0 <- distrib_pdf(parent, 0, pars$orig)
    f0_sub <- if (length(f0) > 1) f0[!is_zero] else f0
    th_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[!is_zero] else x)

    u <- stats::runif(n_pos)
    y[!is_zero] <- distrib_quantile(parent, f0_sub + u * (1 - f0_sub), th_sub)
  }

  y
}

#' @title Zero-Adjusted Discrete Score
#' @name distrib_gradient.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Computes the first derivatives of the hurdle log-mass in closed form. The
#' likelihood SEPARATES, so each parameter is informed by one half of the data
#' only. For the parent's parameters the score is zero at a zero and the
#' truncated score at a positive observation,
#' \deqn{\frac{\partial \ell}{\partial \theta_i}
#'   = \mathbb{I}(y > 0)\left\{s_i(y) + \frac{f(0)}{1 - f(0)}\,s_i(0)\right\},}
#' the second term being the derivative of the normalizing constant
#' \eqn{-\log\{1 - f(0)\}}. For the hurdle parameter,
#' \deqn{\frac{\partial \ell}{\partial \pi}
#'   = \frac{\mathbb{I}(y = 0)}{\pi} - \frac{\mathbb{I}(y > 0)}{1 - \pi},}
#' which is the Bernoulli score of the indicator that the observation is zero.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. `za` rides a logit by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, in
#'   `distrib@params` order.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\pi} the probability of a zero,
#' \eqn{s} the parent's score, \eqn{H} its observed Hessian and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_hessian.ZeroAdjustedDiscreteDistrib()] for the second
#'   order, [distrib_gradient.ZeroInflatedDistrib()], where the two halves do
#'   NOT separate, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#' set.seed(2)
#' y <- distrib_rng(d, 300, theta)
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
#' # A zero tells the parent's parameters nothing at all, the two halves of
#' # the likelihood being separate.
#' unique(g$mu[y == 0])
#'
#' # And the za component is the Bernoulli score of the zero indicator, so it
#' # takes exactly two values.
#' unique(round(g$za, 10))
#' c(1 / 0.4, -1 / 0.6)
S7::method(distrib_gradient, ZeroAdjustedDiscreteDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]

  f0 <- distrib_pdf(parent, 0, pars$orig)
  score_0 <- distrib_gradient(parent, 0, pars$orig)
  correction <- f0 / (1 - f0)

  grad_orig <- distrib_gradient(parent, y, pars$orig)

  res <- list()
  for (nm in names(grad_orig)) {
    term_pos <- grad_orig[[nm]] + correction * score_0[[nm]]
    res[[nm]] <- ifelse(y == 0, 0, term_pos)
  }
  res[[za_name]] <- ifelse(y == 0, 1 / za, -1 / (1 - za))
  res
}

#' @title Zero-Adjusted Discrete Observed Hessian
#' @name distrib_hessian.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Computes the second derivatives of the hurdle log-mass in closed form. The
#' MIXED BLOCKS ARE IDENTICALLY ZERO, at every observation and every parameter
#' value, because the likelihood factorizes into a binary part in \eqn{\pi} and
#' a positive-count part in \eqn{\theta}. The parent block is zero at a zero
#' and, at a positive observation, the parent's own Hessian plus the truncation
#' correction
#' \deqn{H_{\mathrm{corr}} =
#'   \frac{\{1 - f(0)\}\,f''(0) + f'(0)^2}{\{1 - f(0)\}^2},}
#' the second derivative of \eqn{-\log\{1 - f(0)\}}. The hurdle block is
#' \eqn{-1/\pi^2} at a zero and \eqn{-1/(1-\pi)^2} elsewhere.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names]. Every
#'   mixed key holds a vector of exact zeros.
#'
#' @section Notation:
#' \eqn{f} is the parent's mass function, \eqn{\pi} the probability of a zero,
#' \eqn{s} the parent's score, \eqn{H} its observed Hessian and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_gradient.ZeroAdjustedDiscreteDistrib()] for the first
#'   order, [distrib_expected_hessian.ZeroAdjustedDiscreteDistrib()] for the
#'   expectation, [distrib_hessian.ZeroInflatedDistrib()], whose mixed block is
#'   NOT zero, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#' set.seed(2)
#' y <- distrib_rng(d, 300, theta)
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
#' # The mixed block is exactly zero at every observation, where the
#' # zero-inflated wrapper's is not.
#' all(H$mu_za == 0)
#' zi <- zero_inflated(poisson_distrib())
#' unique(distrib_hessian(zi, 0:2, list(mu = 3, zi = 0.4))$mu_zi != 0)
S7::method(distrib_hessian, ZeroAdjustedDiscreteDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]
  n <- length(y)

  f0 <- distrib_pdf(parent, 0, pars$orig)
  grad_0 <- distrib_gradient(parent, 0, pars$orig)
  hess_0_obs <- distrib_hessian(parent, 0, pars$orig)
  h_orig <- distrib_hessian(parent, y, pars$orig)
  denom <- 1 - f0
  pairs <- hess_pairs(names(pars$orig))

  res <- list()

  # Block za-za
  res[[paste0(za_name, "_", za_name)]] <- ifelse(y == 0, -1 / (za^2), -1 / ((1 - za)^2))

  # Mixed blocks are identically 0 (likelihood separation)
  for (nm in names(pars$orig)) {
    res[[paste0(nm, "_", za_name)]] <- rep(0, n)
  }

  # Block theta-theta with truncation correction
  for (nm in names(pairs)) {
    s1 <- grad_0[[pairs[[nm]][1]]]
    s2 <- grad_0[[pairs[[nm]][2]]]
    h_log_0 <- hess_0_obs[[nm]]

    f_prime_1 <- f0 * s1
    f_prime_2 <- f0 * s2
    f_second <- f0 * (h_log_0 + s1 * s2)
    hess_correction <- (denom * f_second + f_prime_1 * f_prime_2) / denom^2

    res[[nm]] <- ifelse(y == 0, 0, h_orig[[nm]] + hess_correction)
  }

  expand_params(res[hess_names(distrib@params)], n)
}

#' @title Zero-Adjusted Discrete Expected Information
#' @name distrib_expected_hessian.ZeroAdjustedDiscreteDistrib
#'
#' @description
#' Computes the expectation of the observed Hessian in closed form. The hurdle
#' block is the Bernoulli information \eqn{-1/\{\pi(1-\pi)\}}, the mixed blocks
#' are exactly zero, and the parent block is \eqn{(1-\pi)} times the truncated
#' parent's own expected Hessian: a positive observation arrives with
#' probability \eqn{1-\pi} and carries the whole of what \eqn{\theta} is
#' estimated from. No component depends on the data, so every returned vector
#' is constant.
#'
#' @details
#' The block diagonal is why a hurdle model is cheap to fit: the two halves can
#' be estimated independently, and the information matrix inverts blockwise.
#' `approx` and `nsim` are accepted so that the signature matches the generic's
#' and neither is read; the expectation here is exact.
#'
#' @param distrib A `ZeroAdjustedDiscreteDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with the parent's parameters followed by `za`.
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
#' \eqn{f} is the parent's mass function, \eqn{\pi} the probability of a zero,
#' \eqn{s} the parent's score, \eqn{H} its observed Hessian and \eqn{\ell} the
#' log-mass of one observation.
#'
#' @seealso [distrib_hessian.ZeroAdjustedDiscreteDistrib()] for the observed
#'   matrix, [fit_distrib()], whose Fisher scoring inverts this, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.4)
#' set.seed(2)
#' y <- distrib_rng(d, 300, theta)
#'
#' EH <- distrib_expected_hessian(d, y, theta)
#' vapply(EH, function(z) z[1], numeric(1))
#'
#' # The hurdle block is the Bernoulli information, and the mixed one is zero.
#' c(reported = EH$za_za[1], bernoulli = -1 / (0.4 * 0.6),
#'   mixed = EH$mu_za[1])
#'
#' # It is what summing the observed Hessian against the mass function gives.
#' sup <- 0:400
#' m <- distrib_pdf(d, sup, theta)
#' Hs <- distrib_hessian(d, sup, theta)
#' rbind(summed = vapply(Hs, function(z) sum(z * m), numeric(1)),
#'       closed = vapply(EH, function(z) z[1], numeric(1)))
S7::method(distrib_expected_hessian, ZeroAdjustedDiscreteDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]
  n <- length(y)

  f0 <- distrib_pdf(parent, 0, pars$orig)
  grad_0 <- distrib_gradient(parent, 0, pars$orig)
  hess_0_obs <- distrib_hessian(parent, 0, pars$orig)
  h_orig_exp <- distrib_expected_hessian(parent, y, pars$orig)
  denom <- 1 - f0
  pairs <- hess_pairs(names(pars$orig))

  res <- list()
  res[[paste0(za_name, "_", za_name)]] <- rep(-1 / (za * (1 - za)), length.out = n)

  for (nm in names(pars$orig)) {
    res[[paste0(nm, "_", za_name)]] <- rep(0, n)
  }

  for (nm in names(pairs)) {
    s1 <- grad_0[[pairs[[nm]][1]]]
    s2 <- grad_0[[pairs[[nm]][2]]]
    h_log_0 <- hess_0_obs[[nm]]

    f_prime_1 <- f0 * s1
    f_prime_2 <- f0 * s2
    f_second <- f0 * (h_log_0 + s1 * s2)
    hess_correction <- (denom * f_second + f_prime_1 * f_prime_2) / denom^2

    E_trunc <- (h_orig_exp[[nm]] - f0 * h_log_0) / denom
    res[[nm]] <- (1 - za) * (E_trunc + hess_correction)
  }

  expand_params(res[hess_names(distrib@params)], n)
}

# ============================== CONTINUOUS ==============================

#' @title Zero-Adjusted Continuous Probability Density Function
#' @name distrib_pdf.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Computes the mixed density
#' \deqn{f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi) f_W(y;\theta) \quad (y \ne 0).}
#' No truncation is needed: a continuous parent has \eqn{P(W = 0) = 0}, so
#' there is no mass to remove before placing the atom and the density is simply
#' scaled by \eqn{1-\pi}.
#'
#' @details
#' The value at zero is a PROBABILITY and the values elsewhere are DENSITIES,
#' so the returned vector mixes two kinds of number. That is what a mixed
#' distribution is, and it is why the object declares
#' [distrib_atoms.ZeroAdjustedContinuousDistrib()]: without that declaration a
#' consumer would integrate the returned function and find \eqn{1-\pi}.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations. Exactly zero gives \eqn{\pi}; any
#'   other value gives \eqn{(1-\pi)} times the parent's density there,
#'   including a point outside the parent's support, where the parent's density
#'   is zero.
#' @param theta A named list with the parent's parameters followed by `za`,
#'   each a numeric vector of length 1 or of the length of `y`.
#' @param log Logical of length 1. When `TRUE` the logarithm is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_atoms.ZeroAdjustedContinuousDistrib()] for the
#'   declaration that makes this a mixed distribution,
#'   [distrib_cdf.ZeroAdjustedContinuousDistrib()] for the distribution
#'   function, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' # A probability at zero, densities elsewhere.
#' c(at_zero = distrib_pdf(d, 0, theta),
#'   at_two = distrib_pdf(d, 2, theta),
#'   scaled_parent = 0.7 * dnorm(2, 1, 2))
#'
#' # The density part integrates to 1 - pi, the atom carrying the rest.
#' integrate(function(z) ifelse(z == 0, 0, distrib_pdf(d, z, theta)),
#'           -Inf, Inf)$value
#'
#' # Which is why the object declares its atom.
#' distrib_atoms(d, theta)
S7::method(distrib_pdf, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix

  log_f_orig <- distrib_pdf(distrib@parent_distrib, y, pars$orig, log = TRUE)
  log_res <- ifelse(y == 0, rep(log(za), length.out = length(y)), log(1 - za) + log_f_orig)

  if (log) log_res else exp(log_res)
}

#' @title Zero-Adjusted Continuous Cumulative Distribution Function
#' @name distrib_cdf.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Computes \eqn{F_Y(q) = (1-\pi)F_W(q;\theta) + \pi\,\mathbb{I}(q \ge 0)}, the
#' parent's own distribution function scaled by \eqn{1-\pi} with a JUMP of size
#' \eqn{\pi} added at zero. That jump is the whole of what makes the law mixed,
#' and its height is exactly the probability
#' [distrib_atoms.ZeroAdjustedContinuousDistrib()] reports.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param q A numeric vector of quantiles. The function is right continuous at
#'   zero, so `q = 0` includes the atom.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm is returned.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities, in \eqn{[0, 1]}.
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_pdf.ZeroAdjustedContinuousDistrib()] for the density,
#'   [distrib_quantile.ZeroAdjustedContinuousDistrib()], which inverts this
#'   across the jump, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' distrib_cdf(d, c(-1, 0, 2), theta)
#'
#' # The jump at zero is the atom's probability.
#' distrib_cdf(d, 0, theta) - distrib_cdf(d, -1e-9, theta)
#'
#' # Away from the atom it is the parent's, scaled and shifted.
#' all.equal(distrib_cdf(d, 2, theta), 0.7 * pnorm(2, 1, 2) + 0.3)
S7::method(distrib_cdf, ZeroAdjustedContinuousDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix

  res <- (1 - za) * distrib_cdf(distrib@parent_distrib, q, pars$orig) + za * (q >= 0)
  res <- pmin(pmax(res, 0), 1)

  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Zero-Adjusted Continuous Quantile Function
#' @name distrib_quantile.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Inverts the mixed distribution function across the jump of size \eqn{\pi} at
#' zero. A probability falling INSIDE the jump, between
#' \eqn{(1-\pi)F_W(0)} and \eqn{(1-\pi)F_W(0) + \pi}, returns exactly `0`;
#' below it the parent's quantile at \eqn{p/(1-\pi)}, and above it the parent's
#' at \eqn{(p - \pi)/(1-\pi)}.
#'
#' @details
#' A whole interval of probabilities therefore maps to the same point, which is
#' what inverting a distribution function with a jump means and is not a
#' failure of the inversion. The interval has width \eqn{\pi}.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param p A numeric vector of probabilities, clamped to \eqn{[0, 1]} after
#'   the `log.p` and `lower.tail` transformations are applied.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles, of the recycled length of `p` and
#'   `theta`.
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_cdf.ZeroAdjustedContinuousDistrib()], which this inverts,
#'   and [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' distrib_quantile(d, c(0.1, 0.3, 0.5, 0.9), theta)
#'
#' # An interval of probabilities of width pi maps to zero, which is what
#' # inverting across a jump means.
#' lo <- distrib_cdf(d, -1e-9, theta)
#' c(just_below = distrib_quantile(d, lo - 1e-6, theta),
#'   inside = distrib_quantile(d, lo + 0.15, theta),
#'   just_above = distrib_quantile(d, lo + 0.3 + 1e-6, theta))
#'
#' # Away from the jump the round trip closes.
#' p <- c(0.05, 0.95)
#' max(abs(distrib_cdf(d, distrib_quantile(d, p, theta), theta) - p))
S7::method(distrib_quantile, ZeroAdjustedContinuousDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  pars <- split_mix_theta(distrib, theta)
  parent <- distrib@parent_distrib

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  p <- pmin(pmax(p, 0), 1)

  all_params <- expand_params(c(list(.p = p), theta))
  p <- all_params$.p
  th_orig <- all_params[names(pars$orig)]
  za <- all_params[[distrib@params[distrib@n_params]]]

  F0_orig <- distrib_cdf(parent, 0, th_orig)
  p_lower <- (1 - za) * F0_orig
  p_upper <- p_lower + za

  q_vals <- numeric(length(p))

  # Left of the jump (possible when the parent has mass below 0)
  idx_left <- (p < p_lower)
  if (any(idx_left)) {
    za_sub <- za[idx_left]
    th_sub <- lapply(th_orig, function(x) if (length(x) > 1) x[idx_left] else x)
    q_vals[idx_left] <- distrib_quantile(parent, p[idx_left] / (1 - za_sub), th_sub)
  }

  # Inside the jump: quantile is 0 (already initialized)

  # Right of the jump
  idx_right <- (p > p_upper)
  if (any(idx_right)) {
    za_sub <- za[idx_right]
    p_trans <- pmin((p[idx_right] - za_sub) / (1 - za_sub), 1)
    th_sub <- lapply(th_orig, function(x) if (length(x) > 1) x[idx_right] else x)
    q_vals[idx_right] <- distrib_quantile(parent, p_trans, th_sub)
  }

  q_vals
}

#' @title Zero-Adjusted Continuous Random Number Generator
#' @name distrib_rng.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Draws zeros with probability \eqn{\pi} and otherwise samples from the parent
#' unchanged. No truncation is involved, unlike the discrete branch: a
#' continuous parent produces an exact zero with probability zero, so a draw
#' from it never collides with the atom.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param n The number of draws, a single non-negative whole number.
#' @param theta A named list with the parent's parameters followed by `za`.
#'
#' @return A numeric vector of length `n`, in which the value `0` appears with
#'   probability \eqn{\pi} and is exact.
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_pdf.ZeroAdjustedContinuousDistrib()] for the law these
#'   are drawn from, [distrib_rng.ZeroAdjustedDiscreteDistrib()], which does
#'   truncate, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' set.seed(3)
#' distrib_rng(d, 8, theta)
#'
#' # A large sample reproduces the atom, and the non-zero draws reproduce the
#' # parent.
#' set.seed(3)
#' big <- distrib_rng(d, 20000, theta)
#' c(sampled_atom = mean(big == 0), parameter = 0.3)
#' round(c(mean = mean(big[big != 0]), sd = sd(big[big != 0])), 2)
S7::method(distrib_rng, ZeroAdjustedContinuousDistrib) <- function(distrib, n, theta) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix

  is_zero <- stats::runif(n) < za
  y <- numeric(n)

  if (any(!is_zero)) {
    th_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[!is_zero] else x)
    y[!is_zero] <- distrib_rng(distrib@parent_distrib, sum(!is_zero), th_sub)
  }
  y
}

#' @title Zero-Adjusted Continuous Score
#' @name distrib_gradient.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Computes the first derivatives of the mixed log-density in closed form. The
#' likelihood separates COMPLETELY here, with no truncation term to correct
#' for: the parent's parameters take the parent's own score at a non-zero
#' observation and zero at the atom, and
#' \deqn{\frac{\partial \ell}{\partial \pi}
#'   = \frac{\mathbb{I}(y = 0)}{\pi} - \frac{\mathbb{I}(y \ne 0)}{1 - \pi},}
#' the Bernoulli score of the indicator that the observation is zero.
#'
#' @details
#' Compare the discrete branch, whose parent block carries the truncation
#' correction \eqn{f(0)s(0)/\{1 - f(0)\}}. There is none here because a
#' continuous parent places no mass at zero, so nothing is removed from it.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, in
#'   `distrib@params` order.
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_hessian.ZeroAdjustedContinuousDistrib()] for the second
#'   order, [distrib_gradient.ZeroAdjustedDiscreteDistrib()], which carries the
#'   truncation term, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#' set.seed(4)
#' y <- distrib_rng(d, 300, theta)
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
#' # Away from the atom it is the parent's own score, with no correction.
#' nz <- which(y != 0)[1]
#' c(mixed = g$mu[nz],
#'   parent = distrib_gradient(gaussian1_distrib(), y[nz],
#'                             theta[c("mu", "sigma")])$mu)
#'
#' # And an observation at the atom says nothing about the parent.
#' unique(g$mu[y == 0])
S7::method(distrib_gradient, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]

  grad_orig <- distrib_gradient(distrib@parent_distrib, y, pars$orig)

  res <- lapply(grad_orig, function(g) ifelse(y == 0, 0, g))
  res[[za_name]] <- ifelse(y == 0, 1 / za, -1 / (1 - za))
  res
}

#' @title Zero-Adjusted Continuous Observed Hessian
#' @name distrib_hessian.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Computes the second derivatives of the mixed log-density in closed form. The
#' MIXED BLOCKS ARE EXACTLY ZERO, the likelihood factorizing into a binary part
#' and a continuous part. The parent block is the parent's own observed Hessian
#' at a non-zero observation and zero at the atom, with no truncation
#' correction; the hurdle block is \eqn{-1/\pi^2} at the atom and
#' \eqn{-1/(1-\pi)^2} elsewhere.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names]. Every
#'   key pairing a parent parameter with `za` holds a vector of exact zeros.
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_gradient.ZeroAdjustedContinuousDistrib()] for the first
#'   order, [distrib_expected_hessian.ZeroAdjustedContinuousDistrib()] for the
#'   expectation, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#' set.seed(4)
#' y <- distrib_rng(d, 300, theta)
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
#' # Every mixed block is exactly zero.
#' c(mu_za = all(H$mu_za == 0), sigma_za = all(H$sigma_za == 0))
S7::method(distrib_hessian, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]
  n <- length(y)

  res <- list()
  res[[paste0(za_name, "_", za_name)]] <- ifelse(y == 0, -1 / (za^2), -1 / ((1 - za)^2))

  for (nm in names(pars$orig)) {
    res[[paste0(nm, "_", za_name)]] <- rep(0, n)
  }

  h_orig <- distrib_hessian(distrib@parent_distrib, y, pars$orig)
  for (nm in names(h_orig)) {
    res[[nm]] <- ifelse(y == 0, 0, h_orig[[nm]])
  }

  expand_params(res[hess_names(distrib@params)], n)
}

#' @title Zero-Adjusted Continuous Expected Information
#' @name distrib_expected_hessian.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Computes the expectation of the observed Hessian in closed form:
#' \deqn{\mathbb{E}[H_{\pi\pi}] = -\frac{1}{\pi(1-\pi)}, \qquad
#'   \mathbb{E}[H_{\theta\pi}] = 0, \qquad
#'   \mathbb{E}[H_{\theta\theta}] = (1-\pi)\,\mathbb{E}[H_W].}
#' The hurdle block is the Bernoulli information, the mixed blocks are exactly
#' zero, and the parent block is the parent's own expectation weighted by the
#' probability \eqn{1-\pi} that an observation reaches it.
#'
#' @details
#' `approx` and `nsim` are accepted so that the signature matches the generic's
#' and neither is read here. They are also NOT forwarded to the parent, so a
#' parent whose own expected Hessian is approximate takes its own defaults.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations, passed to the parent, which uses
#'   its length.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored for the hurdle and mixed blocks. It is also not
#'   forwarded to the parent, which takes its own default. Present so that the
#'   signature matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length `length(y)`, keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{f_W} is the parent's density, \eqn{\pi} the probability of the atom at
#' zero, \eqn{f_Y} the mixed density and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_hessian.ZeroAdjustedContinuousDistrib()] for the observed
#'   matrix, [distrib_expected_hessian.ZeroAdjustedDiscreteDistrib()] for the
#'   discrete branch, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#' set.seed(4)
#' y <- distrib_rng(d, 300, theta)
#'
#' EH <- distrib_expected_hessian(d, y, theta)
#' vapply(EH, function(z) z[1], numeric(1))
#'
#' # The hurdle block is the Bernoulli information and the mixed ones vanish.
#' c(reported = EH$za_za[1], bernoulli = -1 / (0.3 * 0.7),
#'   mixed = EH$mu_za[1])
#'
#' # The parent block is the parent's own, weighted by 1 - pi.
#' c(reported = EH$mu_mu[1], weighted_parent = 0.7 * (-1 / 2^2))
S7::method(distrib_expected_hessian, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]
  n <- length(y)

  res <- list()
  res[[paste0(za_name, "_", za_name)]] <- rep(-1 / (za * (1 - za)), length.out = n)

  for (nm in names(pars$orig)) {
    res[[paste0(nm, "_", za_name)]] <- rep(0, n)
  }

  h_exp_orig <- distrib_expected_hessian(distrib@parent_distrib, y, pars$orig)
  for (nm in names(h_exp_orig)) {
    res[[nm]] <- (1 - za) * h_exp_orig[[nm]]
  }

  expand_params(res[hess_names(distrib@params)], n)
}

#' @title Atoms of a Zero-Adjusted Continuous Distribution
#' @name distrib_atoms.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Declares the single point mass at zero, with probability \eqn{\pi}. This
#' declaration is what marks the object a MIXED distribution: its density
#' integrates to \eqn{1 - \pi}, so a consumer that does not know about the atom
#' will find the missing \eqn{\pi} and call it an error.
#'
#' @details
#' Three things read it. [check_distrib()] adjusts four of its thirteen checks,
#' since the density integrates to \eqn{1-\pi}, the quantile round trip cannot
#' close inside the jump, and a central difference in \eqn{y} across the atom
#' has no derivative to converge to. [expectation()] splits its integral at the
#' atom. And [folded()] REFUSES a parent that declares one, zero being its own
#' preimage under the absolute value while every other point has two.
#'
#' The question is asked of the OBJECT, so the returned probability moves with
#' `theta`.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param theta A named list with the parent's parameters followed by `za`.
#'   Only `za` is read, and only its first element.
#'
#' @return A named list with `y`, the numeric vector `0`, and `p`, the
#'   probability at it.
#'
#' @section Notation:
#' \eqn{\pi} is the probability of the atom at zero.
#'
#' @seealso [distrib_atoms()] for the generic, [check_distrib()] and
#'   [expectation.ZeroAdjustedContinuousDistrib()] for the two consumers, and
#'   [folded()], which rejects a parent that declares one.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' distrib_atoms(d, theta)
#'
#' # A plain continuous parent declares none.
#' distrib_atoms(gaussian1_distrib(), theta[c("mu", "sigma")])
#'
#' # The declaration is what folded() refuses on.
#' try(folded(d))
S7::method(distrib_atoms, ZeroAdjustedContinuousDistrib) <- function(distrib, theta) {
  list(y = 0, p = unname(theta[[distrib@n_params]][1]))
}

#' @title Response Derivative of a Zero-Adjusted Distribution
#'
#' @description
#' Evaluates a response derivative of the PARENT away from the atom and returns
#' `NaN` at it. Away from zero the \eqn{1-\pi} factor is a constant in \eqn{y},
#' so the parent's own derivative is exact and needs no correction; at zero the
#' log-density jumps, \eqn{\log\pi} on one side and
#' \eqn{\log\{(1-\pi)f(y)\}} on the other, and no derivative exists.
#'
#' @details
#' The `NaN` is the point of the function. The finite-difference default
#' inherited from [continuous_distrib()] would straddle the jump and return a
#' NUMBER for it, which is worse than an error: nothing downstream would notice.
#'
#' @param distrib A zero-adjusted distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, the parent's followed by the atom
#'   probability.
#' @param fun The parent's response-derivative function, such as
#'   [distrib_grad_y()] or [distrib_hess_y()]. It is called on the parent with
#'   the non-zero observations and the parent's parameters, subset to those
#'   observations where a parameter varies.
#'
#' @return A numeric vector as long as `y`, `NaN` wherever `y == 0` and the
#'   parent's derivative elsewhere.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{\pi} the probability of the atom and
#' \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_grad_y.ZeroAdjustedContinuousDistrib()] and
#'   [distrib_hess_y.ZeroAdjustedContinuousDistrib()], the two methods it
#'   serves.
#'
#' @keywords internal
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' # NaN at the atom, the parent's derivative elsewhere.
#' distributions7:::za_y_deriv(d, c(-1, 0, 2), theta, distrib_grad_y)
#' distrib_grad_y(gaussian1_distrib(), c(-1, 2), theta[c("mu", "sigma")])
#'
#' # The second order takes the same route.
#' distributions7:::za_y_deriv(d, c(0, 2), theta, distrib_hess_y)
za_y_deriv <- function(distrib, y, theta, fun) {
  pars <- split_mix_theta(distrib, theta)
  out <- rep(NaN, length(y))
  nz <- y != 0
  if (any(nz)) {
    th_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[nz] else x)
    out[nz] <- fun(distrib@parent_distrib, y[nz], th_sub)
  }
  out
}

#' @title Zero-Adjusted Continuous Response Gradient
#' @name distrib_grad_y.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Returns \eqn{\partial\ell/\partial y}, which equals the PARENT's for
#' \eqn{y \ne 0}, the factor \eqn{1-\pi} not depending on \eqn{y}. At the atom
#' the log-density jumps and no derivative exists, so `NaN` is returned. That
#' is deliberate: the numerical fallback would straddle the jump and hand back
#' a number nothing downstream would question.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations. Exactly zero gives `NaN`.
#' @param theta A named list with the parent's parameters followed by `za`.
#'
#' @return A numeric vector as long as `y`.
#'
#' @section Notation:
#' \eqn{\pi} is the probability of the atom at zero and \eqn{\ell} the
#' log-density of one observation.
#'
#' @seealso [distrib_hess_y.ZeroAdjustedContinuousDistrib()] for the second
#'   order, [za_y_deriv()] for the shared body, and [distrib_grad_y()] for the
#'   generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' distrib_grad_y(d, c(-1, 0, 2), theta)
#'
#' # Away from the atom it is the parent's own.
#' distrib_grad_y(gaussian1_distrib(), c(-1, 2), theta[c("mu", "sigma")])
S7::method(distrib_grad_y, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta) {
  za_y_deriv(distrib, y, theta, distrib_grad_y)
}

#' @title Zero-Adjusted Continuous Response Hessian
#' @name distrib_hess_y.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Returns \eqn{\partial^2\ell/\partial y^2}, which equals the PARENT's for
#' \eqn{y \ne 0} and is `NaN` at the atom, for the reason
#' [distrib_grad_y.ZeroAdjustedContinuousDistrib()] gives: the log-density
#' jumps there and no derivative exists.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param y A numeric vector of observations. Exactly zero gives `NaN`.
#' @param theta A named list with the parent's parameters followed by `za`.
#'
#' @return A numeric vector as long as `y`.
#'
#' @section Notation:
#' \eqn{\pi} is the probability of the atom at zero and \eqn{\ell} the
#' log-density of one observation.
#'
#' @seealso [distrib_grad_y.ZeroAdjustedContinuousDistrib()] for the first
#'   order, [za_y_deriv()] for the shared body, and [distrib_hess_y()] for the
#'   generic.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' distrib_hess_y(d, c(-1, 0, 2), theta)
#'
#' # Away from the atom it is the parent's, which for a gaussian is
#' # -1 / sigma^2 everywhere.
#' c(mixed = distrib_hess_y(d, 2, theta), parent = -1 / 2^2)
S7::method(distrib_hess_y, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta) {
  za_y_deriv(distrib, y, theta, distrib_hess_y)
}

#' @title Expectation for Zero-Adjusted Continuous Distributions
#' @name expectation.ZeroAdjustedContinuousDistrib
#'
#' @description
#' Computes \eqn{\mathbb{E}[f(Y)]} by splitting the expectation at the atom,
#' \deqn{\mathbb{E}[f(Y)] = \pi\, f(0) + (1-\pi)\, \mathbb{E}_W[f(W)],}
#' with the second term the PARENT's own expectation. Plain numerical
#' integration over the mixed density would miss the point mass entirely and
#' return \eqn{(1-\pi)\mathbb{E}_W[f]}, silently.
#'
#' @details
#' `f` receives the FULL `theta`, `za` included, even though the integral runs
#' over the parent. The method re-attaches the atom probability inside a
#' wrapper, so a caller's `f` sees the same parameter list at the atom and away
#' from it.
#'
#' @param distrib A `ZeroAdjustedContinuousDistrib` object, from
#'   [zero_adjusted()].
#' @param f A function `f(y, theta, ...)`, receiving the full `theta` including
#'   `za`. It must be vectorized in `y`, as the parent's own [expectation()]
#'   requires.
#' @param theta A named list with the parent's parameters followed by `za`.
#' @param ... Passed to `f` and to the parent's [expectation()].
#'
#' @return A single number, the expectation of `f` under the mixed law.
#'
#' @section Notation:
#' \eqn{\pi} is the probability of the atom at zero, \eqn{W} the parent's
#' variable and \eqn{Y} the zero-adjusted one.
#'
#' @seealso [expectation()] for the generic,
#'   [distrib_atoms.ZeroAdjustedContinuousDistrib()] for the declaration this
#'   rests on, and [base::mean()] and [variance()], the two moments built on
#'   it.
#'
#' @examples
#' d <- zero_adjusted(gaussian1_distrib())
#' theta <- list(mu = 1, sigma = 2, za = 0.3)
#'
#' # The first two moments, against the split written out.
#' c(computed = expectation(d, function(y, theta) y, theta),
#'   theory = 0.7 * 1)
#' c(computed = expectation(d, function(y, theta) y^2, theta),
#'   theory = 0.7 * (1^2 + 2^2))
#'
#' # mean() and variance() are built on the same split.
#' c(mean = mean(d, theta), variance = variance(d, theta))
#' c(theory_var = 0.7 * 5 - (0.7 * 1)^2)
#'
#' # The parent's own expectation misses the atom, and is the wrong answer for
#' # the mixed law by exactly the factor 1 - pi.
#' expectation(gaussian1_distrib(), function(y, theta) y,
#'             theta[c("mu", "sigma")])
S7::method(expectation, ZeroAdjustedContinuousDistrib) <- function(distrib, f, theta, ...) {
  pars <- split_mix_theta(distrib, theta)
  za_name <- distrib@params[distrib@n_params]
  za <- pars$mix

  # Wrap f so that it receives the full theta (with za re-attached), while the
  # integral itself runs over the parent's distribution.
  g <- function(y, theta, .za, ...) {
    full_theta <- c(theta, stats::setNames(list(.za), za_name))
    f(y, full_theta, ...)
  }

  e_parent <- expectation(distrib@parent_distrib, g, pars$orig, .za = za, ...)
  f0 <- f(0, theta, ...)

  za * f0 + (1 - za) * e_parent
}

# --- CONSTRUCTOR WRAPPERS ---

#' @title Zero-Adjusted Distribution Object
#'
#' @description
#' Makes the probability of a zero a parameter of its own, \eqn{\pi}, carried
#' by `za`, and leaves everything else to the parent. What that means depends
#' on the parent's type, and the constructor dispatches on it:
#'
#' - a DISCRETE parent with 0 in its support gives a HURDLE model. The mass the
#'   parent puts at zero is removed, the parent is renormalized over the
#'   positive values, and \eqn{\pi} takes its place.
#' - a CONTINUOUS parent gives a MIXED distribution. Nothing has to be removed,
#'   \eqn{P(Y = 0)} being zero already; a point mass \eqn{\pi} is placed at zero
#'   and the density is scaled by \eqn{1-\pi}.
#'
#' It is the wrapper to reach for when zeros come from their own mechanism, no
#' claim filed or no rainfall, and the parent describes only what happens once
#' that mechanism has been passed. Where the zeros come partly from the parent
#' and partly from a separate source, so that no single zero can be attributed,
#' [zero_inflated()] is the model.
#'
#' @details
#' # A discrete parent: the hurdle
#'
#' \deqn{P(Y = y; \theta, \pi) = \begin{cases} \pi & y = 0 \\
#'   (1 - \pi)\dfrac{f(y; \theta)}{1 - f(0; \theta)} & y > 0. \end{cases}}
#' The division by \eqn{1 - f(0;\theta)} is the TRUNCATION, and it is what
#' separates this from zero-inflation: it makes \eqn{\theta} the parameters of
#' a law on the positive integers instead of the original count process. The
#' log-likelihood separates into a Bernoulli part in \eqn{\pi} and a truncated
#' part in \eqn{\theta}, so every mixed block of the information is exactly
#' zero and the two halves could be fitted separately.
#'
#' # A continuous parent: the mixed law
#'
#' \deqn{f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi) f_W(y;\theta) \quad (y \ne 0).}
#' Here \eqn{f_Y} is a density against Lebesgue measure plus a point mass, so
#' it integrates to \eqn{1 - \pi}; the remainder is the atom, which
#' [distrib_atoms()] reports and [check_distrib()] accounts for. The classical
#' members have a parent on \eqn{(0, \infty)} or \eqn{(0, 1)}, where zero sits
#' on the boundary: the zero-adjusted gamma, inverse gaussian or lognormal for
#' semicontinuous data, and the zero-adjusted beta for proportions.
#'
#' # Choosing between the two wrappers
#'
#' Zero-inflation can only ADD zeros, \eqn{P(Y = 0) = \zeta + (1-\zeta)f(0)}
#' exceeding \eqn{f(0)}; the hurdle replaces \eqn{f(0)} outright and so also
#' covers FEWER zeros than the parent implies. Where both apply they are not
#' nested, and they differ in interpretation more than in fit: zero-inflation
#' keeps \eqn{\theta} as the parameters of the original count process, the
#' hurdle re-reads them as those of a truncated one. Prefer the hurdle when a
#' zero is observable evidence of a distinct decision, and zero-inflation when
#' the two kinds of zero are genuinely indistinguishable.
#'
#' # What the parent must be
#'
#' A discrete parent must have 0 in its support: with \eqn{f(0) = 0} there is
#' no mass to remove. Construction also fails where the result would not be
#' identified:
#'
#' - the parent already models a probability of zero. Zero-truncating a
#'   zero-inflated or zero-adjusted parent cancels its zero parameter out of the
#'   likelihood entirely, leaving an identically zero score;
#' - the support is too small to carry one more parameter. A distribution on
#'   \eqn{k} points has \eqn{k-1} free probabilities, so at least `n_params + 2`
#'   support points are needed. Zero-adjusting a Bernoulli leaves the truncated
#'   part concentrated on \eqn{\{1\}} and `mu` vanishes from the likelihood.
#'
#' A continuous parent supported on the whole line is accepted and gives a
#' spike-at-zero model, but \eqn{y = 0} then no longer identifies its own
#' mechanism: the atom sits where the density is positive. A continuous parent
#' whose support does not reach zero at all, say \eqn{(2, 5)}, is accepted with
#' a WARNING, the atom then being disconnected from the rest of the law, which
#' is legitimate and rarely intended.
#'
#' @section Notation:
#' \eqn{f} is a discrete parent's mass function, \eqn{f_W} a continuous one's
#' density, \eqn{\pi} the probability of a zero, \eqn{\zeta} the inflation
#' probability of the other wrapper, and \eqn{k} the number of support points.
#'
#' @param distrib An object inheriting from `discrete_distrib` with 0 in its
#'   support, or from `continuous_distrib`. An already-wrapped parent and a
#'   discrete support of fewer than `n_params + 2` points are rejected with an
#'   error saying which condition failed.
#' @param link_za The link carrying \eqn{\pi} to the unconstrained scale, a
#'   `linkfunctions7::link` object. Defaults to
#'   [linkfunctions7::logit_link()], which keeps it strictly inside
#'   \eqn{(0, 1)} at every point of the free scale.
#'
#' @return An S7 object of class [ZeroAdjustedDiscreteDistrib] or
#'   [ZeroAdjustedContinuousDistrib], matching the parent's branch. Its
#'   `params` are the parent's followed by `za`, whose bound is \eqn{(0, 1)}
#'   and whose interpretation is `"prob. of zero"`; `distrib_name` is
#'   `"zero-adjusted "` followed by the parent's.
#'
#' @seealso [zero_inflated()] for the mixture counterpart,
#'   [ZeroAdjustedDiscreteDistrib] and [ZeroAdjustedContinuousDistrib] for the
#'   two classes, [distrib_atoms()] for the atom the continuous branch
#'   declares, and [check_distrib()] to validate the result.
#'
#' @examples
#' # Hurdle Poisson: the mass at zero is exactly za, not dpois(0, mu).
#' zap <- zero_adjusted(poisson_distrib())
#' theta <- list(mu = 3, za = 0.3)
#' distrib_pdf(zap, 0:5, theta)
#' c(at_zero = distrib_pdf(zap, 0, theta), parent = dpois(0, 3))
#'
#' # And it can be BELOW the parent's, which zero-inflation cannot reach.
#' c(hurdle = distrib_pdf(zap, 0, list(mu = 3, za = 0.01)),
#'   parent = dpois(0, 3))
#'
#' # The likelihood separates, so every mixed block of the Hessian is zero.
#' set.seed(2)
#' y <- distrib_rng(zap, 300, theta)
#' all(distrib_hessian(zap, y, theta)$mu_za == 0)
#'
#' # Semicontinuous data: a spike at zero and a gamma above it.
#' zag <- zero_adjusted(gamma2_distrib())
#' distrib_atoms(zag, list(mu = 2, sigma2 = 1, za = 0.3))
#'
#' # A fit recovers both halves, each from its own part of the data.
#' set.seed(5)
#' yg <- distrib_rng(zag, 2000, list(mu = 2, sigma2 = 1, za = 0.3))
#' round(coef(fit_distrib(zag, yg)), 3)
#'
#' # Two refusals, each naming the condition that failed.
#' try(zero_adjusted(bernoulli_distrib()))
#' try(zero_adjusted(zero_inflated(poisson_distrib())))
#'
#' @importFrom linkfunctions7 logit_link
#' @export
zero_adjusted <- function(distrib, link_za = logit_link()) {
  if (!S7::S7_inherits(distrib, discrete_distrib) &&
      !S7::S7_inherits(distrib, continuous_distrib)) {
    stop("Input must inherit from 'discrete_distrib' or 'continuous_distrib'.", call. = FALSE)
  }
  check_not_stacked(distrib, "zero_adjusted", "za")

  common <- list(
    parent_distrib = distrib,
    distrib_name = paste0("zero-adjusted ", distrib@distrib_name),
    dimension = distrib@dimension,
    params = c(distrib@params, "za"),
    params_interpretation = c(distrib@params_interpretation, za = "prob. of zero"),
    n_params = distrib@n_params + 1,
    params_bounds = c(distrib@params_bounds, list(za = c(0, 1))),
    link_params = c(distrib@link_params, list(za = link_za)),
    params_smooth = c(param_smoothness(distrib), za = TRUE)
  )

  if (S7::S7_inherits(distrib, discrete_distrib)) {
    if (distrib@bounds[1] > 0) {
      stop(sprintf(paste0(
        "zero_adjusted() requires 0 in the support of '%s', which starts at %g.\n",
        "  A hurdle model removes the mass the parent places at zero and replaces it\n",
        "  by 'za'; with P(Y = 0) = 0 there is nothing to remove or to truncate."
      ), distrib@distrib_name, distrib@bounds[1]), call. = FALSE)
    }
    check_support_is_rich_enough(distrib, "zero_adjusted")
    do.call(ZeroAdjustedDiscreteDistrib, c(common, list(bounds = c(0, distrib@bounds[2]))))
  } else {
    if (distrib@bounds[1] > 0) {
      warning(sprintf(paste0(
        "The support of '%s' starts at %g, so the point mass at zero is disconnected ",
        "from the rest of the distribution. This is well-defined but rarely intended: ",
        "zero-adjustment is normally applied to a density that reaches down to zero."
      ), distrib@distrib_name, distrib@bounds[1]), call. = FALSE)
    }
    do.call(ZeroAdjustedContinuousDistrib, c(common, list(bounds = c(min(0, distrib@bounds[1]), distrib@bounds[2]))))
  }
}
