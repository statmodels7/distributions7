#' @include distrib.R generics.R utility_functions.R numerical_functions.R zero_inflated.R

#' @title S7 Class for Zero-Adjusted Discrete (Hurdle) Distributions
#' @name ZeroAdjustedDiscreteDistrib
#'
#' @description
#' A subclass of \code{discrete_distrib} representing the zero-adjusted (hurdle) version
#' of a wrapped discrete distribution: the probability mass at zero is \emph{replaced} by
#' \eqn{\pi} (parameter \code{za}), and positive values follow the zero-truncated parent.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{discrete_distrib} object.
#' @seealso \code{\link{zero_adjusted}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.ZeroAdjustedDiscreteDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.ZeroAdjustedDiscreteDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.ZeroAdjustedDiscreteDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.ZeroAdjustedDiscreteDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.ZeroAdjustedDiscreteDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.ZeroAdjustedDiscreteDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.ZeroAdjustedDiscreteDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{discrete_distrib}}.
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
#' A subclass of \code{continuous_distrib} representing the zero-adjusted version of a
#' wrapped continuous distribution: a point mass at zero with probability \eqn{\pi}
#' (parameter \code{za}) mixed with the original continuous distribution.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{continuous_distrib} object.
#' @seealso \code{\link{zero_adjusted}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.ZeroAdjustedContinuousDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.ZeroAdjustedContinuousDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.ZeroAdjustedContinuousDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.ZeroAdjustedContinuousDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.ZeroAdjustedContinuousDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.ZeroAdjustedContinuousDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.ZeroAdjustedContinuousDistrib]{distrib_rng()}},
#'   \code{\link[=expectation.ZeroAdjustedContinuousDistrib]{expectation()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
ZeroAdjustedContinuousDistrib <- S7::new_class("ZeroAdjustedContinuousDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib
  )
)

# ============================== DISCRETE (HURDLE) ==============================

#' @title Zero-Adjusted Discrete Probability Mass Function
#' @name distrib_pdf.ZeroAdjustedDiscreteDistrib
#' @description
#' \deqn{P(Y=0) = \pi, \qquad P(Y=y) = (1-\pi)\dfrac{f(y;\theta)}{1-f(0;\theta)} \ \ (y>0)}
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_pdf, ZeroAdjustedDiscreteDistrib) <- function(distrib, y, theta, log = FALSE) {
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
#' @description
#' \deqn{F_{ZA}(q) = \pi + (1-\pi)\dfrac{F(q;\theta) - f(0;\theta)}{1 - f(0;\theta)} \quad (q \ge 0)}
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @seealso \code{\link{zero_adjusted}}
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
#' @description
#' Inverts the hurdle CDF: 0 for \eqn{p \le \pi}, otherwise the parent quantile at
#' \eqn{u(1-f(0)) + f(0)} with \eqn{u = (p-\pi)/(1-\pi)}.
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @seealso \code{\link{zero_adjusted}}
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
#' @description
#' Draws zeros with probability \eqn{\pi} and otherwise samples from the zero-truncated
#' parent via inverse transform sampling.
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @seealso \code{\link{zero_adjusted}}
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

#' @title Zero-Adjusted Discrete Analytical Gradient
#' @name distrib_gradient.ZeroAdjustedDiscreteDistrib
#' @description
#' The hurdle likelihood separates: for the parent's parameters the score at \eqn{y>0} is
#' the parent's score plus the truncation correction \eqn{\dfrac{f(0)}{1-f(0)} S(0)}
#' (and 0 at \eqn{y=0}); for \eqn{\pi}:
#' \deqn{\dfrac{\partial \ell}{\partial \pi} = \mathbb{I}(y=0)\dfrac{1}{\pi} - \mathbb{I}(y>0)\dfrac{1}{1-\pi}}
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{zero_adjusted}}
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

#' @title Zero-Adjusted Discrete Analytical Observed Hessian
#' @name distrib_hessian.ZeroAdjustedDiscreteDistrib
#' @description
#' Observed Hessian of the hurdle model. The mixed blocks are identically zero;
#' the truncation adds the correction
#' \eqn{H_{corr} = \dfrac{(1-f(0))f''(0) + f'(0)^2}{(1-f(0))^2}} for \eqn{y > 0}.
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{zero_adjusted}}
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

#' @title Zero-Adjusted Discrete Analytical Expected Hessian
#' @name distrib_expected_hessian.ZeroAdjustedDiscreteDistrib
#' @description
#' Expected Hessian of the hurdle model:
#' \eqn{E[H_{\pi\pi}] = -\dfrac{1}{\pi(1-\pi)}}, mixed blocks are 0, and
#' \eqn{E[H_{\theta\theta}] = (1-\pi)\left(\dfrac{E[H] - f(0)H(0)}{1-f(0)} + H_{corr}\right)}.
#' @param distrib A \code{ZeroAdjustedDiscreteDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{zero_adjusted}}
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
#' @description
#' \deqn{f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi) f_W(y;\theta) \ \ (y \neq 0)}
#' (mixed density: point mass at 0 plus scaled continuous density).
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_pdf, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta, log = FALSE) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix

  log_f_orig <- distrib_pdf(distrib@parent_distrib, y, pars$orig, log = TRUE)
  log_res <- ifelse(y == 0, rep(log(za), length.out = length(y)), log(1 - za) + log_f_orig)

  if (log) log_res else exp(log_res)
}

#' @title Zero-Adjusted Continuous Cumulative Distribution Function
#' @name distrib_cdf.ZeroAdjustedContinuousDistrib
#' @description
#' \deqn{F_Y(q) = (1-\pi)F_W(q;\theta) + \pi\,\mathbb{I}(q \ge 0)}
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @seealso \code{\link{zero_adjusted}}
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
#' @description
#' Inverts the mixed CDF, handling the jump of size \eqn{\pi} at 0: quantiles falling in
#' the jump are 0; on either side the parent's quantile function is used on the rescaled
#' probability.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @seealso \code{\link{zero_adjusted}}
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
#' @description
#' Draws zeros with probability \eqn{\pi} and otherwise samples from the parent distribution.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @seealso \code{\link{zero_adjusted}}
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

#' @title Zero-Adjusted Continuous Analytical Gradient
#' @name distrib_gradient.ZeroAdjustedContinuousDistrib
#' @description
#' The likelihood separates completely: for the parent's parameters the score is the
#' parent's score at \eqn{y \neq 0} and 0 at \eqn{y=0}; for \eqn{\pi}:
#' \deqn{\dfrac{\partial \ell}{\partial \pi} = \mathbb{I}(y=0)\dfrac{1}{\pi} - \mathbb{I}(y \neq 0)\dfrac{1}{1-\pi}}
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_gradient, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pars <- split_mix_theta(distrib, theta)
  za <- pars$mix
  za_name <- distrib@params[distrib@n_params]

  grad_orig <- distrib_gradient(distrib@parent_distrib, y, pars$orig)

  res <- lapply(grad_orig, function(g) ifelse(y == 0, 0, g))
  res[[za_name]] <- ifelse(y == 0, 1 / za, -1 / (1 - za))
  res
}

#' @title Zero-Adjusted Continuous Analytical Observed Hessian
#' @name distrib_hessian.ZeroAdjustedContinuousDistrib
#' @description
#' Observed Hessian of the zero-adjusted continuous model: the mixed blocks are 0, the
#' \eqn{\pi\pi} block is \eqn{-1/\pi^2} at \eqn{y=0} and \eqn{-1/(1-\pi)^2} otherwise,
#' and the \eqn{\theta\theta} block is the parent's Hessian at \eqn{y \neq 0} (0 at \eqn{y=0}).
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{zero_adjusted}}
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

#' @title Zero-Adjusted Continuous Analytical Expected Hessian
#' @name distrib_expected_hessian.ZeroAdjustedContinuousDistrib
#' @description
#' Expected Hessian: \eqn{E[H_{\pi\pi}] = -\dfrac{1}{\pi(1-\pi)}}, mixed blocks are 0,
#' and \eqn{E[H_{\theta\theta}] = (1-\pi) E[H_W]}.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{zero_adjusted}}
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
#' @description
#' The single point mass at zero, with probability \eqn{\pi}. This is what makes the
#' object a mixed distribution: its density integrates to \eqn{1 - \pi}, not to 1.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A list with \code{y = 0} and \code{p = za}.
#' @seealso \code{\link{zero_adjusted}}, \code{\link{distrib_atoms}}
S7::method(distrib_atoms, ZeroAdjustedContinuousDistrib) <- function(distrib, theta) {
  list(y = 0, p = unname(theta[[distrib@n_params]][1]))
}

#' Response Derivative of a Zero-Adjusted Distribution
#'
#' @description
#' Evaluates a response derivative of the parent away from the atom, and returns
#' \code{NaN} at it.
#'
#' @details
#' The log-density jumps at zero -- \eqn{\log \pi} on one side,
#' \eqn{\log((1-\pi) f(y))} on the other -- so no derivative in \eqn{y} exists
#' there. The finite-difference default inherited from
#' \code{\link{continuous_distrib}} would happily straddle the jump and return a
#' number for it, which is worse than refusing. Away from zero the \eqn{1-\pi}
#' factor is constant in \eqn{y}, so the parent's own derivative is exact and
#' nothing needs correcting.
#'
#' @param distrib A zero-adjusted distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, including the atom probability.
#' @param fun The parent's response-derivative function, e.g.
#'   \code{\link{distrib_grad_y}}.
#'
#' @return A numeric vector as long as \code{y}, \code{NaN} wherever
#'   \code{y == 0}.
#'
#' @seealso \code{\link{zero_adjusted}}
#' @keywords internal
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
#' @description
#' \eqn{\partial \ell / \partial y} equals the parent's for \eqn{y \neq 0}, since the
#' factor \eqn{1-\pi} does not depend on \eqn{y}. At \eqn{y = 0} the log-density jumps
#' to the atom and no derivative exists, so \code{NaN} is returned.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A numeric vector.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_grad_y, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta) {
  za_y_deriv(distrib, y, theta, distrib_grad_y)
}

#' @title Zero-Adjusted Continuous Response Hessian
#' @name distrib_hess_y.ZeroAdjustedContinuousDistrib
#' @description
#' \eqn{\partial^2 \ell / \partial y^2} equals the parent's for \eqn{y \neq 0} and is
#' \code{NaN} at the atom.
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @return A numeric vector.
#' @seealso \code{\link{zero_adjusted}}
S7::method(distrib_hess_y, ZeroAdjustedContinuousDistrib) <- function(distrib, y, theta) {
  za_y_deriv(distrib, y, theta, distrib_hess_y)
}

#' @title Expectation for Zero-Adjusted Continuous Distributions
#' @name expectation.ZeroAdjustedContinuousDistrib
#' @description
#' The zero-adjusted continuous distribution has a point mass at 0 that plain numerical
#' integration would miss. The expectation decomposes exactly as
#' \deqn{E[f(Y)] = \pi f(0) + (1-\pi) E_W[f(W)]}
#' @param distrib A \code{ZeroAdjustedContinuousDistrib} object.
#' @param f A function \code{f(y, theta, ...)} (receives the full theta, including \code{za}).
#' @param theta A list with the parent's parameters followed by \code{za}.
#' @param ... Additional arguments passed to \code{f}.
#' @keywords internal
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

#' Zero-Adjusted Distribution Object
#'
#' @description
#' Creates a zero-adjusted version of an existing distribution: the probability of a
#' zero becomes a parameter of its own, \eqn{\pi} (parameter \code{za}), and everything
#' else is left to the parent. What that means depends on the parent's type, and the
#' constructor dispatches on it:
#' \itemize{
#'   \item \strong{Discrete} (support including 0): a \strong{hurdle} model. The mass the
#'     parent puts at zero is removed, the parent is renormalised over the positive
#'     values, and \eqn{\pi} takes its place.
#'   \item \strong{Continuous}: a \strong{mixed} distribution. Nothing has to be removed,
#'     since \eqn{P(Y = 0) = 0} already; a point mass \eqn{\pi} is placed at zero and the
#'     density is scaled by \eqn{1-\pi}.
#' }
#'
#' Zero-adjustment is the right wrapper when zeros come from their own mechanism ---
#' no claim filed, no purchase made, no rainfall --- and the parent describes only what
#' happens once that mechanism has been passed. When the zeros instead come partly from
#' the parent itself and partly from a separate one, so that no single zero can be
#' attributed, use \code{\link{zero_inflated}}.
#'
#' @param distrib An object inheriting from \code{discrete_distrib} (with 0 in its
#'   support) or from \code{continuous_distrib}.
#' @param link_za A link function object for the zero probability \eqn{\pi}.
#'   Defaults to \code{\link[linkfunctions7]{logit_link}}.
#'
#' @details
#' \strong{Discrete parent (hurdle).}
#' \deqn{
#' P(Y=y; \theta, \pi) =
#' \begin{cases}
#' \pi & y = 0 \\
#' (1 - \pi)\dfrac{f(y; \theta)}{1 - f(0; \theta)} & y > 0
#' \end{cases}
#' }
#' The division by \eqn{1 - f(0;\theta)} is the truncation: it is what distinguishes this
#' from zero-inflation, and what makes \eqn{\theta} the parameters of a law on the
#' positive integers rather than of the original count process. The log-likelihood
#' separates into a Bernoulli part in \eqn{\pi} and a truncated part in \eqn{\theta}, so
#' the mixed blocks of the information matrix are exactly zero and the two halves could
#' in principle be fitted separately.
#'
#' \strong{Continuous parent (mixed).}
#' \deqn{f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi) f_W(y;\theta) \quad (y \neq 0)}
#' Here \eqn{f_Y} is a density with respect to Lebesgue measure plus a point mass at
#' zero, so it integrates to \eqn{1 - \pi}: the remainder is the atom, which
#' \code{\link{distrib_atoms}} reports and \code{\link{check_distrib}} accounts for. The
#' classical members of this family --- zero-adjusted gamma, inverse Gaussian or
#' lognormal for semicontinuous data, zero-adjusted beta for proportions --- all have a
#' parent supported on \eqn{(0, \infty)} or \eqn{(0,1)}, where zero sits on the boundary.
#' A parent supported on the whole line is accepted too and gives a "spike at zero"
#' model, but note that \eqn{y = 0} then no longer identifies its own mechanism: the
#' atom sits where the density is positive.
#'
#' \strong{Choosing between the two wrappers.} Zero-inflation can only add zeros, since
#' \eqn{P(Y=0) = \zeta + (1-\zeta)f(0) > f(0)}; the hurdle replaces \eqn{f(0)} outright
#' and so also covers the case of \emph{fewer} zeros than the parent implies. Where both
#' apply they are not nested, and they differ in interpretation more than in fit:
#' zero-inflation keeps \eqn{\theta} as the parameters of the original count process,
#' while the hurdle re-reads them as those of a truncated one. Prefer the hurdle when a
#' zero is observable evidence of a distinct decision, and zero-inflation when the two
#' kinds of zero are genuinely indistinguishable.
#'
#' \strong{What the parent must be.} A discrete parent must have 0 in its support: with
#' \eqn{f(0) = 0} there is no mass to remove. Construction also fails when the result
#' would not be identified:
#' \itemize{
#'   \item the parent already models a probability of zero --- zero-truncating a
#'     zero-inflated or zero-adjusted distribution cancels its zero parameter out of the
#'     likelihood entirely, leaving an identically zero score;
#'   \item the support is too small to carry one more parameter: a distribution on
#'     \eqn{k} points has \eqn{k-1} free probabilities, so at least \code{n_params + 2}
#'     support points are needed. Zero-adjusting a Bernoulli leaves the truncated part
#'     concentrated on \eqn{\{1\}}, and \code{mu} vanishes from the likelihood.
#' }
#' A continuous parent whose support does not reach zero (say \eqn{(2, 5)}) is accepted
#' with a warning: the atom is then disconnected from the rest of the distribution,
#' which is legitimate but rarely intended.
#'
#' @return An S7 object of class \code{ZeroAdjustedDiscreteDistrib} or
#'   \code{ZeroAdjustedContinuousDistrib}.
#'
#' @examples
#' # Hurdle Poisson: the mass at zero is exactly za, not dpois(0, mu)
#' zap <- zero_adjusted(poisson_distrib())
#' distrib_pdf(zap, 0:5, list(mu = 3, za = 0.3))
#'
#' # Semicontinuous data: a spike at zero and a gamma above it
#' zagamma <- zero_adjusted(gamma_distrib())
#' distrib_atoms(zagamma, list(mu = 2, sigma2 = 1, za = 0.3))
#'
#' # The truncated part of a zero-adjusted Bernoulli has no free parameter
#' try(zero_adjusted(bernoulli_distrib()))
#'
#' @seealso \code{\link{zero_inflated}} for the mixture counterpart,
#'   \code{\link{distrib_atoms}}, \code{\link{check_distrib}}.
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
