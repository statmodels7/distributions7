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

  res <- list()

  # Block za-za
  res[[paste0(za_name, "_", za_name)]] <- ifelse(y == 0, -1 / (za^2), -1 / ((1 - za)^2))

  # Mixed blocks are identically 0 (likelihood separation)
  for (nm in names(pars$orig)) {
    res[[paste0(nm, "_", za_name)]] <- rep(0, n)
  }

  # Block theta-theta with truncation correction
  for (nm in names(h_orig)) {
    parts <- strsplit(nm, "_")[[1]]
    s1 <- grad_0[[parts[1]]]
    s2 <- grad_0[[parts[length(parts)]]]
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

  res <- list()
  res[[paste0(za_name, "_", za_name)]] <- rep(-1 / (za * (1 - za)), length.out = n)

  for (nm in names(pars$orig)) {
    res[[paste0(nm, "_", za_name)]] <- rep(0, n)
  }

  for (nm in names(h_orig_exp)) {
    parts <- strsplit(nm, "_")[[1]]
    s1 <- grad_0[[parts[1]]]
    s2 <- grad_0[[parts[length(parts)]]]
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
#' Creates a zero-adjusted version of an existing distribution object, automatically
#' dispatching on its type:
#' \itemize{
#'   \item \strong{Discrete} (support including 0): a hurdle model, where the probability
#'     at zero is replaced by \eqn{\pi} and positive values follow the zero-truncated parent.
#'   \item \strong{Continuous}: a mixed distribution with a point mass \eqn{\pi} at zero
#'     and the original continuous density scaled by \eqn{1-\pi}.
#' }
#' The mixture probability is exposed as parameter \code{za}.
#'
#' @param distrib An object inheriting from \code{discrete_distrib} or \code{continuous_distrib}.
#' @param link_za A link function object for the zero probability \eqn{\pi}.
#'   Defaults to \code{\link[linkfunctions7]{logit_link}}.
#'
#' @return An S7 object of class \code{ZeroAdjustedDiscreteDistrib} or
#'   \code{ZeroAdjustedContinuousDistrib}.
#'
#' @examples
#' \dontrun{
#' zap <- zero_adjusted(poisson_distrib())
#' distrib_pdf(zap, 0:5, list(mu = 3, za = 0.3))
#'
#' zagamma <- zero_adjusted(gamma_distrib())
#' distrib_rng(zagamma, 10, list(mu = 2, sigma2 = 1, za = 0.3))
#' }
#'
#' @importFrom linkfunctions7 logit_link
#' @export
zero_adjusted <- function(distrib, link_za = logit_link()) {
  if ("za" %in% distrib@params) {
    stop("The parent distribution already has a parameter named 'za'.", call. = FALSE)
  }

  common <- list(
    parent_distrib = distrib,
    distrib_name = paste0("zero-adjusted ", distrib@distrib_name),
    dimension = distrib@dimension,
    params = c(distrib@params, "za"),
    params_interpretation = c(distrib@params_interpretation, za = "prob. of zero"),
    n_params = distrib@n_params + 1,
    params_bounds = c(distrib@params_bounds, list(za = c(0, 1))),
    link_params = c(distrib@link_params, list(za = link_za))
  )

  if (S7::S7_inherits(distrib, discrete_distrib)) {
    if (distrib@bounds[1] > 0) {
      stop("zero_adjusted() requires a discrete distribution with 0 in its support.", call. = FALSE)
    }
    do.call(ZeroAdjustedDiscreteDistrib, c(common, list(bounds = c(0, distrib@bounds[2]))))
  } else if (S7::S7_inherits(distrib, continuous_distrib)) {
    do.call(ZeroAdjustedContinuousDistrib, c(common, list(bounds = c(min(0, distrib@bounds[1]), distrib@bounds[2]))))
  } else {
    stop("Input must inherit from 'discrete_distrib' or 'continuous_distrib'.", call. = FALSE)
  }
}
