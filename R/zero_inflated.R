#' @include distrib.R generics.R utility_functions.R

#' @title S7 Class for Zero-Inflated Distributions
#' @name ZeroInflatedDistrib
#'
#' @description
#' A subclass of \code{discrete_distrib} representing the zero-inflated version of a
#' wrapped discrete distribution: a mixture of a point mass at zero (with probability
#' \eqn{\zeta}) and the original count distribution.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{discrete_distrib} object.
#' @seealso \code{\link{zero_inflated}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.ZeroInflatedDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.ZeroInflatedDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.ZeroInflatedDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.ZeroInflatedDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.ZeroInflatedDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.ZeroInflatedDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.ZeroInflatedDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{discrete_distrib}}.
ZeroInflatedDistrib <- S7::new_class("ZeroInflatedDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib
  )
)

# Internal: split the full theta into the parent's parameters and the mixture one.
split_mix_theta <- function(distrib, theta) {
  n <- distrib@n_params
  list(orig = theta[seq_len(n - 1L)], mix = theta[[n]])
}

#' @title Zero-Inflated Probability Mass Function
#' @name distrib_pdf.ZeroInflatedDistrib
#' @description
#' \deqn{P(Y=y) = \zeta\,\mathbb{I}(y=0) + (1-\zeta) f(y; \theta)}
#' where \eqn{f} is the PMF of the wrapped distribution.
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @param log Logical; if \code{TRUE}, returns the log-probability.
#' @seealso \code{\link{zero_inflated}}
S7::method(distrib_pdf, ZeroInflatedDistrib) <- function(distrib, y, theta, log = FALSE) {
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
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @seealso \code{\link{zero_inflated}}
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
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @seealso \code{\link{zero_inflated}}
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
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @seealso \code{\link{zero_inflated}}
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
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{zero_inflated}}
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
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{zero_inflated}}
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

  res <- list()

  # Block theta-theta: at y=0, w*H(0) + w(1-w)*S(0)S(0)'
  for (nm in names(h_orig)) {
    parts <- strsplit(nm, "_")[[1]]
    p1 <- parts[1]
    p2 <- parts[length(parts)]
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
#' @param distrib A \code{ZeroInflatedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with the parent's parameters followed by \code{zi}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{zero_inflated}}
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

  res <- list()

  # Block theta-theta:
  #   L0 * [w0 H(0) + w0(1-w0) S(0)S(0)'] + (1-zi) * (E[H] - f(0)H(0))
  for (nm in names(h_orig_exp)) {
    parts <- strsplit(nm, "_")[[1]]
    p1 <- parts[1]
    p2 <- parts[length(parts)]
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
#' Creates a zero-inflated version of an existing \strong{discrete} distribution object:
#' a mixture of a point mass at zero (with probability \eqn{\zeta}, parameter \code{zi})
#' and the original count distribution.
#'
#' @param distrib An object inheriting from \code{discrete_distrib} whose support includes 0.
#' @param link_zi A link function object for the zero-inflation probability \eqn{\zeta}.
#'   Defaults to \code{\link[linkfunctions7]{logit_link}}.
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
#' The resulting object supports the full \code{distrib} API: pdf, cdf, quantile, rng,
#' analytical gradient, observed and expected Hessian (all derived from the parent's),
#' plus numerical moments via \code{\link{moment}}.
#'
#' @return An S7 object of class \code{ZeroInflatedDistrib} (inheriting from \code{discrete_distrib}).
#'
#' @examples
#' \dontrun{
#' zip <- zero_inflated(poisson_distrib())
#' distrib_pdf(zip, 0:5, list(mu = 3, zi = 0.2))
#' }
#'
#' @importFrom linkfunctions7 logit_link
#' @export
zero_inflated <- function(distrib, link_zi = logit_link()) {
  if (!S7::S7_inherits(distrib, discrete_distrib)) {
    stop("zero_inflated() requires a discrete distribution ('discrete_distrib').", call. = FALSE)
  }
  if (distrib@bounds[1] > 0) {
    stop("zero_inflated() requires a distribution with 0 in its support.", call. = FALSE)
  }
  if ("zi" %in% distrib@params) {
    stop("The parent distribution already has a parameter named 'zi'.", call. = FALSE)
  }

  ZeroInflatedDistrib(
    parent_distrib = distrib,
    distrib_name = paste0("zero-inflated ", distrib@distrib_name),
    dimension = distrib@dimension,
    bounds = c(min(0, distrib@bounds[1]), distrib@bounds[2]),

    params = c(distrib@params, "zi"),
    params_interpretation = c(distrib@params_interpretation, zi = "prob. of structural zero"),
    n_params = distrib@n_params + 1,

    params_bounds = c(distrib@params_bounds, list(zi = c(0, 1))),
    link_params = c(distrib@link_params, list(zi = link_zi))
  )
}
