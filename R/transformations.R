#' @include distrib.R generics.R

#' @title S7 Class for Variable Transformers
#' @name transformer
#'
#' @description
#' A \code{transformer} object defines the mathematical rules for transforming a random
#' variable \eqn{Y = g(X)}. It is used as input to \code{\link{transformation}}.
#'
#' @param name A character string identifying the transformation.
#' @param trans_fun The forward transformation \eqn{y = g(x)}.
#' @param trans_inv The inverse transformation \eqn{x = g^{-1}(y)}.
#' @param trans_abs_jac The absolute Jacobian of the inverse transformation
#'   \eqn{|J(y)| = |dx/dy|}; must accept a \code{log} argument.
#' @param trans_inv_hessian The second derivative of the inverse transformation \eqn{d^2x/dy^2}.
#' @param grad_log_jac The first derivative of \eqn{\log|J(y)|} with respect to \eqn{y}.
#' @param hess_log_jac The second derivative of \eqn{\log|J(y)|} with respect to \eqn{y}.
#' @param bounds_fun Maps the original support bounds to the transformed ones.
#' @param valid_support Checks whether a support is compatible with the transformation.
#' @param decreasing Logical; \code{TRUE} for monotonically decreasing transformations.
#'
#' @section Methods:
#' No method dispatches on this class: a \code{transformer} is a description of a
#' change of variables, not a distribution. It is consumed by
#' \code{\link{transformation}}, which returns a
#' \code{\link{TransformedDistrib}} carrying the full set of distribution methods.
#'
#' Ready-made transformers: \code{\link{log_transform}}, \code{\link{exp_transform}},
#' \code{\link{sqrt_transform}}, \code{\link{inverse_transform}},
#' \code{\link{power_transform}}, \code{\link{bc_transform}},
#' \code{\link{yj_transform}}, \code{\link{softplus_transform}},
#' \code{\link{asinh_transform}}, \code{\link{logit_transform}},
#' \code{\link{expit_transform}}, \code{\link{affine_transform}}.
#'
#' @seealso \code{\link{transformation}}
#' @export
transformer <- S7::new_class("transformer",
  properties = list(
    name              = S7::class_character,
    trans_fun         = S7::class_function,
    trans_inv         = S7::class_function,
    trans_abs_jac     = S7::class_function,
    trans_inv_hessian = S7::class_function,
    grad_log_jac      = S7::class_function,
    hess_log_jac      = S7::class_function,
    bounds_fun        = S7::class_function,
    valid_support     = S7::class_function,
    decreasing        = S7::class_logical
  )
)

# --- TRANSFORMER CONSTRUCTORS ---

#' Logarithmic Transformation
#'
#' @description Transformer for \eqn{Y = \log(X)}: maps \eqn{(0, \infty)} to the real line.
#' Inverse \eqn{X = e^Y}, Jacobian \eqn{|J| = e^Y}.
#' @return A \code{\link{transformer}} object.
#' @export
log_transform <- function() {
  transformer(
    name = "log",
    valid_support = function(bounds) bounds[1] >= 0,
    bounds_fun = function(bounds) {
      res <- log(bounds)
      res[bounds == 0] <- -Inf
      res
    },
    trans_fun = log,
    trans_inv = exp,
    trans_abs_jac = function(y, log = TRUE) if (log) y else exp(y),
    trans_inv_hessian = function(y) exp(y),
    grad_log_jac = function(y) rep(1, length(y)),
    hess_log_jac = function(y) rep(0, length(y)),
    decreasing = FALSE
  )
}

#' Exponential Transformation
#'
#' @description Transformer for \eqn{Y = e^X}: maps the real line to \eqn{(0, \infty)}.
#' Inverse \eqn{X = \log(Y)}, Jacobian \eqn{|J| = 1/Y}.
#' @return A \code{\link{transformer}} object.
#' @export
exp_transform <- function() {
  transformer(
    name = "exp",
    valid_support = function(bounds) TRUE,
    bounds_fun = exp,
    trans_fun = exp,
    trans_inv = log,
    trans_abs_jac = function(y, log = TRUE) if (log) -log(y) else 1 / y,
    trans_inv_hessian = function(y) -1 / y^2,
    grad_log_jac = function(y) -1 / y,
    hess_log_jac = function(y) 1 / y^2,
    decreasing = FALSE
  )
}

#' Reciprocal (Inverse) Transformation
#'
#' @description Transformer for \eqn{Y = 1/X} (monotonically decreasing on a support not
#' containing 0). Inverse \eqn{X = 1/Y}, Jacobian \eqn{|J| = 1/Y^2}.
#' @return A \code{\link{transformer}} object.
#' @export
inverse_transform <- function() {
  transformer(
    name = "inverse",
    valid_support = function(bounds) !(bounds[1] < 0 && bounds[2] > 0),
    bounds_fun = function(bounds) sort(1 / bounds),
    trans_fun = function(x) 1 / x,
    trans_inv = function(y) 1 / y,
    trans_abs_jac = function(y, log = TRUE) if (log) -2 * log(abs(y)) else 1 / y^2,
    trans_inv_hessian = function(y) 2 / y^3,
    grad_log_jac = function(y) -2 / y,
    hess_log_jac = function(y) 2 / y^2,
    decreasing = TRUE
  )
}

#' Square Root Transformation
#'
#' @description Transformer for \eqn{Y = \sqrt{X}} (requires \eqn{X \ge 0}).
#' Inverse \eqn{X = Y^2}, Jacobian \eqn{|J| = 2Y}.
#' @return A \code{\link{transformer}} object.
#' @export
sqrt_transform <- function() {
  transformer(
    name = "sqrt",
    valid_support = function(bounds) bounds[1] >= 0,
    bounds_fun = sqrt,
    trans_fun = sqrt,
    trans_inv = function(y) y^2,
    trans_abs_jac = function(y, log = TRUE) if (log) log(2) + log(y) else 2 * y,
    trans_inv_hessian = function(y) rep(2, length(y)),
    grad_log_jac = function(y) 1 / y,
    hess_log_jac = function(y) -1 / y^2,
    decreasing = FALSE
  )
}

#' Power Transformation
#'
#' @description Transformer for \eqn{Y = X^p}. Inverse \eqn{X = Y^{1/p}},
#' Jacobian \eqn{|J| = \frac{1}{|p|}|Y|^{1/p - 1}}. Fractional powers require a
#' non-negative support; even integer powers require a support that does not straddle 0.
#' @param p Numeric. The exponent. Defaults to 2.
#' @return A \code{\link{transformer}} object.
#' @export
power_transform <- function(p = 2) {
  transformer(
    name = paste0("power_", p),
    valid_support = function(bounds) {
      is_integer <- (round(p) == p)
      lb <- bounds[1]
      ub <- bounds[2]
      if (!is_integer && lb < 0) {
        FALSE
      } else if (p < 0 && lb <= 0 && ub >= 0) {
        FALSE
      } else if (is_integer && (p %% 2 == 0) && lb < 0 && ub > 0) {
        FALSE
      } else {
        TRUE
      }
    },
    bounds_fun = function(bounds) {
      res <- bounds^p
      is_even_integer <- (round(p) == p && p %% 2 == 0)
      if (p < 0 || (is_even_integer && bounds[2] < 0)) sort(res) else res
    },
    trans_fun = function(x) x^p,
    trans_inv = function(y) sign(y) * abs(y)^(1 / p),
    trans_abs_jac = function(y, log = TRUE) {
      log_J <- -log(abs(p)) + (1 / p - 1) * log(abs(y))
      if (log) log_J else exp(log_J)
    },
    trans_inv_hessian = function(y) (1 / p) * (1 / p - 1) * sign(y) * abs(y)^(1 / p - 2),
    grad_log_jac = function(y) (1 / p - 1) / y,
    hess_log_jac = function(y) -(1 / p - 1) / y^2,
    decreasing = (p < 0)
  )
}

#' Inverse Hyperbolic Sine Transformation
#'
#' @description Transformer for \eqn{Y = \text{asinh}(X)}, a log-like transformation that
#' handles zero and negative values. Inverse \eqn{X = \sinh(Y)}, Jacobian \eqn{|J| = \cosh(Y)}.
#' @return A \code{\link{transformer}} object.
#' @export
asinh_transform <- function() {
  transformer(
    name = "asinh",
    valid_support = function(bounds) TRUE,
    bounds_fun = asinh,
    trans_fun = asinh,
    trans_inv = sinh,
    trans_abs_jac = function(y, log = TRUE) {
      # log(cosh(y)) = |y| + log1p(exp(-2|y|)) - log(2), stable for large |y|
      val <- abs(y) + log1p(exp(-2 * abs(y))) - log(2)
      if (log) val else exp(val)
    },
    trans_inv_hessian = function(y) sinh(y),
    grad_log_jac = function(y) tanh(y),
    hess_log_jac = function(y) 1 - tanh(y)^2,
    decreasing = FALSE
  )
}

#' Box-Cox Transformation
#'
#' @description Transformer for the one-parameter Box-Cox transformation
#' \eqn{Y = (X^\lambda - 1)/\lambda} (with \eqn{Y = \log X} for \eqn{\lambda = 0}).
#' Requires \eqn{X > 0}.
#' @param lambda Numeric. The transformation parameter.
#' @return A \code{\link{transformer}} object.
#' @export
bc_transform <- function(lambda) {
  if (abs(lambda) < 1e-10) {
    l <- log_transform()
    l@name <- "box_cox_0"
    return(l)
  }

  transformer(
    name = paste0("box_cox_", lambda),
    valid_support = function(bounds) bounds[1] >= 0,
    bounds_fun = function(bounds) {
      res <- (bounds^lambda - 1) / lambda
      if (lambda < 0) sort(res) else res
    },
    trans_fun = function(x) (x^lambda - 1) / lambda,
    trans_inv = function(y) {
      base <- lambda * y + 1
      base[base < 0] <- 0
      base^(1 / lambda)
    },
    trans_abs_jac = function(y, log = TRUE) {
      term <- pmax(lambda * y + 1, 1e-16)
      log_J <- ((1 - lambda) / lambda) * log(term)
      if (log) log_J else exp(log_J)
    },
    trans_inv_hessian = function(y) {
      term <- pmax(lambda * y + 1, 0)
      (1 - lambda) * term^((1 - 2 * lambda) / lambda)
    },
    grad_log_jac = function(y) (1 - lambda) / (lambda * y + 1),
    hess_log_jac = function(y) -(lambda * (1 - lambda)) / (lambda * y + 1)^2,
    decreasing = (lambda < 0)
  )
}

#' Yeo-Johnson Transformation
#'
#' @description Transformer for the Yeo-Johnson transformation, extending Box-Cox to
#' negative values: Box-Cox of \eqn{X+1} (parameter \eqn{\lambda}) for \eqn{X \ge 0},
#' and negated Box-Cox of \eqn{|X|+1} (parameter \eqn{2-\lambda}) for \eqn{X < 0}.
#' @param lambda Numeric. The transformation parameter.
#' @return A \code{\link{transformer}} object.
#' @export
yj_transform <- function(lambda) {
  lam2 <- 2 - lambda

  tf <- function(x) {
    y <- numeric(length(x))
    pos <- x >= 0
    if (abs(lambda) < 1e-10) {
      y[pos] <- log1p(x[pos])
    } else {
      y[pos] <- ((x[pos] + 1)^lambda - 1) / lambda
    }
    if (abs(lam2) < 1e-10) {
      y[!pos] <- -log1p(-x[!pos])
    } else {
      y[!pos] <- -((-x[!pos] + 1)^lam2 - 1) / lam2
    }
    y
  }

  transformer(
    name = paste0("yeo_johnson_", lambda),
    valid_support = function(bounds) TRUE,
    bounds_fun = tf,
    trans_fun = tf,
    trans_inv = function(y) {
      x <- numeric(length(y))
      pos <- y >= 0
      if (abs(lambda) < 1e-10) {
        x[pos] <- expm1(y[pos])
      } else {
        x[pos] <- (lambda * y[pos] + 1)^(1 / lambda) - 1
      }
      if (abs(lam2) < 1e-10) {
        x[!pos] <- -expm1(-y[!pos])
      } else {
        x[!pos] <- 1 - (1 - lam2 * y[!pos])^(1 / lam2)
      }
      x
    },
    trans_abs_jac = function(y, log = TRUE) {
      log_J <- numeric(length(y))
      pos <- y >= 0
      if (any(pos)) log_J[pos] <- ((1 - lambda) / lambda) * log(lambda * y[pos] + 1)
      if (any(!pos)) log_J[!pos] <- ((lambda - 1) / lam2) * log(1 - lam2 * y[!pos])
      if (log) log_J else exp(log_J)
    },
    trans_inv_hessian = function(y) {
      h <- numeric(length(y))
      pos <- y >= 0
      if (any(pos)) h[pos] <- (1 - lambda) * (lambda * y[pos] + 1)^((1 - 2 * lambda) / lambda)
      if (any(!pos)) h[!pos] <- (lam2 - 1) * (1 - lam2 * y[!pos])^((1 - 2 * lam2) / lam2)
      h
    },
    grad_log_jac = function(y) {
      g <- numeric(length(y))
      pos <- y >= 0
      if (any(pos)) g[pos] <- (1 - lambda) / (lambda * y[pos] + 1)
      if (any(!pos)) g[!pos] <- (lambda - 1) / (1 - lam2 * y[!pos])
      g
    },
    hess_log_jac = function(y) {
      h <- numeric(length(y))
      pos <- y >= 0
      if (any(pos)) h[pos] <- -(lambda * (1 - lambda)) / (lambda * y[pos] + 1)^2
      if (any(!pos)) h[!pos] <- ((lambda - 1) * lam2) / (1 - lam2 * y[!pos])^2
      h
    },
    decreasing = FALSE
  )
}

#' Affine (Location-Scale) Transformation
#'
#' @description Transformer for \eqn{Y = \text{loc} + \text{scale} \cdot X}.
#' Inverse \eqn{X = (Y - \text{loc})/\text{scale}}, Jacobian \eqn{|J| = 1/|\text{scale}|}.
#' @param loc Numeric. The location shift. Defaults to 0.
#' @param scale Numeric. The scale multiplier (non-zero). Defaults to 1.
#' @return A \code{\link{transformer}} object.
#' @export
affine_transform <- function(loc = 0, scale = 1) {
  if (scale == 0) {
    stop("Scale parameter cannot be equal to zero.", call. = FALSE)
  }

  transformer(
    name = paste0("affine_loc", loc, "_scale", scale),
    valid_support = function(bounds) TRUE,
    bounds_fun = function(bounds) {
      res <- loc + scale * bounds
      if (scale < 0) sort(res) else res
    },
    trans_fun = function(x) loc + scale * x,
    trans_inv = function(y) (y - loc) / scale,
    trans_abs_jac = function(y, log = TRUE) {
      val <- rep(1 / abs(scale), length(y))
      if (log) log(val) else val
    },
    trans_inv_hessian = function(y) rep(0, length(y)),
    grad_log_jac = function(y) rep(0, length(y)),
    hess_log_jac = function(y) rep(0, length(y)),
    decreasing = (scale < 0)
  )
}

#' Logit Transformation
#'
#' @description Transformer for \eqn{Y = \text{logit}(X)}: maps \eqn{(0, 1)} to the real line.
#' Inverse \eqn{X = \text{plogis}(Y)}, Jacobian \eqn{|J| = \text{dlogis}(Y)}.
#' @return A \code{\link{transformer}} object.
#' @export
logit_transform <- function() {
  transformer(
    name = "logit",
    valid_support = function(bounds) bounds[1] >= 0 && bounds[2] <= 1,
    bounds_fun = function(bounds) stats::qlogis(bounds),
    trans_fun = stats::qlogis,
    trans_inv = stats::plogis,
    trans_abs_jac = function(y, log = TRUE) stats::dlogis(y, log = log),
    trans_inv_hessian = function(y) stats::dlogis(y) * (1 - 2 * stats::plogis(y)),
    grad_log_jac = function(y) 1 - 2 * stats::plogis(y),
    hess_log_jac = function(y) -2 * stats::dlogis(y),
    decreasing = FALSE
  )
}

#' Expit (Sigmoid) Transformation
#'
#' @description Transformer for \eqn{Y = \text{plogis}(X)}: maps the real line to \eqn{(0, 1)}.
#' Inverse \eqn{X = \text{logit}(Y)}, Jacobian \eqn{|J| = 1/(Y(1-Y))}.
#' @return A \code{\link{transformer}} object.
#' @export
expit_transform <- function() {
  transformer(
    name = "expit",
    valid_support = function(bounds) TRUE,
    bounds_fun = function(bounds) stats::plogis(bounds),
    trans_fun = stats::plogis,
    trans_inv = stats::qlogis,
    trans_abs_jac = function(y, log = TRUE) {
      val <- -(log(y) + log1p(-y))
      if (log) val else exp(val)
    },
    trans_inv_hessian = function(y) (2 * y - 1) / (y^2 * (1 - y)^2),
    grad_log_jac = function(y) -1 / y + 1 / (1 - y),
    hess_log_jac = function(y) 1 / y^2 + 1 / (1 - y)^2,
    decreasing = FALSE
  )
}

#' Softplus Transformation
#'
#' @description Transformer for the inverse-softplus map \eqn{Y = \frac{1}{a}\log(e^{aX} - 1)},
#' sending \eqn{(0, \infty)} to the real line. Inverse (softplus)
#' \eqn{X = \frac{1}{a}\log(1 + e^{aY})}, Jacobian \eqn{|J| = \text{plogis}(aY)}.
#' @param a Numeric. Positive scale parameter. Defaults to 1.
#' @return A \code{\link{transformer}} object.
#' @export
softplus_transform <- function(a = 1) {
  if (a <= 0) {
    stop("Scale parameter 'a' must be greater than 0.", call. = FALSE)
  }

  transformer(
    name = paste0("softplus_a", round(a, 5)),
    valid_support = function(bounds) bounds[1] >= 0,
    bounds_fun = function(bounds) {
      res <- log(expm1(a * bounds)) / a
      res[bounds == 0] <- -Inf
      res
    },
    trans_fun = function(x) log(expm1(a * x)) / a,
    trans_inv = function(y) pmax(0, y) + log1p(exp(-abs(a * y))) / a,
    trans_abs_jac = function(y, log = TRUE) stats::plogis(a * y, log.p = log),
    trans_inv_hessian = function(y) a * stats::dlogis(a * y),
    grad_log_jac = function(y) a * (1 - stats::plogis(a * y)),
    hess_log_jac = function(y) -(a^2) * stats::dlogis(a * y),
    decreasing = FALSE
  )
}

# --- TRANSFORMED DISTRIBUTION CLASS ---

#' @title S7 Class for Transformed Distributions
#' @name TransformedDistrib
#'
#' @description
#' A subclass of \code{continuous_distrib} representing the distribution of
#' \eqn{Y = g(X)}, where \eqn{X} follows a wrapped continuous distribution and \eqn{g}
#' is a bijective \code{\link{transformer}}.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{continuous_distrib} object.
#' @param transformer The \code{\link{transformer}} defining \eqn{g}.
#' @seealso \code{\link{transformation}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.TransformedDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.TransformedDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.TransformedDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.TransformedDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.TransformedDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.TransformedDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.TransformedDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
TransformedDistrib <- S7::new_class("TransformedDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib,
    transformer = transformer
  )
)

#' @title Transformed Probability Density Function
#' @name distrib_pdf.TransformedDistrib
#' @description
#' Change of variables: \eqn{f_Y(y) = f_X(g^{-1}(y)) \cdot |J(y)|}. Computed in log space;
#' singular log-densities are clamped to avoid \code{Inf - Inf} during integration.
#' @param distrib A \code{TransformedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @seealso \code{\link{transformation}}
S7::method(distrib_pdf, TransformedDistrib) <- function(distrib, y, theta, log = FALSE) {
  tr <- distrib@transformer
  log_pdf_x <- distrib_pdf(distrib@parent_distrib, tr@trans_inv(y), theta, log = TRUE)
  log_J <- tr@trans_abs_jac(y, log = TRUE)

  sing <- is.infinite(log_pdf_x) & log_pdf_x > 0
  if (any(sing)) log_pdf_x[sing] <- log(.Machine$double.xmax)

  log_val <- log_pdf_x + log_J
  # A zero parent density wins over an infinite Jacobian (-Inf + Inf would be NaN)
  log_val[is.infinite(log_pdf_x) & log_pdf_x < 0] <- -Inf

  if (log) log_val else exp(log_val)
}

#' @title Transformed Cumulative Distribution Function
#' @name distrib_cdf.TransformedDistrib
#' @description
#' \eqn{F_Y(q) = F_X(g^{-1}(q))}, with tails swapped for decreasing transformations.
#' @param distrib A \code{TransformedDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @seealso \code{\link{transformation}}
S7::method(distrib_cdf, TransformedDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  tr <- distrib@transformer
  if (tr@decreasing) lower.tail <- !lower.tail
  distrib_cdf(distrib@parent_distrib, tr@trans_inv(q), theta, lower.tail = lower.tail, log.p = log.p)
}

#' @title Transformed Quantile Function
#' @name distrib_quantile.TransformedDistrib
#' @description
#' \eqn{Q_Y(p) = g(Q_X(p))}, with tails swapped for decreasing transformations.
#' @param distrib A \code{TransformedDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @seealso \code{\link{transformation}}
S7::method(distrib_quantile, TransformedDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  tr <- distrib@transformer
  if (tr@decreasing) lower.tail <- !lower.tail
  tr@trans_fun(distrib_quantile(distrib@parent_distrib, p, theta, lower.tail = lower.tail, log.p = log.p))
}

#' @title Transformed Random Number Generator
#' @name distrib_rng.TransformedDistrib
#' @description Draws from the parent distribution and applies \eqn{g}.
#' @param distrib A \code{TransformedDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list of the parent's parameters.
#' @seealso \code{\link{transformation}}
S7::method(distrib_rng, TransformedDistrib) <- function(distrib, n, theta) {
  distrib@transformer@trans_fun(distrib_rng(distrib@parent_distrib, n, theta))
}

#' @title Transformed Analytical Gradient
#' @name distrib_gradient.TransformedDistrib
#' @description
#' The Jacobian does not depend on the parameters, so the score of the transformed model
#' equals the parent's score evaluated at \eqn{x = g^{-1}(y)}.
#' @param distrib A \code{TransformedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{transformation}}
S7::method(distrib_gradient, TransformedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  distrib_gradient(distrib@parent_distrib, distrib@transformer@trans_inv(y), theta)
}

#' @title Transformed Analytical Observed Hessian
#' @name distrib_hessian.TransformedDistrib
#' @description
#' The observed Hessian of the transformed model equals the parent's observed Hessian
#' evaluated at \eqn{x = g^{-1}(y)}.
#' @param distrib A \code{TransformedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{transformation}}
S7::method(distrib_hessian, TransformedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  distrib_hessian(distrib@parent_distrib, distrib@transformer@trans_inv(y), theta)
}

#' @title Transformed Analytical Expected Hessian
#' @name distrib_expected_hessian.TransformedDistrib
#' @description
#' Since \eqn{\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) + \log|J(y)|} and the Jacobian
#' does not depend on \eqn{\theta}, the expected Hessian of the transformed model is
#' \emph{exactly} the parent's expected Hessian (the expectation is just re-parameterized
#' by the change of variables). No Monte Carlo approximation is needed.
#' @param distrib A \code{TransformedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list of the parent's parameters.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{transformation}}
S7::method(distrib_expected_hessian, TransformedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  distrib_expected_hessian(distrib@parent_distrib, distrib@transformer@trans_inv(y), theta)
}

# --- CONSTRUCTOR WRAPPER ---

#' Apply a Variable Transformation to a Distribution Object
#'
#' @description
#' Creates a new distribution object for \eqn{Y = g(X)}, where \eqn{X} follows an existing
#' \strong{continuous} distribution and \eqn{g} is a bijective transformation described by
#' a \code{\link{transformer}} object.
#'
#' @param distrib An object inheriting from \code{continuous_distrib}.
#' @param transformer A \code{\link{transformer}} object (e.g. \code{\link{log_transform}()},
#'   \code{\link{bc_transform}(0.5)}, \code{\link{affine_transform}(1, 2)}).
#'
#' @details
#' The density follows the change-of-variables formula
#' \deqn{f_Y(y) = f_X(g^{-1}(y)) \cdot \left|\dfrac{d}{dy} g^{-1}(y)\right|}
#' CDF, quantiles and RNG are obtained by mapping through \eqn{g} (with tails swapped for
#' decreasing transformations). Since \eqn{g} does not depend on the parameters, the
#' score, observed Hessian and expected Hessian coincide with the parent's, evaluated at
#' \eqn{x = g^{-1}(y)}. Moments are available numerically via \code{\link{moment}}.
#'
#' @return An S7 object of class \code{TransformedDistrib} (inheriting from \code{continuous_distrib}).
#'
#' @examples
#' \dontrun{
#' # A lognormal built by transformation, equal to lognormal_distrib()
#' logn <- transformation(gaussian_distrib(), exp_transform())
#' distrib_pdf(logn, 2, list(mu = 0, sigma = 1))
#' dlnorm(2, 0, 1)
#' }
#'
#' @seealso \code{\link{transformer}}, \code{\link{log_transform}}, \code{\link{exp_transform}},
#'   \code{\link{affine_transform}}, \code{\link{bc_transform}}, \code{\link{yj_transform}}
#' @export
transformation <- function(distrib, transformer) {
  if (!S7::S7_inherits(distrib, continuous_distrib)) {
    stop("transformation() currently supports only continuous distributions.", call. = FALSE)
  }
  if (!S7::S7_inherits(transformer, distributions7::transformer)) {
    stop("Argument 'transformer' must be a 'transformer' object.", call. = FALSE)
  }
  if (!transformer@valid_support(distrib@bounds)) {
    stop(sprintf(
      "The '%s' transformation is not valid for the support of the '%s' distribution.",
      transformer@name, distrib@distrib_name
    ), call. = FALSE)
  }

  interpretation <- paste0(distrib@params_interpretation, " (", distrib@distrib_name, " scale)")
  names(interpretation) <- distrib@params

  TransformedDistrib(
    parent_distrib = distrib,
    transformer = transformer,
    distrib_name = paste0(transformer@name, "(", distrib@distrib_name, ")"),
    dimension = distrib@dimension,
    bounds = transformer@bounds_fun(distrib@bounds),

    params = distrib@params,
    params_interpretation = interpretation,
    n_params = distrib@n_params,

    params_bounds = distrib@params_bounds,
    link_params = distrib@link_params,
    # g does not depend on theta, so a kink in the log-likelihood survives the
    # change of variables untouched: affine(laplace) is still non-smooth in mu.
    params_smooth = param_smoothness(distrib)
  )
}
