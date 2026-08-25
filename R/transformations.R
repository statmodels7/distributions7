#' @include distrib.R generics.R
NULL

#' @title S7 Class for Variable Transformers
#' @name transformer
#'
#' @description
#' A `transformer` records the mathematical rules of a change of variables
#' \eqn{Y = g(X)}: the map, its inverse, the Jacobian of that inverse and the
#' two derivatives of the Jacobian's logarithm, together with the bookkeeping
#' that says which supports the map applies to and what it does to them. It
#' carries no distribution and answers no distribution generic; it is consumed
#' by [transformation()], which returns a [TransformedDistrib].
#'
#' @details
#' # What a map must be to be one
#'
#' \eqn{g} has to be BIJECTIVE on the parent's support. A map that is two to
#' one, such as the absolute value, has no inverse to carry a density through
#' and cannot be a transformer at all; [folded()] handles that case by adding
#' the two preimages instead.
#'
#' # Writing your own
#'
#' The ten properties are all required, and `trans_abs_jac` must accept a `log`
#' argument, because [distrib_pdf.TransformedDistrib()] works on the log scale
#' throughout. `valid_support` takes the parent's `bounds` and answers whether
#' the map applies; `bounds_fun` takes the same and returns the transformed
#' support. Nothing validates the derivatives, so
#' [check_distrib()] on the resulting distribution is what catches a
#' transcription error in them.
#'
#' @param name A single string identifying the transformation. It composes the
#'   result's `distrib_name` as `"name(parent)"` unless `transformation()` is
#'   given a `new_name`.
#' @param trans_fun The forward map \eqn{y = g(x)}, a function of one numeric
#'   vector returning one of the same length.
#' @param trans_inv The inverse map \eqn{x = g^{-1}(y)}, likewise.
#' @param trans_abs_jac The absolute Jacobian of the inverse,
#'   \eqn{\lvert J(y)\rvert = \lvert dx/dy\rvert}. It MUST accept a `log`
#'   argument and return the logarithm when it is `TRUE`.
#' @param trans_inv_hessian The second derivative of the inverse map,
#'   \eqn{d^2x/dy^2}.
#' @param grad_log_jac The first derivative of \eqn{\log\lvert J(y)\rvert} with
#'   respect to \eqn{y}.
#' @param hess_log_jac The second derivative of \eqn{\log\lvert J(y)\rvert}
#'   with respect to \eqn{y}.
#' @param bounds_fun A function of the parent's `bounds` returning the
#'   transformed support, a numeric vector of length 2.
#' @param valid_support A function of the parent's `bounds` returning a single
#'   logical: whether the map applies there. It is what
#'   `transformation(gaussian1_distrib(), log_transform())` fails on.
#' @param decreasing Logical of length 1. `TRUE` for a monotonically decreasing
#'   map, which swaps the tails in the distribution and quantile functions.
#'
#' @return An S7 object of class `transformer`, carrying the ten properties
#'   above and nothing else.
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{J} the Jacobian of its inverse, \eqn{X}
#' the parent's variable and \eqn{Y = g(X)} the transformed one.
#'
#' @section Methods:
#' None. A `transformer` is a description of a change of variables, not a
#' distribution. It is consumed by [transformation()], which returns a
#' [TransformedDistrib] carrying the full set of distribution methods.
#'
#' The twelve ready-made transformers: [log_transform()], [exp_transform()],
#' [sqrt_transform()], [inverse_transform()], [power_transform()],
#' [bc_transform()], [yj_transform()], [softplus_transform()],
#' [asinh_transform()], [logit_transform()], [expit_transform()] and
#' [affine_transform()].
#'
#' @seealso [transformation()], which consumes one, [TransformedDistrib] for
#'   what it produces, and [folded()] for a map that cannot be a transformer.
#'
#' @examples
#' tr <- log_transform()
#' S7::S7_inherits(tr, transformer)
#' tr@name
#' S7::prop_names(tr)
#'
#' # The map, its inverse and the Jacobian of that inverse.
#' c(forward = tr@trans_fun(exp(1)), inverse = tr@trans_inv(1),
#'   log_jacobian = tr@trans_abs_jac(1, log = TRUE))
#'
#' # It says which supports it applies to, and what it does to them.
#' c(on_positives = tr@valid_support(c(0, Inf)),
#'   on_the_line = tr@valid_support(c(-Inf, Inf)))
#' tr@bounds_fun(c(0, Inf))
#'
#' # And a decreasing map declares itself, which swaps the tails downstream.
#' vapply(list(log_transform(), inverse_transform(), expit_transform()),
#'        function(z) z@decreasing, TRUE)
#'
#' # A transformer acts only through transformation().
#' distrib_pdf(transformation(gamma2_distrib(), tr), 0,
#'             list(mu = 2, sigma2 = 1))
#'
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gamma2_distrib(), log_transform())
#' distrib_pdf(d, 0, list(mu = 2, sigma2 = 1))
#'
#' @seealso [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#'
#' @seealso [log_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gamma2_distrib(), inverse_transform())
#' distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gamma2_distrib(), sqrt_transform())
#' distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gamma2_distrib(), power_transform(p = 2))
#' distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gaussian1_distrib(), asinh_transform())
#' distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gamma2_distrib(), bc_transform(lambda = 0.5))
#' distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gaussian1_distrib(), yj_transform(lambda = 0.5))
#' distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gaussian1_distrib(), affine_transform(loc = 1, scale = 2))
#' distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(beta1_distrib(), logit_transform())
#' distrib_pdf(d, 0, list(mu = 0.4, phi = 5))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gaussian1_distrib(), expit_transform())
#' distrib_pdf(d, 0.5, list(mu = 0, sigma = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [softplus_transform()], [asinh_transform()], [logit_transform()], [affine_transform()], [transformation()]
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
#' @return A [transformer()] object.
#' @examples
#' d <- transformation(gamma2_distrib(), softplus_transform(a = 1))
#' distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#'
#' @seealso [log_transform()], [exp_transform()], [sqrt_transform()], [inverse_transform()], [power_transform()], [bc_transform()], [yj_transform()], [asinh_transform()], [logit_transform()], [expit_transform()], [affine_transform()], [transformation()]
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
#' The S7 class of the distribution of \eqn{Y = g(X)}, where \eqn{X} follows a
#' wrapped continuous distribution and \eqn{g} is a BIJECTIVE [transformer()].
#' Its density is the change-of-variables formula
#' \deqn{f_Y(y) = f_X(g^{-1}(y))\,\lvert J(y)\rvert,
#'   \qquad J(y) = \frac{d}{dy}g^{-1}(y),}
#' and it carries exactly the parent's parameters: the transformation adds none
#' and removes none.
#'
#' @details
#' # Why the derivatives are the parent's
#'
#' \eqn{g} does not depend on \eqn{\theta}, so the Jacobian is a constant in
#' the parameters and
#' \eqn{\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert}
#' differentiates in \eqn{\theta} to the parent's own derivative at
#' \eqn{x = g^{-1}(y)}. The score, the observed Hessian and the EXPECTED
#' Hessian are therefore the parent's, exactly, with no Monte Carlo anywhere.
#' The response derivatives are not: those pick the Jacobian up and come from
#' the base class.
#'
#' Build one with [transformation()], which checks that the parent is
#' continuous and that the transformer is valid on its support. This page
#' documents the raw S7 constructor, which checks neither.
#'
#' @param parent_distrib The wrapped `continuous_distrib` object.
#' @param transformer The [transformer()] defining \eqn{g}.
#' @inheritParams distrib
#'
#' @return An S7 object of class `TransformedDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. It carries `parent_distrib` and
#'   `transformer` beside the parent's properties. For an object built by
#'   [transformation()], `params`, `params_bounds` and `link_params` are the
#'   parent's unchanged; `bounds` is the transformer's image of the parent's;
#'   `params_interpretation` is the parent's with the parent's name appended in
#'   brackets, since a mean on the transformed scale is not the parent's mean;
#'   and `distrib_name` is `"g(parent)"` or whatever `new_name` said.
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{J} the Jacobian of its inverse,
#' \eqn{f_X} the parent's density, \eqn{f_Y} the transformed one and
#' \eqn{\ell} a log-density.
#'
#' @seealso [transformation()] to build one, [transformer()] for the map,
#'   [folded()] for a map that is NOT bijective and so cannot be a transformer,
#'   and [fixed()] for the wrapper that removes parameters.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.TransformedDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.TransformedDistrib],
#'   [`distrib_gradient()`][distrib_gradient.TransformedDistrib],
#'   [`distrib_hessian()`][distrib_hessian.TransformedDistrib],
#'   [`distrib_pdf()`][distrib_pdf.TransformedDistrib],
#'   [`distrib_quantile()`][distrib_quantile.TransformedDistrib],
#'   [`distrib_rng()`][distrib_rng.TransformedDistrib]
#'
#' The third and fourth derivatives come from the shared wrapper machinery in
#' `wrapper_derivatives.R`; everything else is inherited from
#' [continuous_distrib()].
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The parameters are the parent's, and the support is the image of its own.
#' d@params
#' d@bounds
#' d@distrib_name
#'
#' # The exponential of a gaussian is the lognormal, whose second parameter is
#' # the variance rather than the standard deviation.
#' y <- c(0.5, 1, 3)
#' all.equal(distrib_pdf(d, y, theta),
#'           distrib_pdf(lognormal1_distrib(), y,
#'                       list(mu = 0.5, sigma2 = 0.8^2)))
#'
#' # The interpretation says which scale a parameter lives on.
#' d@params_interpretation
TransformedDistrib <- S7::new_class("TransformedDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib,
    transformer = transformer
  )
)

#' @title Transformed Probability Density Function
#' @name distrib_pdf.TransformedDistrib
#'
#' @description
#' Computes the change-of-variables density
#' \deqn{f_Y(y) = f_X(g^{-1}(y))\,\lvert J(y)\rvert,}
#' entirely on the LOG scale and exponentiated at the end, so that a density
#' and a Jacobian which would each overflow or underflow on their own combine
#' as a sum of representable numbers.
#'
#' @details
#' Two singular cases are handled explicitly. A parent log-density that is
#' `+Inf`, which a density with a pole produces, is clamped to
#' `log(.Machine$double.xmax)`, so that a quadrature over the transformed
#' density meets a large number instead of `Inf - Inf`. A parent density of
#' exactly zero WINS over an infinite Jacobian, so the result is `-Inf`. The
#' point is outside the transformed support and `-Inf` says so, where the
#' unguarded sum would give `NaN`.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param y A numeric vector of observations on the TRANSFORMED scale. A point
#'   outside the transformed support gives `0`, the parent's density at its
#'   preimage being zero there.
#' @param theta A named list of the PARENT's parameters, unchanged by the
#'   transformation.
#' @param log Logical of length 1. When `TRUE` the log-density is returned,
#'   which is the quantity actually computed. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{J} the Jacobian of its inverse and
#' \eqn{f_X}, \eqn{f_Y} the parent's and the transformed density.
#'
#' @seealso [distrib_cdf.TransformedDistrib()] for the distribution function,
#'   [transformation()] for the family, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' y <- c(0.5, 1, 3)
#'
#' distrib_pdf(d, y, theta)
#'
#' # The exponential of a gaussian is the lognormal.
#' all.equal(distrib_pdf(d, y, theta), dlnorm(y, 0.5, 0.8))
#'
#' # It is a density: it integrates to one over the transformed support.
#' integrate(function(z) distrib_pdf(d, z, theta), 0, Inf)$value
#'
#' # Written out, the parent's density at the preimage times the Jacobian.
#' tr <- d@transformer
#' all.equal(distrib_pdf(d, y, theta),
#'           dnorm(tr@trans_inv(y), 0.5, 0.8) * tr@trans_abs_jac(y))
S7::method(distrib_pdf, TransformedDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
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
#'
#' @description
#' Computes \eqn{F_Y(q) = F_X(g^{-1}(q))} by evaluating the parent's own
#' distribution function at the preimage. For a DECREASING transformation the
#' tails swap: `lower.tail` is inverted before the parent is called, because
#' \eqn{Y \le q} is \eqn{X \ge g^{-1}(q)} there.
#'
#' @details
#' `lower.tail` and `log.p` are passed THROUGH to the parent and not applied
#' afterwards, so a parent with an accurate upper tail or an accurate
#' log-probability keeps that accuracy here.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param q A numeric vector of quantiles on the transformed scale.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}. For a decreasing transformation this
#'   is inverted before reaching the parent.
#' @param log.p Logical of length 1. When `TRUE` the logarithm is returned,
#'   computed by the parent. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities, in \eqn{[0, 1]}.
#'
#' @section Notation:
#' \eqn{g} is the transformation and \eqn{F_X}, \eqn{F_Y} the parent's and the
#' transformed distribution function.
#'
#' @seealso [distrib_pdf.TransformedDistrib()] for the density,
#'   [distrib_quantile.TransformedDistrib()], which inverts this, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' q <- c(0.5, 1, 3)
#'
#' distrib_cdf(d, q, theta)
#' all.equal(distrib_cdf(d, q, theta), plnorm(q, 0.5, 0.8))
#'
#' # A decreasing transformation swaps the tails: the reciprocal of a gamma.
#' ig <- transformation(gamma1_distrib(), inverse_transform())
#' th2 <- list(mu = 2, phi = 0.3)
#' ig@transformer@decreasing
#' c(transformed = distrib_cdf(ig, 0.5, th2),
#'   parent_upper = distrib_cdf(gamma1_distrib(), 2, th2, lower.tail = FALSE))
S7::method(distrib_cdf, TransformedDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  tr <- distrib@transformer
  if (tr@decreasing) lower.tail <- !lower.tail
  distrib_cdf(distrib@parent_distrib, tr@trans_inv(q), theta, lower.tail = lower.tail, log.p = log.p)
}

#' @title Transformed Quantile Function
#' @name distrib_quantile.TransformedDistrib
#'
#' @description
#' Computes \eqn{Q_Y(p) = g(Q_X(p))} by evaluating the parent's own quantile
#' function and mapping the result forward. For a DECREASING transformation the
#' tails swap, so `lower.tail` is inverted before the parent is called and the
#' \eqn{p}-th quantile of \eqn{Y} is \eqn{g} of the \eqn{(1-p)}-th of \eqn{X}.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param p A numeric vector of probabilities.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}. For a decreasing transformation this is inverted before
#'   reaching the parent.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm
#'   and passed as such to the parent. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles on the transformed scale.
#'
#' @section Notation:
#' \eqn{g} is the transformation and \eqn{Q_X}, \eqn{Q_Y} the parent's and the
#' transformed quantile function.
#'
#' @seealso [distrib_cdf.TransformedDistrib()], which this inverts, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' p <- c(0.1, 0.5, 0.9)
#'
#' distrib_quantile(d, p, theta)
#' all.equal(distrib_quantile(d, p, theta), qlnorm(p, 0.5, 0.8))
#'
#' # The round trip closes, the transformation being a bijection.
#' max(abs(distrib_cdf(d, distrib_quantile(d, p, theta), theta) - p))
#'
#' # Under a decreasing transformation the p-th quantile is g of the
#' # (1 - p)-th of the parent.
#' ig <- transformation(gamma1_distrib(), inverse_transform())
#' th2 <- list(mu = 2, phi = 0.3)
#' c(transformed = distrib_quantile(ig, 0.9, th2),
#'   mapped = 1 / distrib_quantile(gamma1_distrib(), 0.1, th2))
S7::method(distrib_quantile, TransformedDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  tr <- distrib@transformer
  if (tr@decreasing) lower.tail <- !lower.tail
  tr@trans_fun(distrib_quantile(distrib@parent_distrib, p, theta, lower.tail = lower.tail, log.p = log.p))
}

#' @title Transformed Random Number Generator
#' @name distrib_rng.TransformedDistrib
#'
#' @description
#' Draws `n` values from the parent and applies \eqn{g} to them. That is the
#' DEFINITION of \eqn{Y = g(X)}, not an approximation of it, so the draws are
#' exact whatever route the parent's own generator takes and consume exactly
#' what the parent consumes from R's stream.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param n The number of draws, a single non-negative whole number.
#' @param theta A named list of the parent's parameters.
#'
#' @return A numeric vector of length `n`, on the transformed scale.
#'
#' @seealso [distrib_pdf.TransformedDistrib()] for the density these are drawn
#'   from, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#'
#' set.seed(1)
#' distrib_rng(d, 5, theta)
#'
#' # It is the parent's draw mapped forward, from the same seed.
#' set.seed(1)
#' exp(distrib_rng(gaussian1_distrib(), 5, theta))
#'
#' # And a large sample reproduces the transformed distribution function.
#' set.seed(2)
#' big <- distrib_rng(d, 20000, theta)
#' c(sampled = mean(big < 2), exact = distrib_cdf(d, 2, theta))
S7::method(distrib_rng, TransformedDistrib) <- function(distrib, n, theta) {
  distrib@transformer@trans_fun(distrib_rng(distrib@parent_distrib, n, theta))
}

#' @title Transformed Score
#' @name distrib_gradient.TransformedDistrib
#'
#' @description
#' Returns the parent's score evaluated at \eqn{x = g^{-1}(y)}. The Jacobian
#' does not depend on the parameters, so the score of the transformed model IS
#' the parent's; the change of variables moves where the score is read, not
#' what it is.
#'
#' @details
#' The transformation does not depend on \eqn{\theta}, so
#' \eqn{\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert}
#' has a second term that is a CONSTANT in the parameters and differentiates
#' away. Every derivative in \eqn{\theta} is therefore the parent's own,
#' evaluated at \eqn{x = g^{-1}(y)}, and nothing is recomputed. The response
#' derivatives are the exception and come from the base class, the Jacobian
#' depending on \eqn{y}.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param y A numeric vector of observations on the transformed scale.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The links are the parent's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, in the parent's
#'   order.
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{J} the Jacobian of its inverse and
#' \eqn{\ell} a log-density.
#'
#' @seealso [distrib_hessian.TransformedDistrib()] for the second order,
#'   [distrib_grad_y()], which is NOT the parent's, and [distrib_gradient()]
#'   for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' set.seed(2)
#' y <- distrib_rng(d, 40, theta)
#'
#' g <- distrib_gradient(d, y, theta)
#' vapply(g, sum, numeric(1))
#'
#' # It is literally the parent's score at the preimage.
#' identical(g, distrib_gradient(gaussian1_distrib(), log(y), theta))
#'
#' # And it agrees with a numerical derivative of the TRANSFORMED
#' # log-likelihood, the Jacobian having differentiated away.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
S7::method(distrib_gradient, TransformedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  distrib_gradient(distrib@parent_distrib, distrib@transformer@trans_inv(y), theta)
}

#' @title Transformed Observed Hessian
#' @name distrib_hessian.TransformedDistrib
#'
#' @description
#' Returns the parent's observed Hessian evaluated at \eqn{x = g^{-1}(y)}, for
#' the reason the score has: the Jacobian is a constant in \eqn{\theta} and
#' differentiates away twice over.
#'
#' @details
#' The transformation does not depend on \eqn{\theta}, so
#' \eqn{\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert}
#' has a second term that is a CONSTANT in the parameters and differentiates
#' away. Every derivative in \eqn{\theta} is therefore the parent's own,
#' evaluated at \eqn{x = g^{-1}(y)}, and nothing is recomputed. The response
#' derivatives are the exception and come from the base class, the Jacobian
#' depending on \eqn{y}.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param y A numeric vector of observations on the transformed scale.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{J} the Jacobian of its inverse and
#' \eqn{\ell} a log-density.
#'
#' @seealso [distrib_gradient.TransformedDistrib()] for the first order,
#'   [distrib_expected_hessian.TransformedDistrib()] for the expectation, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' set.seed(2)
#' y <- distrib_rng(d, 40, theta)
#'
#' H <- distrib_hessian(d, y, theta)
#' vapply(H, sum, numeric(1))
#'
#' # The parent's, at the preimage.
#' identical(H, distrib_hessian(gaussian1_distrib(), log(y), theta))
#'
#' # Against a numerical Hessian of the transformed log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' Hn <- numDeriv::hessian(ll, unlist(theta))
#' ref <- vapply(distributions7:::hess_pairs(d@params),
#'               function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
#'               numeric(1))
#' max(abs(vapply(H, sum, numeric(1)) - ref))
S7::method(distrib_hessian, TransformedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  distrib_hessian(distrib@parent_distrib, distrib@transformer@trans_inv(y), theta)
}

#' @title Transformed Expected Information
#' @name distrib_expected_hessian.TransformedDistrib
#'
#' @description
#' Returns the parent's expected Hessian, EXACTLY. Since
#' \eqn{\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert}
#' and the Jacobian does not depend on \eqn{\theta}, the expectation of the
#' transformed model's second derivative is the parent's expectation
#' reparametrized by the change of variables, which is the same number. No
#' approximation runs and no Monte Carlo is needed, whatever `approx` says.
#'
#' @details
#' This is the strongest of the three delegations: the score and the observed
#' Hessian are the parent's AT A POINT, and this one is the parent's as a
#' whole. A family that is a transformation of a family with a closed-form
#' information therefore has one too, at no cost.
#'
#' @param distrib A `TransformedDistrib` object, from [transformation()].
#' @param y A numeric vector of observations on the transformed scale. Read
#'   only through the parent, which uses its length.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored: the delegation is exact. Present so that the
#'   signature matches the generic's. Note that it is not forwarded, so a
#'   parent whose OWN expected Hessian is approximate takes its own default.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{J} the Jacobian of its inverse and
#' \eqn{\ell} a log-density.
#'
#' @seealso [distrib_hessian.TransformedDistrib()] for the observed matrix,
#'   [fit_distrib()], whose Fisher scoring inverts this, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- transformation(gaussian1_distrib(), exp_transform())
#' theta <- list(mu = 0.5, sigma = 0.8)
#' set.seed(2)
#' y <- distrib_rng(d, 40, theta)
#'
#' EH <- distrib_expected_hessian(d, y, theta)
#' vapply(EH, function(z) z[1], numeric(1))
#'
#' # It is the parent's, at the preimage.
#' identical(EH, distrib_expected_hessian(gaussian1_distrib(), log(y), theta))
#'
#' # So a lognormal built this way inherits the gaussian's closed form, whose
#' # mean block is -1 / sigma^2.
#' c(reported = EH$mu_mu[1], theory = -1 / 0.8^2)
S7::method(distrib_expected_hessian, TransformedDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  distrib_expected_hessian(distrib@parent_distrib, distrib@transformer@trans_inv(y), theta)
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Apply a Variable Transformation to a Distribution Object
#'
#' @description
#' Builds the distribution of \eqn{Y = g(X)}, where \eqn{X} follows an existing
#' CONTINUOUS distribution and \eqn{g} is a bijective transformation described
#' by a [transformer()]. The result carries exactly the parent's parameters:
#' the transformation adds none and removes none, and the density is the
#' change-of-variables formula
#' \deqn{f_Y(y) = f_X(g^{-1}(y)) \cdot \left\lvert\frac{d}{dy} g^{-1}(y)\right\rvert.}
#'
#' @details
#' # What comes from where
#'
#' The distribution function, the quantile function and the generator are
#' obtained by mapping through \eqn{g}, with the tails swapped for a decreasing
#' transformation. Since \eqn{g} does not depend on the parameters, the score,
#' the observed Hessian and the expected Hessian COINCIDE with the parent's,
#' evaluated at \eqn{x = g^{-1}(y)}; nothing is recomputed and no accuracy is
#' lost. The moments are available numerically through [moment()].
#'
#' # What is rejected
#'
#' A discrete parent, since a change of variables needs a density; a
#' transformer that is not valid on the parent's support, which the transformer
#' itself decides through its `valid_support`; and a `new_name` that is not a
#' single non-empty string. A map that is not injective, such as the absolute
#' value, cannot be a transformer at all and is [folded()] instead.
#'
#' # Naming the result
#'
#' Several standard families are a transformation of one already here, and
#' `new_name` lets the result carry the name it is known by instead of the
#' recipe that produced it: the reciprocal of a gamma is an inverse gamma, the
#' exponential of a logistic a log-logistic, the exponential of an exponential
#' a Pareto. Only the printed name changes, and nothing about the distribution
#' depends on it.
#'
#' @section Notation:
#' \eqn{g} is the transformation, \eqn{f_X} the parent's density, \eqn{f_Y} the
#' transformed one and \eqn{\theta} the parameters shared by both.
#'
#' @param distrib An object inheriting from `continuous_distrib`. A discrete
#'   one is rejected.
#' @param transformer A [transformer()] object, such as [log_transform()],
#'   `bc_transform(0.5)` or `affine_transform(1, 2)`. It is rejected where its
#'   own `valid_support` says it does not apply to `distrib@bounds`.
#' @param new_name A single non-empty string naming the result, or `NULL`, the
#'   default, to compose the parent's name with the transformer's as
#'   `"g(parent)"`.
#'
#' @return An S7 object of class [TransformedDistrib], carrying
#'   `parent_distrib` and `transformer`. Its `params`, `params_bounds` and
#'   `link_params` are the parent's unchanged; `bounds` is the transformer's
#'   image of the parent's; `params_interpretation` names the parent's scale in
#'   brackets; and `params_smooth` is the parent's, a kink in the parameters
#'   surviving a change of variables in the response untouched.
#'
#' @seealso [transformer()] for the map and the twelve ready-made ones,
#'   [TransformedDistrib] for the class, [folded()] for a map that is not
#'   injective, and [fixed()], [truncated()], [zero_inflated()] and
#'   [zero_adjusted()] for the other wrappers.
#'
#' @examples
#' # A lognormal built by transformation, equal to lognormal1_distrib() with
#' # its second parameter read as a variance.
#' logn <- transformation(gaussian1_distrib(), exp_transform())
#' logn@params
#' c(built = distrib_pdf(logn, 2, list(mu = 0, sigma = 1)), dlnorm = dlnorm(2))
#'
#' # The reciprocal of a gamma is an inverse gamma, and can say so.
#' ig <- transformation(gamma2_distrib(), inverse_transform(),
#'                      new_name = "inverse gamma")
#' ig@distrib_name
#'
#' # The derivatives are the parent's at the preimage, so a fit of the
#' # transformed family is the parent's fit of the transformed data.
#' set.seed(1)
#' y <- distrib_rng(logn, 500, list(mu = 0, sigma = 1))
#' rbind(transformed = coef(fit_distrib(logn, y)),
#'       parent_on_logs = coef(fit_distrib(gaussian1_distrib(), log(y))))
#'
#' # Three refusals, each naming the condition that failed.
#' try(transformation(poisson_distrib(), log_transform()))
#' try(transformation(gaussian1_distrib(), log_transform()))
#' try(transformation(gaussian1_distrib(), exp_transform(), new_name = ""))
#'
#' @export
transformation <- function(distrib, transformer, new_name = NULL) {
  if (!S7::S7_inherits(distrib, continuous_distrib)) {
    stop("transformation() currently supports only continuous distributions.", call. = FALSE)
  }
  if (!S7::S7_inherits(transformer, distributions7::transformer)) {
    stop("Argument 'transformer' must be a 'transformer' object.", call. = FALSE)
  }
  if (!is.null(new_name) &&
      (!is.character(new_name) || length(new_name) != 1L || is.na(new_name) ||
       !nzchar(new_name))) {
    stop("'new_name' must be NULL or a single non-empty string.", call. = FALSE)
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
    distrib_name = if (is.null(new_name)) {
      paste0(transformer@name, "(", distrib@distrib_name, ")")
    } else {
      new_name
    },
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
