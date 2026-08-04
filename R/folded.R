#' @include distrib.R generics.R
NULL

# The helpers this file calls -- log_deriv(), bell_f_ratio(), parent_ell(),
# memo_ratio(), assemble_deriv() -- live in wrapper_derivatives.R and are
# reached at run time, so they need no @include. fold_deriv_k() is the one
# exception, being called at the top level to register a method, and it is
# therefore defined here rather than beside its siblings: an @include on
# wrapper_derivatives.R would pull that file ahead of the wrapper classes it
# registers methods on, and the package would not load.

# Folding at zero is not a change of variable. The map y -> |y| is two to one,
# so it has no inverse and cannot be a transformer(): what the density does is
# add the two preimages,
#   L(x; theta) = f(x; theta) + f(-x; theta),   x >= 0,
# which leaves the parameters where they were and adds no parameter of its own,
# exactly as truncation does.
#
# Every derivative follows from the machinery of wrapper_derivatives.R without
# an addition. With w = f(x)/L,
#   d^B L / L = w (d^B f(x) / f(x)) + (1 - w) (d^B f(-x) / f(-x)),
# both terms being the complete Bell quantity the parent already supplies, one
# at +x and one at -x; and l = log L is then the moment-to-cumulant relation
# over those ratios. At first order this reads w s(x) + (1 - w) s(-x), the
# score of a two-component mixture, which is what a fold is.

#' @title S7 Class for Folded Distributions
#' @name FoldedDistrib
#'
#' @description
#' A subclass of \code{continuous_distrib} representing the distribution of
#' \eqn{|Y|} when \eqn{Y} follows the wrapped distribution.
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{continuous_distrib} object.
#' @return An object of class \code{FoldedDistrib}.
#' @seealso \code{\link{folded}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.FoldedDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.FoldedDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.FoldedDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_gradient.FoldedDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_grad_y.FoldedDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_hessian.FoldedDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_hess_y.FoldedDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_pdf.FoldedDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_rng.FoldedDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}, whose
#' numerical quantile inverts the exact folded distribution function.
FoldedDistrib <- S7::new_class("FoldedDistrib",
  parent = continuous_distrib,
  properties = list(parent_distrib = distrib)
)


#' The Two Preimages of a Folded Point
#'
#' @description
#' The parent's density at \eqn{+x} and at \eqn{-x}, their sum, and the weight
#' \eqn{w = f(x)/L} the first carries.
#'
#' @details
#' Every method of \code{\link{FoldedDistrib}} needs the same four quantities,
#' and computing them once keeps the parent's density from being evaluated
#' twice per call. Points outside the folded support contribute nothing and are
#' returned with a zero density rather than being dropped, so that the result
#' aligns with the input.
#'
#' @param parent The wrapped distribution.
#' @param x A numeric vector of folded observations.
#' @param theta A named list of the parent's parameters.
#'
#' @return A list with \code{fp}, \code{fm}, \code{L} and \code{w}.
#'
#' @seealso \code{\link{folded}}
#'
#' @keywords internal
fold_parts <- function(parent, x, theta) {
  fp <- distrib_pdf(parent, x, theta)
  fm <- distrib_pdf(parent, -x, theta)
  fp[!is.finite(fp)] <- 0
  fm[!is.finite(fm)] <- 0
  L <- fp + fm
  list(fp = fp, fm = fm, L = L, w = ifelse(L > 0, fp / L, 0.5))
}


#' The Block Ratios of a Folded Density
#'
#' @description
#' Returns a memoised function giving \eqn{d^B L / L} for any block, which is
#' what \code{\link{log_deriv}} consumes.
#'
#' @details
#' The ratio is the parent's complete Bell polynomial at each preimage,
#' weighted by which preimage the point came from:
#' \eqn{w\,(d^B f(x)/f(x)) + (1-w)\,(d^B f(-x)/f(-x))}. Both parent
#' evaluations are fetched once and the result memoised, since a partition sum
#' at fourth order asks for the same blocks repeatedly.
#'
#' @param parent The wrapped distribution.
#' @param x A numeric vector of folded observations.
#' @param theta A named list of the parent's parameters.
#' @param order The highest order needed, 1 to 4.
#' @param params The parent's parameter names, in declaration order.
#' @param w The weight of the positive preimage, from
#'   \code{\link{fold_parts}}.
#'
#' @return A function of one block, returning that ratio's vector.
#'
#' @seealso \code{\link{folded}}, \code{\link{fold_parts}}
#'
#' @keywords internal
fold_ratio <- function(parent, x, theta, order, params, w) {
  ell_p <- parent_ell(parent, x, theta, order, params)
  ell_m <- parent_ell(parent, -x, theta, order, params)
  memo_ratio(function(block) {
    w * bell_f_ratio(block, ell_p) + (1 - w) * bell_f_ratio(block, ell_m)
  }, params)
}


#' Derivatives of a Folded Distribution
#'
#' @description
#' Builds the order-\code{k} derivative method for \code{\link{folded}}.
#'
#' @details
#' The ratios handed to \code{\link{log_deriv}} are
#' \eqn{d^B L / L = w\,(d^B f(x)/f(x)) + (1-w)\,(d^B f(-x)/f(-x))} with
#' \eqn{w = f(x)/L}, each term a complete Bell polynomial in the parent's own
#' log-derivatives evaluated at one of the two preimages. Folding adds no
#' parameter, so every index is an index of the parent's.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @seealso \code{\link{folded}}, \code{\link{fold_ratio}}
#'
#' @keywords internal
fold_deriv_k <- function(order) {
  function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
    if (expected) {
      return(expected_derivative(distrib, y, theta, order = order,
                                 approx = match.arg(approx), nsim = nsim))
    }
    parent <- distrib@parent_distrib
    params <- distrib@params
    p <- fold_parts(parent, y, theta)
    ratio <- fold_ratio(parent, y, theta, order, params, p$w)
    assemble_deriv(distrib, order, function(idx) log_deriv(idx, ratio))
  }
}


# --- S7 METHODS IMPLEMENTATION ---

#' @title Folded Density
#' @name distrib_pdf.FoldedDistrib
#' @description
#' \deqn{L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge 0}
#' the two preimages of \eqn{x} under the absolute value added together.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{folded}}
S7::method(distrib_pdf, FoldedDistrib) <- function(distrib, y, theta, log = FALSE) {
  p <- fold_parts(distrib@parent_distrib, y, theta)
  out <- ifelse(y < 0, 0, p$L)
  if (log) log(out) else out
}

#' @title Folded Distribution Function
#' @name distrib_cdf.FoldedDistrib
#' @description
#' \deqn{P(|Y| \le q) = F(q) - F(-q)}
#' @param distrib A \code{FoldedDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(|Y| \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{folded}}
S7::method(distrib_cdf, FoldedDistrib) <- function(distrib, q, theta,
                                                   lower.tail = TRUE,
                                                   log.p = FALSE) {
  parent <- distrib@parent_distrib
  p <- distrib_cdf(parent, q, theta) - distrib_cdf(parent, -q, theta)
  p[q < 0] <- 0
  p <- pmin(pmax(p, 0), 1)
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

#' @title Folded Random Generation
#' @name distrib_rng.FoldedDistrib
#' @description Draws from the parent and takes the absolute value, which is
#'   the definition rather than an approximation of it.
#' @param distrib A \code{FoldedDistrib} object.
#' @param n The number of draws.
#' @param theta A named list of the parent's parameters.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{folded}}
S7::method(distrib_rng, FoldedDistrib) <- function(distrib, n, theta) {
  abs(distrib_rng(distrib@parent_distrib, n, theta))
}

#' @title Folded Analytical Gradient
#' @name distrib_gradient.FoldedDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial \theta_i}
#'       = w\, s_i(x) + (1-w)\, s_i(-x), \qquad w = \dfrac{f(x)}{L(x)}}
#' the score of a two-component mixture, the components being the two
#' preimages.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list, one component per parameter.
#' @seealso \code{\link{folded}}
S7::method(distrib_gradient, FoldedDistrib) <- function(distrib, y, theta,
                                                        scale = c("parameter", "link"), ...) {
  parent <- distrib@parent_distrib
  params <- distrib@params
  p <- fold_parts(parent, y, theta)
  gp <- distrib_gradient(parent, y, theta)
  gm <- distrib_gradient(parent, -y, theta)
  stats::setNames(
    lapply(params, function(nm) p$w * gp[[nm]] + (1 - p$w) * gm[[nm]]),
    params
  )
}

#' @title Folded Analytical Observed Hessian
#' @name distrib_hessian.FoldedDistrib
#' @description
#' The moment-to-cumulant relation applied to the ratios
#' \eqn{d^B L / L = w\,(d^B f(x)/f(x)) + (1-w)\,(d^B f(-x)/f(-x))}, which at
#' second order is the familiar mixture Hessian.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list, one component per pair of parameters.
#' @seealso \code{\link{folded}}
S7::method(distrib_hessian, FoldedDistrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"), ...) {
  parent <- distrib@parent_distrib
  params <- distrib@params
  p <- fold_parts(parent, y, theta)
  ratio <- fold_ratio(parent, y, theta, 2L, params, p$w)
  nm <- hess_names(params)
  pairs <- hess_pairs(params)
  stats::setNames(
    lapply(seq_along(nm), function(k) {
      log_deriv(c(pairs[[k]][1], pairs[[k]][2]), ratio)
    }), nm
  )
}

#' @title Folded Analytical Third-Order Derivatives
#' @name distrib_deriv3.FoldedDistrib
#' @description Third-order derivatives from the same partition sums as the
#'   Hessian; see \code{\link{distrib_hessian.FoldedDistrib}}.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expectation is approximated.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx How the expectation is approximated when requested.
#' @param nsim Monte Carlo sample size.
#' @param ... Unused.
#' @return A named list of third-order components.
#' @seealso \code{\link{folded}}
S7::method(distrib_deriv3, FoldedDistrib) <- fold_deriv_k(3L)

#' @title Folded Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.FoldedDistrib
#' @description Fourth-order derivatives from the same partition sums as the
#'   Hessian; see \code{\link{distrib_hessian.FoldedDistrib}}.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param expected Logical; if \code{TRUE}, the expectation is approximated.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx How the expectation is approximated when requested.
#' @param nsim Monte Carlo sample size.
#' @param ... Unused.
#' @return A named list of fourth-order components.
#' @seealso \code{\link{folded}}
S7::method(distrib_deriv4, FoldedDistrib) <- fold_deriv_k(4L)

#' @title Folded Response Gradient
#' @name distrib_grad_y.FoldedDistrib
#' @description
#' \deqn{\dfrac{\partial \ell}{\partial x} = w\, g(x) - (1-w)\, g(-x)}
#' with \eqn{g = \partial \log f/\partial y}. The minus sign is the chain rule
#' through the reflected preimage.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{folded}}
S7::method(distrib_grad_y, FoldedDistrib) <- function(distrib, y, theta, ...) {
  parent <- distrib@parent_distrib
  p <- fold_parts(parent, y, theta)
  p$w * distrib_grad_y(parent, y, theta) -
    (1 - p$w) * distrib_grad_y(parent, -y, theta)
}

#' @title Folded Response Hessian
#' @name distrib_hess_y.FoldedDistrib
#' @description
#' \deqn{\dfrac{\partial^2 \ell}{\partial x^2}
#'       = w\left(h(x) + g(x)^2\right) + (1-w)\left(h(-x) + g(-x)^2\right)
#'         - \left(w g(x) - (1-w) g(-x)\right)^2}
#' with \eqn{g} and \eqn{h} the parent's first and second response
#' derivatives.
#' @param distrib A \code{FoldedDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{folded}}
S7::method(distrib_hess_y, FoldedDistrib) <- function(distrib, y, theta, ...) {
  parent <- distrib@parent_distrib
  p <- fold_parts(parent, y, theta)
  gp <- distrib_grad_y(parent, y, theta)
  gm <- distrib_grad_y(parent, -y, theta)
  hp <- distrib_hess_y(parent, y, theta)
  hm <- distrib_hess_y(parent, -y, theta)
  first <- p$w * gp - (1 - p$w) * gm
  p$w * (hp + gp^2) + (1 - p$w) * (hm + gm^2) - first^2
}

#' Does a Distribution Declare Atoms
#'
#' @description
#' Answers whether a distribution registers \code{\link{distrib_atoms}} for
#' itself rather than inheriting the base class's empty answer.
#'
#' @details
#' The question is asked of the class rather than of a value, since a
#' constructor has no \code{theta} to evaluate the generic at, and the answer
#' is structural in any case. The class a method was registered on is what
#' settles it: \code{identical()} on the method object would not, S7 wrapping
#' it.
#'
#' The argument is named \code{parent} deliberately. The base class of this
#' package is called \code{distrib}, so an argument of that name would shadow
#' it and the comparison below would test the class against the distribution
#' object instead -- answering \code{TRUE} for every family, which is exactly
#' what it did before this helper existed.
#'
#' @param parent A distribution object.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{folded}}, \code{\link{distrib_atoms}}
#'
#' @keywords internal
declares_atoms <- function(parent) {
  m <- tryCatch(S7::method(distrib_atoms, S7::S7_class(parent)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  reg <- tryCatch(attr(m, "signature")[[1]], error = function(e) NULL)
  !is.null(reg) && !is_class(reg, distrib)
}

# --- CONSTRUCTOR WRAPPER ---

#' Fold a Distribution at Zero
#'
#' @description
#' Wraps a continuous distribution into the distribution of the absolute value
#' of its variable.
#'
#' @param distrib A \code{continuous_distrib} object whose support reaches
#'   below zero.
#'
#' @details
#' Folding is not a change of variable. The map \eqn{y \mapsto |y|} is two to
#' one and has no inverse, so it cannot be a
#' \code{\link{transformer}}: instead of carrying a density through a Jacobian
#' it adds the two preimages,
#' \deqn{L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge 0,}
#' with distribution function \eqn{F(x) - F(-x)}. No parameter is added and
#' none is removed, as with \code{\link{truncated}}.
#'
#' Every derivative comes from the wrapper machinery unchanged. Writing
#' \eqn{w = f(x)/L(x)} for the weight of the positive preimage, the ratios that
#' machinery consumes are
#' \deqn{\dfrac{d^B L}{L} = w \dfrac{d^B f(x)}{f(x)}
#'       + (1-w) \dfrac{d^B f(-x)}{f(-x)},}
#' each term a complete Bell polynomial in the parent's own log-derivatives,
#' one evaluated at \eqn{+x} and one at \eqn{-x}; \eqn{\log L} then follows
#' from the moment-to-cumulant relation. At first order this is
#' \eqn{w s(x) + (1-w) s(-x)}, the score of a two-component mixture, which is
#' what a fold is.
#'
#' The parent must reach below zero. A distribution already supported on the
#' non-negative half line folds to itself, and returning it unchanged would
#' hide a mistaken call rather than report it; the same check makes
#' \code{folded()} of a folded distribution an error.
#'
#' A parent with an atom is refused as well. The point zero is its own
#' preimage while every other point has two, so an atom at zero would be
#' counted twice by the sum above and an atom elsewhere would be moved onto its
#' reflection.
#'
#' \strong{The half-normal} is \code{fixed(folded(gaussian_distrib()), mu = 0)},
#' and the folded normal proper is \code{folded(gaussian_distrib())}.
#'
#' \strong{The sign of a symmetric parent's location is not identified.} When
#' the parent is symmetric about its location, \eqn{f(-x; \mu) = f(x; -\mu)},
#' so the two terms of \eqn{L} merely swap and the likelihood is
#' \emph{exactly} even in \eqn{\mu}: a fit converges to \eqn{+\hat\mu} or
#' \eqn{-\hat\mu} according to where it started, at the same maximised value
#' to every digit. This is a property of the model rather than of the fitting,
#' and it is not refused, the folded normal being a standard family; what is
#' estimable is \eqn{|\mu|} together with the remaining parameters. Holding
#' the location at zero removes the question and gives the half-normal. A
#' parent that is not symmetric about its location, such as
#' \code{\link{skewnormal_distrib}}, has no such invariance and its sign is
#' identified.
#'
#' @return An S7 object of class \code{\link{FoldedDistrib}}.
#'
#' @seealso \code{\link{truncated}}, \code{\link{fixed}},
#'   \code{\link{transformation}}
#'
#' @examples
#' d <- folded(gaussian_distrib())
#' theta <- list(mu = 0.5, sigma = 1)
#' distrib_pdf(d, c(0, 0.5, 2), theta)
#'
#' # the half-normal: a folded gaussian with its location held at zero
#' hn <- fixed(folded(gaussian_distrib()), mu = 0)
#' hn@params
#' distrib_pdf(hn, c(0.5, 1), list(sigma = 2))
#'
#' @export
folded <- function(distrib) {
  if (!S7::S7_inherits(distrib, continuous_distrib)) {
    stop("folded() supports continuous distributions only.", call. = FALSE)
  }
  if (distrib@bounds[1] >= 0) {
    stop(sprintf(paste0(
      "'%s' is supported on [%s, %s], which the absolute value leaves alone.\n",
      "  folded() would return the same distribution, so the call is a mistake\n",
      "  rather than a no-op and is reported as one."
    ), distrib@distrib_name, format(distrib@bounds[1]),
    format(distrib@bounds[2])), call. = FALSE)
  }
  if (declares_atoms(distrib)) {
    stop(sprintf(paste0(
      "'%s' carries an atom, and folding would misplace it: zero is its own\n",
      "  preimage while every other point has two, so an atom at zero would be\n",
      "  counted twice and one elsewhere moved onto its reflection."
    ), distrib@distrib_name), call. = FALSE)
  }

  FoldedDistrib(
    parent_distrib = distrib,
    distrib_name = sprintf("folded %s", distrib@distrib_name),
    dimension = distrib@dimension,
    bounds = c(0, max(abs(distrib@bounds))),
    params = distrib@params,
    params_interpretation = distrib@params_interpretation,
    n_params = distrib@n_params,
    params_bounds = distrib@params_bounds,
    link_params = distrib@link_params,
    params_smooth = param_smoothness(distrib)
  )
}
