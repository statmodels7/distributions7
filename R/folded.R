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
#' The S7 class of the distribution of \eqn{|Y|} when \eqn{Y} follows the
#' wrapped continuous distribution, with density
#' \deqn{L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge 0,}
#' the two preimages of \eqn{x} added together. It inherits from
#' `continuous_distrib` and carries the SAME parameters as its parent: folding
#' adds none and removes none.
#'
#' Build one with [folded()], which checks that the parent reaches below zero
#' and carries no atom. This page documents the raw S7 constructor, which
#' validates neither.
#'
#' @param parent_distrib The wrapped `continuous_distrib` object.
#' @param ... The properties of the parent `distrib` class, listed under Value.
#'
#' @return An S7 object of class `FoldedDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. It carries `parent_distrib`
#'   beside the parent's `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params` and
#'   `params_smooth`. For an object built by [folded()] the parameters are the
#'   wrapped family's, `bounds` becomes `c(0, max(abs(parent bounds)))`, and
#'   `distrib_name` is `"folded "` followed by the parent's.
#'
#' @seealso [folded()] to build one, [truncated()] for the other wrapper that
#'   adds no parameter, and [fixed()], which with `mu = 0` turns a folded
#'   gaussian into the half-normal.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.FoldedDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.FoldedDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.FoldedDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.FoldedDistrib],
#'   [`distrib_gradient()`][distrib_gradient.FoldedDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.FoldedDistrib],
#'   [`distrib_hessian()`][distrib_hessian.FoldedDistrib],
#'   [`distrib_pdf()`][distrib_pdf.FoldedDistrib],
#'   [`distrib_rng()`][distrib_rng.FoldedDistrib]
#'
#' Everything else is inherited from [continuous_distrib()], whose numerical
#' quantile inverts the exact folded distribution function above.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The parameters are the parent's, unchanged, and the support is the
#' # non-negative half line.
#' d@params
#' d@bounds
#' d@distrib_name
#'
#' # The parent is kept, so a method can reach it at both preimages.
#' d@parent_distrib@distrib_name
FoldedDistrib <- S7::new_class("FoldedDistrib",
  parent = continuous_distrib,
  properties = list(parent_distrib = distrib)
)


#' @title The Two Preimages of a Folded Point
#'
#' @description
#' Computes the parent's density at \eqn{+x} and at \eqn{-x}, their sum, and
#' the weight \eqn{w = f(x)/L(x)} of the positive preimage. Every method of
#' [FoldedDistrib] needs the same four quantities, so they are formed once per
#' call: the density is \eqn{L}, the score is \eqn{w s(x) + (1-w) s(-x)}, and
#' every higher-order ratio is a \eqn{w}-weighted average of the parent's own.
#'
#' @details
#' At \eqn{x = 0} the two preimages coincide, so \eqn{w = 1/2} whatever the
#' parameters are and the folded density is exactly twice the parent's. Far out
#' in the tail of a parent centered above zero, \eqn{w} approaches one and the
#' fold becomes invisible.
#'
#' @param parent The wrapped `continuous_distrib` object.
#' @param x A numeric vector of points at which to evaluate. Negative values
#'   are not screened out here; the calling method zeroes them.
#' @param theta A named list of the parent's parameters.
#'
#' @return A named list of four numeric vectors of the recycled length: `fp`
#'   and `fm`, the parent's density at \eqn{+x} and \eqn{-x}; `L`, their sum;
#'   and `w`, the ratio `fp / L`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{w} the weight
#' of the positive preimage and \eqn{s} the parent's score.
#'
#' @seealso [fold_ratio()] for the higher-order quantities built on these, and
#'   [distrib_pdf.FoldedDistrib()] for the first consumer.
#'
#' @examples
#' g <- gaussian1_distrib()
#' theta <- list(mu = 0.5, sigma = 1.2)
#' p <- distributions7:::fold_parts(g, c(0, 1, 3), theta)
#' str(p)
#'
#' # L is the sum of the two preimages, which is the folded density.
#' all.equal(p$L, dnorm(c(0, 1, 3), 0.5, 1.2) + dnorm(-c(0, 1, 3), 0.5, 1.2))
#'
#' # At zero the two preimages coincide, so the weight is exactly one half.
#' p$w[1]
#'
#' # And far into the tail the fold becomes invisible: w approaches one.
#' distributions7:::fold_parts(g, c(3, 6, 12), theta)$w
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


#' @title The Block Ratios of a Folded Density
#'
#' @description
#' Returns a memoized function giving \eqn{d^B L / L} for any block \eqn{B} of
#' parameter names. That is what [log_deriv()] consumes to assemble a
#' derivative of \eqn{\log L} of any order. The ratio is the parent's complete
#' Bell polynomial evaluated at each preimage and averaged with the weight
#' \eqn{w},
#' \deqn{\frac{d^B L}{L} = w \frac{d^B f(x)}{f(x)}
#'   + (1-w) \frac{d^B f(-x)}{f(-x)}.}
#'
#' @details
#' The memoization matters because a partition sum asks for the same block many
#' times: at order four the same singleton block appears in most of the
#' partitions. Each distinct block costs two evaluations of the parent's
#' derivatives, one at \eqn{+x} and one at \eqn{-x}, and no more.
#'
#' @param parent The wrapped `continuous_distrib` object.
#' @param x A numeric vector of points.
#' @param theta A named list of the parent's parameters.
#' @param order The highest derivative order the caller will ask for, a single
#'   whole number from 1 to 4. Blocks longer than this are not prepared.
#' @param params The parameter names, as `distrib@params`.
#' @param w The weight of the positive preimage, from [fold_parts()].
#'
#' @return A function of one character vector, the block, returning a numeric
#'   vector of the recycled length of `x`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{B} a multiset
#' of parameter names, \eqn{d^B} the corresponding partial derivative and
#' \eqn{w} the weight of the positive preimage.
#'
#' @seealso [fold_parts()] for the weight, [log_deriv()] for the partition sum
#'   this feeds, and [bell_f_ratio()] for the parent's own Bell quantity.
#'
#' @examples
#' g <- gaussian1_distrib()
#' theta <- list(mu = 0.5, sigma = 1.2)
#' x <- c(0.5, 1, 3)
#' w <- distributions7:::fold_parts(g, x, theta)$w
#' r <- distributions7:::fold_ratio(g, x, theta, 2L, c("mu", "sigma"), w)
#'
#' # A singleton block is the weighted average of the parent's own score.
#' sp <- distrib_gradient(g, x, theta)$mu
#' sm <- distrib_gradient(g, -x, theta)$mu
#' all.equal(r("mu"), w * sp + (1 - w) * sm)
#'
#' # Which at first order IS the folded score, log L having no other term.
#' all.equal(r("mu"), distrib_gradient(folded(g), x, theta)$mu)
#'
#' # A two-element block is the second-order ratio, not the second derivative
#' # of the logarithm: the two differ by the square of the first.
#' c(ratio = r(c("mu", "mu"))[2],
#'   log_second = distrib_hessian(folded(g), x, theta)$mu_mu[2],
#'   difference = r(c("mu", "mu"))[2] - r("mu")[2]^2)
#'
#' @keywords internal
fold_ratio <- function(parent, x, theta, order, params, w) {
  ell_p <- parent_ell(parent, x, theta, order, params)
  ell_m <- parent_ell(parent, -x, theta, order, params)
  memo_ratio(function(block) {
    w * bell_f_ratio(block, ell_p) + (1 - w) * bell_f_ratio(block, ell_m)
  }, params)
}


#' @title Build a Folded Derivative Method of a Given Order
#'
#' @description
#' Returns the method [folded()] registers for `distrib_deriv3()` or
#' `distrib_deriv4()`. Both orders run the same three steps, so the body is
#' written once and the order is closed over: form the weight with
#' [fold_parts()], build the block ratios with [fold_ratio()], and hand them to
#' [log_deriv()] for every component of the order.
#'
#' @details
#' The ratios are \eqn{d^B L/L = w\,(d^B f(x)/f(x)) + (1-w)\,(d^B f(-x)/f(-x))},
#' and [log_deriv()] turns them into derivatives of \eqn{\log L} by the
#' moment-to-cumulant relation. Nothing about the fold is written out at third
#' or fourth order; the same two functions serve every order the parent
#' supplies.
#'
#' @param order The derivative order, `3L` or `4L`.
#'
#' @return A function with the signature of `distrib_deriv3()`, suitable for
#'   `S7::method(...) <- `.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{w} the weight
#' of the positive preimage and \eqn{B} a multiset of parameter names.
#'
#' @seealso [distrib_deriv3.FoldedDistrib()] and
#'   [distrib_deriv4.FoldedDistrib()], the two methods it builds, and
#'   [log_deriv()] for the partition sum.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' set.seed(1)
#' y <- distrib_rng(d, 6, theta)
#'
#' m3 <- distributions7:::fold_deriv_k(3L)
#'
#' # It builds the method the class registers, so the two agree.
#' all.equal(m3(d, y, theta), distrib_deriv3(d, y, theta))
#'
#' # And the order it was built at fixes the component count.
#' c(order3 = length(m3(d, y, theta)),
#'   order4 = length(distributions7:::fold_deriv_k(4L)(d, y, theta)))
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
#'
#' @description
#' Computes the folded density
#' \deqn{L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge 0,}
#' the parent's density at the two preimages of \eqn{x} under the absolute
#' value, added together. Below zero the folded variable places no mass, so the
#' density is exactly `0` and its logarithm `-Inf`.
#'
#' @details
#' At \eqn{x = 0} the two preimages coincide and the folded density is twice
#' the parent's. The addition is what separates a fold from a change of
#' variable: the map is two to one, so there is no Jacobian to divide by, and
#' [transformer()] cannot express it.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations. Negative values return `0`
#'   without an error, the folded support being \eqn{[0, \infty)}.
#' @param theta A named list of the PARENT's parameters, each a numeric vector
#'   of length 1 or of the length of `y`. Folding adds and removes nothing.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   The logarithm is taken of the sum, not inside the parent, so it underflows
#'   where the sum does. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one and \eqn{\theta} the
#' parameters shared by both.
#'
#' @seealso [distrib_cdf.FoldedDistrib()] for the distribution function,
#'   [distrib_gradient.FoldedDistrib()] for the score, [folded()] for the
#'   family, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#'
#' distrib_pdf(d, c(0, 0.5, 2), theta)
#'
#' # Which is the parent's density at the two preimages, added.
#' y <- c(0, 0.5, 2)
#' all.equal(distrib_pdf(d, y, theta),
#'           dnorm(y, 0.5, 1.2) + dnorm(-y, 0.5, 1.2))
#'
#' # At zero the preimages coincide, so the density is twice the parent's.
#' c(folded = distrib_pdf(d, 0, theta), twice = 2 * dnorm(0, 0.5, 1.2))
#'
#' # Below zero there is no mass.
#' distrib_pdf(d, c(-1, -0.1), theta)
#'
#' # And it integrates to one over the half line.
#' integrate(function(z) distrib_pdf(d, z, theta), 0, Inf)$value
S7::method(distrib_pdf, FoldedDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  p <- fold_parts(distrib@parent_distrib, y, theta)
  out <- ifelse(y < 0, 0, p$L)
  if (log) log(out) else out
}

#' @title Folded Distribution Function
#' @name distrib_cdf.FoldedDistrib
#'
#' @description
#' Computes \eqn{P(|Y| \le q) = F(q) - F(-q)} from the parent's own
#' distribution function, exactly and with no quadrature. It is the difference
#' of two calls on the parent, clamped to \eqn{[0, 1]} against rounding and set
#' to zero below the support.
#'
#' @details
#' `lower.tail = FALSE` and `log.p = TRUE` are formed from the computed
#' probability, as \eqn{1 - p} and \eqn{\log p}, not by a separate route
#' through the parent. Far into the upper tail that difference cancels, so a
#' caller who needs the survival function of a folded distribution to many
#' digits should expect the loss the subtraction implies.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param q A numeric vector of quantiles. Values below zero give probability
#'   `0`.
#' @param theta A named list of the parent's parameters.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(|Y| \le q)}; when `FALSE` they are
#'   \eqn{P(|Y| > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm is returned.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities, in \eqn{[0, 1]}.
#'
#' @section Notation:
#' \eqn{F} is the parent's distribution function and \eqn{Y} the parent's
#' variable.
#'
#' @seealso [distrib_pdf.FoldedDistrib()] for the density,
#'   [distrib_quantile()], which the class inherits and which inverts this
#'   exactly, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#'
#' q <- c(0, 1, 3)
#' distrib_cdf(d, q, theta)
#'
#' # Which is the difference of the parent's own distribution function.
#' all.equal(distrib_cdf(d, q, theta),
#'           pnorm(q, 0.5, 1.2) - pnorm(-q, 0.5, 1.2))
#'
#' # Both tails and the logarithm.
#' distrib_cdf(d, 1, theta, lower.tail = FALSE)
#' distrib_cdf(d, 1, theta, log.p = TRUE)
#'
#' # The inherited quantile inverts it, so the round trip closes.
#' p <- c(0.1, 0.5, 0.9)
#' max(abs(distrib_cdf(d, distrib_quantile(d, p, theta), theta) - p))
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
#'
#' @description
#' Draws `n` values from the parent and takes the absolute value. That is the
#' DEFINITION of the folded variable, not an approximation of it, so the draws
#' are exact whatever route the parent's own generator takes. They consume
#' whatever the parent consumes from R's stream, and are reproducible under
#' [base::set.seed()].
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param n The number of draws, a single non-negative whole number.
#' @param theta A named list of the parent's parameters.
#'
#' @return A numeric vector of length `n`, every value non-negative.
#'
#' @seealso [distrib_pdf.FoldedDistrib()] for the density these are drawn
#'   from, [folded()] for the family, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#'
#' set.seed(1)
#' distrib_rng(d, 5, theta)
#'
#' # It is the parent's draw with the sign removed, from the same seed.
#' set.seed(1)
#' abs(distrib_rng(gaussian1_distrib(), 5, theta))
#'
#' # And a large sample reproduces the folded density.
#' set.seed(2)
#' big <- distrib_rng(d, 20000, theta)
#' c(sampled = mean(big < 1),
#'   exact = distrib_cdf(d, 1, theta))
S7::method(distrib_rng, FoldedDistrib) <- function(distrib, n, theta) {
  abs(distrib_rng(distrib@parent_distrib, n, theta))
}

#' @title Folded Score
#' @name distrib_gradient.FoldedDistrib
#'
#' @description
#' Computes the first derivatives of the folded log-density,
#' \deqn{\frac{\partial \ell}{\partial \theta_i}
#'       = w\, s_i(x) + (1-w)\, s_i(-x), \qquad w = \frac{f(x)}{L(x)},}
#' the score of a TWO-COMPONENT MIXTURE whose components are the two preimages.
#' It needs nothing of the fold beyond the weight: \eqn{s} is the parent's own
#' score, evaluated at \eqn{+x} and at \eqn{-x}.
#'
#' @details
#' The weight is where the folded score departs from the parent's. Near zero
#' \eqn{w} is one half and the two preimages contribute equally; far into the
#' tail of a parent centered above zero \eqn{w} approaches one and the folded
#' score approaches the parent's.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The links are the parent's, so the two scales
#'   differ exactly where the parent's do.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, in
#'   `distrib@params` order, which is the parent's order.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{s} the
#' parent's score, \eqn{w} the weight of the positive preimage and \eqn{\ell}
#' the folded log-density.
#'
#' @seealso [distrib_hessian.FoldedDistrib()] for the second order,
#'   [fold_parts()] for the weight, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' set.seed(2)
#' y <- distrib_rng(d, 30, theta)
#'
#' g <- distrib_gradient(d, y, theta)
#' vapply(g, sum, numeric(1))
#'
#' # Against a numerical derivative of the folded log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#'
#' # It is the mixture score, written out from the parent's own.
#' p <- distributions7:::fold_parts(gaussian1_distrib(), y, theta)
#' sp <- distrib_gradient(gaussian1_distrib(), y, theta)$mu
#' sm <- distrib_gradient(gaussian1_distrib(), -y, theta)$mu
#' all.equal(g$mu, p$w * sp + (1 - p$w) * sm)
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

#' @title Folded Observed Hessian
#' @name distrib_hessian.FoldedDistrib
#'
#' @description
#' Computes the second derivatives of the folded log-density by the
#' moment-to-cumulant relation applied to the block ratios
#' \deqn{\frac{d^B L}{L} = w \frac{d^B f(x)}{f(x)}
#'   + (1-w) \frac{d^B f(-x)}{f(-x)},}
#' which at second order is the familiar mixture Hessian: the weighted average
#' of the two components' second-order ratios, less the square of the weighted
#' average of their first-order ones.
#'
#' @details
#' No formula is written out for the fold. [fold_ratio()] supplies the ratios
#' and [log_deriv()] turns them into derivatives of \eqn{\log L}, the same two
#' functions the third and fourth orders use. Each ratio is a complete Bell
#' polynomial in the parent's own log-derivatives, so a parent with closed
#' forms gives the folded family closed forms.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list of the parent's parameters.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{B} a multiset
#' of parameter names, \eqn{w} the weight of the positive preimage and
#' \eqn{\ell} the folded log-density.
#'
#' @seealso [distrib_gradient.FoldedDistrib()] for the first order,
#'   [distrib_deriv3.FoldedDistrib()] for the third, [fold_ratio()] for the
#'   ratios, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' set.seed(2)
#' y <- distrib_rng(d, 30, theta)
#'
#' H <- distrib_hessian(d, y, theta)
#' vapply(H, sum, numeric(1))
#'
#' # Against a numerical Hessian of the folded log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' Hn <- numDeriv::hessian(ll, unlist(theta))
#' ref <- vapply(distributions7:::hess_pairs(d@params),
#'               function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
#'               numeric(1))
#' max(abs(vapply(H, sum, numeric(1)) - ref))
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

#' @title Folded Third Derivatives
#' @name distrib_deriv3.FoldedDistrib
#'
#' @description
#' Computes every third derivative of the folded log-density from the SAME
#' partition sums the Hessian uses, one order up. The block ratios are
#' \deqn{\frac{d^B L}{L} = w \frac{d^B f(x)}{f(x)}
#'   + (1-w) \frac{d^B f(-x)}{f(-x)},}
#' each term a complete Bell polynomial in the parent's own log-derivatives,
#' and [log_deriv()] turns them into derivatives of \eqn{\log L} by the
#' moment-to-cumulant relation. A parent with closed forms to third order gives
#' the folded family closed forms to third order.
#'
#' @details
#' With `expected = TRUE` the expectation goes to [expected_derivative()],
#' which reads `approx` and `nsim`; there is no closed-form expectation for a
#' fold in general.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list of the parent's parameters.
#' @param expected Logical of length 1. When `TRUE` the expectation is
#'   approximated by [expected_derivative()]. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx The strategy [expected_derivative()] uses when
#'   `expected = TRUE`: one of `"integrate"` (the default), `"bartlett"`,
#'   `"mc"` or `"opg"`. Ignored otherwise.
#' @param nsim The Monte Carlo sample size when `approx = "mc"`. Defaults to
#'   `10000`. Ignored otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, keyed and ordered as
#'   `deriv_names(distrib@params, 3)`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{B} a multiset
#' of parameter names and \eqn{w} the weight of the positive preimage.
#'
#' @seealso [distrib_deriv4.FoldedDistrib()] for the next order,
#'   [distrib_hessian.FoldedDistrib()] for the second, [fold_deriv_k()] for the
#'   shared body, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' set.seed(2)
#' y <- distrib_rng(d, 30, theta)
#'
#' d3 <- distrib_deriv3(d, y, theta)
#' names(d3)
#'
#' # Against one stencil on the analytic Hessian.
#' h <- 1e-4
#' tp <- theta; tp$sigma <- tp$sigma + h
#' tm <- theta; tm$sigma <- tm$sigma - h
#' c(exact = sum(d3[["mu_mu_sigma"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["mu_mu"]]) -
#'              sum(distrib_hessian(d, y, tm)[["mu_mu"]])) / (2 * h))
S7::method(distrib_deriv3, FoldedDistrib) <- fold_deriv_k(3L)

#' @title Folded Fourth Derivatives
#' @name distrib_deriv4.FoldedDistrib
#'
#' @description
#' Computes every fourth derivative of the folded log-density from the same
#' partition sums as [distrib_deriv3.FoldedDistrib()], one order up: the block
#' ratios \eqn{d^B L/L = w\,(d^B f(x)/f(x)) + (1-w)\,(d^B f(-x)/f(-x))} handed
#' to [log_deriv()]. A parent with closed forms to fourth order gives the
#' folded family closed forms to fourth order.
#'
#' @details
#' With `expected = TRUE` the expectation goes to [expected_derivative()],
#' which reads `approx` and `nsim`.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list of the parent's parameters.
#' @param expected Logical of length 1. When `TRUE` the expectation is
#'   approximated by [expected_derivative()]. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx The strategy [expected_derivative()] uses when
#'   `expected = TRUE`. Ignored otherwise.
#' @param nsim The Monte Carlo sample size when `approx = "mc"`. Defaults to
#'   `10000`. Ignored otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, keyed and ordered as
#'   `deriv_names(distrib@params, 4)`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{L} the folded one, \eqn{B} a multiset
#' of parameter names and \eqn{w} the weight of the positive preimage.
#'
#' @seealso [distrib_deriv3.FoldedDistrib()] for the order below,
#'   [fold_deriv_k()] for the shared body, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' set.seed(2)
#' y <- distrib_rng(d, 30, theta)
#'
#' d4 <- distrib_deriv4(d, y, theta)
#' length(d4)
#'
#' # Against one stencil on the analytic third order.
#' h <- 1e-4
#' tp <- theta; tp$sigma <- tp$sigma + h
#' tm <- theta; tm$sigma <- tm$sigma - h
#' c(exact = sum(d4[["mu_mu_sigma_sigma"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["mu_mu_sigma"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["mu_mu_sigma"]])) / (2 * h))
S7::method(distrib_deriv4, FoldedDistrib) <- fold_deriv_k(4L)

#' @title Folded Response Gradient
#' @name distrib_grad_y.FoldedDistrib
#'
#' @description
#' Computes the derivative of the folded log-density in the response,
#' \deqn{\frac{\partial \ell}{\partial x} = w\, g(x) - (1-w)\, g(-x),}
#' with \eqn{g = \partial\log f/\partial y} the parent's own response gradient
#' and \eqn{w} the weight of the positive preimage. The MINUS sign is the chain
#' rule through the reflected preimage: moving \eqn{x} up moves \eqn{-x} down.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list of the parent's parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{g} its response gradient, \eqn{w} the
#' weight of the positive preimage and \eqn{\ell} the folded log-density.
#'
#' @seealso [distrib_hess_y.FoldedDistrib()] for the second order,
#'   [distrib_gradient.FoldedDistrib()] for the derivatives in the parameters,
#'   and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' y <- c(0.2, 1, 3)
#'
#' distrib_grad_y(d, y, theta)
#'
#' # Against a numerical derivative of the folded log-density.
#' vapply(y, function(v)
#'   numDeriv::grad(function(z) distrib_pdf(d, z, theta, log = TRUE), v),
#'   numeric(1))
#'
#' # The minus sign written out: the reflected preimage pushes the other way.
#' g0 <- gaussian1_distrib()
#' w <- distributions7:::fold_parts(g0, y, theta)$w
#' all.equal(distrib_grad_y(d, y, theta),
#'           w * distrib_grad_y(g0, y, theta) -
#'             (1 - w) * distrib_grad_y(g0, -y, theta))
S7::method(distrib_grad_y, FoldedDistrib) <- function(distrib, y, theta, ...) {
  parent <- distrib@parent_distrib
  p <- fold_parts(parent, y, theta)
  p$w * distrib_grad_y(parent, y, theta) -
    (1 - p$w) * distrib_grad_y(parent, -y, theta)
}

#' @title Folded Response Hessian
#' @name distrib_hess_y.FoldedDistrib
#'
#' @description
#' Computes the second derivative of the folded log-density in the response,
#' \deqn{\frac{\partial^2 \ell}{\partial x^2}
#'       = w\left(h(x) + g(x)^2\right) + (1-w)\left(h(-x) + g(-x)^2\right)
#'         - \left(w g(x) - (1-w) g(-x)\right)^2,}
#' with \eqn{g} and \eqn{h} the parent's first and second response derivatives.
#' The last term is the square of [distrib_grad_y.FoldedDistrib()]'s value, so
#' the expression is the mixture's second-order ratio less the square of its
#' first: the same moment-to-cumulant step the parameter derivatives take,
#' applied in the response.
#'
#' @param distrib A `FoldedDistrib` object, from [folded()].
#' @param y A numeric vector of observations, non-negative.
#' @param theta A named list of the parent's parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the recycled length of `y` and `theta`.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{g} and \eqn{h} its first and second
#' response derivatives, \eqn{w} the weight of the positive preimage and
#' \eqn{\ell} the folded log-density.
#'
#' @seealso [distrib_grad_y.FoldedDistrib()] for the first order,
#'   [distrib_hessian.FoldedDistrib()] for the parameter derivatives, and
#'   [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1.2)
#' y <- c(0.2, 1, 3)
#'
#' distrib_hess_y(d, y, theta)
#'
#' # Against a numerical second derivative.
#' vapply(y, function(v)
#'   numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
#'                     v)[1, 1],
#'   numeric(1))
#'
#' # Far from zero the fold is invisible and the curvature is the parent's,
#' # which for a gaussian is -1 / sigma^2.
#' c(folded = distrib_hess_y(d, 8, theta), parent = -1 / 1.2^2)
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

#' @title Does a Distribution Declare Atoms
#'
#' @description
#' Answers whether a distribution registers [distrib_atoms()] for itself, which
#' is how [folded()] decides that a parent carries a point mass and must be
#' rejected. The question is asked of the CLASS the method is registered on,
#' not of a value: a distribution with an atom declares one at every parameter
#' setting, and a value taken at one setting could be empty by accident.
#'
#' @details
#' The argument is named `parent` deliberately. The base class of this package
#' is called `distrib`, and an argument of that name would shadow it: the
#' comparison meant for the base class would then be against the object. That
#' defect has been met before in this package and is avoided by naming.
#'
#' @param parent A `distrib` object.
#'
#' @return `TRUE` when [distrib_atoms()] is registered on a class strictly
#'   below `distrib`, `FALSE` when the method comes from the base class or is
#'   absent.
#'
#' @seealso [distrib_atoms()] for the generic, [folded()], which consults this,
#'   and [zero_adjusted()], which produces a parent it rejects.
#'
#' @examples
#' # A plain gaussian declares none.
#' distributions7:::declares_atoms(gaussian1_distrib())
#'
#' # A zero-adjusted continuous parent is a mixed distribution and does.
#' distributions7:::declares_atoms(zero_adjusted(gaussian1_distrib()))
#'
#' # Which is why folded() rejects the second by name.
#' try(folded(zero_adjusted(gaussian1_distrib())))
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

#' @title Fold a Distribution at Zero
#'
#' @description
#' Wraps a continuous distribution into the distribution of the absolute value
#' of its variable. The result has density
#' \deqn{L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge 0,}
#' distribution function \eqn{F(x) - F(-x)}, and exactly the parent's
#' parameters: folding adds none and removes none, as [truncated()] does not.
#' The half-normal is `fixed(folded(gaussian1_distrib()), mu = 0)`, and the
#' folded normal proper is `folded(gaussian1_distrib())`.
#'
#' @details
#' # A fold is not a change of variable
#'
#' The map \eqn{y \mapsto |y|} is TWO TO ONE and has no inverse, so it cannot
#' be a [transformer()]. Instead of carrying a density through a Jacobian it
#' ADDS the two preimages. Every point above zero has two of them and zero has
#' one, which is why the constructor is careful about atoms.
#'
#' # Where the derivatives come from
#'
#' Every one comes from the wrapper machinery unchanged. Writing
#' \eqn{w = f(x)/L(x)} for the weight of the positive preimage, the ratios that
#' machinery consumes are
#' \deqn{\frac{d^B L}{L} = w \frac{d^B f(x)}{f(x)}
#'       + (1-w) \frac{d^B f(-x)}{f(-x)},}
#' each term a complete Bell polynomial in the parent's own log-derivatives,
#' one at \eqn{+x} and one at \eqn{-x}; \eqn{\log L} then follows by the
#' moment-to-cumulant relation. At first order this reads
#' \eqn{w s(x) + (1-w) s(-x)}, the score of a two-component mixture, which is
#' what a fold is. A parent analytic to fourth order gives a folded family
#' analytic to fourth order.
#'
#' # What is rejected
#'
#' A parent that does not reach below zero folds to itself, so the call would
#' be a no-op; returning it unchanged would hide a mistaken call rather than
#' report it, and the same check makes `folded()` of a folded distribution an
#' error. A parent with an atom is rejected too: zero is its own preimage while
#' every other point has two, so an atom at zero would be counted twice by the
#' sum above and one elsewhere would be moved onto its reflection. A discrete
#' parent is rejected outright.
#'
#' # The sign of a symmetric parent's location is not identified
#'
#' When the parent is symmetric about its location, \eqn{f(-x; \mu) = f(x;
#' -\mu)}, so the two terms of \eqn{L} merely swap and the likelihood is
#' EXACTLY even in \eqn{\mu}. A fit converges to \eqn{+\hat\mu} or
#' \eqn{-\hat\mu} according to where it started, at the same maximized value to
#' every digit: measured on 400 draws at \eqn{\mu = 1.5}, two fits started at
#' \eqn{\pm 1} reach \eqn{\pm 1.50494} with the same log-likelihood of
#' \eqn{-511.2516}. This is a property of the model, not of the fitting, and it
#' is not rejected, the folded normal being a standard family; what is
#' estimable is \eqn{|\mu|} together with the remaining parameters. Holding the
#' location at zero removes the question and gives the half-normal. A parent
#' that is not symmetric about its location, such as [skewnormal1_distrib()],
#' has no such invariance and its sign is identified.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{F} its distribution function, \eqn{L}
#' the folded density, \eqn{w} the weight of the positive preimage, \eqn{s} the
#' parent's score and \eqn{B} a multiset of parameter names.
#'
#' @param distrib A `continuous_distrib` object whose support reaches below
#'   zero and which declares no atom. A discrete distribution, a distribution
#'   supported on \eqn{[0, \infty)}, an already folded distribution, and a
#'   zero-adjusted continuous one are each rejected with an error saying which
#'   condition failed.
#'
#' @return An S7 object of class [FoldedDistrib] carrying `parent_distrib`. Its
#'   `params`, `params_bounds` and `link_params` are the parent's unchanged;
#'   `bounds` is `c(0, max(abs(parent bounds)))`; `distrib_name` is `"folded "`
#'   followed by the parent's; and `params_smooth` is the parent's smoothness.
#'
#' @seealso [truncated()], the other wrapper that adds no parameter,
#'   [fixed()], which gives the half-normal, [transformation()] for a map that
#'   IS one to one, and [FoldedDistrib] for the class.
#'
#' @examples
#' d <- folded(gaussian1_distrib())
#' theta <- list(mu = 0.5, sigma = 1)
#' distrib_pdf(d, c(0, 0.5, 2), theta)
#'
#' # The half-normal: a folded gaussian with its location held at zero, whose
#' # density is twice the gaussian's and whose mean is sigma sqrt(2 / pi).
#' hn <- fixed(folded(gaussian1_distrib()), mu = 0)
#' hn@params
#' all.equal(distrib_pdf(hn, c(0.5, 1), list(sigma = 2)),
#'           2 * dnorm(c(0.5, 1), 0, 2))
#' c(mean = mean(hn, list(sigma = 2)), theory = 2 * sqrt(2 / pi))
#'
#' # The sign of a symmetric parent's location is not identified: two fits
#' # started either side reach the same maximum at opposite signs.
#' set.seed(3)
#' z <- distrib_rng(d, 400, list(mu = 1.5, sigma = 1))
#' f1 <- fit_distrib(d, z, start = list(mu = 1, sigma = 1))
#' f2 <- fit_distrib(d, z, start = list(mu = -1, sigma = 1))
#' rbind(from_plus = c(coef(f1), logLik = as.numeric(logLik(f1))),
#'       from_minus = c(coef(f2), logLik = as.numeric(logLik(f2))))
#'
#' # A parent that is not symmetric about its location has no such invariance.
#' sn <- folded(skewnormal1_distrib())
#' set.seed(4)
#' zs <- distrib_rng(sn, 300, list(mu = 1, sigma = 1, alpha = 3))
#' lf <- function(m)
#'   sum(distrib_pdf(sn, zs, list(mu = m, sigma = 1, alpha = 3), log = TRUE))
#' c(at_plus_1 = lf(1), at_minus_1 = lf(-1))
#'
#' # Three rejections, each naming the condition that failed.
#' try(folded(gamma1_distrib()))
#' try(folded(folded(gaussian1_distrib())))
#' try(folded(poisson_distrib()))
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
