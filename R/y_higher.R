#' @include y_derivatives.R generics.R
#' @include gaussian1_distrib.R gaussian2_distrib.R gaussian3_distrib.R
#' @include cauchy_distrib.R logistic_distrib.R laplace_distrib.R
#' @include laplace2_distrib.R enet_distrib.R pseudohuber_distrib.R
#' @include student_t1_distrib.R
#' @include skewnormal1_distrib.R skewnormal2_distrib.R skewt_distrib.R
#' @include gumbel_distrib.R
#' @include gamma1_distrib.R gamma2_distrib.R chisq_distrib.R
#' @include exponential_distrib.R beta1_distrib.R beta2_distrib.R
#' @include weibull1_distrib.R gengamma1_distrib.R
#' @include invgauss1_distrib.R invgauss2_distrib.R lognormal1_distrib.R
#' @include gpd_distrib.R vonmises1_distrib.R vonmises2_distrib.R
#' @include reparametrize.R
NULL

# ===========================================================================
# Third and fourth derivatives of the log-density with respect to the
# response.
#
# For a family in which the response enters only as y - mu there is nothing
# to derive: differentiating in y and in the location are the same operation
# up to a sign,
#
#   d^k l / dy^k = (-1)^k  d^k l / dmu^k,
#
# so every location family inherits these orders from the parameter
# derivatives it already has. What is left is the families on a half line,
# where the two are unrelated, and those take one stencil of the requested
# order applied to the log-density -- never a difference of a difference.
# ===========================================================================

#' Numerical Third and Fourth Response Derivatives
#'
#' @description
#' One central stencil of the requested order applied to
#' \code{distrib_pdf(..., log = TRUE)}.
#'
#' @details
#' The stencil reaches two steps either side, so the step is clamped to half
#' of what \code{\link{fd_steps_y}} allows: a Gamma observation near zero
#' would otherwise be differentiated at a point outside the support. The
#' relative step is \eqn{\varepsilon^{1/(k+2)}}, which balances the
#' \eqn{h^{2}} truncation against the \eqn{\varepsilon/h^{k}} rounding.
#'
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 3 or 4.
#' @param h_rel Relative finite-difference step.
#'
#' @return A numeric vector the length of \code{y}.
#'
#' @seealso \code{\link{numerical_hess_y}}
#'
#' @examples
#' numerical_deriv_y(gaussian1_distrib(), c(-1, 0, 1),
#'                   list(mu = 0, sigma = 1), order = 3)
#'
#' @export
numerical_deriv_y <- function(distrib, y, theta, order,
                              h_rel = .Machine$double.eps^(1 / (order + 2))) {
  h <- fd_steps_y(y, distrib@bounds, h_rel) / 2
  ld <- function(s) distrib_pdf(distrib, y + s * h, theta, log = TRUE)
  if (order == 3L) {
    (-0.5 * ld(-2) + ld(-1) - ld(1) + 0.5 * ld(2)) / h^3
  } else {
    (ld(-2) - 4 * ld(-1) + 6 * ld(0) - 4 * ld(1) + ld(2)) / h^4
  }
}

#' Third and Fourth Derivatives With Respect to the Response
#'
#' @description
#' \eqn{\partial^{3}\ell/\partial y^{3}} and
#' \eqn{\partial^{4}\ell/\partial y^{4}}, completing the sequence begun by
#' \code{\link{distrib_grad_y}} and \code{\link{distrib_hess_y}}.
#'
#' @details
#' A family whose response enters only as \eqn{y - \mu} gets these from the
#' derivatives it already has in the location, since
#' \eqn{\partial^{k}\ell/\partial y^{k} = (-1)^{k}\,
#' \partial^{k}\ell/\partial\mu^{k}}; the others take one stencil of the
#' requested order on the log-density.
#'
#' As with the orders below, a discrete family has no such derivative and the
#' generic rejects rather than returning a difference across the lattice.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param ... Passed to methods.
#'
#' @return A numeric vector the length of \code{y}.
#'
#' @seealso \code{\link{distrib_hess_y}}, \code{\link{numerical_deriv_y}}
#'
#' @examples
#' distrib_deriv3_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#'
#' @export
distrib_deriv3_y <- S7::new_generic(
  "distrib_deriv3_y", "distrib",
  function(distrib, y, theta, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @rdname distrib_deriv3_y
#'
#' @examples
#' distrib_deriv4_y(logistic_distrib(), 0.5, list(mu = 0, sigma = 1))
#'
#' @export
distrib_deriv4_y <- S7::new_generic(
  "distrib_deriv4_y", "distrib",
  function(distrib, y, theta, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  })

#' @title Default Third and Fourth Response Derivatives
#' @name distrib_deriv3_y.continuous_distrib
#' @description One stencil of the order asked for, on the log-density.
#' @param distrib A \code{continuous_distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(distrib_deriv3_y, continuous_distrib) <- function(distrib, y, theta,
                                                             ...) {
  numerical_deriv_y(distrib, y, theta, 3L)
}

#' @rdname distrib_deriv3_y.continuous_distrib
#' @name distrib_deriv4_y.continuous_distrib
S7::method(distrib_deriv4_y, continuous_distrib) <- function(distrib, y, theta,
                                                             ...) {
  numerical_deriv_y(distrib, y, theta, 4L)
}

#' The Response Derivative of a Location Family
#'
#' @description
#' Builds the order-\code{k} response derivative of a family whose response
#' enters only as \eqn{y - \mu}, as \eqn{(-1)^{k}} times the pure derivative
#' in the location.
#'
#' @details
#' The identity is exact and needs no formula of its own, which is the point:
#' the location derivative is already written, often as a compiled kernel, and
#' the response derivative is the same number with a sign. It also fixes the
#' two orders below, where \eqn{\partial\ell/\partial y = -\partial\ell/
#' \partial\mu} and \eqn{\partial^{2}\ell/\partial y^{2} =
#' \partial^{2}\ell/\partial\mu^{2}}, and the tests check the new orders
#' against exactly that.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @keywords internal
loc_deriv_y_k <- function(order) {
  key <- paste(rep("mu", order), collapse = "_")
  sgn <- (-1)^order
  function(distrib, y, theta, ...) {
    d <- if (order == 3L) distrib_deriv3(distrib, y, theta) else
                          distrib_deriv4(distrib, y, theta)
    rep_len(sgn * d[[key]], length(y))
  }
}

# The families whose response enters only as y - mu. A family on a half line
# is not among them: there the location and the response are unrelated
# directions and the fallback stencil is what applies.
for (.cls in list(Gaussian1Distrib, Gaussian2Distrib, Gaussian3Distrib,
                  CauchyDistrib, LogisticDistrib, LaplaceDistrib,
                  Laplace2Distrib, EnetDistrib, PseudoHuberDistrib,
                  StudentT1Distrib, SkewNormal1Distrib,
                  SkewNormal2Distrib, SkewTDistrib, GumbelDistrib)) {
  S7::method(distrib_deriv3_y, .cls) <- loc_deriv_y_k(3L)
  S7::method(distrib_deriv4_y, .cls) <- loc_deriv_y_k(4L)
}
rm(.cls)


# --- the families whose response is not a pure location ---------------------
#
# Every one of them already carries a closed first and second response
# derivative, and the third is the same elementary function differentiated once
# more: the log-density is a sum of terms in log(y), log(1-y), a power of y, a
# logarithm of an affine function of y, or a cosine. The terms are written once
# below and each family is a sum of them, which is shorter than four formulas
# per family and leaves one place for each rule to be wrong in.

#' The k-th Response Derivative of an Elementary Term
#'
#' @description
#' \code{dy_log} differentiates \eqn{c\log y}, \code{dy_log1m} differentiates
#' \eqn{c\log(1-y)}, \code{dy_pow} differentiates \eqn{c\,y^{p}},
#' \code{dy_logaff} differentiates \eqn{c\log(a + by)} and \code{dy_cos}
#' differentiates \eqn{c\cos(y-m)}.
#'
#' @details
#' A term linear in \eqn{y} is \code{dy_pow} at \eqn{p = 1} and vanishes from
#' the second order, so it needs no case of its own.
#'
#' @param c The coefficient.
#' @param y The response.
#' @param k The derivative order.
#' @param p The exponent, for \code{dy_pow}.
#' @param a,b The affine coefficients, for \code{dy_logaff}.
#' @param m The location, for \code{dy_cos}.
#'
#' @return A numeric vector the length of \code{y}.
#'
#' @keywords internal
dy_log <- function(c, y, k) c * (-1)^(k - 1L) * factorial(k - 1L) / y^k

#' @rdname dy_log
#' @keywords internal
dy_log1m <- function(c, y, k) -c * factorial(k - 1L) / (1 - y)^k

#' @rdname dy_log
#' @keywords internal
dy_pow <- function(c, p, y, k) c * prod(p - seq_len(k) + 1) * y^(p - k)

#' @rdname dy_log
#' @keywords internal
dy_logaff <- function(c, a, b, y, k) {
  c * (-1)^(k - 1L) * factorial(k - 1L) * b^k / (a + b * y)^k
}

#' @rdname dy_log
#' @keywords internal
dy_cos <- function(c, m, y, k) c * cos(y - m + k * pi / 2)

#' The k-th Response Derivative of a Function of the Log Response
#'
#' @description
#' Evaluates \eqn{\partial^{k} g(\log y)/\partial y^{k}} from the derivatives
#' of \eqn{g}, by
#' \eqn{y^{-k}\sum_{j} s(k, j)\, g^{(j)}(\log y)} with \eqn{s(k, j)} the signed
#' Stirling numbers of the first kind.
#'
#' @param gd A list whose \eqn{j}-th element is \eqn{g^{(j)}(\log y)}.
#' @param y The response.
#' @param k The derivative order, 1 to 4.
#'
#' @return A numeric vector the length of \code{y}.
#'
#' @keywords internal
dy_of_log <- function(gd, y, k) {
  s <- list(c(1), c(-1, 1), c(2, -3, 1), c(-6, 11, -6, 1))[[k]]
  acc <- 0
  for (j in seq_len(k)) acc <- acc + s[j] * gd[[j]]
  acc / y^k
}

#' Register the Third and Fourth Response Derivatives of a Family
#'
#' @description
#' Turns a function of \code{(distrib, y, theta, k)} into the two methods, so
#' that a family writes its rule once instead of twice.
#'
#' @param cls The S7 class.
#' @param f The rule.
#'
#' @return Invisibly \code{NULL}; called for the registration.
#'
#' @keywords internal
register_dy_k <- function(cls, f) {
  S7::method(distrib_deriv3_y, cls) <- function(distrib, y, theta, ...) {
    rep_len(f(distrib, y, theta, 3L), length(y))
  }
  S7::method(distrib_deriv4_y, cls) <- function(distrib, y, theta, ...) {
    rep_len(f(distrib, y, theta, 4L), length(y))
  }
  invisible(NULL)
}

# gamma and its special cases: (a - 1) log y - rate * y, so only the
# logarithm survives past the first order
register_dy_k(Gamma1Distrib, function(distrib, y, theta, k) {
  dy_log(gamma1_shape_rate(theta)$shape - 1, y, k)
})
register_dy_k(Gamma2Distrib, function(distrib, y, theta, k) {
  dy_log(theta[[1]]^2 / theta[[2]] - 1, y, k)
})
register_dy_k(ChisqDistrib, function(distrib, y, theta, k) {
  dy_log(theta[[1]] / 2 - 1, y, k)
})
# the exponential has unit shape, so the log-density is linear in the response
# and everything above the first order is exactly zero
register_dy_k(ExponentialDistrib, function(distrib, y, theta, k) {
  0 * y
})

# beta: (alpha - 1) log y + (beta - 1) log(1 - y)
register_dy_k(Beta1Distrib, function(distrib, y, theta, k) {
  a <- theta[[1]] * theta[[2]]
  b <- (1 - theta[[1]]) * theta[[2]]
  dy_log(a - 1, y, k) + dy_log1m(b - 1, y, k)
})
register_dy_k(Beta2Distrib, function(distrib, y, theta, k) {
  dy_log(theta[[1]] - 1, y, k) + dy_log1m(theta[[2]] - 1, y, k)
})

# weibull: (sigma - 1) log y - (y/mu)^sigma
register_dy_k(Weibull1Distrib, function(distrib, y, theta, k) {
  sigma <- theta[[2]]
  dy_log(sigma - 1, y, k) + dy_pow(-theta[[1]]^-sigma, sigma, y, k)
})

# generalized gamma: (d - 1) log y - (y/a)^p
register_dy_k(GenGamma1Distrib, function(distrib, y, theta, k) {
  p <- theta[[3]]
  dy_log(theta[[2]] - 1, y, k) + dy_pow(-theta[[1]]^-p, p, y, k)
})

# inverse gaussian: -1.5 log y - y/(2 phi mu^2) - 1/(2 phi y)
register_dy_k(InvGauss1Distrib, function(distrib, y, theta, k) {
  dy_log(-1.5, y, k) + dy_pow(-1 / (2 * theta[[2]]), -1, y, k)
})
register_dy_k(InvGauss2Distrib, function(distrib, y, theta, k) {
  dy_log(-1.5, y, k) + dy_pow(-theta[[2]] / 2, -1, y, k)
})

# lognormal: g(log y) with g(u) = -u - (u - mu)^2 / (2 sigma^2)
register_dy_k(Lognormal1Distrib, function(distrib, y, theta, k) {
  s2 <- theta[[2]]
  u <- log(y)
  gd <- list(-1 - (u - theta[[1]]) / s2, rep(-1 / s2, length(y)),
             rep(0, length(y)), rep(0, length(y)))
  dy_of_log(gd, y, k)
})

# generalized Pareto: -(1 + 1/xi) log(1 + xi y / sigma). The coefficient is
# written as xi^k + xi^(k-1) rather than (1 + 1/xi) xi^k, which is the same
# number and stays finite as xi goes to zero, where the family is exponential
# and every order above the first vanishes.
register_dy_k(GPDDistrib, function(distrib, y, theta, k) {
  sigma <- theta[[1]]
  xi <- theta[[2]]
  cf <- -(xi^k + xi^(k - 1L)) / sigma^k
  cf * (-1)^(k - 1L) * factorial(k - 1L) / (1 + xi * y / sigma)^k
})

# von Mises: kappa cos(y - mu), whose derivatives are the same four functions
# in rotation
register_dy_k(VonMises1Distrib, function(distrib, y, theta, k) {
  dy_cos(theta[[2]], theta[[1]], y, k)
})
register_dy_k(VonMises2Distrib, function(distrib, y, theta, k) {
  dy_cos(vm2_parts(theta)$kappa, theta[[1]], y, k)
})

# a reparametrization acts on the parameters and the derivative is taken in the
# response, so the two do not interact and the parent's answer is the answer
S7::method(distrib_deriv3_y, ReparamContinuousDistrib) <-
  function(distrib, y, theta, ...) {
    distrib_deriv3_y(distrib@parent_distrib, y, reparam_theta(distrib, theta), ...)
  }
S7::method(distrib_deriv4_y, ReparamContinuousDistrib) <-
  function(distrib, y, theta, ...) {
    distrib_deriv4_y(distrib@parent_distrib, y, reparam_theta(distrib, theta), ...)
  }
