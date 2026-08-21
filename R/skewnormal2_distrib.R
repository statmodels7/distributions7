#' @include distrib.R generics.R skewnormal1_distrib.R reparametrize.R moments.R
NULL

# The skew normal in Azzalini's CENTERED parametrization: the mean, the standard
# deviation and the skewness itself, rather than the location, the scale and
# the shape.
#
# This is a family of its own and not a reparametrize() of skewnormal1, for two
# reasons that both come from the map passing through a cube root. It carries a
# sign, which a jet cannot take; and its derivative is unbounded as the
# skewness goes to zero, so what makes the parametrization worth having is a
# CANCELLATION between terms that individually diverge. Measured at gamma1 =
# 1e-4, d alpha / d gamma1 is 258 while the variance of the score in gamma1 is
# 0.158, the same value it has at gamma1 = 0.05: the divergence cancels.

#' The Constant Behind the Centered Parametrization
#'
#' @description
#' \eqn{b = \sqrt{2/\pi}}, which is \eqn{\mathbb{E}[|Z|]} for a standard normal
#' and appears in every quantity of the skew normal's first moment.
#'
#' @return A single number.
#'
#' @seealso \code{\link{skewnormal2_distrib}}
#'
#' @keywords internal
sn_b <- function() sqrt(2 / pi)

#' The Largest Skewness a Skew Normal Can Reach
#'
#' @description
#' The supremum of \eqn{|\gamma_1|} over the family, attained only in the limit
#' \eqn{\alpha \to \pm\infty}: \eqn{(4-\pi)/2 \cdot (b/\sqrt{1-b^2})^3}, about
#' 0.9952717.
#'
#' @details
#' A skewness beyond it belongs to no skew normal, which is why the constructor
#' bounds the parameter there rather than letting the map return a
#' \code{NaN} several frames down. It is also the reason the skew \eqn{t}
#' exists.
#'
#' @return A single number.
#'
#' @seealso \code{\link{skewnormal2_distrib}}
#'
#' @keywords internal
sn_max_skew <- function() {
  b <- sn_b()
  (4 - pi) / 2 * (b / sqrt(1 - b^2))^3
}

#' From the Centered Parameters to the Direct Ones
#'
#' @description
#' The map \eqn{(\mu, \sigma, \gamma_1) \mapsto (\xi, \omega, \alpha)} that
#' takes the mean, the standard deviation and the skewness to the location, the
#' scale and the shape.
#'
#' @details
#' With \eqn{b = \sqrt{2/\pi}},
#' \deqn{c = \mathrm{sign}(\gamma_1)
#'           \left(\dfrac{2|\gamma_1|}{4-\pi}\right)^{1/3}, \qquad
#'       \mu_z = \dfrac{c}{\sqrt{1+c^2}}, \qquad
#'       \delta = \dfrac{\mu_z}{b}, \qquad
#'       \alpha = \dfrac{\delta}{\sqrt{1-\delta^2}},}
#' and then \eqn{\omega = \sigma/\sqrt{1-\mu_z^2}} and
#' \eqn{\xi = \mu - \omega\mu_z}.
#'
#' The function is written once and used twice: on plain numbers for the
#' density, and on jets for the derivatives. The sign is taken from the plain
#' value before any jet is seeded, which is what a jet cannot do for itself and
#' what makes this a family rather than a \code{\link{reparametrize}}.
#'
#' @param mu,sigma,gamma1 The centered parameters, numbers or jets.
#' @param s The sign of \eqn{\gamma_1}, taken from its plain value.
#'
#' @return A named list with \code{mu}, \code{sigma} and \code{alpha}, the
#'   parent's parameters.
#'
#' @seealso \code{\link{skewnormal2_distrib}}
#'
#' @keywords internal
sn_cp_to_dp <- function(mu, sigma, gamma1, s) {
  b <- sn_b()
  # |gamma1| as s * gamma1: away from zero the sign is locally constant, so
  # this is exact and carries the right derivatives when the argument is a jet.
  cc <- s * (2 * (s * gamma1) / (4 - pi))^(1 / 3)
  muz <- cc / sqrt(1 + cc^2)
  del <- muz / b
  om <- sigma / sqrt(1 - muz^2)
  list(mu = mu - om * muz,
       sigma = om,
       alpha = del / sqrt(1 - del^2))
}

#' @title S7 Class for the Skew Normal in Its Centered Parametrization
#' @name SkewNormal2Distrib
#'
#' @description A subclass of \code{continuous_distrib} for the skew normal
#'   written in its mean, standard deviation and skewness.
#' @inheritParams distrib
#' @return An object of class \code{SkewNormal2Distrib}.
#' @seealso \code{\link{skewnormal2_distrib}}, \code{\link{skewnormal1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.SkewNormal2Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.SkewNormal2Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.SkewNormal2Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.SkewNormal2Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.SkewNormal2Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.SkewNormal2Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.SkewNormal2Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.SkewNormal2Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.SkewNormal2Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
SkewNormal2Distrib <- S7::new_class("SkewNormal2Distrib",
                                    parent = continuous_distrib)

#' The Direct Parameters a Centered Triple Implies
#'
#' @description
#' Runs \code{\link{sn_cp_to_dp}} on plain numbers, which every probability
#' function needs before delegating to \code{\link{skewnormal1_distrib}}.
#'
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#'
#' @return A named list of the direct parameters.
#'
#' @seealso \code{\link{skewnormal2_distrib}}
#'
#' @keywords internal
sn2_theta <- function(theta) {
  g <- theta[[3]]
  s <- ifelse(g >= 0, 1, -1)
  sn_cp_to_dp(theta[[1]], theta[[2]], g, s)
}

#' Derivatives of the Skew Normal in Its Centered Parametrization
#'
#' @description
#' The parent's derivatives carried into the centered coordinates by the
#' partition sum of \code{\link{chain_derivatives}}.
#'
#' @param distrib A \code{\link{SkewNormal2Distrib}} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param order The derivative order.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#'
#' @details
#' The map to the direct parametrization runs through
#' \eqn{r = \sqrt[3]{2\gamma_1/(4-\pi)}}, whose derivative grows like
#' \eqn{\gamma_1^{-2/3}}, so at zero skewness the map is not differentiable
#' and the chain rule is asked for a quantity that does not exist. The first
#' derivatives of the log-density survive the limit -- the map's factor
#' cancels and they approach a finite value from both sides -- but the second
#' ones diverge at that rate, which is a property of the CENTERED
#' parametrization and not of the family. It is rejected here, where the map
#' is used and the reason can be named, rather than left to reach a
#' comparison against \code{NA} several frames further on.
#'
#' @return A named list of component vectors.
#'
#' @seealso \code{\link{skewnormal2_distrib}}
#'
#' @keywords internal
sn2_chain <- function(distrib, y, theta, order, expected = FALSE) {
  theta <- align_theta(distrib, theta)
  if (any(theta[[3L]] == 0)) {
    stop(paste0(
      "The centered parametrization has no derivatives at zero skewness:\n",
      "  the map to the direct parameters runs through the cube root of\n",
      "  gamma1, whose derivative is unbounded there. The first derivatives\n",
      "  of the log-density have a finite limit and the second ones grow\n",
      "  like gamma1^(-2/3), so the point is excluded rather than\n",
      "  approximated. skewnormal1_distrib() carries the same family in the\n",
      "  direct parametrization, whose derivatives at alpha = 0 are ordinary\n",
      "  numbers."), call. = FALSE)
  }
  chain_derivatives(
    parent = skewnormal1_distrib(),
    y = y,
    th_par = sn2_theta(theta),
    maps = md_skewnormal2(theta[1:3]),
    new_params = distrib@params,
    order = order,
    expected = expected
  )
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Skew Normal Density in the Centered Parametrization
#' @name distrib_pdf.SkewNormal2Distrib
#' @description The skew normal density at the implied direct parameters.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_pdf, SkewNormal2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(skewnormal1_distrib(), y, sn2_theta(theta), log = log)
}

#' @title Skew Normal Distribution Function in the Centered Parametrization
#' @name distrib_cdf.SkewNormal2Distrib
#' @description The skew normal distribution function, through Owen's T.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns log-probabilities.
#' @param ... Passed on.
#' @return A numeric vector.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_cdf, SkewNormal2Distrib) <- function(distrib, q, theta,
                                                         lower.tail = TRUE,
                                                         log.p = FALSE, ...) {
  distrib_cdf(skewnormal1_distrib(), q, sn2_theta(theta),
              lower.tail = lower.tail, log.p = log.p, ...)
}

#' @title Skew Normal Quantile Function in the Centered Parametrization
#' @name distrib_quantile.SkewNormal2Distrib
#' @description The parent's quantile function at the implied parameters.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is a log-probability.
#' @param ... Passed on.
#' @return A numeric vector.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_quantile, SkewNormal2Distrib) <- function(distrib, p, theta,
                                                              lower.tail = TRUE,
                                                              log.p = FALSE, ...) {
  distrib_quantile(skewnormal1_distrib(), p, sn2_theta(theta),
                   lower.tail = lower.tail, log.p = log.p, ...)
}

#' @title Skew Normal Random Generation in the Centered Parametrization
#' @name distrib_rng.SkewNormal2Distrib
#' @description The parent's generator at the implied parameters.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param n The number of draws.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @return A numeric vector.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_rng, SkewNormal2Distrib) <- function(distrib, n, theta) {
  distrib_rng(skewnormal1_distrib(), n, sn2_theta(theta))
}

#' @title Skew Normal Gradient in the Centered Parametrization
#' @name distrib_gradient.SkewNormal2Distrib
#' @description
#' The parent's score carried by the Jacobian of the centered-to-direct map. The
#' components in \eqn{\gamma_1} stay of order one however small \eqn{\gamma_1}
#' is, although the Jacobian itself grows without bound: the divergent parts
#' cancel, which is what the centered parametrization is for.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_gradient, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"), ...) {
  sn2_chain(distrib, y, theta, 1L)
}

#' @title Skew Normal Observed Hessian in the Centered Parametrization
#' @name distrib_hessian.SkewNormal2Distrib
#' @description The second-order chain rule through the same map.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_hessian, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  sn2_chain(distrib, y, theta, 2L)[hess_names(distrib@params)]
}

#' @title Skew Normal Expected Hessian in the Centered Parametrization
#' @name distrib_expected_hessian.SkewNormal2Distrib
#' @description
#' The parent's expected information carried by the same congruence. It is
#' \strong{non-singular at zero skewness}, which the direct parametrization's
#' is not: there the score for \eqn{\alpha} is proportional to the score for
#' the location and the information loses a rank.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Passed to the parent.
#' @param nsim Passed to the parent.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_expected_hessian, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                                      scale = c("parameter", "link"),
                                                                      approx = c("bartlett", "integrate", "mc", "opg"),
                                                                      nsim = 10000, ...) {
  sn2_chain(distrib, y, theta, 2L, expected = TRUE)[hess_names(distrib@params)]
}

#' @title The Centered Skew Normal Does Not Write Its Expected Information Out
#' @name expected_hessian_exact.SkewNormal2Distrib
#' @description
#' The method above is the CHAIN onto \code{\link{skewnormal1_distrib}}, whose
#' expected information is the base class's quadrature, so the registration
#' says where the arithmetic is assembled and not that it is closed form.
#' @details
#' Measured at 100 observations it costs 5220 ms against the parent's 2230 --
#' more than what it chains onto, the chain being paid on top -- where the
#' families that do write it out answer in a median of 0.183 ms. Reported as
#' exact, it made \code{\link{fit_distrib}} reject a legitimate
#' \code{fisher_scoring(approx = )} here with a message that was untrue.
#' @param x A \code{SkewNormal2Distrib} object.
#' @param ... Unused.
#' @return \code{FALSE}.
#' @seealso \code{\link{expected_hessian_exact}}
#' @keywords internal
S7::method(expected_hessian_exact, SkewNormal2Distrib) <- function(x, ...) {
  expected_hessian_exact(skewnormal1_distrib())
}

#' @title Skew Normal Third-Order Derivatives in the Centered Parametrization
#' @name distrib_deriv3.SkewNormal2Distrib
#' @description The partition sum at order three.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Passed to the parent.
#' @param nsim Passed to the parent.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_deriv3, SkewNormal2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  sn2_chain(distrib, y, theta, 3L, expected = expected)
}

#' @title Skew Normal Fourth-Order Derivatives in the Centered Parametrization
#' @name distrib_deriv4.SkewNormal2Distrib
#' @description The partition sum at order four.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Passed to the parent.
#' @param nsim Passed to the parent.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_deriv4, SkewNormal2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  sn2_chain(distrib, y, theta, 4L, expected = expected)
}

#' @title Skew Normal Response Derivatives in the Centered Parametrization
#' @name distrib_grad_y.SkewNormal2Distrib
#' @description The parent's, unchanged: the coordinates change, the response
#'   does not.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param ... Passed on.
#' @return A numeric vector.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_grad_y, SkewNormal2Distrib) <- function(distrib, y, theta, ...) {
  distrib_grad_y(skewnormal1_distrib(), y, sn2_theta(theta), ...)
}

#' @title Skew Normal Second Response Derivative in the Centered Parametrization
#' @name distrib_hess_y.SkewNormal2Distrib
#' @description The parent's, unchanged.
#' @param distrib A \code{SkewNormal2Distrib} object.
#' @param y The response.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param ... Passed on.
#' @return A numeric vector.
#' @seealso \code{\link{skewnormal2_distrib}}
S7::method(distrib_hess_y, SkewNormal2Distrib) <- function(distrib, y, theta, ...) {
  distrib_hess_y(skewnormal1_distrib(), y, sn2_theta(theta), ...)
}

#' @title Mean of the Skew Normal in the Centered Parametrization
#' @name mean.SkewNormal2Distrib
#' @description \eqn{\mu}, the parameter itself: that is what centered means.
#' @param x A \code{SkewNormal2Distrib} object.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(mean, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 3L, 0) + theta[[1]]
}

#' @title Variance of the Skew Normal in the Centered Parametrization
#' @name variance.SkewNormal2Distrib
#' @description \eqn{\sigma^2}, the square of the parameter.
#' @param x A \code{SkewNormal2Distrib} object.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(variance, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 3L, 0)
}

#' @title Skewness of the Skew Normal in the Centered Parametrization
#' @name skewness.SkewNormal2Distrib
#' @description \eqn{\gamma_1}, the parameter itself.
#' @param x A \code{SkewNormal2Distrib} object.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(skewness, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 3L, 0) + theta[[3]]
}

#' @title Kurtosis of the Skew Normal in the Centered Parametrization
#' @name kurtosis.SkewNormal2Distrib
#' @description
#' The parent's, at the implied direct parameters: the centered parametrization
#' fixes the first three moments and leaves the fourth to follow.
#' @param x A \code{SkewNormal2Distrib} object.
#' @param theta A list with \code{mu}, \code{sigma} and \code{gamma1}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(kurtosis, SkewNormal2Distrib) <- function(x, theta, ...) {
  kurtosis(skewnormal1_distrib(), sn2_theta(align_theta(x, theta)))
}


#' Skew Normal Distribution in Its Centered Parametrization
#'
#' @description
#' Creates a skew normal distribution object parametrized by its mean, its
#' standard deviation and its skewness.
#'
#' @details
#' The direct parametrization of \code{\link{skewnormal1_distrib}} carries a
#' location, a scale and a shape, none of which is a moment. Here all three
#' parameters are moments, and the family is Azzalini's centered
#' parametrization.
#'
#' \strong{Why this is not a reparametrize().} The map passes through
#' \eqn{c = \mathrm{sign}(\gamma_1)(2|\gamma_1|/(4-\pi))^{1/3}}, and two things
#' follow. It carries a sign, which a jet cannot take of itself, so the sign is
#' read off the plain value before any jet is seeded. And
#' \eqn{\partial\alpha/\partial\gamma_1} grows without bound as
#' \eqn{\gamma_1 \to 0}: measured, 3.9 at \eqn{\gamma_1 = 0.5} and 258 at
#' \eqn{10^{-4}}. The value of the parametrization is that the score
#' does \strong{not} follow it, the divergent contributions canceling, so the
#' variance of the score in \eqn{\gamma_1} is 0.158 at \eqn{\gamma_1 = 0.05}
#' and 0.158 again at \eqn{0.01}.
#'
#' \strong{What that costs in arithmetic.} The cancellation is between terms of
#' size proportional to the Jacobian, so the significant digits lost grow like
#' the logarithm of it: negligible over the range a fit visits, and severe only
#' within a few multiples of \eqn{10^{-8}} of exact symmetry. The
#' \strong{expected information}, unlike the direct parametrization's, is
#' non-singular at zero skewness, which is the property the parametrization
#' exists for.
#'
#' \strong{The bound on the skewness.} A skew normal cannot reach
#' \eqn{|\gamma_1| > 0.9952717}, whatever \eqn{\alpha} is, so the parameter is
#' bounded there. That ceiling is the reason the skew \eqn{t} exists.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{2}{\omega}\,\phi\!\left(\frac{y-\xi}{\omega}\right)\Phi\!\left(\alpha\,\frac{y-\xi}{\omega}\right), \qquad (\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}}
#'
#' @param link_mu Link function for the mean. Defaults to the identity.
#' @param link_sigma Link function for the standard deviation. Defaults to the
#'   log.
#' @param link_gamma1 Link function for the skewness. Defaults to a link
#'   bounded to the reachable interval.
#'
#' @return An S7 object of class \code{\link{SkewNormal2Distrib}}.
#'
#' @seealso \code{\link{skewnormal1_distrib}}, \code{\link{skewt_distrib}}
#'
#' @examples
#' d <- skewnormal2_distrib()
#' theta <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # all three parameters are moments, which is what centered means
#' c(mean = mean(d, theta), sd = sqrt(variance(d, theta)),
#'   skewness = skewness(d, theta))
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2014). \emph{The Skew-Normal and Related
#' Families}. Cambridge University Press. The centred parametrization is
#' section 3.1.4.
#'
#' @export
skewnormal2_distrib <- function(link_mu = identity_link(),
                                link_sigma = log_link(),
                                link_gamma1 = bounded_link(
                                  lwr = -sn_max_skew(), upr = sn_max_skew()
                                )) {
  g <- sn_max_skew()
  SkewNormal2Distrib(
    distrib_name = "skew normal2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma", "gamma1"),
    params_interpretation = c(mu = "mean", sigma = "standard deviation",
                              gamma1 = "skewness"),
    n_params = 3,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf),
                         gamma1 = c(-g, g)),
    link_params = list(mu = link_mu, sigma = link_sigma,
                       gamma1 = link_gamma1)
  )
}
