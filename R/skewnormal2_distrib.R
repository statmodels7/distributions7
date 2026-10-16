#' @include distrib.R generics.R skewnormal1_distrib.R reparametrize.R moments.R dexpected_hessian.R
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

#' @title The Constant Behind the Centered Parametrization
#'
#' @description
#' Returns \eqn{b = \sqrt{2/\pi} \approx 0.7978846}, which is
#' \eqn{E[|Z|]} for a standard Gaussian \eqn{Z}. It is the one constant that
#' the skew normal's first moment introduces, and it appears in every quantity of
#' the centered parametrization: the mean is \eqn{\xi + \omega b\delta}, the
#' variance \eqn{\omega^2(1-b^2\delta^2)}, and the largest reachable skewness
#' is built from \eqn{b/\sqrt{1-b^2}}.
#'
#' @return A single number.
#'
#' @seealso [sn_max_skew()], which is written in terms of it, and
#'   [skewnormal2_distrib()] for the parametrization.
#'
#' @examples
#' distributions7:::sn_b()
#'
#' # It is the mean of the absolute value of a standard Gaussian.
#' set.seed(1)
#' c(constant = distributions7:::sn_b(), sample = mean(abs(rnorm(1e6))))
#'
#' @keywords internal
sn_b <- function() sqrt(2 / pi)

#' @title The Largest Skewness That a Skew Normal Can Reach
#'
#' @description
#' Returns \eqn{\sup|\gamma_1|} over the skew normal family,
#' \deqn{\dfrac{4-\pi}{2}\left(\dfrac{b}{\sqrt{1-b^2}}\right)^3
#'       \approx 0.9952717,}
#' with \eqn{b = \sqrt{2/\pi}}. The bound is attained only in the limit
#' \eqn{\alpha \to \pm\infty}, where the family degenerates to a half-normal.
#'
#' @details
#' A skewness beyond it belongs to no skew normal, so
#' [skewnormal2_distrib()] bounds `gamma1` there and gives it a
#' [linkfunctions7::bounded_link()] over the open interval. Without the bound
#' the centered-to-direct map would return `NaN` without any indication of the
#' cause.
#'
#' The ceiling is approached slowly: numerically, the skewness is 0.9556 at
#' \eqn{\alpha = 10} and 0.99527 at \eqn{\alpha = 10^4}. Data whose absolute
#' skewness exceeds this value need [skewt_distrib()].
#'
#' @return A single number.
#'
#' @seealso [sn_b()] for the constant from which it is built,
#'   [skewnormal2_distrib()] for the parameter that it bounds, and
#'   [skewt_distrib()] for the family that goes further.
#'
#' @examples
#' distributions7:::sn_max_skew()
#'
#' # The direct parametrization approaches it and does not pass it.
#' d1 <- skewnormal1_distrib()
#' vapply(c(10, 1e4, 1e8),
#'        function(a) skewness(d1, list(mu = 0, sigma = 1, alpha = a)), 0)
#'
#' # It is where the centered parametrization's link is bounded.
#' skewnormal2_distrib()@params_bounds$gamma1
#'
#' @keywords internal
sn_max_skew <- function() {
  b <- sn_b()
  (4 - pi) / 2 * (b / sqrt(1 - b^2))^3
}

#' @title One Minus the Squared Skewness Parameter
#'
#' @description
#' Returns \eqn{1 - q}, with
#' \eqn{q = (\lvert\gamma_1\rvert/\gamma_{\max})^{2/3}} and \eqn{\gamma_{\max}}
#' the ceiling of [sn_max_skew()]. With \eqn{c} and \eqn{\delta} as in
#' [sn_cp_to_dp()], \eqn{1 - q = (1 + c^2)(1 - \delta^2)}, and the shape of the
#' direct parametrization is \eqn{\alpha = c/(b\sqrt{1-q})}, so this is the
#' quantity by which it divides. It is evaluated as
#' \eqn{-\mathrm{expm1}\{(2/3)\log(\lvert\gamma_1\rvert/\gamma_{\max})\}}, which
#' involves no cancellation.
#'
#' @details
#' The identity is exact rather than an expansion near the ceiling. With
#' \eqn{c = (2\lvert\gamma_1\rvert/(4-\pi))^{1/3}} the map forms
#' \eqn{\mu_z = c/\sqrt{1+c^2}} and \eqn{\delta = \mu_z/b}, and \eqn{c^2}
#' equals \eqn{(\lvert\gamma_1\rvert/\gamma_{\max})^{2/3}\,b^2/(1-b^2)}, so the
#' ratio to the ceiling carries the whole of the degeneracy.
#'
#' The two forms that this function replaces differ from it only at the top of
#' the range, and both subtract two nearly equal numbers there.
#' `1 - (mu_z/b)^2` reaches exactly zero at the largest skewness that
#' [linkfunctions7::bounded_link()] can produce, and the difference
#' \eqn{b^2 + (b^2-1)c^2} reaches \eqn{-1.11\times 10^{-16}}. At that skewness, one unit in the last place
#' inside the bound, this function returns \eqn{1.48\times 10^{-16}}, and the
#' shape that follows from it is \eqn{1.36\times 10^{8}} and finite.
#'
#' @param gamma1 The skewness, a numeric vector. The function does not
#'   validate it: the result is positive strictly inside
#'   \eqn{(-\gamma_{\max}, \gamma_{\max})}, zero at either bound and negative
#'   outside them.
#' @param s The sign of `gamma1`, \eqn{\pm 1}, taken by the caller from its
#'   plain value.
#'
#' @return A numeric vector of the length of `gamma1`.
#'
#' @seealso [sn_cp_to_dp()] and [md_skewnormal2()], the two places that divide
#'   by it, and [sn_max_skew()] for the ceiling against which it is written.
#'
#' @examples
#' g <- distributions7:::sn_max_skew() * c(0.5, 0.9, 0.999, 1 - 1e-15)
#' distributions7:::sn_one_minus_delta2(g, 1)
#'
#' # The expression it replaces has already lost the last of those.
#' b <- distributions7:::sn_b()
#' cc <- (2 * g / (4 - pi))^(1 / 3)
#' 1 - (cc / sqrt(1 + cc^2) / b)^2
#'
#' @keywords internal
sn_one_minus_delta2 <- function(gamma1, s) {
  -expm1((2 / 3) * log(s * gamma1 / sn_max_skew()))
}

#' @title From the Centered Parameters to the Direct Ones
#'
#' @description
#' Maps \eqn{(\mu, \sigma, \gamma_1)}, the mean, the standard deviation and the
#' skewness, to \eqn{(\xi, \omega, \alpha)}, the location, the scale and the
#' shape that [skewnormal1_distrib()] takes. Every probability function of the
#' centered family calls it and then delegates to the direct one.
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
#' The caller supplies the sign, so the body reads
#' \eqn{s\,(2s\gamma_1/(4-\pi))^{1/3}} with no `abs()` in it. Away from zero
#' the sign is locally constant, so the expression equals the cube-root term
#' of the map. The derivatives of the map are computed in closed form by
#' [md_skewnormal2()].
#'
#' @param mu,sigma,gamma1 The centered parameters: the mean, the standard
#'   deviation and the skewness, each a numeric vector. `sigma` must be
#'   positive and `gamma1` must lie strictly inside
#'   \eqn{(-0.9952717, 0.9952717)}; the function does not validate them.
#' @param s The sign of `gamma1`, \eqn{\pm 1}, taken by the caller from its
#'   plain value.
#'
#' @return A named list with `mu`, `sigma` and `alpha`, the direct parameters,
#'   each of the length of the recycled inputs. The names are the parent's, so
#'   the result can be passed to [skewnormal1_distrib()]'s methods as they
#'   stand; `mu` there is the location \eqn{\xi} and `sigma` the scale
#'   \eqn{\omega}.
#'
#' @seealso [sn2_theta()], which supplies the sign and calls this;
#'   [md_skewnormal2()] for the derivative tables of the same map; and
#'   [skewnormal2_distrib()] for the family.
#'
#' @examples
#' # The mean and the standard deviation are not the location and the scale.
#' distributions7:::sn_cp_to_dp(0, 1, 0.5, 1)
#'
#' # At zero skewness the map is the identity on the first two, and the shape
#' # is zero: the Gaussian sits at the same point in both parametrizations.
#' distributions7:::sn_cp_to_dp(3, 2, 0, 1)
#'
#' # Round trip: the direct parameters reproduce the centered moments.
#' dp <- distributions7:::sn_cp_to_dp(0, 1, 0.5, 1)
#' d1 <- skewnormal1_distrib()
#' c(mean = mean(d1, dp), sd = sqrt(variance(d1, dp)), skew = skewness(d1, dp))
#'
#' @keywords internal
sn_cp_to_dp <- function(mu, sigma, gamma1, s) {
  b <- sn_b()
  # |gamma1| as s * gamma1: away from zero the sign is locally constant, so
  # this is exact and carries the right derivatives when the argument is a jet.
  cc <- s * (2 * (s * gamma1) / (4 - pi))^(1 / 3)
  # The three expressions are md_skewnormal2()'s, so the value and its
  # derivatives are read off one algebraic form. The scale is a product
  # rather than a quotient because 1 + cc^2 and 1 - muz^2 are reciprocal,
  # and the shape divides by sn_one_minus_delta2(), which is the only place
  # the parametrization degenerates.
  list(mu = mu - sigma * cc,
       sigma = sigma * sqrt(1 + cc^2),
       alpha = cc / (b * sqrt(sn_one_minus_delta2(gamma1, s))))
}

#' @title Skew Normal Distribution Class, Centered Parametrization
#' @name SkewNormal2Distrib
#'
#' @description
#' This is the S7 class of the skew normal written in its first three moments:
#' the mean \eqn{\mu}, the standard deviation \eqn{\sigma} and the skewness
#' \eqn{\gamma_1}. It is the same law as [SkewNormal1Distrib], reached through
#' the map of [sn_cp_to_dp()]; the two classes differ only in which three
#' numbers name a member of the family.
#'
#' The expected information of this parametrization is non-singular at
#' \eqn{\gamma_1 = 0}, where the expected information of the direct
#' parametrization loses a rank; at
#' \eqn{\sigma = 1} its eigenvalues tend to 2, 1 and \eqn{1/6} as the skewness
#' goes to zero.
#'
#' Build one with [skewnormal2_distrib()], which supplies the three link
#' functions and bounds the skewness at \eqn{\pm 0.9952717}. This page
#' documents the raw S7 constructor, which takes the parent's properties; the
#' validator of `distrib` checks them, and the constructor adds no check of
#' its own.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `SkewNormal2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [skewnormal2_distrib()] the properties hold `"skew normal2"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma", "gamma1")`, the
#'   interpretations `c(mu = "mean", sigma = "standard deviation", gamma1 =
#'   "skewness")`, `3`, and the domains \eqn{(-\infty,\infty)},
#'   \eqn{(0,\infty)} and \eqn{(-0.9952717, 0.9952717)}.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.SkewNormal2Distrib],
#'   [`distrib_cdf()`][distrib_cdf.SkewNormal2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.SkewNormal2Distrib] and
#'   [`distrib_rng()`][distrib_rng.SkewNormal2Distrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.SkewNormal2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.SkewNormal2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.SkewNormal2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.SkewNormal2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.SkewNormal2Distrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.SkewNormal2Distrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.SkewNormal2Distrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.SkewNormal2Distrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.SkewNormal2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewNormal2Distrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.SkewNormal2Distrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.SkewNormal2Distrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.SkewNormal2Distrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.SkewNormal2Distrib],
#'   [`distrib_cross3_y()`][distrib_cross3_y.SkewNormal2Distrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.SkewNormal2Distrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.SkewNormal2Distrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.SkewNormal2Distrib],
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.SkewNormal2Distrib],
#'   [`distrib_deriv3_cdf()`][distrib_deriv3_cdf.SkewNormal2Distrib] and
#'   [`distrib_deriv4_cdf()`][distrib_deriv4_cdf.SkewNormal2Distrib];
#' - the moments: [`mean()`][mean.SkewNormal2Distrib],
#'   [`variance()`][variance.SkewNormal2Distrib],
#'   [`skewness()`][skewness.SkewNormal2Distrib] and
#'   [`kurtosis()`][kurtosis.SkewNormal2Distrib];
#' - the fitting helpers: [expected_hessian_costly()].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @section The point at zero skewness:
#' The map runs through a cube root, so \eqn{\partial\alpha/\partial\gamma_1}
#' is unbounded at \eqn{\gamma_1 = 0}. The score has a finite limit there and
#' is returned, and so is the expected information. A derivative of the
#' log-density with \eqn{k \ge 2} differentiations in \eqn{\gamma_1} diverges
#' like \eqn{\gamma_1^{4/3-k}} (\eqn{\gamma_1^{-2/3}} for \eqn{k = 2}). The
#' methods of order two and more, and the derivatives of the expected
#' information, signal an error at \eqn{\gamma_1 = 0} exactly. The density and the distribution function equal
#' the Gaussian's there.
#'
#' @seealso [skewnormal2_distrib()] to build one;
#'   [skewnormal1_distrib()] for the direct parametrization;
#'   [sn_cp_to_dp()] for the map between them.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # The skewness is bounded, which the direct parametrization's shape is not.
#' d@params_bounds$gamma1
#'
#' # All three parameters are moments, which is why it is called "centered".
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' c(mean = mean(d, th), sd = sqrt(variance(d, th)), skewness = skewness(d, th))
SkewNormal2Distrib <- S7::new_class("SkewNormal2Distrib",
                                    parent = continuous_distrib)

#' @title The Direct Parameters That a Centered Triple Implies
#'
#' @description
#' Takes the sign of the skewness off its plain value and runs
#' [sn_cp_to_dp()], returning the location, scale and shape that
#' [skewnormal1_distrib()] takes. Every probability function of the centered
#' family calls this and then delegates.
#'
#' @param theta A list with `mu`, `sigma` and `gamma1`, in that order, each a
#'   numeric vector. It is read positionally, so it must already be aligned;
#'   the methods that call it have been through [align_theta()].
#'
#' @return A named list with `mu`, `sigma` and `alpha`: the location
#'   \eqn{\xi}, the scale \eqn{\omega} and the shape \eqn{\alpha}. The names
#'   are the parent's, so the result passes straight into
#'   [skewnormal1_distrib()]'s methods.
#'
#' @seealso [sn_cp_to_dp()] for the map itself.
#'
#' @examples
#' distributions7:::sn2_theta(list(mu = 0, sigma = 1, gamma1 = 0.5))
#'
#' # A negative skewness reflects the shape and moves the location the other way.
#' distributions7:::sn2_theta(list(mu = 0, sigma = 1, gamma1 = -0.5))
#'
#' @keywords internal
sn2_theta <- function(theta) {
  g <- theta[[3]]
  s <- ifelse(g >= 0, 1, -1)
  dp <- sn_cp_to_dp(theta[[1]], theta[[2]], g, s)
  sn2_reject_unmappable(dp, theta)
  dp
}

#' @title Reject Centered Parameters That the Map Cannot Carry
#'
#' @description
#' Checks that the direct parameters produced by [sn_cp_to_dp()] are usable,
#' and signals an error in the centered family's own terms when they are not.
#'
#' @details
#' Every probability function of [skewnormal2_distrib()] evaluates the parent
#' at the mapped parameters, and the parent validates what it is handed
#' against its own domains. Without this check a caller who wrote `gamma1`
#' reads an error about `alpha`, a parameter that the model does not have, and
#' about `"skew normal1"`, a family that the call does not name. The message
#' reports the centered parameter that is responsible and the value that it
#' holds.
#'
#' The check covers the failures of the delegation itself.
#' [skewnormal2_distrib()] bounds `gamma1`, and every generic checks the
#' parameters against their domains before it calls the method, so a
#' skewness outside its domain is reported before the map runs. The remaining
#' case is a `sigma` and a `gamma1`, each inside its own domain, whose implied
#' scale or location is not finite.
#'
#' @param dp The direct parameters, as returned by [sn_cp_to_dp()].
#' @param theta The centered parameters that the caller supplied, ordered as
#'   `c("mu", "sigma", "gamma1")`.
#'
#' @return Invisibly `NULL`; raises an error naming `"skew normal2"` if any
#'   mapped parameter is not finite or the mapped scale is not positive.
#'
#' @seealso [sn2_theta()], which calls it, and [sn_cp_to_dp()] for the map.
#'
#' @examples
#' # A skewness at the ceiling is reported in the centered family's own terms.
#' th <- list(mu = 0, sigma = 1, gamma1 = distributions7:::sn_max_skew())
#' try(distributions7:::sn2_theta(th))
#'
#' # A scale and a skewness that each lie inside their own domains can still
#' # imply a scale that is not finite.
#' th2 <- list(mu = 0, sigma = 1.1e308, gamma1 = 0.99)
#' try(distributions7:::sn2_theta(th2))
#'
#' @keywords internal
sn2_reject_unmappable <- function(dp, theta) {
  if (all(is.finite(dp$mu)) && all(is.finite(dp$sigma)) &&
      all(dp$sigma > 0) && all(is.finite(dp$alpha))) {
    return(invisible(NULL))
  }

  n <- max(lengths(dp))
  mu <- rep_len(dp$mu, n)
  sg <- rep_len(dp$sigma, n)
  al <- rep_len(dp$alpha, n)
  i <- which(!(is.finite(mu) & is.finite(sg) & sg > 0 & is.finite(al)))[1L]
  at <- function(k) format(rep_len(theta[[k]], n)[i])

  cause <- if (!is.finite(al[i])) {
    paste0("  'gamma1' = ", at(3L), " is at the ceiling ",
           format(sn_max_skew()), " of the skewness that a skew normal\n",
           "  can carry, where the shape to which it maps is not finite")
  } else if (!is.finite(sg[i]) || sg[i] <= 0) {
    paste0("  'sigma' = ", at(2L), " and 'gamma1' = ", at(3L),
           " imply a scale that is not finite")
  } else {
    paste0("  'mu' = ", at(1L), ", 'sigma' = ", at(2L),
           " and 'gamma1' = ", at(3L),
           " imply a location that is not finite")
  }
  stop("Invalid parameter value(s) for the 'skew normal2' distribution:\n",
       cause, ".", call. = FALSE)
}

#' @title Skew Normal Density in the Centered Parametrization
#' @name distrib_pdf.SkewNormal2Distrib
#'
#' @description
#' Computes the skew normal density at the direct parameters that the centered
#' triple implies. With
#' \eqn{(\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)} from
#' [sn_cp_to_dp()] and \eqn{z = (y-\xi)/\omega},
#' \deqn{f(y; \mu, \sigma, \gamma_1)
#'       = \dfrac{2}{\omega}\,\phi(z)\,\Phi(\alpha z).}
#' The density is the same function of \eqn{y} as the one that
#' [distrib_pdf.SkewNormal1Distrib()] computes; only the three numbers that
#' name it differ.
#'
#' Unlike the derivatives, the density is defined at \eqn{\gamma_1 = 0} and
#' equals the Gaussian's there: the map itself is continuous through zero, and
#' only its derivative is not.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param y A numeric vector of observations, anywhere on the real line.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `y`. `sigma` must be
#'   strictly positive and `gamma1` must lie in the open interval
#'   \eqn{(-0.9952717, 0.9952717)}; the generic signals an error for a value
#'   outside it before the map is evaluated.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma), length(gamma1))`.
#'
#' @section Notation:
#' \eqn{\mu}, \eqn{\sigma} and \eqn{\gamma_1} are the mean, the standard
#' deviation and the skewness; \eqn{\xi}, \eqn{\omega} and \eqn{\alpha} are the
#' location, scale and shape that they imply; \eqn{\phi} and \eqn{\Phi} are the
#' standard Gaussian density and distribution function.
#'
#' @seealso [distrib_pdf.SkewNormal1Distrib()] for the same density in the
#'   direct parametrization, [sn_cp_to_dp()] for the map, and [distrib_pdf()]
#'   for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' d1 <- skewnormal1_distrib()
#' y <- c(-1, 0.3, 1.7)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # The same law is reached through the map.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(d1, y, distributions7:::sn2_theta(th)))
#'
#' # It integrates to one, and its first moment is the parameter mu.
#' c(mass = integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value,
#'   mean = integrate(function(v) v * distrib_pdf(d, v, th), -Inf, Inf)$value)
#'
#' # At zero skewness the density is the Gaussian's, where the derivatives
#' # are not defined.
#' all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, gamma1 = 0)), dnorm(y))
S7::method(distrib_pdf, SkewNormal2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(skewnormal1_distrib(), y, sn2_theta(theta), log = log)
}

#' @title Skew Normal Distribution Function in the Centered Parametrization
#' @name distrib_cdf.SkewNormal2Distrib
#'
#' @description
#' Computes the skew normal distribution function at the direct parameters that
#' the centered triple implies, through Azzalini's Owen's T identity
#' \eqn{F(q) = \Phi(z) - 2T(z, \alpha)} with \eqn{z = (q-\xi)/\omega}. The
#' arithmetic is that of [distrib_cdf.SkewNormal1Distrib()]; this method
#' supplies the translated parameters.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `q`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Passed to [distrib_cdf.SkewNormal1Distrib()].
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms with `log.p = TRUE`.
#'
#' @section Notation:
#' \eqn{\gamma_1} is the skewness, \eqn{(\xi, \omega, \alpha)} are the implied
#' location, scale and shape, and \eqn{T} is Owen's T function.
#'
#' @seealso [distrib_cdf.SkewNormal1Distrib()] for the identity,
#'   [distrib_quantile.SkewNormal2Distrib()] for its inverse, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' q <- c(-2, -0.5, 0.5, 2)
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # The distribution function is compared with a direct quadrature of the
#' # density.
#' rbind(owen = distrib_cdf(d, q, th),
#'       quadrature = vapply(q, function(u)
#'         integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#'
#' # At zero skewness it is the Gaussian's.
#' all.equal(distrib_cdf(d, q, list(mu = 0, sigma = 1, gamma1 = 0)), pnorm(q))
#'
#' # A positive skewness puts more than half the mass below the mean.
#' distrib_cdf(d, 0, th)
S7::method(distrib_cdf, SkewNormal2Distrib) <- function(distrib, q, theta,
                                                         lower.tail = TRUE,
                                                         log.p = FALSE, ...) {
  distrib_cdf(skewnormal1_distrib(), q, sn2_theta(theta),
              lower.tail = lower.tail, log.p = log.p, ...)
}

#' @title Skew Normal Quantile Function in the Centered Parametrization
#' @name distrib_quantile.SkewNormal2Distrib
#'
#' @description
#' Computes the quantiles of the skew normal at the direct parameters that the
#' centered triple implies. The skew normal has no closed-form quantile
#' function, so the value comes from [continuous_distrib()]'s root finding on
#' the distribution function, reached through
#' [distrib_quantile.continuous_distrib()].
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms when `log.p = TRUE`.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of the length of `p`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Passed to [distrib_quantile.continuous_distrib()], which does not
#'   use them.
#'
#' @return A numeric vector of quantiles, of the length of the recycled inputs.
#'
#' @seealso [distrib_cdf.SkewNormal2Distrib()], which it inverts,
#'   [distrib_rng.SkewNormal2Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # A round trip through the distribution function recovers the probabilities.
#' p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' q <- distrib_quantile(d, p, th)
#' rbind(quantile = q, back = distrib_cdf(d, q, th))
#'
#' # A right-skewed density has its median below its mean.
#' c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
S7::method(distrib_quantile, SkewNormal2Distrib) <- function(distrib, p, theta,
                                                              lower.tail = TRUE,
                                                              log.p = FALSE, ...) {
  distrib_quantile(skewnormal1_distrib(), p, sn2_theta(theta),
                   lower.tail = lower.tail, log.p = log.p, ...)
}

#' @title Skew Normal Random Generation in the Centered Parametrization
#' @name distrib_rng.SkewNormal2Distrib
#'
#' @description
#' Draws from the skew normal at the direct parameters that the centered triple
#' implies, through [distrib_rng.SkewNormal1Distrib()]'s stochastic
#' representation. The draws are exact and cost two `rnorm` calls, whatever the
#' skewness is.
#'
#' @param distrib A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, each a
#'   numeric vector of length 1 or of length `n`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws, whose first three sample moments
#'   estimate `mu`, `sigma^2` and `gamma1`.
#'
#' @seealso [distrib_rng.SkewNormal1Distrib()] for the representation,
#'   [mean.SkewNormal2Distrib()] for the moments that the draws reproduce, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 3, sigma = 2, gamma1 = 0.6)
#'
#' set.seed(2)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # The three parameters are the three sample moments.
#' rbind(sample = c(mean(x), sd(x), mean((x - mean(x))^3) / sd(x)^3),
#'       parameter = c(th$mu, th$sigma, th$gamma1))
S7::method(distrib_rng, SkewNormal2Distrib) <- function(distrib, n, theta, ...) {
  distrib_rng(skewnormal1_distrib(), n, sn2_theta(theta))
}

#' @title The Skewness Below Which the Expected Information Is a Series
#' @name sn2_ge
#' @description Returns the bound on \eqn{|\gamma_1|} below which
#'   [distrib_expected_hessian.SkewNormal2Distrib()] and its derivatives come
#'   from the series in \eqn{r = (\gamma_1/c)^{1/3}} rather than from
#'   quadrature. It matches `SN2_GE` in `src/pt_skewnormal2.h`.
#' @return A single number.
#' @keywords internal
sn2_ge <- function() 3e-3

#' @title Reject the Zero Skewness Where a Derivative Diverges
#' @name sn2_reject_zero
#' @description Signals an error when any skewness is exactly zero. The
#'   log-density's series in \eqn{r} has nonzero \eqn{r^4} and \eqn{r^5} terms,
#'   so every observed derivative of order two or more in \eqn{\gamma_1}, and
#'   the derivatives of the expected information in \eqn{\gamma_1}, grow like a
#'   negative power of \eqn{\gamma_1} and are infinite at zero.
#' @param theta An aligned parameter list.
#' @param what A short description of the quantity, for the message.
#' @return `NULL`, invisibly, when every skewness differs from zero.
#' @keywords internal
sn2_reject_zero <- function(theta, what) {
  if (any(theta[[3L]] == 0)) {
    stop(paste0(
      "The ", what, " of the centered skew normal is infinite at zero\n",
      "  skewness: observation by observation the log-density and the\n",
      "  distribution function carry a term in gamma1^(4/3), whose second\n",
      "  derivative grows like gamma1^(-2/3). The score, the expected\n",
      "  information and the first derivatives of the distribution function\n",
      "  are finite there.\n",
      "  skewnormal1_distrib() carries the same family in the direct\n",
      "  parametrization, whose derivatives at alpha = 0 are ordinary numbers."),
      call. = FALSE)
  }
  invisible(NULL)
}

#' @title Skew Normal Derivatives in the Centered Parametrization
#' @name distrib_gradient.SkewNormal2Distrib
#' @aliases distrib_hessian.SkewNormal2Distrib distrib_deriv3.SkewNormal2Distrib
#'   distrib_deriv4.SkewNormal2Distrib distrib_deriv5.SkewNormal2Distrib
#'
#' @description
#' Return the derivatives of the log-density in \eqn{(\mu, \sigma, \gamma_1)}
#' of orders one to five. Each order is computed by its own compiled kernel.
#'
#' @details
#' With \eqn{w = (y - \mu)/\sigma}, \eqn{c = (4 - \pi)/2} and
#' \eqn{r = \mathrm{sign}(\gamma_1)(|\gamma_1|/c)^{1/3}}, the log-density is
#' \deqn{\ell = \tfrac12\log\frac{2}{\pi} - \log\sigma - \tfrac12\log(1 + r^2) - \tfrac12 z^2 + \log\Phi(x),
#'   \quad z = \frac{w + r}{\sqrt{1 + r^2}}, \quad
#'   x = \frac{r\,z}{\sqrt{b^2 - (1 - b^2) r^2}},}
#' with \eqn{b = \sqrt{2/\pi}}, and a derivative in \eqn{\gamma_1} is
#' \eqn{(3 c r^2)^{-1}\partial_r}. Every component is a combination of the
#' derivatives of \eqn{F = \ell + \log\sigma} in \eqn{w} and \eqn{\gamma_1},
#' derived offline. Their closed forms in \eqn{x}, \eqn{w}, \eqn{r} and
#' \eqn{\phi(x)/\Phi(x)} cancel terms of order \eqn{r^{-2k}} as
#' \eqn{\gamma_1 \to 0}, so where \eqn{|x| < 0.4} and \eqn{|r| < 0.4} the
#' kernels take the derivatives in \eqn{\gamma_1} from the series
#' \eqn{F = \sum_{n \ge 3} F_n(w) r^n}, whose coefficients are polynomials in
#' \eqn{w} computed offline at 60 digits. Far below the step of the skewing
#' function, where \eqn{x < -1} and outside the region of the series, the
#' closed forms in \eqn{\phi(x)/\Phi(x)} cancel as well, and the kernels
#' compose the derivatives instead: those of \eqn{1/\sqrt{b^2 - (1-b^2)r^2}},
#' \eqn{(1+r^2)^{-1/2}} and \eqn{\log(1+r^2)} by Faa di Bruno in \eqn{r}, then
#' in \eqn{\gamma_1} through \eqn{r(\gamma_1)}, and the outer functions
#' \eqn{-u^2/2} and \eqn{\log\Phi} with \eqn{z} and \eqn{x}, which are linear
#' in \eqn{w}; the derivatives of \eqn{\log\Phi} come from the cumulant route
#' of [distrib_deriv3.SkewNormal1Distrib()], which does not cancel there.
#'
#' The series has no \eqn{r} and no \eqn{r^2} term, so the score is finite at
#' \eqn{\gamma_1 = 0}; its \eqn{r^4} and \eqn{r^5} terms do not vanish, so the
#' derivatives of order two or more in \eqn{\gamma_1} diverge there and the
#' methods of order two and more signal an error at zero skewness.
#'
#' @param distrib A `SkewNormal2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`.
#' @param scale `"parameter"` or `"link"`. At orders one to four the generic
#'   applies the link scale after dispatch; the fifth-order method applies it
#'   itself, through [deriv5_scale()].
#' @param expected Logical; for orders three and four, whether the expected
#'   derivative is returned, by [expected_derivative()].
#' @param approx,nsim Passed to [expected_derivative()] when `expected` is
#'   `TRUE`.
#' @param ... Unused.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A named list with one numeric vector per component: 3, 6, 10, 15
#'   and 21 components at orders one to five.
#'
#' @seealso [distrib_expected_hessian.SkewNormal2Distrib()].
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 0.2, sigma = 1.3, gamma1 = 0.4)
#' y <- c(-1.7, 0.3, 2.4)
#' distrib_gradient(d, y, th)
#' # the score is finite at zero skewness
#' distrib_gradient(d, y, list(mu = 0.2, sigma = 1.3, gamma1 = 0))$gamma1
S7::method(distrib_gradient, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"),
                                                             ..., threads = 1L) {
  skewnormal2_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

S7::method(distrib_hessian, SkewNormal2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ..., threads = 1L) {
  sn2_reject_zero(theta, "observed Hessian")
  skewnormal2_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

S7::method(distrib_deriv3, SkewNormal2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
    threads = 1L) {
  sn2_reject_zero(theta, "third derivative")
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    skewnormal2_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

S7::method(distrib_deriv4, SkewNormal2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
    threads = 1L) {
  sn2_reject_zero(theta, "fourth derivative")
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    skewnormal2_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  }
}

S7::method(distrib_deriv5, SkewNormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  sn2_reject_zero(theta, "fifth derivative")
  deriv5_scale(distrib, y, theta,
               skewnormal2_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads),
               match.arg(scale))
}


#' @title Skew Normal Expected Information in the Centered Parametrization
#' @name distrib_expected_hessian.SkewNormal2Distrib
#' @aliases distrib_dexpected_hessian.SkewNormal2Distrib
#'   distrib_d2expected_hessian.SkewNormal2Distrib
#'
#' @description
#' Return the expected information and its first two derivatives in the
#' parameters. For \eqn{|\gamma_1| <} [sn2_ge()] they come from their series in
#' \eqn{r = (\gamma_1/c)^{1/3}}; elsewhere they come from the quadrature of
#' [loc_scale_expected()] over the family's own observed derivatives, because
#' the family is location-scale in \eqn{(\mu, \sigma)} at fixed
#' \eqn{\gamma_1}. The quadrature is centered at the step of the skewing
#' function, \eqn{w = -r} ([loc_scale_center()]), which at a skewness near the
#' ceiling is too sharp for nodes clustered at zero.
#'
#' @details
#' The series is the series in \eqn{r} of the observed components integrated
#' term by term against the series of the density, with gaussian moments,
#' computed offline at 60 digits. It is asymptotic rather than convergent,
#' and below the bound its terms fall under \eqn{10^{-20}} before they turn.
#' The information is finite at \eqn{\gamma_1 = 0}, where it is
#' \eqn{\mathrm{diag}(1/\sigma^2, 2/\sigma^2, 1/6)} with the sign of a Hessian,
#' but it is
#' not analytic in \eqn{\gamma_1} there: \eqn{E[\ell_{\gamma_1\gamma_1}]}
#' carries a term in \eqn{\gamma_1^{2/3}} and \eqn{E[\ell_{\mu\gamma_1}]} one
#' in \eqn{\gamma_1^{4/3}}. Its derivatives in \eqn{\gamma_1} are therefore
#' infinite at zero skewness, where the two derivative methods signal an
#' error.
#'
#' @param distrib A `SkewNormal2Distrib` object.
#' @param y A numeric vector of observations, read for its length.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Accepted for the generics' signatures and unused.
#' @param ... Unused.
#' @param threads A single positive integer, passed to the kernels behind the
#'   quadrature.
#'
#' @return A named list keyed as [hess_names()], [dexpected_names()] or
#'   [d2expected_names()].
#'
#' @seealso [loc_scale_expected()], [sn2_ge()].
#'
#' @examples
#' d <- skewnormal2_distrib()
#' # The information is finite at zero skewness, where that of the direct
#' # parametrization is singular.
#' distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, gamma1 = 0))
S7::method(distrib_expected_hessian, SkewNormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  sn2_expected_parts(distrib, theta, 0L, length(y), threads)
}

S7::method(distrib_dexpected_hessian, SkewNormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  sn2_reject_zero(theta, "derivative of the expected information")
  dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                     function(k) sn2_expected_parts(distrib, theta, k, length(y),
                                                    threads))
}

S7::method(distrib_d2expected_hessian, SkewNormal2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
    threads = 1L) {
  sn2_reject_zero(theta, "second derivative of the expected information")
  dexpected_analytic(distrib, y, theta, match.arg(scale), 2L, threads,
                     function(k) sn2_expected_parts(distrib, theta, k, length(y),
                                                    threads))
}

#' @title The Expected Information of the Centered Skew Normal, by Region
#' @name sn2_expected_parts
#' @description Returns the expected information (`k = 0`) or its derivatives
#'   (`k = 1, 2`), from the series kernels where \eqn{|\gamma_1| <}
#'   [sn2_ge()] and from [loc_scale_expected()] elsewhere. The region is chosen
#'   for each observation from its own skewness.
#' @param distrib A `SkewNormal2Distrib` object.
#' @param theta The parameters.
#' @param k The order: 0, 1 or 2.
#' @param n The number of observations.
#' @param threads Passed to [loc_scale_expected()].
#' @return A named list of numeric vectors of length `n`.
#' @keywords internal
sn2_expected_parts <- function(distrib, theta, k, n, threads = 1L) {
  theta <- align_theta(distrib, theta)
  sg <- rep_len(theta[[2L]], n)
  g <- rep_len(theta[[3L]], n)
  small <- abs(g) < sn2_ge()
  ser_fn <- switch(k + 1L, skewnormal2_expected_series_cpp,
                   skewnormal2_dexpected1_series_cpp,
                   skewnormal2_dexpected2_series_cpp)
  if (all(small)) return(ser_fn(sg, g))
  if (!any(small)) return(loc_scale_expected(distrib, theta, k, n, threads))
  out_q <- loc_scale_expected(
    distrib, lapply(theta, function(v) rep_len(v, n)[!small]), k, sum(!small), threads)
  out_s <- ser_fn(sg[small], g[small])
  out <- lapply(names(out_q), function(nm) {
    v <- numeric(n)
    v[!small] <- out_q[[nm]]
    if (!is.null(out_s[[nm]])) v[small] <- out_s[[nm]]
    v
  })
  names(out) <- names(out_q)
  out
}

S7::method(expected_hessian_costly, SkewNormal2Distrib) <- function(x, ...) TRUE


#' @title Skew Normal Derivatives in the Response, Centered Parametrization
#' @name distrib_grad_y.SkewNormal2Distrib
#' @aliases distrib_hess_y.SkewNormal2Distrib distrib_cross2_y.SkewNormal2Distrib
#'   distrib_grad_y_hess.SkewNormal2Distrib distrib_hess_y_hess.SkewNormal2Distrib
#'
#' @description
#' Return the first and second derivatives of the log-density in the response,
#' and the mixed derivatives of orders one and two in the response and one and
#' two in \eqn{(\mu, \sigma, \gamma_1)}. Each derivative comes from its own
#' compiled kernel, written in the same way as the parameter derivatives (see
#' [distrib_gradient.SkewNormal2Distrib()]).
#'
#' @param distrib A `SkewNormal2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `gamma1`.
#' @param scale `"parameter"` or `"link"`, for the mixed derivatives.
#' @param ... Unused.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Defaults to `1L`.
#'
#' @return A numeric vector for the derivatives in the response; a named list
#'   for the mixed derivatives.
#'
#' @seealso [distrib_gradient.SkewNormal2Distrib()].
#'
#' @examples
#' d <- skewnormal2_distrib()
#' distrib_grad_y(d, c(-1.7, 0.3), list(mu = 0.2, sigma = 1.3, gamma1 = 0.4))
S7::method(distrib_grad_y, SkewNormal2Distrib) <- function(distrib, y, theta, ...,
                                                           threads = 1L) {
  skewnormal2_dy1_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)$y
}

S7::method(distrib_hess_y, SkewNormal2Distrib) <- function(distrib, y, theta, ...,
                                                           threads = 1L) {
  skewnormal2_dy2_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)$y
}


#' @title Mean of the Skew Normal in the Centered Parametrization
#' @name mean.SkewNormal2Distrib
#'
#' @description
#' Returns \eqn{\mu}, the first parameter, which in this parametrization is the
#' mean by construction. The method adds [moment_const()] to it, which recycles
#' the value to the length that the three parameters imply, so a `theta` whose
#' components vary by observation yields one mean per observation.
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of the length that the recycled
#'   parameters imply.
#'
#' @seealso [variance.SkewNormal2Distrib()] and
#'   [skewness.SkewNormal2Distrib()], the other two parameters read back;
#'   [kurtosis.SkewNormal2Distrib()], which is not a parameter; and
#'   [mean.SkewNormal1Distrib()] for the same quantity in the direct
#'   parametrization, where it is not.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' mean(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#'
#' # The result has one value per observation when a parameter varies.
#' mean(d, list(mu = c(0, 3, 7), sigma = 2, gamma1 = 0.6))
#'
#' @keywords internal
S7::method(mean, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 3L, 0) + theta[[1]]
}

#' @title Variance of the Skew Normal in the Centered Parametrization
#' @name variance.SkewNormal2Distrib
#'
#' @description
#' Returns \eqn{\sigma^2}, the square of the second parameter, which in this
#' parametrization is the variance by construction. [moment_const()] recycles
#' the result to the length that the three parameters imply.
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, of the length that the recycled
#'   parameters imply.
#'
#' @seealso [mean.SkewNormal2Distrib()] and [skewness.SkewNormal2Distrib()],
#'   and [variance.SkewNormal1Distrib()], where the scale is not the standard
#'   deviation.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' variance(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#'
#' # It agrees with a quadrature of the second central moment.
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' c(parameter = variance(d, th),
#'   quadrature = integrate(function(v) v^2 * distrib_pdf(d, v, th),
#'                          -Inf, Inf)$value)
#'
#' @keywords internal
S7::method(variance, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 3L, 0)
}

#' @title Skewness of the Skew Normal in the Centered Parametrization
#' @name skewness.SkewNormal2Distrib
#'
#' @description
#' Returns \eqn{\gamma_1}, the third parameter, which in this parametrization
#' is the standardized third central moment by construction. [moment_const()]
#' recycles the result to the length that the three parameters imply.
#'
#' The value cannot leave \eqn{(-0.9952717, 0.9952717)}, because the
#' constructor bounds the parameter at the supremum that the family reaches;
#' see [sn_max_skew()].
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of skewnesses, of the length that the recycled
#'   parameters imply.
#'
#' @seealso [sn_max_skew()] for the bound,
#'   [skewness.SkewNormal1Distrib()] for the same quantity computed from a
#'   shape, and [kurtosis.SkewNormal2Distrib()], which is not free.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' skewness(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#'
#' # It agrees with a quadrature of the standardized third central moment.
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#' c(parameter = skewness(d, th),
#'   quadrature = integrate(function(v) v^3 * distrib_pdf(d, v, th),
#'                          -Inf, Inf)$value)
#'
#' @keywords internal
S7::method(skewness, SkewNormal2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  moment_const(theta, 3L, 0) + theta[[3]]
}

#' @title Kurtosis of the Skew Normal in the Centered Parametrization
#' @name kurtosis.SkewNormal2Distrib
#'
#' @description
#' Returns the excess kurtosis, computed from
#' [kurtosis.SkewNormal1Distrib()] at the implied direct parameters. It is the
#' one moment that this parametrization does not name: once the first three
#' moments are fixed, all three parameters are used up, so the fourth follows
#' from them.
#'
#' The consequence for use is that a skew normal cannot match an arbitrary
#' first four moments. At \eqn{\gamma_1 = 0.5} the excess kurtosis is 0.347,
#' and it is determined by the skewness alone.
#'
#' @param x A `SkewNormal2Distrib` object, from [skewnormal2_distrib()].
#' @param theta A named list with components `mu`, `sigma` and `gamma1`, in any
#'   order; it is aligned here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of excess kurtoses, of the length that the recycled
#'   parameters imply.
#'
#' @seealso [skewness.SkewNormal2Distrib()], which does name a parameter;
#'   [kurtosis.SkewNormal1Distrib()] for the closed form; and
#'   [skewt_distrib()], whose degrees of freedom give the fourth moment a
#'   parameter of its own.
#'
#' @examples
#' d <- skewnormal2_distrib()
#'
#' # The kurtosis depends on the skewness only: it does not change when the
#' # mean and the standard deviation change.
#' c(a = kurtosis(d, list(mu = 0, sigma = 1, gamma1 = 0.5)),
#'   b = kurtosis(d, list(mu = 9, sigma = 4, gamma1 = 0.5)))
#'
#' # It rises with the skewness, and is zero at symmetry.
#' vapply(c(0.001, 0.3, 0.6, 0.9),
#'        function(g) kurtosis(d, list(mu = 0, sigma = 1, gamma1 = g)), 0)
#'
#' @keywords internal
S7::method(kurtosis, SkewNormal2Distrib) <- function(x, theta, ...) {
  kurtosis(skewnormal1_distrib(), sn2_theta(align_theta(x, theta)))
}


#' @title Skew Normal Distribution Object, Centered Parametrization
#'
#' @description
#' Builds a skew normal distribution object parametrized by its mean
#' \eqn{\mu}, its standard deviation \eqn{\sigma} and its skewness
#' \eqn{\gamma_1}. It is the same family as [skewnormal1_distrib()], named by
#' three moments instead of by a location, a scale and a shape.
#'
#' The parametrization is Azzalini's. Its expected information stays
#' non-singular at \eqn{\gamma_1 = 0}, where the expected information of the
#' direct parametrization loses a rank. The map between the two runs through a
#' cube root, so the derivative of the shape in \eqn{\gamma_1} is unbounded at
#' zero skewness and the log-density has no derivative of order two or more in
#' \eqn{\gamma_1} there.
#'
#' @param link_mu A `linkfunctions7` link object for the mean, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_sigma A link object for the standard deviation, which must be
#'   strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_gamma1 A link object for the skewness, which must lie strictly
#'   inside \eqn{(-0.9952717, 0.9952717)}. Defaults to
#'   [linkfunctions7::bounded_link()] over exactly that interval, so any real
#'   linear predictor maps to a skewness that the family can reach.
#'
#' @details
#' # Relation to reparametrize()
#'
#' The map passes through
#' \eqn{c = \mathrm{sign}(\gamma_1)(2|\gamma_1|/(4-\pi))^{1/3}}, which carries
#' a sign. For this reason [sn2_theta()] takes the sign from the plain value of
#' `gamma1` and passes it to [sn_cp_to_dp()] as an argument, and the body of
#' the map contains no `abs()`. The derivatives of the log-density in
#' the centered parameters are written out per order in compiled kernels (see
#' [distrib_gradient.SkewNormal2Distrib()]); near \eqn{\gamma_1 = 0}, where
#' their closed forms cancel, they come from the series of the log-density in
#' \eqn{\gamma_1^{1/3}}, and far below the step of the skewing function from a
#' composition that does not cancel.
#'
#' # Behavior near zero skewness
#'
#' \eqn{\partial\alpha/\partial\gamma_1} grows without bound as
#' \eqn{\gamma_1 \to 0}, like \eqn{\gamma_1^{-2/3}}. The score does not follow
#' it: the divergent
#' contributions cancel, so the information in \eqn{\gamma_1} tends to
#' \eqn{1/6} and the whole matrix stays positive definite at symmetry.
#'
#' The observed curvature does diverge, at the rate that the cube root sets:
#' \eqn{\gamma_1^{-2/3}}. The derivatives of order two and more are
#' therefore rejected at \eqn{\gamma_1 = 0} exactly, with a message naming the
#' cause; the score and the expected information are returned there. The density,
#' the distribution function, the quantile function, the generator and both
#' response derivatives are fine there and equal the Gaussian's.
#'
#' For \eqn{|\gamma_1|} below [sn2_ge()] the expected information is computed
#' from series in \eqn{\gamma_1^{1/3}}, which contain no cancellation. A
#' problem that is symmetric by construction is better posed in
#' [skewnormal1_distrib()], since the observed derivatives of order two and
#' more in the centered parametrization are rejected at \eqn{\gamma_1 = 0}.
#'
#' # The bound on the skewness
#'
#' A skew normal cannot reach \eqn{|\gamma_1| > 0.9952717} whatever its shape,
#' so `gamma1` is bounded there and carries a bounded link by default. A model
#' that needs a larger skewness can use [skewt_distrib()].
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (-\infty, \infty)}
#' - \eqn{\sigma \in (0, \infty)}
#' - \eqn{\gamma_1 \in (-0.9952717, 0.9952717)}
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{2}{\omega}\,\phi\!\left(\frac{y-\xi}{\omega}\right)\Phi\!\left(\alpha\,\frac{y-\xi}{\omega}\right), \qquad (\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)}
#' on \eqn{y \in \mathbb{R}}, with \eqn{\mathrm{DP}} the map of
#' [sn_cp_to_dp()].
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2},
#'       \qquad \gamma_1(Y) = \gamma_1.}
#'
#' @return An S7 object of class [SkewNormal2Distrib], inheriting from
#'   `continuous_distrib`. Its `params` are `c("mu", "sigma", "gamma1")`, its
#'   `bounds` `c(-Inf, Inf)`, and its `link_params` the three links given here.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
#' Families*. Cambridge University Press. The centered parametrization is
#' section 3.1.4.
#'
#' @seealso [skewnormal1_distrib()] for the direct parametrization,
#'   [sn_cp_to_dp()] for the map between them,
#'   [skewt_distrib()] for a family that reaches a larger skewness, and
#'   [SkewNormal2Distrib] for the class and its method list.
#'
#' @examples
#' d <- skewnormal2_distrib()
#' th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
#'
#' # All three parameters are moments, which is why it is called "centered".
#' c(mean = mean(d, th), sd = sqrt(variance(d, th)),
#'   skewness = skewness(d, th))
#'
#' # The fourth moment is not free: it follows from the skewness.
#' vapply(c(0.001, 0.5, 0.9),
#'        function(g) kurtosis(d, list(mu = 0, sigma = 1, gamma1 = g)), 0)
#'
#' # The information stays invertible at symmetry, where the direct
#' # parametrization's does not.
#' smallest_eigenvalue <- function(dd, p) {
#'   e <- distrib_expected_hessian(dd, 0, p)
#'   M <- matrix(c(e[[1]], e[[4]], e[[5]],
#'                 e[[4]], e[[2]], e[[6]],
#'                 e[[5]], e[[6]], e[[3]]), 3, 3)   # hess_names() order
#'   min(eigen(-M, only.values = TRUE)$values)
#' }
#' c(centered = smallest_eigenvalue(d, list(mu = 0, sigma = 1, gamma1 = 1e-6)),
#'   direct = smallest_eigenvalue(skewnormal1_distrib(),
#'                                list(mu = 0, sigma = 1, alpha = 0)))
#'
#' # A fit recovers all three moments.
#' set.seed(11)
#' x <- distrib_rng(d, 4000, list(mu = 3, sigma = 2, gamma1 = 0.6))
#' coef(fit_distrib(d, x))
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
